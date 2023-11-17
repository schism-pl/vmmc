use crate::particle::{self, IsParticle, Particle, ParticleId, VParticle};
use crate::patchy_discs::{PatchyDiscParams, PatchyDiscsPotential};
use crate::position::Position;
use crate::simbox::SimBox;
use crate::stats::RunStats;
use crate::{params::DIMENSION, position::DimVec};
use anyhow::{anyhow, Result};
use rand::rngs::SmallRng;
use rand::Rng;
use rand_distr::num_traits::Zero;
use std::collections::{HashSet, VecDeque};
use std::f32::INFINITY;

// NOTE: I'm just assuming things aren't isotropic
// NOTE: assuming not is_repulsive
// NOTE: assuming no non-pairwise forces
// NOTE: I'm eliding the nmoving < cutoff at callsites since I check to make sure we don't recheck neighbors

// TODO: documentation!
// TODO: equation map

// TODO: in debug mode, show that average pair energy decreases (increases?) over time

#[derive(Clone, Copy, Debug)]
pub enum MoveDir {
    Forward,
    Backward,
}

impl From<MoveDir> for f64 {
    fn from(movedir: MoveDir) -> Self {
        match movedir {
            MoveDir::Forward => 1.0,
            MoveDir::Backward => -1.0,
        }
    }
}
#[derive(Debug)]
pub struct VirtualMoves {
    // particle moving, delta position, delta orientation
    inner: Vec<(ParticleId, VParticle)>,
}

impl VirtualMoves {
    fn new() -> Self {
        Self { inner: Vec::new() }
    }

    fn push(&mut self, id: ParticleId, p: VParticle) {
        self.inner.push((id, p))
    }
}

// particle move
#[derive(Debug)]
pub struct ProposedMove {
    seed_id: ParticleId,
    step_size: f64,
    is_rotation: bool,
    vec: DimVec, // unit vec
}

impl ProposedMove {
    pub fn new(seed_id: ParticleId, step_size: f64, is_rotation: bool, vec: DimVec) -> Self {
        Self {
            seed_id,
            step_size,
            is_rotation,
            vec,
        }
    }

    pub fn step_factor(&self, dir: MoveDir) -> f64 {
        f64::from(dir) * self.step_size
    }

    pub fn is_rotation(&self) -> bool {
        self.is_rotation
    }

    pub fn vec(&self) -> DimVec {
        self.vec
    }

    pub fn scaled_vec(&self, dir: MoveDir) -> DimVec {
        self.vec.scalar_mul(self.step_factor(dir))
    }
}

pub struct VmmcParams {
    prob_translate: f64,
    max_trial_translation: f64,
    max_trial_rotation: f64,
    reference_radius: f64,
}

impl VmmcParams {
    pub fn new(
        prob_translate: f64,
        max_trial_translation: f64,
        max_trial_rotation: f64,
        reference_radius: f64,
    ) -> Self {
        Self {
            prob_translate,
            max_trial_translation,
            max_trial_rotation,
            reference_radius,
        }
    }
}

pub struct Vmmc {
    particles: Vec<Particle>,
    simbox: SimBox,
    potential: PatchyDiscsPotential,
    params: VmmcParams,
}
impl Vmmc {
    pub fn new(
        particles: Vec<Particle>,
        simbox: SimBox,
        params: VmmcParams,
        pd_params: PatchyDiscParams,
    ) -> Self {
        let potential = PatchyDiscsPotential::new(pd_params);
        Self {
            particles,
            simbox,
            potential,
            params,
        }
    }

    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    pub fn simbox(&self) -> &SimBox {
        &self.simbox
    }

    pub fn potential(&self) -> &PatchyDiscsPotential {
        &self.potential
    }

    pub fn particle(&self, p_id: ParticleId) -> &Particle {
        &self.particles[p_id as usize]
    }

    // get energy
    pub fn get_particle_energy(&self, p: &Particle) -> f64 {
        let mut energy = 0.0;
        let interactions = self
            .potential
            .determine_interactions(&self.simbox, &self.particles, p);
        for &neighbor_id in interactions.iter() {
            let neighbor = self.particle(neighbor_id);
            energy += self
                .potential
                .compute_pair_energy(self.simbox(), p, neighbor);
        }
        energy
    }

    pub fn get_average_energy(&self) -> f64 {
        let mut total_energy = 0.0;
        for p in self.particles() {
            total_energy += self.get_particle_energy(p);
        }
        // divide by an extra 2 so we don't double count bonds (we should only count p0->p1 but not also p1->p0)
        total_energy / (self.particles.len() as f64 * 2.0)
    }

    // TODO: differentially test?
    fn compute_stokes_radius(&self, vmoves: &VirtualMoves, mov: &ProposedMove) -> f64 {
        let mut stokes_radius = 0.0;
        // Calculate center of mass of the moving cluster (translations only).
        let center_of_mass = self.center_of_mass(vmoves);

        for (_, vp) in vmoves.inner.iter() {
            let delta = if !mov.is_rotation {
                self.simbox.sep_in_box(vp.orig_pos(), center_of_mass)
            } else {
                self.simbox
                    .sep_in_box(vp.orig_pos(), self.particle(mov.seed_id).pos())
            };

            stokes_radius += delta.cross_prod_magnitude_sqd(mov.vec);
        }

        let scale_factor = self.params.reference_radius
            / (self.params.reference_radius + (stokes_radius / vmoves.inner.len() as f64).sqrt());
        if mov.is_rotation {
            scale_factor * scale_factor * scale_factor
        } else {
            scale_factor
        }
    }

    // get center of mass for vmov set
    // TODO: boundary conditions?
    fn center_of_mass(&self, vmoves: &VirtualMoves) -> Position {
        let mut com = Position::zeroes();
        for (_, vp) in vmoves.inner.iter() {
            com += vp.orig_pos(); //particle.pos();
        }
        com.scalar_div(vmoves.inner.len() as f64)
    }

    // returns difference in position between final and original
    fn calculate_motion(
        &self,
        particle: &Particle,
        mov: &ProposedMove,
        seed: &Particle,
        dir: MoveDir,
    ) -> VParticle {
        let mut final_p = particle.pos();
        let mut final_or = particle.or();

        if !mov.is_rotation() {
            final_p += mov.scaled_vec(dir);
        } else {
            let rel_pos = particle.pos() - seed.pos();
            final_p += rel_pos.rotated_by(mov.vec(), mov.step_factor(dir));

            final_or += final_or.rotated_by(mov.vec(), mov.step_factor(dir));
        }

        final_p = self.simbox.map_pos_into_box(final_p);
        VParticle::new(particle.pos(), particle.or(), final_p, final_or)
    }

    fn choose_random_move(&self, rng: &mut SmallRng) -> ProposedMove {
        // 1. Choose a particle that will lead the move
        let seed_id = self.choose_seed(rng);
        // 2. Choose a direction (unit vector) for the move
        let rand_vec = Position::unit_vector(rng);
        // 3. Choose a move type (translation or rotation)
        let is_rotation = rng.gen::<f64>() >= self.params.prob_translate;
        // 4. Choose a size for the move
        let step_size = if is_rotation {
            // Rotation
            let r: f64 = rng.gen();
            self.params.max_trial_rotation * r.powf(1.0 / (DIMENSION as f64))
        } else {
            // Translate
            // Scale step-size to uniformly sample unit sphere/circle.
            let r: f64 = rng.gen();
            // random number between (-1.0 and 1.0) * max_trial_translation
            self.params.max_trial_translation * (2.0 * r - 1.0)
        };

        ProposedMove::new(seed_id, step_size, is_rotation, rand_vec)
    }

    // Note: each particle gets to attempt to link to every neighbor, even if neighbor
    // has tried to be recruited before
    fn recruit_cluster(&self, rng: &mut SmallRng, mov: &ProposedMove) -> Result<VirtualMoves> {
        let seed = self.particle(mov.seed_id);
        // particles in the cluster who still need to make their linking pass
        let mut worklist = VecDeque::from([mov.seed_id]);
        // particles who have made their linking pass
        let mut seen = HashSet::new();
        let mut vmoves = VirtualMoves::new();

        while !worklist.is_empty() {
            // A new particle tries to link to its neighbors
            let id = worklist.pop_front().unwrap();
            let particle = self.particle(id);
            seen.insert(id);

            let final_p = self.calculate_motion(particle, mov, seed, MoveDir::Forward);

            // this particle moves
            vmoves.push(id, final_p.clone());

            // particle / seen / final_p
            let reverse_p = self.calculate_motion(particle, mov, seed, MoveDir::Backward);
            let interactions =
                self.potential
                    .determine_interactions(&self.simbox, &self.particles, particle);

            for &neighbor_id in interactions.iter() {
                let neighbor = self.particle(neighbor_id);

                let (link_weight, reverse_link_weight) =
                    self.compute_link_weights(particle, neighbor, &final_p, &reverse_p);

                log::debug!(
                    "p{:?} is attempting to link p{:?} with weights fwd:{:?} back:{:?}",
                    id,
                    neighbor_id,
                    link_weight,
                    reverse_link_weight
                );

                if rng.gen::<f64>() <= link_weight {
                    if rng.gen::<f64>() <= reverse_link_weight / link_weight {
                        // The neighbor has been linked.
                        // Add it to the move and attempt to recruit its neighbors
                        if !seen.contains(&neighbor_id) {
                            worklist.push_back(neighbor_id);
                        }
                    } else {
                        log::debug!(
                            "P{:?} got frustrated. Link weight = {:?} Reverse = {:?}",
                            id,
                            link_weight,
                            reverse_link_weight
                        );
                        return Err(anyhow!("Particle got frustrated"));
                    }
                }
            }
        }
        Ok(vmoves)
    }

    // equation 6 from the paper
    fn compute_link_weights(
        &self,
        particle: &Particle,
        interacting_p: &Particle,
        final_p: &VParticle,
        reverse_p: &VParticle,
    ) -> (f64, f64) {
        let init_energy = self
            .potential
            .compute_pair_energy(&self.simbox, particle, interacting_p);
        let final_energy = self
            .potential
            .compute_pair_energy(&self.simbox, final_p, interacting_p);
        let reverse_energy =
            self.potential
                .compute_pair_energy(&self.simbox, reverse_p, interacting_p);

        let link_weight = 0.0_f64.max(1.0 - (init_energy - final_energy).exp());
        let reverse_link_weight = 0.0_f64.max(1.0 - (init_energy - reverse_energy).exp());
        // if final_energy.is_infinite() || reverse_energy.is_infinite() {
        //     println!("init = {:?} final = {:?} reverse = {:?} lw = {:?} rlw = {:?} rlw/lw={:?} inf/inf = {:?}", init_energy, final_energy, reverse_energy, link_weight, reverse_link_weight, reverse_link_weight/link_weight, f64::INFINITY / f64::INFINITY);
        // }

        // if reverse_energy.is_infinite(){
        //     println!("init = {:?} final reverse = {:?} lw = {:?} rlw = {:?}", init_energy, reverse_energy, link_weight, reverse_link_weight);
        // }
        // assert!(init_energy.is_finite());
        // assert!(final_energy.is_finite());
        // assert!(reverse_energy.is_finite());

        debug_assert!(
            link_weight.is_finite() && !link_weight.is_nan() && !link_weight.is_subnormal()
        );
        debug_assert!(
            reverse_link_weight.is_finite()
                && !reverse_link_weight.is_nan()
                && !reverse_link_weight.is_subnormal()
        );
        (link_weight, reverse_link_weight)
    }

    fn choose_seed(&self, rng: &mut SmallRng) -> ParticleId {
        rng.gen_range(0..self.particles.len() as u16)
    }

    fn commit_moves(&mut self, vmoves: &VirtualMoves) {
        for (p_id, vp) in vmoves.inner.iter() {
            let p = &mut self.particles[*p_id as usize];
            self.simbox.move_particle_tenancy(*p_id, p.pos(), vp.pos());
            p.update_pos(vp.pos());
            p.update_or(vp.or());
        }
    }

    fn revert_moves(&mut self, vmoves: &VirtualMoves) {
        for (p_id, vp) in vmoves.inner.iter() {
            let p = &mut self.particles[*p_id as usize];
            self.simbox
                .move_particle_tenancy(*p_id, p.pos(), vp.orig_pos());
            p.update_pos(vp.orig_pos());
            p.update_or(vp.orig_or());
        }
    }

    fn attempt_commit(
        &mut self,
        rng: &mut SmallRng,
        mov: &ProposedMove,
        vmoves: &VirtualMoves,
    ) -> Result<()> {
        // approximate Stokes scaling factor
        let scale_factor = if vmoves.inner.len() > 1 {
            self.compute_stokes_radius(vmoves, mov)
        } else {
            1.0
        };
        // Stokes drag rejection
        if rng.gen::<f64>() > scale_factor {
            return Err(anyhow!("Stokes drag rejection"));
        }

        self.commit_moves(vmoves);

        // check for overlaps
        for (p_id, _) in vmoves.inner.iter() {
            let p = self.particle(*p_id);
            if self.overlaps(p) {
                self.revert_moves(vmoves);
                return Err(anyhow!("Overlapping particle"));
            }
        }

        Ok(())
    }

    // check if a particle overlaps any other particles
    fn overlaps(&self, p: &Particle) -> bool {
        for neighbor_id in self.simbox.get_neighbors(p) {
            let neighbor = self.particle(neighbor_id);
            if neighbor == p {
                continue;
            }
            // dist < 1.0 (hard sphere radius)
            if self.simbox.sep_in_box(p.pos(), neighbor.pos()).norm() < 1.0 {
                return true;
            }
        }
        false
    }

    fn sim_has_overlaps(&self) -> bool {
        self.particles.iter().any(|p| self.overlaps(p))
    }

    // run a bunch of checks on the vmmc object to make sure it looks good
    // meant to be used in debug_assert!
    pub fn well_formed(&self) -> bool {
        // check for particle overlaps?
        for p in self.particles.iter() {
            // check that position is in the box
            if !self.simbox.pos_in_box(p.pos()) {
                panic!("Pos not in box!: {:?}", p);
            }

            // check that magnitude of orientation is about 1.0
            // this is because it should be a unit vector
            // allow for a bit of deviation due to rounding errors
            if (p.or().norm() - 1.0).abs() >= 0.0001 {
                panic!("Orientation not normalized!: {:?} off by ", p);
            }
        }
        // check that distance between all particles is at least 1.0
        // (assuming radius is 0.5) (hard spheres are not allowed to overlap)
        if self.sim_has_overlaps() {
            panic!("Sim has overlaps!");
        }

        // check that tenency array is synced with particle positions
        for p in self.particles.iter() {
            // println!("{:?} {:?} {:?} {:?}", p.id(), self.simbox.get_cell_id(p.pos()), self.simbox.get_cell(p.pos()), p.pos() );
            if !self.simbox.in_cell(p) {
                panic!("Particle not in correct location in tenancy array");
            }
        }

        // check that every particle is in tenency array exactly once
        let mut seen = HashSet::new();
        for cell in self.simbox.cells().iter() {
            for tenant in cell.iter() {
                if *tenant != ParticleId::MAX {
                    if seen.contains(tenant) {
                        panic!("Duplicate particle in tenancy array");
                    }
                    seen.insert(*tenant);
                }
            }
        }
        if seen.len() != self.particles.len() {
            panic!("Wrong number of particles in tenancy array");
        }
        for p in self.particles.iter() {
            if !seen.contains(&p.id()) {
                panic!("Particle missing from tenancy array");
            }
        }

        true
    }

    fn step(&mut self, rng: &mut SmallRng, stats: &mut RunStats) -> Result<()> {
        debug_assert!(self.well_formed());
        stats.record_attempt();
        let mov = self.choose_random_move(rng);
        debug_assert!((mov.vec().norm() - 1.0).abs() < 0.00001);
        log::debug!("Chose a random move: {:?}", mov);
        let virtual_moves = self.recruit_cluster(rng, &mov)?;
        log::debug!("Found a promising set of moves: {:?}", virtual_moves);
        if self.attempt_commit(rng, &mov, &virtual_moves).is_ok() {
            stats.record_accept(mov.is_rotation(), virtual_moves.inner.len());
        }

        Ok(())
    }

    pub fn step_n(&mut self, n: usize, rng: &mut SmallRng) {
        let mut run_stats = RunStats::new(self.particles.len());
        for idx in 0..n {
            log::info!("Successful moves: {:?}/{:?}", run_stats.num_accepts(), idx);
            let _ = self.step(rng, &mut run_stats);
        }
        // println!("{:?}", run_stats);
    }
}
