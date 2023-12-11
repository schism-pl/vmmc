use crate::consts::PARTICLE_RADIUS;
use crate::particle::{IsParticle, Particle, ParticleId, Particles, VParticle};
use crate::patchy_discs::PatchyDiscsPotential;
use crate::position::Position;
use crate::simbox::SimBox;
use crate::stats::RunStats;
use crate::{consts::DIMENSION, position::DimVec};
use anyhow::{anyhow, Result};
use rand::rngs::SmallRng;
use rand::Rng;
use std::collections::{HashSet, VecDeque};

// NOTE: I'm just assuming things aren't isotropic
// NOTE: assuming not is_repulsive
// NOTE: assuming no non-pairwise forces

// TODO: documentation!
// TODO: equation map

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
    cluster_cutoff: usize,
    vec: DimVec, // unit vec
}

impl ProposedMove {
    pub fn new(
        seed_id: ParticleId,
        step_size: f64,
        is_rotation: bool,
        cluster_cutoff: usize,
        vec: DimVec,
    ) -> Self {
        Self {
            seed_id,
            step_size,
            is_rotation,
            cluster_cutoff,
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
    max_translation: f64,
    max_rotation: f64,
}

impl VmmcParams {
    pub fn new(prob_translate: f64, max_translation: f64, max_rotation: f64) -> Self {
        Self {
            prob_translate,
            max_translation,
            max_rotation,
        }
    }
}

pub struct Vmmc {
    simbox: SimBox,
    potential: PatchyDiscsPotential,
    params: VmmcParams,
}
impl Vmmc {
    pub fn new(simbox: SimBox, params: VmmcParams, interaction_energy: f64) -> Self {
        let potential = PatchyDiscsPotential::new(simbox.shapes(), interaction_energy);
        Self {
            simbox,
            potential,
            params,
        }
    }

    pub fn particles(&self) -> &Particles {
        self.simbox.particles()
    }

    pub fn particles_mut(&mut self) -> &mut Particles {
        self.simbox.particles_mut()
    }

    pub fn simbox(&self) -> &SimBox {
        &self.simbox
    }

    pub fn simbox_mut(&mut self) -> &mut SimBox {
        &mut self.simbox
    }

    pub fn potential(&self) -> &PatchyDiscsPotential {
        &self.potential
    }

    pub fn particle(&self, p_id: ParticleId) -> &Particle {
        self.simbox.particle(p_id)
    }

    pub fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(&self, p1: &P1, p2: &P2) -> f64 {
        self.potential.compute_pair_energy(&self.simbox, p1, p2)
    }

    pub fn determine_interactions(&self, p: &Particle) -> Vec<ParticleId> {
        self.potential.determine_interactions(&self.simbox, p)
    }

    // get energy
    pub fn get_particle_energy(&self, p: &Particle) -> f64 {
        let mut energy = 0.0;
        let interactions = self.determine_interactions(p);
        for &neighbor_id in interactions.iter() {
            let neighbor = self.particle(neighbor_id);
            energy += self.compute_pair_energy(p, neighbor);
        }
        energy
    }

    pub fn get_average_energy(&self) -> f64 {
        let mut total_energy = 0.0;
        for p in self.particles().iter() {
            total_energy += self.get_particle_energy(p);
        }
        // divide by an extra 2 so we don't double count bonds (we should only count p0->p1 but not also p1->p0)
        total_energy / (self.particles().num_particles() as f64 * 2.0)
    }

    // TODO: what is this?
    fn compute_stokes_radius(&self, vmoves: &VirtualMoves, mov: &ProposedMove) -> f64 {
        let mut stokes_radius = 0.0;

        for (_, vp) in vmoves.inner.iter() {
            let delta = if !mov.is_rotation {
                let center_of_mass = self.center_of_mass(vmoves);
                self.simbox.sep_in_box(vp.orig_pos(), center_of_mass)
            } else {
                self.simbox
                    .sep_in_box(vp.orig_pos(), self.particle(mov.seed_id).pos())
            };

            stokes_radius += delta.cross_prod_magnitude_sqd(mov.vec);
        }
        stokes_radius /= vmoves.inner.len() as f64;
        stokes_radius = stokes_radius.sqrt();

        let scale_factor = PARTICLE_RADIUS / (PARTICLE_RADIUS + stokes_radius);
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
    // TODO: is this function right? it reads wierd
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
        VParticle::new(
            particle.pos(),
            particle.or(),
            final_p,
            final_or,
            particle.shape_id(),
        )
    }

    fn choose_random_move(&self, rng: &mut SmallRng) -> ProposedMove {
        // 1. Choose a particle that will lead the move
        let seed_id = self.choose_random_p_id(rng);
        // 2. Choose a direction (unit vector) for the move
        let rand_vec = Position::unit_vector(rng);
        // 3. Choose a move type (translation or rotation)
        let is_rotation = rng.gen::<f64>() >= self.params.prob_translate;
        // 4. Pick a cluster size cutoff
        // I don't know why we use this distribution, but it has something to do with particle choice fairness
        let cluster_cutoff = (1.0 / rng.gen::<f64>()).floor() as usize;
        // 5. Choose a size for the move
        let step_size = if is_rotation {
            // Rotation
            let r: f64 = rng.gen();
            self.params.max_rotation * r.powf(1.0 / (DIMENSION as f64))
        } else {
            // Translate
            // Scale step-size to uniformly sample unit sphere/circle.
            let r: f64 = rng.gen();
            // random number between (-1.0 and 1.0) * max_trial_translation
            self.params.max_translation * (2.0 * r - 1.0)
        };

        ProposedMove::new(seed_id, step_size, is_rotation, cluster_cutoff, rand_vec)
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

        let mut cluster_size = 1;

        while !worklist.is_empty() {
            if cluster_size > mov.cluster_cutoff {
                return Err(anyhow!("Cluster size surpassed chosen cutoff"));
            }
            cluster_size += 1;
            // A new particle tries to link to its neighbors
            let id = worklist.pop_front().unwrap();

            // in case we get duplicate entries in worklist
            if seen.contains(&id) {
                continue;
            }
            let particle = self.particle(id);

            seen.insert(id);

            let final_p = self.calculate_motion(particle, mov, seed, MoveDir::Forward);

            // this particle moves
            vmoves.push(id, final_p.clone());

            let reverse_p = self.calculate_motion(particle, mov, seed, MoveDir::Backward);
            let interactions = self.determine_interactions(particle);

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
        let init_energy = self.compute_pair_energy(particle, interacting_p);
        let final_energy = self.compute_pair_energy(final_p, interacting_p);
        let reverse_energy = self.compute_pair_energy(reverse_p, interacting_p);

        let link_weight = 0.0_f64.max(1.0 - (init_energy - final_energy).exp());
        let reverse_link_weight = 0.0_f64.max(1.0 - (init_energy - reverse_energy).exp());

        (link_weight, reverse_link_weight)
    }

    pub fn choose_random_p_id(&self, rng: &mut SmallRng) -> ParticleId {
        assert_ne!(self.particles().num_particles(), 0);
        loop {
            let p_id = rng.gen_range(0..self.particles().num_particles() as u16);
            if self.particles().is_active_particle(p_id) {
                return p_id;
            }
        }
    }

    fn commit_moves(&mut self, vmoves: &VirtualMoves) {
        for (p_id, vp) in vmoves.inner.iter() {
            let old_pos = self.simbox.particles().particle(*p_id).pos();
            self.simbox.move_particle_tenancy(*p_id, old_pos, vp.pos());
            let p = self.simbox.particles_mut().particle_mut(*p_id);
            p.update_pos(vp.pos());
            p.update_or(vp.or());
        }
    }

    fn revert_moves(&mut self, vmoves: &VirtualMoves) {
        for (p_id, vp) in vmoves.inner.iter() {
            let old_pos = self.simbox.particles().particle(*p_id).pos();
            self.simbox
                .move_particle_tenancy(*p_id, old_pos, vp.orig_pos());
            let p = &mut self.simbox.particles_mut().particle_mut(*p_id);
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
        self.particles().iter().any(|p| self.overlaps(p))
    }

    // run a bunch of checks on the vmmc object to make sure it looks good
    // meant to be used in debug_assert!
    pub fn well_formed(&self) -> bool {
        // check for particle overlaps?
        for p in self.particles().iter() {
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
        for p in self.particles().iter() {
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
        if seen.len() != self.particles().num_particles() {
            panic!("Wrong number of particles in tenancy array");
        }
        for p in self.particles().iter() {
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
        // if virtual_moves.inner.len() >= self.particles().num_particles() {
        //     println!(
        //         "{:?}",
        //         virtual_moves
        //             .inner
        //             .iter()
        //             .map(|(a, _)| *a)
        //             .collect::<Vec<ParticleId>>()
        //     );
        //     panic!("oh no");
        // }
        if self.attempt_commit(rng, &mov, &virtual_moves).is_ok() {
            stats.record_accept();
        }

        Ok(())
    }

    pub fn step_n(&mut self, n: usize, rng: &mut SmallRng) -> RunStats {
        let mut run_stats = RunStats::new();
        for idx in 0..n {
            log::info!("Successful moves: {:?}/{:?}", run_stats.num_accepts(), idx);
            let _ = self.step(rng, &mut run_stats);
        }
        run_stats
    }

    pub fn set_interaction_energy(&mut self, interaction_energy: f64) {
        self.potential.set_interaction_energy(interaction_energy)
    }

    // pub fn step_n(&mut self, n: usize, rng: &mut SmallRng) -> RunStats {
    //     let mut run_stats = RunStats::new(self.particles().len());
    //     for idx in 0..n {
    //         log::info!("Successful moves: {:?}/{:?}", run_stats.num_accepts(), idx);
    //         let _ = self.step(rng, &mut run_stats);
    //     }
    //     run_stats
    // }
}
