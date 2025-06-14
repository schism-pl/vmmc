use crate::consts::{MAX_ROTATION, MAX_TRANSLATION, PARTICLE_RADIUS, PROB_TRANSLATE};
use crate::particle::{IsParticle, Particle, ParticleId, Particles, VParticle};
use crate::position::DimVec;
use crate::position::{random_unit_vec, Position};
use crate::potentials::patchy_discs::PatchyDiscsPotential;
use crate::simbox::SimBox;
use crate::stats::RunStats;
use crate::Prng;
use anyhow::{anyhow, Result};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::collections::{HashSet, VecDeque};
use std::mem::size_of_val;

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
        self.vec.scale_by(self.step_factor(dir))
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Vmmc {
    simbox: SimBox,
    potential: PatchyDiscsPotential,
    dynamic_particle_count: bool,
}
impl Vmmc {
    pub fn new(simbox: SimBox, interaction_energy: f64, dynamic_particle_count: bool) -> Self {
        let potential = PatchyDiscsPotential::new(interaction_energy);
        Self {
            simbox,
            potential,
            dynamic_particle_count,
        }
    }

    pub fn particles(&self) -> &Particles {
        self.simbox.particles()
    }

    pub fn dynamic_particle_count(&self) -> bool {
        self.dynamic_particle_count
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

    pub fn particle_mut(&mut self, p_id: ParticleId) -> &mut Particle {
        self.simbox.particle_mut(p_id)
    }

    pub fn needed_mem(&self) -> usize {
        let particles_mem = self.particles().needed_mem();
        let tenancy_mem = size_of_val(self.simbox.cells());
        let shapes_mem = size_of_val(self.simbox.shapes());
        particles_mem + tenancy_mem + shapes_mem
    }

    pub fn max_needed_mem(&self) -> usize {
        let particles_mem = self.particles().max_needed_mem();
        let tenancy_mem = size_of_val(self.simbox.cells());
        let shapes_mem = size_of_val(self.simbox.shapes());
        particles_mem + tenancy_mem + shapes_mem
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

            stokes_radius += delta.cross_prod_sqd(mov.vec);
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
        com.div_by(vmoves.inner.len() as f64)
    }

    // returns difference in position between final and original
    // Invariant: not overlapping before implies not overlapping after
    // Invariant, distance between root and position should be same before and after
    // TODO: simplify
    fn calculate_motion(
        &self,
        particle: &Particle,
        mov: &ProposedMove,
        seed: &Particle,
        roots: &Option<(VParticle, VParticle)>,
        dir: MoveDir,
    ) -> VParticle {
        let mut final_p = particle.pos();
        let mut final_or = particle.or();

        if !mov.is_rotation() {
            final_p += mov.scaled_vec(dir);
        } else {
            if let Some((fwd_root, rev_root)) = roots {
                // Rotating (not the seed)
                assert_eq!(fwd_root.orig_pos(), rev_root.orig_pos());
                let displacement = self
                    .simbox()
                    .sep_in_box(particle.pos(), fwd_root.orig_pos());
                let theta = mov.step_factor(dir);
                let delta = displacement.rotated_by(theta); // magnitiude should be same before and after?
                                                            // println!("displacement = {:?} delta = {:?} theta = {:?}", displacement, delta, theta);
                                                            // assert_eq!(delta.l2_norm(), displacement.l2_norm());
                final_or = final_or.rotated_by(theta);
                if let MoveDir::Forward = dir {
                    final_p = fwd_root.pos() + delta;
                } else {
                    final_p = rev_root.pos() + delta;
                }
            } else {
                // We are rotating the seed!
                // No update to position, it just turns
                assert_eq!(seed.id(), particle.id());
                // let displacement = self.simbox().sep_in_box(particle.pos(), seed.pos());
                let theta = mov.step_factor(dir);
                // final_p += displacement.rotated_by(theta);
                final_or = final_or.rotated_by(theta);
            }
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

    fn choose_random_move(&self, rng: &mut Prng) -> ProposedMove {
        // 1. Choose a particle that will lead the move
        let seed_id = self.choose_random_p_id(rng);
        // 2. Choose a direction (unit vector) for the move (taken from maxwell-boltzman distribution)
        let rand_vec = random_unit_vec(rng);
        // 3. Choose a move type (translation or rotation)
        let is_rotation = rng.random::<f64>() >= PROB_TRANSLATE;
        // 4. Pick a cluster size cutoff
        // I don't know why we use this distribution, but it has something to do with particle choice fairness
        let cluster_cutoff = (1.0 / rng.random::<f64>()).floor() as usize;
        // 5. Choose a size for the move
        let step_size = if is_rotation {
            // Rotation
            let r: f64 = rng.random();
            MAX_ROTATION * (2.0 * r - 1.0)
            // MAX_ROTATION * r.powf(0.5)
        } else {
            // Translate
            // Scale step-size to uniformly sample unit sphere/circle.
            let r: f64 = rng.random();
            // random number between (-1.0 and 1.0) * max_translation
            MAX_TRANSLATION * (2.0 * r - 1.0)
        };

        ProposedMove::new(seed_id, step_size, is_rotation, cluster_cutoff, rand_vec)
    }

    // Note: each particle gets to attempt to link to every neighbor, even if neighbor
    // has tried to be recruited before
    fn recruit_cluster(&self, rng: &mut Prng, mov: &ProposedMove) -> Result<VirtualMoves> {
        let seed = self.particle(mov.seed_id);
        // particles in the cluster who still need to make their linking pass
        let mut worklist = VecDeque::from([(mov.seed_id, None)]);
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
            let (id, roots) = worklist.pop_front().unwrap();

            // in case we get duplicate entries in worklist
            if seen.contains(&id) {
                continue;
            }
            let particle = self.particle(id);

            seen.insert(id);

            // let final_p = self.calculate_motion(particle, mov, seed, MoveDir::Forward);
            let final_p = self.calculate_motion(particle, mov, seed, &roots, MoveDir::Forward);
            // println!("final_p = {:?}", final_p);
            // println!("final_p2 = {:?}", final_p2);
            // assert_eq!(final_p, final_p2);

            // this particle moves
            vmoves.push(id, final_p.clone());

            // let reverse_p = self.calculate_motion(particle, mov, seed, MoveDir::Backward);
            let reverse_p = self.calculate_motion(particle, mov, seed, &roots, MoveDir::Backward);

            let interactions = self.determine_interactions(particle);

            for &neighbor_id in interactions.iter() {
                let neighbor = self.particle(neighbor_id);

                let (link_weight, reverse_link_weight) =
                    self.compute_link_weights(particle, neighbor, &final_p, &reverse_p);

                if rng.random::<f64>() <= link_weight {
                    if rng.random::<f64>() <= reverse_link_weight / link_weight {
                        // The neighbor has been linked.
                        // Add it to the move and attempt to recruit its neighbors
                        if !seen.contains(&neighbor_id) {
                            worklist.push_back((
                                neighbor_id,
                                Some((final_p.clone(), reverse_p.clone())),
                            ));
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

    pub fn choose_random_p_id(&self, rng: &mut Prng) -> ParticleId {
        assert_ne!(self.particles().num_particles(), 0);
        loop {
            let p_id = rng.random_range(0..self.particles().particle_watermark() as u16);
            if self.particles().is_active_particle(p_id) {
                return p_id;
            }
        }
    }

    /// Commits a set of virtual moves to simulator state
    /// Returns Ok if commit was accepted, Err if there was an overlap
    /// The simulator will remain in a valid state regardless of whether this function returns Ok or Err.
    fn commit_moves(&mut self, vmoves: &VirtualMoves) -> Result<()> {
        // 1) remove all moved particles from the tenancy array
        for (p_id, _) in vmoves.inner.iter() {
            let curr_cell = self.simbox.get_cell_id(self.particle(*p_id).pos());
            self.simbox.delete_p_from_cell(*p_id, curr_cell);
        }

        // 2) check for overlaps while inserting particles into their new positions
        // We need to insert as we go to detect overlaps between moving particles
        let mut inserted = Vec::new();
        for (p_id, vp) in vmoves.inner.iter() {
            if self.simbox.would_overlap(*p_id, vp.pos()) {
                // Overlap detected - remove all previously inserted particles
                for &(restore_id, restore_cell) in inserted.iter() {
                    self.simbox.delete_p_from_cell(restore_id, restore_cell);
                }
                // Restore all particles to their original cells
                for (p_id, vp) in vmoves.inner.iter() {
                    let orig_cell = self.simbox.get_cell_id(vp.orig_pos());
                    self.simbox.insert_p_into_cell(*p_id, orig_cell);
                }
                return Err(anyhow!("Overlap"));
            }

            // No overlap - insert particle into new cell
            let new_cell = self.simbox.get_cell_id(vp.pos());
            self.simbox.insert_p_into_cell(*p_id, new_cell);
            inserted.push((*p_id, new_cell));
        }

        // 3) No overlaps found - update particle positions
        for (p_id, vp) in vmoves.inner.iter() {
            let p = self.simbox.particles_mut().particle_mut(*p_id);
            p.update_pos(vp.pos());
            p.update_or(vp.or());
        }
        Ok(())
    }

    fn attempt_commit(
        &mut self,
        rng: &mut Prng,
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
        if rng.random::<f64>() > scale_factor {
            return Err(anyhow!("Stokes drag rejection"));
        }

        self.commit_moves(vmoves)?;
        Ok(())
    }

    fn sim_has_overlaps(&self) -> bool {
        self.particles().iter().any(|p| self.simbox().overlaps(p))
    }

    // vmoves_well_formed
    // Translation: delta is same for all
    // Rotation: all moves stay the same distance from seed

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
            if (p.or().l2_norm() - 1.0).abs() >= 0.0001 {
                panic!("Orientation not normalized!: {:?} off by ", p);
            }
        }

        // check for overlaps
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

    pub fn step(&mut self, rng: &mut Prng, stats: &mut RunStats) -> Result<()> {
        // 0. If there are no particles, just skip
        if self.particles().num_particles() == 0 {
            return Ok(());
        }
        // debug_assert!(self.well_formed());
        stats.record_attempt();
        let mov = self.choose_random_move(rng);
        debug_assert!((mov.vec().l2_norm() - 1.0).abs() < 0.00001);
        log::debug!("Chose a random move: {:?}", mov);
        let virtual_moves = self.recruit_cluster(rng, &mov)?;
        log::debug!("Found a promising set of moves: {:?}", virtual_moves);
        if self.attempt_commit(rng, &mov, &virtual_moves).is_ok() {
            stats.record_accept();
        }
        // if self.sim_has_overlaps(){
        //     println!("Overlap on move: {:?}", mov);
        //     println!("Step factor = {:?}", mov.step_factor(MoveDir::Forward));
        //     for (p_id,vmov) in virtual_moves.inner.iter() {
        //         println!("{}: ({} -> {}) ({} -> {})", p_id, vmov.orig_pos(), vmov.pos(), vmov.orig_or(), vmov.or());
        //     }
        //     assert!(false);
        // }
        Ok(())
    }

    pub fn step_n(&mut self, n: usize, rng: &mut Prng) -> RunStats {
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
}
