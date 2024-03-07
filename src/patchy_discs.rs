use crate::consts::PARTICLE_DIAMETER;
use crate::particle::{IsParticle, Particle, ParticleId};
use crate::simbox::SimBox;

// Note: Pairwise potentials are just a filter map

pub struct PatchyDiscsPotential {
    interaction_energy: f64,
}

impl PatchyDiscsPotential {
    pub fn new(interaction_energy: f64) -> Self {
        Self { interaction_energy }
    }

    pub fn interaction_energy(&self) -> f64 {
        self.interaction_energy
    }

    pub fn set_interaction_energy(&mut self, interaction_energy: f64) {
        self.interaction_energy = interaction_energy
    }

    pub fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64 {
        let m0 = simbox.morphology(particle0);
        let m1 = simbox.morphology(particle1);

        let p0 = particle0.pos();
        let p1 = particle1.pos();

        // norm_sqd < 1.0 => norm < 1.0
        // TODO: explicity calc diff here (d0 in steve's code)
        let dist_sqd = simbox.sep_in_box(p0, p1).l2_norm_sqd();

        // overlap!
        if dist_sqd < PARTICLE_DIAMETER {
            return f64::INFINITY;
        }

        // TODO: change for particle overlap (to max)
        // if dist_sqd.sqrt() > PARTICLE_DIAMETER + m0.max_patch_radius() + m1.max_patch_radius() {
        //     return 0.0;
        // }
        if dist_sqd.sqrt() > PARTICLE_DIAMETER + m0.max_patch_radius().max(m1.max_patch_radius()) {
            return 0.0;
        }

        // check all pairs of patches
        for (p_idx0, patch0) in m0.patches().iter().enumerate() {
            // Compute position of patch i on first disc.
            let pc0 = simbox.patch_center(particle0, p_idx0);
            for (p_idx1, patch1) in m1.patches().iter().enumerate() {
                // if these patches aren't compatible, skip
                if patch0.chemtype() != patch1.chemtype() {
                    continue;
                }

                let pc1 = simbox.patch_center(particle1, p_idx1);
                let patch_dist = simbox.sep_in_box(pc0, pc1).l2_norm();

                // TODO: change to min for patch overlap
                if patch_dist < patch0.radius().max(patch1.radius()) {
                    // theres no way for more than 2 patches to interact between 2 particles
                    // TODO: write out exact conditions for this and assert it
                    // appendix E of https://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.011044
                    return -self.interaction_energy;
                }
            }
        }

        0.0
    }

    pub fn determine_interactions(&self, simbox: &SimBox, p: &Particle) -> Vec<ParticleId> {
        // log::debug!("Determining interactions for p{:?}", p.id());
        let mut interactions = Vec::new();
        for neighbor_id in simbox.get_neighbors(p.pos()) {
            let neighbor = simbox.particle(neighbor_id); //&simbox.particles()[neighbor_id as usize];
            if neighbor == p {
                continue;
            }
            let energy = self.compute_pair_energy(simbox, p, neighbor);

            // particles interact!
            if energy < 0.0 {
                let m = simbox.morphology(p);
                assert_ne!(interactions.len(), m.patches().len());
                interactions.push(neighbor_id);
            }
        }

        interactions
    }
}

pub struct AngularPatchyDiscsPotential {
    interaction_energy: f64,
}

impl AngularPatchyDiscsPotential {
    pub fn new(interaction_energy: f64) -> Self {
        Self { interaction_energy }
    }

    pub fn interaction_energy(&self) -> f64 {
        self.interaction_energy
    }

    pub fn set_interaction_energy(&mut self, interaction_energy: f64) {
        self.interaction_energy = interaction_energy
    }

    pub fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64 {
        let m0 = simbox.morphology(particle0);
        let m1 = simbox.morphology(particle1);

        let p0 = particle0.pos();
        let p1 = particle1.pos();

        let diff = simbox.sep_in_box(p0, p1);

        // norm_sqd < 1.0 => norm < 1.0
        let dist = diff.l2_norm_sqd().sqrt();

        // overlap!
        if dist < PARTICLE_DIAMETER {
            return f64::INFINITY;
        }

        if dist > PARTICLE_DIAMETER + m0.max_patch_radius().max(m1.max_patch_radius()) {
            return 0.0;
        }

        let or0 = particle0.or();
        let or1 = particle1.or();

        // calculate angle between difference vector and particle orientationt to find closest patch
        let angle0 = diff.dot_prod(or0).acos();
        let angle1 = diff.dot_prod(or1).acos();

        if let Some(patch0) = m0.closest_patch(angle0) {
            if let Some(patch1) = m1.closest_patch(angle1) {
                if dist <= PARTICLE_DIAMETER + patch0.radius().max(patch1.radius()) {
                    return -self.interaction_energy;
                }
            }
        }
        0.0
    }

    pub fn determine_interactions(&self, simbox: &SimBox, p: &Particle) -> Vec<ParticleId> {
        // log::debug!("Determining interactions for p{:?}", p.id());
        let mut interactions = Vec::new();
        for neighbor_id in simbox.get_neighbors(p.pos()) {
            let neighbor = simbox.particle(neighbor_id); //&simbox.particles()[neighbor_id as usize];
            if neighbor == p {
                continue;
            }
            let energy = self.compute_pair_energy(simbox, p, neighbor);

            // particles interact!
            if energy < 0.0 {
                let m = simbox.morphology(p);
                assert_ne!(interactions.len(), m.patches().len());
                interactions.push(neighbor_id);
            }
        }

        interactions
    }
}
