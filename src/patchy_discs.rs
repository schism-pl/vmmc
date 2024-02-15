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
        let dist_sqd = simbox.sep_in_box(p0, p1).norm_sqd();

        // overlap!
        if dist_sqd < PARTICLE_DIAMETER {
            return f64::INFINITY;
        }

        if dist_sqd.sqrt() > PARTICLE_DIAMETER + m0.max_patch_radius() + m1.max_patch_radius() {
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
                let patch_dist = simbox.sep_in_box(pc0, pc1).norm();

                if patch_dist < patch0.radius() + patch1.radius() {
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
