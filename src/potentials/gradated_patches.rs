use crate::consts::NC_DIAG_LEN;
use crate::morphology::CoreShape;
// use crate::consts::PARTICLE_DIAMETER;
use crate::particle::IsParticle;
use crate::simbox::SimBox;

use super::Potential;

// Note: Pairwise potentials are just a filter map

pub struct GcPotential {
    interaction_energy: f64,
}

impl GcPotential {
    pub fn new(interaction_energy: f64) -> Self {
        Self { interaction_energy }
    }

    pub fn nc_compute_pair_energy<P1: IsParticle, P2: IsParticle>(
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
        // if dist < PARTICLE_DIAMETER {
        //     return f64::INFINITY;
        // }
        if simbox.nc_overlaps_with(particle0, particle1) {
            return f64::INFINITY;
        }

        let max_interaction_dist = m0.max_patch_radius().max(m1.max_patch_radius());

        // can't possibly interact!
        if dist > NC_DIAG_LEN + max_interaction_dist {
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

                if patch_dist < patch0.radius() + patch1.radius() {
                    // theres no way for more than 2 patches to interact between 2 particles
                    // TODO: write out exact conditions for this and assert it
                    // appendix E of https://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.011044
                    let scaling_factor =
                        (max_interaction_dist - (dist - NC_DIAG_LEN)) / max_interaction_dist;
                    return -self.interaction_energy * scaling_factor;
                }
            }
        }

        0.0
    }

    pub fn hs_compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        _simbox: &SimBox,
        _particle0: &P1,
        _particle1: &P2,
    ) -> f64 {
        unimplemented!()
        // let m0 = simbox.morphology(particle0);
        // let m1 = simbox.morphology(particle1);
        // let p0 = particle0.pos();
        // let p1 = particle1.pos();

        // let diff = simbox.sep_in_box(p0, p1);

        // // norm_sqd < 1.0 => norm < 1.0
        // let dist = diff.l2_norm_sqd().sqrt();

        // // overlap!
        // if dist < PARTICLE_DIAMETER {
        //     return f64::INFINITY;
        // }

        // let max_interaction_dist = m0.max_patch_radius() + m1.max_patch_radius();

        // // can't possibly interact!
        // if dist > PARTICLE_DIAMETER + max_interaction_dist {
        //     return 0.0;
        // }

        // // check all pairs of patches
        // for (p_idx0, patch0) in m0.patches().iter().enumerate() {
        //     // Compute position of patch i on first disc.
        //     let pc0 = simbox.patch_center(particle0, p_idx0);
        //     for (p_idx1, patch1) in m1.patches().iter().enumerate() {
        //         // if these patches aren't compatible, skip
        //         if patch0.chemtype() != patch1.chemtype() {
        //             continue;
        //         }

        //         let pc1 = simbox.patch_center(particle1, p_idx1);
        //         let patch_dist = simbox.sep_in_box(pc0, pc1).l2_norm();

        //         if patch_dist < patch0.radius() + patch1.radius() {
        //             // theres no way for more than 2 patches to interact between 2 particles
        //             // TODO: write out exact conditions for this and assert it
        //             // appendix E of https://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.011044
        //             if self.repulsed {
        //                 let scaling_factor = (max_interaction_dist - (dist - PARTICLE_DIAMETER))
        //                     / max_interaction_dist;
        //                 return -self.interaction_energy * scaling_factor;
        //             } else {
        //                 return -self.interaction_energy;
        //             }
        //         }
        //     }
        // }

        // 0.0
    }
}

impl Potential for GcPotential {
    fn interaction_energy(&mut self) -> &mut f64 {
        &mut self.interaction_energy
    }

    fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64 {
        match simbox.morphology(particle0).shape() {
            CoreShape::Circle => self.hs_compute_pair_energy(simbox, particle0, particle1),
            CoreShape::Square => self.nc_compute_pair_energy(simbox, particle0, particle1),
        }
    }
}
