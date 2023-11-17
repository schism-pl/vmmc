use crate::consts::DIMENSION;
use crate::particle::{IsParticle, Particle, ParticleId};
use crate::position::{Orientation, Position};
use crate::simbox::SimBox;
use std::f64::consts::PI;

// Note: Pairwise potentials are just a filter map

pub struct PatchyDiscParams {
    num_patches: usize,
    interaction_energy: f64,
    interaction_range: f64,
    sqd_cutoff_distance: f64,
    sqd_cutoff_max: f64,
}

impl PatchyDiscParams {
    pub fn new(num_patches: usize, interaction_energy: f64, interaction_range: f64) -> Self {
        Self {
            num_patches,
            interaction_energy,
            interaction_range,
            // sqd distance that particles interact
            sqd_cutoff_distance: interaction_range * interaction_range,
            // upper bound on distance that particles can be before they won't interact
            sqd_cutoff_max: (1.0 + interaction_range) * (1.0 + interaction_range),
        }
    }
}

pub struct PatchyDiscsPotential {
    sin_theta: Vec<f64>,
    cos_theta: Vec<f64>,
    params: PatchyDiscParams,
}

impl PatchyDiscsPotential {
    pub fn new(params: PatchyDiscParams) -> Self {
        let patch_seperation: f64 = 2.0 * PI / (params.num_patches as f64);
        let mut sin_theta: Vec<_> = Vec::new();
        let mut cos_theta: Vec<_> = Vec::new();
        for idx in 0..params.num_patches {
            let theta = (idx as f64) * patch_seperation;
            sin_theta.push(theta.sin());
            cos_theta.push(theta.cos());
        }
        Self {
            sin_theta,
            cos_theta,
            params,
        }
    }

    pub fn pos_on_disc(
        &self,
        simbox: &SimBox,
        idx: usize,
        p: Position,
        or: Orientation,
    ) -> Position {
        let x = p.x() + 0.5 * (or.x() * self.cos_theta[idx] - or.y() * self.sin_theta[idx]);
        let y = p.y() + 0.5 * (or.x() * self.sin_theta[idx] + or.y() * self.cos_theta[idx]);
        simbox.map_pos_into_box(Position::new([x, y])) // simbox map into pos
    }

    pub fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64 {
        debug_assert!(DIMENSION == 2); // patchy disks is for 2d

        let p0 = particle0.pos();
        let or0 = particle0.or();
        let p1 = particle1.pos();
        let or1 = particle1.or();

        // norm_sqd < 1.0 => norm < 1.0
        let dist_sqd = simbox.sep_in_box(p0, p1).norm_sqd();

        if dist_sqd < 1.0 {
            return f64::INFINITY;
        }

        // theres no way they can interact
        if dist_sqd > self.params.sqd_cutoff_max {
            return 0.0;
        }

        // assert!(dist >= 1.0); // particles should not overlap

        // check all pairs of patches
        for idx in 0..self.params.num_patches {
            // Compute position of patch i on first disc.
            let new_p0 = self.pos_on_disc(simbox, idx, p0, or0);
            for jdx in 0..self.params.num_patches {
                let new_p1 = self.pos_on_disc(simbox, jdx, p1, or1);
                let sqd_dist = simbox.sep_in_box(new_p0, new_p1).norm_sqd();

                if sqd_dist < self.params.sqd_cutoff_distance {
                    // theres no way for more than 2 patches to interaact between 2 particles
                    // TODO: write out exact conditions for this and assert it
                    return -self.params.interaction_energy;
                }
            }
        }

        0.0
    }

    pub fn determine_interactions(
        &self,
        simbox: &SimBox,
        particles: &[Particle],
        p: &Particle,
    ) -> Vec<ParticleId> {
        log::debug!("Determining interactions for p{:?}", p.id());
        let mut interactions = Vec::new();
        for neighbor_id in simbox.get_neighbors(p) {
            let neighbor = &particles[neighbor_id as usize];
            if neighbor == p {
                continue;
            }
            let energy = self.compute_pair_energy(simbox, p, neighbor);
            log::debug!(
                "p{:?} tries to interact with p{:?} with energy={:?}",
                p.id(),
                neighbor_id,
                energy
            );

            // particles interact!
            if energy < 0.0 {
                if interactions.len() == self.params.num_patches {
                    println!("{:?} {:?}", p, interactions);
                    println!("--- {:?}", neighbor);
                    for &interacting_p_id in interactions.iter() {
                        println!("--- {:?}", &particles[interacting_p_id as usize]);
                    }
                }
                assert_ne!(interactions.len(), self.params.num_patches);
                interactions.push(neighbor_id);
            }
        }

        interactions
    }
}
