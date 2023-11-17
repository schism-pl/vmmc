use crate::params::DIMENSION;
use crate::particle::{IsParticle, Particle, ParticleId};
use crate::position::{Orientation, Position};
use crate::simbox::SimBox;
use std::f64::consts::PI;

// Note: Pairwise potentials are just a filter map

pub struct PatchyDiscParams {
    max_interactions: usize,
    interaction_energy: f64,
    sqd_cutoff_distance: f64,
}

impl PatchyDiscParams {
    pub fn new(max_interactions: usize, interaction_energy: f64, sqd_cutoff_distance: f64) -> Self {
        Self {
            max_interactions,
            interaction_energy,
            sqd_cutoff_distance,
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
        let patch_seperation: f64 = 2.0 * PI / (params.max_interactions as f64);
        let mut sin_theta: Vec<_> = Vec::new();
        let mut cos_theta: Vec<_> = Vec::new();
        for idx in 0..params.max_interactions {
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

    pub fn compute_pair_energy<
        P1: IsParticle + std::fmt::Debug,
        P2: IsParticle + std::fmt::Debug,
    >(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64 {
        assert!(DIMENSION == 2); // patchy disks is for 2d

        let p0 = particle0.pos();
        let or0 = particle0.or();
        let p1 = particle1.pos();
        let or1 = particle1.or();

        // we could make this norm_sqd if that made things faster?
        let dist = simbox.sep_in_box(p0, p1).norm_sqd();

        if dist < 1.0 {
            return f64::INFINITY;
        }
        // assert!(dist >= 1.0); // particles should not overlap

        // check all pairs of patches
        let mut energy = 0.0;
        for idx in 0..self.params.max_interactions {
            // Compute position of patch i on first disc.
            let new_p0 = self.pos_on_disc(&simbox, idx, p0, or0);
            for jdx in 0..self.params.max_interactions {
                let new_p1 = self.pos_on_disc(&simbox, jdx, p1, or1);
                let sqd_dist = simbox.sep_in_box(new_p0, new_p1).norm_sqd();

                if sqd_dist < self.params.sqd_cutoff_distance {
                    energy = self.params.interaction_energy;
                }
            }
        }

        -energy
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
                if interactions.len() == self.params.max_interactions {
                    println!("{:?} {:?}", p, interactions);
                    println!("--- {:?}", neighbor);
                    for &interacting_p_id in interactions.iter() {
                        println!("--- {:?}", &particles[interacting_p_id as usize]);
                    }
                }
                assert_ne!(interactions.len(), self.params.max_interactions);
                interactions.push(neighbor_id);
            }
        }

        interactions
    }
}
