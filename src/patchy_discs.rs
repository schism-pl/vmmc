use crate::consts::PARTICLE_RADIUS;
use crate::morphology::Morphology;
use crate::particle::{IsParticle, Particle, ParticleId, ShapeId};
use crate::position::{Orientation, Position};
use crate::simbox::SimBox;
use std::f64::consts::PI;

// Note: Pairwise potentials are just a filter map

pub struct PatchyDiscsPotential {
    // maps shape_id -> theta -> sin(theta) for each patch angle
    sin_thetas: Vec<Vec<f64>>,
    // maps shape_id -> theta -> cos(theta) for each patch angle
    cos_thetas: Vec<Vec<f64>>,
    interaction_energy: f64,
}

impl PatchyDiscsPotential {
    pub fn new(shapes: &[Morphology], interaction_energy: f64) -> Self {
        let mut sin_thetas: Vec<_> = Vec::new();
        let mut cos_thetas: Vec<_> = Vec::new();
        for shape in shapes {
            let mut sin_theta: Vec<_> = Vec::new();
            let mut cos_theta: Vec<_> = Vec::new();
            for patch in shape.patches() {
                let theta = patch.theta() * PI / 180.0; // convert degrees to radians
                sin_theta.push(theta.sin());
                cos_theta.push(theta.cos());
            }
            sin_thetas.push(sin_theta);
            cos_thetas.push(cos_theta);
        }

        println!(
            "{:?} {:?} {:?} {:?}",
            shapes, interaction_energy, sin_thetas, cos_thetas
        );
        Self {
            sin_thetas,
            cos_thetas,
            interaction_energy,
        }
    }

    pub fn interaction_energy(&self) -> f64 {
        self.interaction_energy
    }

    pub fn set_interaction_energy(&mut self, interaction_energy: f64) {
        self.interaction_energy = interaction_energy
    }

    pub fn pos_on_disc(
        &self,
        simbox: &SimBox,
        p_idx: usize,
        p: Position,
        or: Orientation,
        shape_id: ShapeId,
    ) -> Position {
        let sin_theta = self.sin_thetas[shape_id as usize][p_idx];
        let cos_theta = self.cos_thetas[shape_id as usize][p_idx];
        let x = p.x() + PARTICLE_RADIUS * (or.x() * cos_theta - or.y() * sin_theta);
        let y = p.y() + PARTICLE_RADIUS * (or.x() * sin_theta + or.y() * cos_theta);
        simbox.map_pos_into_box(Position::new([x, y])) // simbox map into pos
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
        let or0 = particle0.or();
        let p1 = particle1.pos();
        let or1 = particle1.or();

        // norm_sqd < 1.0 => norm < 1.0
        let dist_sqd = simbox.sep_in_box(p0, p1).norm_sqd();

        if dist_sqd < 1.0 {
            return f64::INFINITY;
        }

        // theres no way they can interact
        if dist_sqd > m0.sqd_cutoff_max().max(m1.sqd_cutoff_max()) {
            return 0.0;
        }

        // check all pairs of patches
        for (p_idx0, patch0) in m0.patches().iter().enumerate() {
            // Compute position of patch i on first disc.
            let new_p0 = self.pos_on_disc(simbox, p_idx0, p0, or0, particle0.shape_id());
            for (p_idx1, patch1) in m1.patches().iter().enumerate() {
                // if these patches aren't compatible, skip
                if patch0.color() != patch1.color() {
                    continue;
                }

                let new_p1 = self.pos_on_disc(simbox, p_idx1, p1, or1, particle1.shape_id());
                let sqd_dist = simbox.sep_in_box(new_p0, new_p1).norm_sqd();

                if sqd_dist < patch0.radius_sqd().max(patch1.radius_sqd()) {
                    // theres no way for more than 2 patches to interact between 2 particles
                    // TODO: write out exact conditions for this and assert it
                    return -self.interaction_energy;
                }
            }
        }

        0.0
    }

    pub fn determine_interactions(&self, simbox: &SimBox, p: &Particle) -> Vec<ParticleId> {
        // log::debug!("Determining interactions for p{:?}", p.id());
        let mut interactions = Vec::new();
        for neighbor_id in simbox.get_neighbors(p) {
            let neighbor = &simbox.particles()[neighbor_id as usize];
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
