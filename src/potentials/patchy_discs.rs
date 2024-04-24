use crate::consts::PARTICLE_DIAMETER;
use crate::particle::{IsParticle, Particle, ParticleId};
use crate::position::Orientation;
use crate::simbox::SimBox;
use std::f64::consts::PI;

// fn is_radian(theta: f64) -> bool {
//     (theta.is_normal() || theta.is_zero()) && theta.is_sign_positive() && theta <= (2.0*PI)
// }

// Note: Pairwise potentials are just a filter map

fn calc_angle(or: Orientation, other: Orientation) -> f64 {
    // This clamp seems to be necessary purely for numerical error reasons
    let dot = other.dot_prod(or).clamp(-1.0, 1.0); // clamp is here for numerical reasons. TODO: remove
    let cross = other.cross_prod(or);

    if cross >= 0.0 {
        dot.acos()
    } else {
        2.0 * PI - dot.acos()
    }
}

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

        // vector from p1 -> p0
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

        // calculate angle between difference vector and particle orientation to find closest patch
        let normalized_diff = diff.div_by(dist);
        let angle0 = calc_angle(normalized_diff, -or0); // TODO: why is this negative?
        let angle1 = calc_angle(normalized_diff, or1);
        // if !(angle0.is_normal() || angle0.is_zero()) || !(angle1.is_normal() || angle1.is_zero()) {
        //     println!("angle0 = {angle0}, angle1 = {angle1}, normalized_diff = {normalized_diff} or0 = {or0} or1 = {or1}");
        //     assert!(angle0.is_normal() || angle0.is_zero());
        //     assert!(angle1.is_normal() || angle1.is_zero());
        // }

        // println!("p0 angle = {angle0} p1 angle = {angle1}");

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
                // if interactions.len() == m.patches().len(){
                //     println!("{:?}", p);

                //     for interacting_p in interactions.iter() {
                //         println!("{:?}", simbox.particle(*interacting_p));
                //     }
                //     println!("{:?}", simbox.particle(neighbor_id));
                // }
                assert_ne!(interactions.len(), m.patches().len());
                interactions.push(neighbor_id);
            }
        }

        interactions
    }
}
