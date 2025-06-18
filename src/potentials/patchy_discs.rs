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
    sigma: f64
}

impl PatchyDiscsPotential {
    pub fn new(interaction_energy: f64, sigma: f64) -> Self {
        Self { interaction_energy , sigma}
    }

    pub fn interaction_energy(&self) -> f64 {
        self.interaction_energy
    }

    pub fn sigma(&self) -> f64 {
        self.sigma
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
        //below is the former code which only calculated the closest patches

        //to get the patch location, we need the angle of the nanoparticle, the angle of the patch, and
        //then multiply by the radius?
        if let Some(p0idx) = m0.closest_patch_index(angle0) {
            let patch0 = &m0.patches()[p0idx];
            if let Some(p1idx) = m1.closest_patch_index(angle1) {
                let patch1 = &m0.patches()[p1idx];
                if patch0.chemtype() != patch1.chemtype() {
                    return 0.0;
                }
                if patch0.repulsive() != patch1.repulsive() {
                    return 0.0;
                }
                if dist <= PARTICLE_DIAMETER + patch0.radius().max(patch1.radius()) {
                    //println!("patch accepted pD:{}",pD);
                    if patch0.repulsive() {
                        //return self.interaction_energy;
                    } else {
                        return -self.interaction_energy;
                    }
                } 
            }
        }
        //this new code will instead loop and calculate all patches
        //I don't understand how dist is calculated, but it is also the distance between particles,
        //not patches, so i will find this in a new way.
        //this code doesnt work yet, but is for attractive potential, rather than touching
        0.0
    }

    pub fn direct_interaction<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> bool {
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
            return false;
        }

        if dist > PARTICLE_DIAMETER + m0.max_patch_radius().max(m1.max_patch_radius()) {
            return false;
        }

        let or0 = particle0.or();
        let or1 = particle1.or();

        // calculate angle between difference vector and particle orientation to find closest patch
        let normalized_diff = diff.div_by(dist);
        let angle0 = calc_angle(normalized_diff, -or0); // TODO: why is this negative?
        let angle1 = calc_angle(normalized_diff, or1);
        if let Some(p0idx) = m0.closest_patch_index(angle0) {
            let patch0 = &m0.patches()[p0idx];
            if let Some(p1idx) = m1.closest_patch_index(angle1) {
                let patch1 = &m0.patches()[p1idx];
                if patch0.chemtype() != patch1.chemtype() {
                    return false;
                }
                if patch0.repulsive() != patch1.repulsive() {
                    return false;
                }
                if dist <= PARTICLE_DIAMETER + patch0.radius().max(patch1.radius()) {
                    return true;
                }
            }
        }


        return false
    }


    pub fn determine_interactions(&self, simbox: &SimBox, p: &Particle) -> Vec<ParticleId> {
        let mut interactions = Vec::new();
        for neighbor_id in simbox.get_neighbors(p.pos()) {
            let neighbor = simbox.particle(neighbor_id); //&simbox.particles()[neighbor_id as usize];
            if neighbor == p {
                continue;
            }
            if self.direct_interaction(simbox, p, neighbor) {
                // if interactions.len() == m.patches().len(){
                //     println!("{:?}", p);
                //     for interacting_p in interactions.iter() {
                //         println!("{:?}", simbox.particle(*interacting_p));
                //     }
                //     println!("{:?}", simbox.particle(neighbor_id));
                // }
                //if(interactions.len() != m.patches().len())
                //{
                //    println!("i: {}, p: {}",interactions.len(),m.patches().len()); 
                //}else{
                //    println!("interactions match");
                //}
                //assert_ne!(interactions.len(), m.patches().len());
                //println!("interaction len: {}",interactions.len());
                interactions.push(neighbor_id);//i am removing this stuff, its limited to the
                //number of patches in another file and causing memory errors. i dont see the point
                //of this.
            }
        }

        interactions
    }
}
