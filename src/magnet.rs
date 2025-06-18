use crate::consts::PARTICLE_DIAMETER;
use crate::particle::{IsParticle, Particle, ParticleId};
use crate::position::Orientation;
use crate::simbox::SimBox;
use std::f64::consts::PI;
use std::f64;

fn magnitude(vec: &Vec<f64>) -> f64 {
    vec.iter().map(|&x| x.powi(2)).sum::<f64>().sqrt()
}


// fn is_radian(theta: f64) -> bool {
//     (theta.is_normal() || theta.is_zero()) && theta.is_sign_positive() && theta <= (2.0*PI)
// }

// Note: Pairwise potentials are just a filter map

pub struct DipoleHamiltonian {
    H_vec: Vec<f64>,
    m_magnitude: f64,
    m_vec: Vec<f64>,
    H_magnitude: f64,
}

impl DipoleHamiltonian {
    //pub fn new(H_magnitude: f64,H_vec: Vec<f64>,m_magnitude: f64,m_vec: Vec<f64>) -> Self {
    //    Self { H_magnitude,H_vec,m_magnitude,m_vec }
    //}
    pub fn new(H_vec: Vec<f64>,m_vec: Vec<f64>) -> Self {
        let H_magnitude = (magnitude(&H_vec));
        let m_magnitude = (magnitude(&m_vec));
        Self {H_vec,m_magnitude,m_vec,H_magnitude}
    }

    pub fn m_magnitude(&self) -> f64 {
        self.m_magnitude
    }

    pub fn H_vec(&self) -> Vec<f64> {
        self.H_vec.clone()
    }

    pub fn m_vec(&self) -> Vec<f64> {
        self.m_vec.clone()
    }

    pub fn H_magnitude(&self) -> f64 {
        self.H_magnitude
    }

    //pub fn set_H_magnitude(&self, H_magnitude: f64) -> f64 {
    //    self.H_magnitude = H_magnitude
    //}

    pub fn set_m_magnitude(&mut self, m_magnitude: f64) {
        self.m_magnitude = m_magnitude
    }

    pub fn set_H_vec(&mut self, H_vec: Vec<f64>) {
        self.H_vec = H_vec
    }

    //pub fn m_vec(&self, m_vec: Vec<f64>) -> Vec<f64> {
    //    self.m_vec = m_vec
    //}



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

        let mu0 = 1.0;//magnetic permeability constant should be 1.25663706127(20)×10−6 N⋅A−2, but ill make it 1 in the code for simplicity.
        
        let m_vec0 = vec![self.m_vec[0]+m0.dstr()[0],self.m_vec[1]+m0.dstr()[1],self.m_vec[2]+m0.dstr()[2]];
        let m_vec1 = vec![self.m_vec[0]+m1.dstr()[0],self.m_vec[1]+m1.dstr()[1],self.m_vec[2]+m1.dstr()[2]];

        let mag0 = magnitude(&m_vec0);
        let mag1 = magnitude(&m_vec1);

        let dotProduct = (diff.x()*m_vec0[0]+diff.y()*m_vec0[1]+0.0*m_vec0[2]); 
        let cosTheta = dotProduct/(mag0*dist);

        let E = (1.0-3.0*cosTheta*cosTheta)*(mu0/(4.0*PI))*((mag0*mag1)/(dist*dist*dist));//the radius between particles and the dipole moment term are perpendicular so cos_theta is 0 leading to this reduced equation
        return E;
    }

    pub fn field_energy<P1: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
    ) -> f64 {
        let m0 = simbox.morphology(particle0);
        let m_vec0 = [self.m_vec[0]+m0.dstr()[0],self.m_vec[1]+m0.dstr()[1],self.m_vec[2]+m0.dstr()[2]];
        let p0 = particle0.pos(); 
        let dot_product = (self.H_vec[0]*m_vec0[0]+self.H_vec[1]*m_vec0[1]+self.H_vec[2]*m_vec0[2]);
        //equation is E = m_magnitude * H_magnitude * cos_theta, but cos_theta = dot_product /
        //m_magnitude * H_magnitude so E = dot_product
        return dot_product;
    
    }

}
