use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Patch {
    radius: f64,  // radius of patch (in units of particle diameter)
    theta: f64,   // angle in degrees
    chemtype: u8, // patches must have same chemtype to be compatible
}

impl Patch {
    pub fn new(radius: f64, theta: f64, chemtype: u8) -> Self {
        Self {
            radius,
            theta,
            chemtype,
        }
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    pub fn chemtype(&self) -> u8 {
        self.chemtype
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Morphology {
    patches: Vec<Patch>,
    max_radius: f64,
    // maps theta -> sin(theta) for each patch angle
    sin_theta: Vec<f64>,
    // maps theta -> cos(theta) for each patch angle
    cos_theta: Vec<f64>,
}

impl Morphology {
    pub fn new(patches: Vec<Patch>) -> Self {
        let max_radius = patches
            .iter()
            .map(|p| p.radius)
            .fold(f64::MIN, |a, b| a.max(b));

        let mut sin_theta: Vec<_> = Vec::new();
        let mut cos_theta: Vec<_> = Vec::new();
        for patch in &patches {
            let theta = patch.theta() * PI / 180.0; // convert degrees to radians
            sin_theta.push(theta.sin());
            cos_theta.push(theta.cos());
        }

        Self {
            patches,
            max_radius,
            sin_theta,
            cos_theta,
        }
    }

    pub fn patches(&self) -> &[Patch] {
        &self.patches
    }

    pub fn max_patch_radius(&self) -> f64 {
        self.max_radius
    }

    pub fn sin_theta(&self, patch_idx: usize) -> f64 {
        self.sin_theta[patch_idx]
    }

    pub fn cos_theta(&self, patch_idx: usize) -> f64 {
        self.cos_theta[patch_idx]
    }

    pub fn regular_3patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 120.0, 0);
        let p2 = Patch::new(radius, 240.0, 0);
        Self::new(vec![p0, p1, p2])
    }

    pub fn regular_4patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 90.0, 0);
        let p2 = Patch::new(radius, 180.0, 0);
        let p3 = Patch::new(radius, 270.0, 0);
        Self::new(vec![p0, p1, p2, p3])
    }

    pub fn regular_6patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 60.0, 0);
        let p2 = Patch::new(radius, 120.0, 0);
        let p3 = Patch::new(radius, 180.0, 0);
        let p4 = Patch::new(radius, 240.0, 0);
        let p5 = Patch::new(radius, 300.0, 0);
        Self::new(vec![p0, p1, p2, p3, p4, p5])
    }
}
