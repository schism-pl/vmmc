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
}

impl Morphology {
    pub fn new(patches: Vec<Patch>) -> Self {
        let max_radius = patches
            .iter()
            .map(|p| p.radius)
            .fold(f64::MIN, |a, b| a.max(b));

        Self {
            patches,
            max_radius,
        }
    }

    pub fn patches(&self) -> &[Patch] {
        &self.patches
    }

    pub fn max_patch_radius(&self) -> f64 {
        self.max_radius
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
