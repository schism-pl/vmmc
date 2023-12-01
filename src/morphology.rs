pub struct Patch {
    radius: f64,
    theta: f64,      // angle in degrees
    color: u8,       // rename color?
    radius_sqd: f64, // cache radius squared
}

impl Patch {
    pub fn new(radius: f64, theta: f64, color: u8) -> Self {
        Self {
            radius,
            theta,
            color,
            radius_sqd: radius * radius,
        }
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn radius_sqd(&self) -> f64 {
        self.radius_sqd
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    pub fn color(&self) -> u8 {
        self.color
    }
}

pub struct Morphology {
    patches: Vec<Patch>,
    // max distance that this particle can interact with another
    // used for optimization
    sqd_cutoff_max: f64,
}

impl Morphology {
    pub fn new(patches: Vec<Patch>) -> Self {
        // max of all squared cutoff distances
        let sqd_cutoff_max = patches
            .iter()
            .map(|p| (1.0 + p.radius) * (1.0 + p.radius))
            .fold(f64::MIN, |a, b| a.min(b));
        Self {
            patches,
            sqd_cutoff_max,
        }
    }

    pub fn patches(&self) -> &[Patch] {
        &self.patches
    }

    pub fn sqd_cutoff_max(&self) -> f64 {
        self.sqd_cutoff_max
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
