pub struct Patch {
    radius: f64,
    offset: f64, // (degrees ???)
    color: u8,   // rename color?
}

impl Patch {
    pub fn new(radius: f64, offset: f64, color: u8) -> Self {
        Self {
            radius,
            offset,
            color,
        }
    }
}

pub struct Morphology {
    patches: Vec<Patch>,
}

impl Morphology {
    pub fn new(patches: Vec<Patch>) -> Self {
        Self { patches }
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
