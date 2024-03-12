use rand::rngs::SmallRng;
use rand::Rng;
use rand_distr::{Distribution, Normal};

use fixed::types::I32F32;
use xyzvec::XYVec;

pub type Num = I32F32;
pub type DimVec = XYVec<Num>;
pub type Orientation = DimVec;
pub type PosDifference = DimVec;
pub type Position = DimVec;

#[macro_export]
macro_rules! num {
    // `()` indicates that the macro takes no argument.
    ($e: expr) => {
        // The macro will expand into the contents of this block.
        (crate::types::Num::from_num($e))
    };
}

// Random num between 0 and 1
pub fn rand_num(rng: &mut SmallRng) -> Num {
    num!(rng.gen::<f64>())
}

pub fn rand_num_range(rng: &mut SmallRng, lower: f64, upper: f64) -> Num {
    num!(rng.gen_range(lower..upper))
}

// rng.gen_range(self.min_x()..self.max_x())

pub fn random_dimvec(rng: &mut SmallRng) -> DimVec {
    let normal = Normal::new(0.0, 1.0).unwrap();
    let x = normal.sample(rng);
    let y = normal.sample(rng);
    DimVec::new([num!(x), num!(y)])
}

pub fn random_unit_vec(rng: &mut SmallRng) -> DimVec {
    let rand_pos = random_dimvec(rng);
    let norm = rand_pos.l2_norm();
    rand_pos.div_by(norm)
}
