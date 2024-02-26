use rand::rngs::SmallRng;
use rand_distr::{Distribution, Normal};

use xyzvec::XYVec;

pub type DimVec = XYVec<f64>;
pub type Orientation = DimVec;
pub type PosDifference = DimVec;
pub type Position = DimVec;

pub fn random_dimvec(rng: &mut SmallRng) -> DimVec {
    let normal = Normal::new(0.0, 1.0).unwrap();
    let x = normal.sample(rng);
    let y = normal.sample(rng);
    DimVec::new([x, y])
}

pub fn random_unit_vec(rng: &mut SmallRng) -> DimVec {
    let rand_pos = random_dimvec(rng);
    let norm = rand_pos.l2_norm();
    rand_pos.div_by(norm)
}
