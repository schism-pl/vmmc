use crate::Prng;
use rand::Rng;
use xyzvec::XYVec;

pub type DimVec = XYVec<f64>;
pub type Orientation = DimVec;
pub type PosDifference = DimVec;
pub type Position = DimVec;

pub fn random_dimvec(rng: &mut Prng) -> DimVec {
    // let normal = Normal::new(0.0, 1.0).unwrap();
    let x: f64 = rng.random_range(-1.0..1.0);
    let y: f64 = rng.random_range(-1.0..1.0);
    DimVec::new([x, y])
}

pub fn random_unit_vec(rng: &mut Prng) -> DimVec {
    let rand_pos = random_dimvec(rng);
    let norm = rand_pos.l2_norm();
    rand_pos.div_by(norm)
}
