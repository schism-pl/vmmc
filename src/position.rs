use rand::rngs::SmallRng;
use rand_distr::{Distribution, Normal};

use xyzvec::XYVec;

// use crate::consts::DIMENSION;

// #[derive(Clone, PartialEq, Copy)]
// pub struct DimVec {
//     inner: [f64; DIMENSION],
// }

pub type DimVec = XYVec<f64>;
pub type Orientation = DimVec;
pub type PosDifference = DimVec;
pub type Position = DimVec;

pub fn random_dimvec(rng: &mut SmallRng) -> DimVec {
    let normal = Normal::new(0.0, 1.0).unwrap();
    // let mut r = Position::zeroes();
    let x = normal.sample(rng);
    let y = normal.sample(rng);
    DimVec::new([x, y])
}

pub fn random_unit_vec(rng: &mut SmallRng) -> DimVec {
    let rand_pos = random_dimvec(rng);
    let norm = rand_pos.l2_norm();
    rand_pos.div_by(norm)
}

//     // pick dimvec from normal distribution
//     pub fn random(rng: &mut SmallRng) -> Self {
//         let normal = Normal::new(0.0, 1.0).unwrap();
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = normal.sample(rng);
//         }
//         r
//     }

//     // normalized radom dimvec
//     pub fn rand_unit_vector(rng: &mut SmallRng) -> Self {
//         let rand_pos = DimVec::random(rng);
//         let norm = rand_pos.norm();
//         rand_pos.scalar_div(norm)
//     }

// impl DimVec {
//     pub fn new(inner: [f64; DIMENSION]) -> Self {
//         Self { inner }
//     }

//     pub fn zeroes() -> Self {
//         Position {
//             inner: [0.0; DIMENSION],
//         }
//     }

//     pub fn x(&self) -> f64 {
//         self.inner[0]
//     }

//     pub fn y(&self) -> f64 {
//         self.inner[1]
//     }

//     pub fn z(&self) -> f64 {
//         self.inner[2]
//     }

//     pub fn scalar_mul(&self, d: f64) -> Position {
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = self.inner[idx] * d;
//         }
//         r
//     }

//     pub fn scalar_div(&self, d: f64) -> Position {
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = self.inner[idx] / d;
//         }
//         r
//     }

//     // pick dimvec from normal distribution
//     pub fn random(rng: &mut SmallRng) -> Self {
//         let normal = Normal::new(0.0, 1.0).unwrap();
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = normal.sample(rng);
//         }
//         r
//     }

//     // normalized radom dimvec
//     pub fn rand_unit_vector(rng: &mut SmallRng) -> Self {
//         let rand_pos = DimVec::random(rng);
//         let norm = rand_pos.norm();
//         rand_pos.scalar_div(norm)
//     }

//     pub fn norm(&self) -> f64 {
//         let mut r = 0.0;
//         for idx in 0..DIMENSION {
//             r += self.inner[idx] * self.inner[idx];
//         }
//         r.sqrt()
//     }

//     pub fn norm_sqd(&self) -> f64 {
//         let mut r = 0.0;
//         for idx in 0..DIMENSION {
//             r += self.inner[idx] * self.inner[idx];
//         }
//         r
//     }

//     pub fn cross_prod(&self, other: Position) -> f64 {
//         self.x() * other.y() - self.y() * other.x()
//     }

//     pub fn dot_prod(&self, other: Position) -> f64 {
//         self.x() * other.x() + self.y() * other.y()
//     }

//     pub fn cross_prod_magnitude_sqd(&self, other: Position) -> f64 {
//         if DIMENSION == 2 {
//             let a1: f64 = self.x() * other.y() - self.y() * other.x();
//             a1 * a1
//         } else if DIMENSION == 3 {
//             let a1: f64 = self.x() * other.y() - self.y() * other.x();
//             let a2: f64 = self.y() * other.z() - self.z() * other.y();
//             let a3: f64 = self.z() * other.x() - self.x() * other.z();

//             a1 * a1 + a2 * a2 + a3 * a3
//         } else {
//             panic!("Dimension is not 2 or 3")
//         }
//     }

//     pub fn shifted_by(&self, x: f64, y: f64) -> Self {
//         let new_x = self.x() + x;
//         let new_y = self.y() + y;
//         Self {
//             inner: [new_x, new_y],
//         }
//     }

//     // TODO: what are the units of theta? (update docs if its not radians)
//     pub fn rotated_by(&self, theta: f64) -> Self {
//         let c = theta.cos();
//         let s = theta.sin();

//         let x = (self.x() * c - self.y() * s) - self.x();
//         let y = self.x() * s + self.y() * c - self.y();
//         Self::new([x, y])
//     }

//     // reenable if we do 3d
//     // fn rotated_by_3d(&self, _other: Self, _theta: f64) -> Self {
//     //     panic!("Reenable if we do 3d");
//     //     // let c  = theta.cos();
//     //     // let s = theta.sin();

//     //     // let cross_prod = self.cross_prod(other);
//     //     // let x = ((self.x() - other.x()*cross_prod))*(c - 1.0) + (other.z()*self.y() - other.y()*self.z())*s;
//     //     // let y = ((self.y() - other.y()*cross_prod))*(c - 1.0) + (other.x()*self.z() - other.z()*self.x())*s;
//     //     // let z = ((self.z() - other.z()*cross_prod))*(c - 1.0) + (other.y()*self.x() - other.x()*self.y())*s;
//     //     // Self::new([x,y,z])
//     // }

//     // pub fn rotated_by(&self, other: Self, theta: f64) -> Self {
//     //     match DIMENSION {
//     //         2 => self.rotated_by_2d(theta),
//     //         3 => self.rotated_by_3d(other, theta),
//     //         _ => panic!("Dimension is not 2 or 3"),
//     //     }
//     // }

//     pub fn div_u32(self, other: [u32; 2]) -> Self {
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = self.inner[idx] / other[idx] as f64;
//         }
//         r
//     }
// }

// impl Add for Position {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = self.inner[idx] + other.inner[idx];
//         }
//         r
//     }
// }

// impl AddAssign for Position {
//     fn add_assign(&mut self, other: Position) {
//         for idx in 0..DIMENSION {
//             self.inner[idx] += other.inner[idx];
//         }
//     }
// }

// impl Sub for Position {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         let mut r = Position::zeroes();
//         for idx in 0..DIMENSION {
//             r.inner[idx] = self.inner[idx] - other.inner[idx];
//         }
//         r
//     }
// }

// impl SubAssign for Position {
//     fn sub_assign(&mut self, other: Position) {
//         for idx in 0..DIMENSION {
//             self.inner[idx] -= other.inner[idx];
//         }
//     }
// }

// impl fmt::Debug for DimVec {
//     fn fmt(&self, f: &mut Formatter) -> fmt::Result {
//         write!(f, "({:.3}, {:.3})", self.x(), self.y())
//     }
// }
