use std::{
    fmt::{self, Formatter},
    ops::{Add, AddAssign, Div, Sub, SubAssign},
};

use rand::rngs::ThreadRng;
use rand_distr::{Distribution, Normal};

use crate::params::DIMENSION;

#[derive(Clone, PartialEq, Copy)]
pub struct DimVec {
    inner: [f64; DIMENSION],
}

pub type Orientation = DimVec;
pub type PosDifference = DimVec;
pub type Position = DimVec;

impl DimVec {
    pub fn new(inner: [f64; DIMENSION]) -> Self {
        Self { inner }
    }

    pub fn zeroes() -> Self {
        Position {
            inner: [0.0; DIMENSION],
        }
    }

    pub fn x(&self) -> f64 {
        self.inner[0]
    }

    pub fn y(&self) -> f64 {
        self.inner[1]
    }

    pub fn z(&self) -> f64 {
        self.inner[2]
    }

    pub fn scalar_mul(&self, d: f64) -> Position {
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = self.inner[idx] * d;
        }
        r
    }

    pub fn scalar_div(&self, d: f64) -> Position {
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = self.inner[idx] / d;
        }
        r
    }

    pub fn random(rng: &mut ThreadRng) -> Self {
        let normal = Normal::new(0.0, 1.0).unwrap();
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = normal.sample(rng);
        }
        r
    }

    // picks a random position then does a scalar div by maginitude
    pub fn unit_vector(rng: &mut ThreadRng) -> Self {
        let rand_pos = Position::random(rng);
        let norm = rand_pos.norm();
        rand_pos.scalar_div(norm)
    }

    pub fn norm(&self) -> f64 {
        let mut r = 0.0;
        for idx in 0..DIMENSION {
            r += self.inner[idx] * self.inner[idx];
        }
        r.sqrt()
    }

    pub fn norm_sqd(&self) -> f64 {
        let mut r = 0.0;
        for idx in 0..DIMENSION {
            r += self.inner[idx] * self.inner[idx];
        }
        r
    }

    // pub fn floor(&self) -> [usize; DIMENSION] {
    //     let mut r = [0;DIMENSION];
    //     for idx in 0..DIMENSION {
    //         r[idx] = self.inner[idx].floor() as usize;
    //     }
    //     r
    // }

    pub fn cross_prod_magnitude_sqd(&self, other: Position) -> f64 {
        if DIMENSION == 2 {
            let a1: f64 = self.x() * other.y() - self.y() * other.x();
            a1 * a1
        } else if DIMENSION == 3 {
            let a1: f64 = self.x() * other.y() - self.y() * other.x();
            let a2: f64 = self.y() * other.z() - self.z() * other.y();
            let a3: f64 = self.z() * other.x() - self.x() * other.z();

            a1 * a1 + a2 * a2 + a3 * a3
        } else {
            panic!("Dimension is not 2 or 3")
        }
    }

    pub fn shifted_by(&self, x: f64, y: f64) -> Self {
        let new_x = self.x() + x;
        let new_y = self.y() + y;
        Self {
            inner: [new_x, new_y],
        }
    }

    fn rotated_by_2d(&self, theta: f64) -> Self {
        let c = theta.cos();
        let s = theta.sin();

        let x = (self.x() * c - self.y() * s) - self.x();
        let y = self.x() * s + self.y() * c - self.y();
        Self::new([x, y])
    }

    // reenable if we do 3d
    fn rotated_by_3d(&self, _other: Self, _theta: f64) -> Self {
        panic!("Reenable if we do 3d");
        // let c  = theta.cos();
        // let s = theta.sin();

        // let cross_prod = self.cross_prod(other);
        // let x = ((self.x() - other.x()*cross_prod))*(c - 1.0) + (other.z()*self.y() - other.y()*self.z())*s;
        // let y = ((self.y() - other.y()*cross_prod))*(c - 1.0) + (other.x()*self.z() - other.z()*self.x())*s;
        // let z = ((self.z() - other.z()*cross_prod))*(c - 1.0) + (other.y()*self.x() - other.x()*self.y())*s;
        // Self::new([x,y,z])
    }

    pub fn rotated_by(&self, other: Self, theta: f64) -> Self {
        match DIMENSION {
            2 => self.rotated_by_2d(theta),
            3 => self.rotated_by_3d(other, theta),
            _ => panic!("Dimension is not 2 or 3"),
        }
    }

    pub fn div_u32(self, other: [u32; 2]) -> Self {
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = self.inner[idx] / other[idx] as f64;
        }
        r
    }
}

impl Add for Position {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = self.inner[idx] + other.inner[idx];
        }
        r
    }
}

impl AddAssign for Position {
    fn add_assign(&mut self, other: Position) {
        for idx in 0..DIMENSION {
            self.inner[idx] += other.inner[idx];
        }
    }
}

impl Sub for Position {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut r = Position::zeroes();
        for idx in 0..DIMENSION {
            r.inner[idx] = self.inner[idx] - other.inner[idx];
        }
        r
    }
}

impl SubAssign for Position {
    fn sub_assign(&mut self, other: Position) {
        for idx in 0..DIMENSION {
            self.inner[idx] -= other.inner[idx];
        }
    }
}

impl fmt::Debug for DimVec {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({:.3}, {:.3})", self.x(), self.y())
    }
}
