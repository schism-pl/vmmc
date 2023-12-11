use std::ops::{Add, AddAssign};

#[derive(Debug)]
pub struct RunStats {
    num_attempts: usize,
    num_accepts: usize,
    // cluster_translations: Vec<usize>,
    // cluster_rotations: Vec<usize>,
}

impl RunStats {
    pub fn new() -> Self {
        let num_attempts = 0;
        let num_accepts = 0;
        // let cluster_translations = vec![0; num_particles];
        // let cluster_rotations = vec![0; num_particles];
        RunStats {
            num_attempts,
            num_accepts,
            // cluster_translations,
            // cluster_rotations,
        }
    }

    pub fn record_attempt(&mut self) {
        self.num_attempts += 1;
    }

    pub fn record_accept(&mut self) {
        self.num_accepts += 1;
        // if cluster_size >= self.cluster_rotations.len() {
        //     println!(
        //         "cluster size = {:?}, total_particles = {:?}",
        //         cluster_size,
        //         self.cluster_rotations.len()
        //     );
        //     panic!("more particles in cluster than total particles?")
        // }
        // if is_rotation {
        //     self.cluster_rotations[cluster_size] += 1;
        // } else {
        //     self.cluster_translations[cluster_size] += 1;
        // }
    }

    pub fn num_attempts(&self) -> usize {
        self.num_attempts
    }

    pub fn num_accepts(&self) -> usize {
        self.num_accepts
    }
}

impl Add for RunStats {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            num_attempts: self.num_attempts() + other.num_attempts(),
            num_accepts: self.num_accepts() + other.num_accepts(),
        }
    }
}
