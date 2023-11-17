#[derive(Debug)]
pub struct RunStats {
    num_attempts: usize,
    num_accepts: usize,
    cluster_translations: Vec<usize>,
    cluster_rotations: Vec<usize>,
}

impl RunStats {
    pub fn new(num_particles: usize) -> Self {
        let num_attempts = 0;
        let num_accepts = 0;
        let cluster_translations = vec![0; num_particles];
        let cluster_rotations = vec![0; num_particles];
        RunStats {
            num_attempts,
            num_accepts,
            cluster_translations,
            cluster_rotations,
        }
    }

    pub fn record_attempt(&mut self) {
        self.num_attempts += 1;
    }

    pub fn record_accept(&mut self, is_rotation: bool, cluster_size: usize) {
        self.num_accepts += 1;
        assert!(cluster_size < self.cluster_rotations.len());
        if is_rotation {
            self.cluster_rotations[cluster_size] += 1;
        } else {
            self.cluster_translations[cluster_size] += 1;
        }
    }

    pub fn num_attempts(&self) -> usize {
        self.num_accepts
    }

    pub fn num_accepts(&self) -> usize {
        self.num_accepts
    }
}
