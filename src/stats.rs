use std::ops::Add;

#[derive(Debug)]
pub struct RunStats {
    num_attempts: usize,
    num_accepts: usize,
}

impl RunStats {
    pub fn new() -> Self {
        let num_attempts = 0;
        let num_accepts = 0;
        RunStats {
            num_attempts,
            num_accepts,
        }
    }

    pub fn record_attempt(&mut self) {
        self.num_attempts += 1;
    }

    pub fn record_accept(&mut self) {
        self.num_accepts += 1;
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

impl Default for RunStats {
    fn default() -> Self {
        Self::new()
    }
}
