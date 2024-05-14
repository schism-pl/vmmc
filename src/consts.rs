use std::f64::consts::SQRT_2;

/// Changing consts in this file will probably break the simulation
/// please don't do that
// TODO: explain why each of these is hard-coded
pub const MAX_PARTICLES_PER_CELL: usize = 4;
// pub const PARTICLE_DIAMETER: f64 = 1.0;
// pub const PARTICLE_RADIUS: f64 = 0.5;
pub const NC_SIDE_LEN: f64 = 1.0;
pub const NC_HALF_SIDE_LEN: f64 = NC_SIDE_LEN / 2.0;
pub const NC_DIAG_LEN: f64 = 1.0 * SQRT_2;
pub const NC_HALF_DIAG_LEN: f64 = NC_DIAG_LEN / 2.0;

pub const MAX_PARTICLES: usize = 10000;
pub const MAX_INITIAL_PACKING_FRACTION: f64 = 0.5;

// less important, there just isn't a reason for users to mess with this

pub const PROB_TRANSLATE: f64 = 0.5;
pub const MAX_TRANSLATION: f64 = 0.3;
pub const MAX_ROTATION: f64 = 0.2;
