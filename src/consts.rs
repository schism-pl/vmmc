use crate::types::Num;

/// Changing consts in this file will probably break the simulation
/// please don't do that
// TODO: explain why each of these is hard-coded
pub const MAX_PARTICLES_PER_CELL: usize = 4;
pub const PARTICLE_DIAMETER: Num = 1.0;
pub const PARTICLE_RADIUS: Num = 0.5;
pub const MAX_PARTICLES: usize = 10000;
pub const MAX_INITIAL_PACKING_FRACTION: Num = 0.5;

// less important, there just isn't a reason for users to mess with this

pub const PROB_TRANSLATE: Num = 0.5;
pub const MAX_TRANSLATION: Num = 0.3;
pub const MAX_ROTATION: Num = 0.2;
