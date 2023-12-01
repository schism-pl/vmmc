use rand::rngs::SmallRng;

use crate::position::{Orientation, Position};

pub type ParticleId = u16;
pub type ShapeId = u16;

pub trait IsParticle {
    fn pos(&self) -> Position;
    fn or(&self) -> Orientation;
}

// pub struct Patch {
//     radius:
//     color:
//     offset:
// }

// pub struct Morphology {
//     patches: Vec<Patch>,
// }

#[derive(PartialEq, Debug)]
pub struct Particle {
    id: ParticleId,
    shape_id: ShapeId,
    pos: Position,
    or: Orientation,
}

impl IsParticle for Particle {
    fn pos(&self) -> Position {
        self.pos
    }

    fn or(&self) -> Orientation {
        self.or
    }
}

impl Particle {
    pub fn new(id: ParticleId, pos: Position, or: Orientation) -> Self {
        Self {
            id,
            pos,
            shape_id: 0,
            or,
        }
    }

    pub fn id(&self) -> ParticleId {
        self.id
    }

    pub fn update_pos(&mut self, new_pos: Position) {
        self.pos = new_pos;
    }

    pub fn update_or(&mut self, new_or: Orientation) {
        self.or = new_or;
    }

    pub fn rotate(&mut self, other: Orientation, theta: f64) {
        self.or = self.or.rotated_by(other, theta)
    }

    pub fn random(rng: &mut SmallRng, id: ParticleId) -> Self {
        Particle::new(id, Position::random(rng), Orientation::random(rng))
    }
}

#[derive(PartialEq, Debug, Clone)]
pub struct VParticle {
    orig_pos: Position,
    orig_or: Orientation,
    pos: Position,
    or: Orientation,
}

impl VParticle {
    pub fn new(orig_pos: Position, orig_or: Orientation, pos: Position, or: Orientation) -> Self {
        Self {
            orig_pos,
            orig_or,
            pos,
            or,
        }
    }

    pub fn orig_pos(&self) -> Position {
        self.orig_pos
    }

    pub fn orig_or(&self) -> Orientation {
        self.orig_or
    }
}

impl IsParticle for VParticle {
    fn pos(&self) -> Position {
        self.pos
    }

    fn or(&self) -> Orientation {
        self.or
    }
}
