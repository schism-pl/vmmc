use serde::{Deserialize, Serialize};

use crate::consts::MAX_PARTICLES;
use crate::position::{Orientation, Position};
use std::mem::{size_of, size_of_val};

pub type ParticleId = u16;
pub type ShapeId = u16;

pub trait IsParticle {
    fn pos(&self) -> Position;
    fn or(&self) -> Orientation;
    fn shape_id(&self) -> ShapeId;
}
// Clone is implemented to enable quickcheck
#[derive(PartialEq, Debug, Clone, Serialize, Deserialize)]
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

    fn shape_id(&self) -> ShapeId {
        self.shape_id
    }
}

impl Particle {
    pub fn new(id: ParticleId, pos: Position, or: Orientation, shape_id: ShapeId) -> Self {
        Self {
            id,
            pos,
            shape_id,
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

    // pub fn random(rng: &mut Prng, id: ParticleId, shape_id: ShapeId) -> Self {
    //     Particle::new(id, random_dimvec(rng), random_unit_vec(rng), shape_id)
    // }
}

#[derive(PartialEq, Debug, Clone)]
pub struct VParticle {
    orig_pos: Position,
    orig_or: Orientation,
    pos: Position,
    or: Orientation,
    shape_id: ShapeId,
}

impl VParticle {
    pub fn new(
        orig_pos: Position,
        orig_or: Orientation,
        pos: Position,
        or: Orientation,
        shape_id: ShapeId,
    ) -> Self {
        Self {
            orig_pos,
            orig_or,
            pos,
            or,
            shape_id,
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

    fn shape_id(&self) -> ShapeId {
        self.shape_id
    }
}

/// Collection of particles
/// Need somewhat complicated data structure to efficiently insert/remove particles
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Particles {
    particles: Vec<Option<Particle>>,
    reserve: Vec<ParticleId>,
    num_particles: usize,
}

impl Particles {
    pub fn new(particles: Vec<Particle>) -> Self {
        let num_particles = particles.len();
        Self {
            particles: particles.into_iter().map(Some).collect(),
            reserve: Vec::new(),
            num_particles,
        }
    }

    pub fn particle(&self, p_id: ParticleId) -> &Particle {
        match &self.particles[p_id as usize] {
            Some(p) => p,
            None => panic!("That particle doesn't exist!"),
        }
    }

    pub fn is_active_particle(&self, p_id: ParticleId) -> bool {
        self.particles[p_id as usize].is_some()
    }

    pub fn particle_mut(&mut self, p_id: ParticleId) -> &mut Particle {
        match &mut self.particles[p_id as usize] {
            Some(p) => p,
            None => panic!("That particle doesn't exist!"),
        }
    }

    pub fn num_particles(&self) -> usize {
        self.num_particles
    }

    // Highest particle ID ever assigned. Used to select random particles
    pub fn particle_watermark(&self) -> usize {
        self.particles.len()
    }

    pub fn needed_mem(&self) -> usize {
        size_of::<Option<Particle>>() * self.num_particles + size_of_val(&self.reserve)
    }

    pub fn max_needed_mem(&self) -> usize {
        size_of::<Option<Particle>>() * MAX_PARTICLES
    }

    pub fn insert(&mut self, p: Particle) {
        self.num_particles += 1;
        let p_id = p.id();
        if (p.id() as usize) < self.particles.len() {
            assert!(self.particles[p_id as usize].is_none());
            self.particles[p_id as usize] = Some(p);
        } else {
            self.particles.push(Some(p))
        }
    }

    pub fn get_unused_p_id(&mut self) -> ParticleId {
        match self.reserve.pop() {
            Some(p_id) => p_id,
            None => self.particles.len() as ParticleId,
        }
    }

    pub fn push_unused_p_id(&mut self, p_id: ParticleId) {
        self.reserve.push(p_id);
    }

    pub fn delete(&mut self, p_id: ParticleId) {
        self.num_particles -= 1;
        assert!(self.particles[p_id as usize].is_some());
        self.reserve.push(p_id);
        self.particles[p_id as usize] = None;
    }

    pub fn iter(&self) -> ParticleIterator {
        ParticleIterator::new(&self.particles)
    }
}

pub struct ParticleIterator<'a> {
    particles: &'a [Option<Particle>],
    index: usize,
}

impl<'a> ParticleIterator<'a> {
    fn new(particles: &'a [Option<Particle>]) -> Self {
        Self {
            particles,
            index: 0,
        }
    }
}

impl<'a> Iterator for ParticleIterator<'a> {
    type Item = &'a Particle;
    /// scan forward until we find a particle or run out of particles
    fn next(&mut self) -> Option<Self::Item> {
        // let mut particle = None;
        while self.index < self.particles.len() {
            let particle = match &self.particles[self.index] {
                Some(p) => Some(p),
                None => None,
            };
            self.index += 1;
            if particle.is_some() {
                return particle;
            }
        }
        None
    }
}
