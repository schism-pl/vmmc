use rand::{rngs::SmallRng, Rng};

use crate::consts::{MAX_PARTICLES_PER_CELL, PARTICLE_DIAMETER, PARTICLE_RADIUS};
use crate::morphology::Morphology;
use crate::particle::Particles;
use crate::position::random_unit_vec;
use crate::{
    particle::{IsParticle, Particle, ParticleId},
    position::{DimVec, Position},
};

type CellId = usize;
type Cell = [ParticleId; MAX_PARTICLES_PER_CELL];

// cells_per_axis * cells_per_axis
type CellGrid = Vec<Cell>;

fn map_into_range(p: f64, lower: f64, upper: f64) -> f64 {
    if p < lower {
        p + (upper - lower)
    } else if p >= upper {
        p - (upper - lower)
    } else {
        p
    }
}

fn map_id_into_range(p: i32, lower: i32, upper: i32) -> i32 {
    if p < lower {
        p + (upper - lower)
    } else if p >= upper {
        p - (upper - lower)
    } else {
        p
    }
}

// assumed to be periodic in every dimension
// dimension[i] = cell_dimensions[i] * cells_per_axis[i]
// Note: cell_dimensions * cells_per_axis must = dimensions

// Clone is implemented to enable quickcheck
#[derive(Clone, Debug)]
pub struct SimBox {
    dimensions: DimVec,
    cells_per_axis: [usize; 2],
    cell_dimensions: DimVec,
    particles: Particles,
    shapes: Vec<Morphology>, // maps shape_id to morphology
    cells: CellGrid,
    tenants: Vec<u8>,
}

// x_idx * cells_per_axis[1] * MAX_PARTICLES_PER_CELL + y_idx * MAX_PARTICLES_PER_CELL
// ParticleIdx::MAX = empty
impl SimBox {
    pub fn new(
        dimensions: DimVec,
        cells_per_axis: [usize; 2],
        cell_dimensions: DimVec,
        particles: Particles,
        shapes: Vec<Morphology>,
    ) -> Self {
        let cells = vec![[ParticleId::MAX; 4]; cells_per_axis[0] * cells_per_axis[1]];
        let tenants = vec![0_u8; cells.len()];
        Self {
            dimensions,
            cells_per_axis,
            cell_dimensions,
            particles,
            shapes,
            cells,
            tenants,
        }
    }

    /// if particle `p_id` were at position `pos`, would it cause any overlaps?
    pub fn would_overlap(&self, p_id: ParticleId, pos: Position) -> bool {
        for other_id in self.get_neighbors(pos) {
            let other = self.particle(other_id);
            if other_id == p_id {
                continue;
            }
            if self.sep_in_box(pos, other.pos()).l2_norm() < PARTICLE_DIAMETER {
                return true;
            }
        }
        false
    }

    // check if a particle overlaps any other particles
    /// must work even if p has invalid tenancy info due to its use in commit_moves
    pub fn overlaps(&self, p: &Particle) -> bool {
        self.would_overlap(p.id(), p.pos())
    }

    // uniform distribution of morphologies present in `shapes`
    pub fn new_with_randomized_particles(
        dimensions: DimVec,
        max_interaction_range: f64,
        num_particles: usize,
        shapes: Vec<Morphology>,
        rng: &mut SmallRng,
    ) -> Self {
        let cells_x_axis = (dimensions.x() / max_interaction_range).floor();
        let cells_y_axis = (dimensions.y() / max_interaction_range).floor();

        let cell_width = dimensions.x() / cells_x_axis;
        let cell_height = dimensions.y() / cells_y_axis;

        assert!(cell_width >= max_interaction_range);
        assert!(cell_height >= max_interaction_range);

        let cells_per_axis = [cells_x_axis as usize, cells_y_axis as usize];

        let cell_dimensions = DimVec::new([cell_width, cell_height]);

        // We initilize dummy simbox with no particles
        let mut simbox = Self::new(
            dimensions,
            cells_per_axis,
            cell_dimensions,
            Particles::new(Vec::new()),
            shapes,
        );

        for _ in 0..num_particles {
            let p = simbox.new_random_particle(rng);
            simbox.insert_particle(p);
        }
        simbox
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn particle(&self, p_id: ParticleId) -> &Particle {
        self.particles.particle(p_id)
    }

    pub fn particle_mut(&mut self, p_id: ParticleId) -> &mut Particle {
        self.particles.particle_mut(p_id)
    }

    pub fn particles(&self) -> &Particles {
        &self.particles
    }

    pub fn particles_mut(&mut self) -> &mut Particles {
        &mut self.particles
    }

    pub fn shapes(&self) -> &[Morphology] {
        &self.shapes
    }

    pub fn morphology<P: IsParticle>(&self, p: &P) -> &Morphology {
        &self.shapes[p.shape_id() as usize]
    }

    pub fn volume(&self) -> f64 {
        self.dimensions.x() * self.dimensions.y()
    }

    pub fn insert_particle(&mut self, p: Particle) {
        let cell_id = self.get_cell_id(p.pos());
        self.insert_p_into_cell(p.id(), cell_id);
        self.particles.insert(p);
    }

    pub fn remove_particle(&mut self, p_id: ParticleId) {
        let p = self.particle(p_id);
        let cell_id = self.get_cell_id(p.pos());
        self.delete_p_from_cell(p_id, cell_id);
        self.particles.delete(p_id)
    }

    // TODO: clean up this accursed function
    // Note: removes a particle from reserved particle ids even if particle isn't inserted
    pub fn new_random_particle(&mut self, rng: &mut SmallRng) -> Particle {
        // check if a position overlaps any existing
        fn pos_has_overlap(simbox: &SimBox, pos: &Position) -> bool {
            for other in simbox.particles.iter() {
                if simbox.sep_in_box(*pos, other.pos()).l2_norm() < PARTICLE_DIAMETER {
                    return true;
                }
            }
            false
        }

        let p_id = self.particles.get_unused_p_id();
        let shape_id = rng.gen_range(0..self.shapes.len()) as u16; // choose a uniform random shape
        let mut pos = self.random_pos(rng);
        while pos_has_overlap(self, &pos) {
            pos = self.random_pos(rng);
        }
        let or = random_unit_vec(rng);
        Particle::new(p_id, pos, or, shape_id)
    }

    pub fn get_cell_id(&self, pos: Position) -> CellId {
        // adjust so all indexes are positive
        // indexes start at bottom left and scan bottom to top
        let x = pos.x() + (self.dimensions.x() / 2.0);
        let y = pos.y() + (self.dimensions.y() / 2.0);
        assert!(x != self.dimensions.x());
        assert!(y != self.dimensions.y());

        let x_idx = (x / self.cell_dimensions.x()).floor() as usize;
        let y_idx = (y / self.cell_dimensions.y()).floor() as usize;
        x_idx * self.cells_per_axis[1] + y_idx
    }

    pub fn get_cell(&self, pos: Position) -> &Cell {
        let cell_id = self.get_cell_id(pos);
        &self.cells[cell_id]
    }

    pub fn get_cell_mut(&mut self, pos: Position) -> &mut Cell {
        let cell_id = self.get_cell_id(pos);
        &mut self.cells[cell_id]
    }

    pub fn in_cell(&self, p: &Particle) -> bool {
        let cell = self.get_cell(p.pos());
        cell.iter().any(|c| c == &p.id())
    }

    pub fn get_neighbor(&self, p: &Particle, x: f64, y: f64) -> &Cell {
        let neighbor_pos = p.pos().translated_by(x, y);
        let neighbor_pos = self.map_pos_into_box(neighbor_pos);
        self.get_cell(neighbor_pos)
    }

    fn delete_p_from_cell(&mut self, p_id: ParticleId, cell_id: CellId) {
        let cell = &mut self.cells[cell_id];
        for idx in 0..MAX_PARTICLES_PER_CELL {
            if cell[idx] == p_id {
                cell[idx] = ParticleId::MAX;
                self.tenants[cell_id] -= 1;
                return;
            }
        }
        panic!("Tried to delete particle from wrong cell")
    }

    fn insert_p_into_cell(&mut self, p_id: ParticleId, cell_id: CellId) {
        let cell = &mut self.cells[cell_id];
        for idx in 0..MAX_PARTICLES_PER_CELL {
            if cell[idx] == ParticleId::MAX {
                cell[idx] = p_id;
                self.tenants[cell_id] += 1;
                return;
            }
        }
        panic!("Tried to insert particle into full cell")
    }

    pub fn remove_particle_tenancy(&mut self, p_id: ParticleId) {
        let pos = self.particle(p_id).pos();
        let cell_id = self.get_cell_id(pos);
        self.delete_p_from_cell(p_id, cell_id);
    }

    pub fn insert_particle_tenancy(&mut self, p_id: ParticleId) {
        let pos = self.particle(p_id).pos();
        let cell_id = self.get_cell_id(pos);
        self.insert_p_into_cell(p_id, cell_id);
    }

    // TODO: make pretty
    // TODO: optimize since we know x_off and y_off will only ever be 1 or -1?
    pub fn get_neighbor_id(&self, id: CellId, x_off: i32, y_off: i32) -> CellId {
        let x = id / self.cells_per_axis[1];
        let y = id - x * self.cells_per_axis[1];
        let x = x as i32;
        let y = y as i32;

        let new_x = map_id_into_range(x + x_off, 0, self.cells_per_axis[0] as i32);
        let new_y = map_id_into_range(y + y_off, 0, self.cells_per_axis[1] as i32);

        debug_assert!(new_x >= 0 && new_x < self.cells_per_axis[0] as i32);
        debug_assert!(new_y >= 0 && new_y < self.cells_per_axis[1] as i32);
        let new_x = new_x as CellId;
        let new_y = new_y as CellId;

        new_x * self.cells_per_axis[1] + new_y
    }

    // ugly function, but surprisingly effective
    pub fn get_neighbors(&self, pos: Position) -> Vec<ParticleId> {
        let mut r: Vec<ParticleId> = Vec::new();

        let center_id = self.get_cell_id(pos);

        let up_left_id = self.get_neighbor_id(center_id, -1, 1);
        let up_id = self.get_neighbor_id(center_id, 0, 1);
        let up_right_id = self.get_neighbor_id(center_id, 1, 1);
        let left_id = self.get_neighbor_id(center_id, -1, 0);
        let right_id = self.get_neighbor_id(center_id, 1, 0);
        let bottom_left_id = self.get_neighbor_id(center_id, -1, -1);
        let bottom_id = self.get_neighbor_id(center_id, 0, -1);
        let bottom_right_id = self.get_neighbor_id(center_id, 1, -1);

        if self.tenants[up_left_id] != 0 {
            r.extend(self.cells[up_left_id]);
        }
        if self.tenants[up_id] != 0 {
            r.extend(self.cells[up_id]);
        }
        if self.tenants[up_right_id] != 0 {
            r.extend(self.cells[up_right_id]);
        }

        if self.tenants[left_id] != 0 {
            r.extend(self.cells[left_id]);
        }
        if self.tenants[center_id] != 0 {
            r.extend(self.cells[center_id]);
        }
        if self.tenants[right_id] != 0 {
            r.extend(self.cells[right_id]);
        }

        if self.tenants[bottom_left_id] != 0 {
            r.extend(self.cells[bottom_left_id]);
        }
        if self.tenants[bottom_id] != 0 {
            r.extend(self.cells[bottom_id]);
        }
        if self.tenants[bottom_right_id] != 0 {
            r.extend(self.cells[bottom_right_id]);
        }

        // TODO: remove this iter cloned thing
        let result = r
            .iter()
            .filter(|&&p_id| p_id != ParticleId::MAX)
            .cloned()
            .collect();

        result
    }

    pub fn min_x(&self) -> f64 {
        -0.5 * self.dimensions.x()
    }
    pub fn max_x(&self) -> f64 {
        0.5 * self.dimensions.x()
    }
    pub fn min_y(&self) -> f64 {
        -0.5 * self.dimensions.y()
    }
    pub fn max_y(&self) -> f64 {
        0.5 * self.dimensions.y()
    }

    pub fn cell_dimensions(&self) -> DimVec {
        self.cell_dimensions
    }

    pub fn cells_per_axis(&self) -> [usize; 2] {
        self.cells_per_axis
    }

    pub fn pos_in_box(&self, pos: Position) -> bool {
        pos.x() >= self.min_x()
            && pos.x() < self.max_x()
            && pos.y() >= self.min_y()
            && pos.y() < self.max_y()
    }

    // TODO: patch_center and send to simbox
    // sin_thetas/cos_thetas should be part of morphology
    // function then goes to simbox
    pub fn patch_center<P: IsParticle>(&self, p: &P, patch_idx: usize) -> Position {
        let or = p.or();
        let sin_theta = self.morphology(p).sin_theta(patch_idx);
        let cos_theta = self.morphology(p).cos_theta(patch_idx);
        let x = p.pos().x() + PARTICLE_RADIUS * (or.x() * cos_theta - or.y() * sin_theta);
        let y = p.pos().y() + PARTICLE_RADIUS * (or.x() * sin_theta + or.y() * cos_theta);
        self.map_pos_into_box(Position::new([x, y])) // simbox map into pos
    }

    pub fn map_pos_into_box(&self, pos: Position) -> Position {
        let x = map_into_range(pos.x(), self.min_x(), self.max_x());
        let y = map_into_range(pos.y(), self.min_y(), self.max_y());
        // debug_assert!(self.pos_in_box(pos));
        Position::new([x, y])
    }

    pub fn sep_in_box(&self, p0: Position, p1: Position) -> DimVec {
        self.map_pos_into_box(p0 - p1)
    }

    pub fn random_pos(&self, rng: &mut SmallRng) -> Position {
        let x = rng.gen_range(self.min_x()..self.max_x());
        let y = rng.gen_range(self.min_y()..self.max_y());
        Position::new([x, y])
    }
}
