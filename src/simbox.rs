use rand::{rngs::SmallRng, Rng};

use crate::consts::{DIMENSION, MAX_PARTICLES_PER_CELL};
use crate::position::Orientation;
use crate::{
    particle::{IsParticle, Particle, ParticleId},
    position::{DimVec, Position},
};

type CellId = usize;
type Cell = [ParticleId; MAX_PARTICLES_PER_CELL];

// cells_per_axis * cells_per_axis
type CellGrid = Vec<Cell>; //

// TODO: test boundary conditions extensively
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
pub struct SimBox {
    dimensions: DimVec,
    cells_per_axis: [usize; DIMENSION],
    cell_dimensions: DimVec,
    particles: Vec<Particle>,
    cells: CellGrid,
    tenants: Vec<u8>,
}

// x_idx * cells_per_axis[1] * MAX_PARTICLES_PER_CELL + y_idx * MAX_PARTICLES_PER_CELL
// ParticleIdx::MAX = empty
impl SimBox {
    pub fn new(
        dimensions: [f64; DIMENSION],
        cells_per_axis: [usize; DIMENSION],
        cell_dimensions: [f64; DIMENSION],
        particles: Vec<Particle>,
    ) -> Self {
        // assert_eq!(dimensions = cells_per_axis * cell_dimensions)
        // let particles = Vec::new(); // TODO: initalize!
        // let tenants = CellGrid::new(cells_per_axis, particles);
        let cells = vec![[ParticleId::MAX; 4]; cells_per_axis[0] * cells_per_axis[1]];
        let tenants = vec![0_u8; cells.len()];
        let mut r = Self {
            dimensions: DimVec::new(dimensions),
            cells_per_axis,
            cell_dimensions: DimVec::new(cell_dimensions),
            particles,
            cells,
            tenants,
        };

        for p_id in 0..r.particles.len() {
            // for p in r.particles.iter() {
            let cell_id = r.get_cell_id(r.particles[p_id].pos());
            r.insert_p_into_cell(p_id as u16, cell_id);
        }
        r
    }

    pub fn new_empty(
        dimensions: [f64; DIMENSION],
        cells_per_axis: [usize; DIMENSION],
        cell_dimensions: [f64; DIMENSION],
    ) -> Self {
        Self::new(dimensions, cells_per_axis, cell_dimensions, Vec::new())
    }

    // TODO: dedup with other impl of randomized_particles and overlaps
    // check if a particle overlaps any other particles
    fn overlaps(&self, pos: &Position, particles: &[Particle]) -> bool {
        for other in particles.iter() {
            // dist < 1.0 (hard sphere radius)
            if self.sep_in_box(*pos, other.pos()).norm() < 1.0 {
                return true;
            }
        }
        false
    }

    // TODO: dedup with other impl of randomized_particles and overlaps
    pub fn new_with_randomized_particles(
        dimensions: [f64; DIMENSION],
        cells_per_axis: [usize; DIMENSION],
        cell_dimensions: [f64; DIMENSION],
        n: usize,
        rng: &mut SmallRng,
    ) -> Self {
        // We initilize dummy simbox with no particles
        let simbox = Self::new_empty(dimensions, cells_per_axis, cell_dimensions);

        let mut particles = Vec::new();
        for idx in 0..n {
            let mut pos = simbox.random_pos(rng);
            while simbox.overlaps(&pos, &particles) {
                pos = simbox.random_pos(rng);
            }
            let or = Orientation::unit_vector(rng);
            let particle = Particle::new(idx as u16, pos, or);
            particles.push(particle);
        }
        // generate the real simbox
        SimBox::new(dimensions, cells_per_axis, cell_dimensions, particles)
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    pub fn particles_mut(&mut self) -> &mut [Particle] {
        &mut self.particles
    }

    pub fn get_cell_id(&self, pos: Position) -> CellId {
        // adjust so all indexes are positive
        // indexes start at bottom left and scan bottom to top
        let x_idx =
            ((pos.x() + (self.dimensions.x() / 2.0)) / self.cell_dimensions.x()).floor() as usize;
        let y_idx =
            ((pos.y() + (self.dimensions.y() / 2.0)) / self.cell_dimensions.y()).floor() as usize;
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
        let neighbor_pos = p.pos().shifted_by(x, y);
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

    pub fn move_particle_tenancy(
        &mut self,
        p_id: ParticleId,
        old_pos: Position,
        new_pos: Position,
    ) {
        let old_cell_id = self.get_cell_id(old_pos);
        let new_cell_id = self.get_cell_id(new_pos);
        if old_cell_id != new_cell_id {
            // 1. remove particleId from previous location
            self.delete_p_from_cell(p_id, old_cell_id);
            // 2. add particle to its new location
            self.insert_p_into_cell(p_id, new_cell_id);
        }
    }

    // TODO: make pretty
    // TODO: optimize since we know x_off and y_off will only ever be 1 or -1?
    pub fn get_neighbor_id(&self, id: CellId, x_off: i32, y_off: i32) -> CellId {
        // let num_cells = self.cell_dimensions[0] * self.cell_dimensions[1];

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

    // TODO: improve by using a neighbors data structure?
    pub fn get_neighbors(&self, p: &Particle) -> Vec<ParticleId> {
        let mut r: Vec<ParticleId> = Vec::new();

        let center_id = self.get_cell_id(p.pos());

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
    pub fn min_z(&self) -> f64 {
        -0.5 * self.dimensions.z()
    }
    pub fn max_z(&self) -> f64 {
        0.5 * self.dimensions.z()
    }

    pub fn pos_in_box(&self, pos: Position) -> bool {
        if DIMENSION == 2 {
            pos.x() >= self.min_x()
                && pos.x() < self.max_x()
                && pos.y() >= self.min_y()
                && pos.y() < self.max_y()
        } else if DIMENSION == 3 {
            pos.x() >= self.min_x()
                && pos.x() < self.max_x()
                && pos.y() >= self.min_y()
                && pos.y() < self.max_y()
                && pos.z() >= self.min_z()
                && pos.z() < self.max_z()
        } else {
            panic!("Dimension is not 2 or 3")
        }
    }

    pub fn map_pos_into_box(&self, pos: Position) -> Position {
        if DIMENSION == 2 {
            let x = map_into_range(pos.x(), self.min_x(), self.max_x());
            let y = map_into_range(pos.y(), self.min_y(), self.max_y());
            // debug_assert!(self.pos_in_box(pos));
            Position::new([x, y])
        } else if DIMENSION == 3 {
            let _x = map_into_range(pos.x(), self.min_x(), self.max_x());
            let _y = map_into_range(pos.y(), self.min_y(), self.max_y());
            let _z = map_into_range(pos.z(), self.min_z(), self.max_z());
            panic!("DIMENSION IS SET TO 2")
            // Position::new([x,y,z])
        } else {
            panic!("Dimension is not 2 or 3")
        }
    }

    pub fn sep_in_box(&self, p0: Position, p1: Position) -> DimVec {
        assert_eq!(DIMENSION, 2);
        self.map_pos_into_box(p0 - p1)
    }

    pub fn random_pos(&self, rng: &mut SmallRng) -> Position {
        let x = rng.gen_range(self.min_x()..self.max_x());
        let y = rng.gen_range(self.min_y()..self.max_y());
        Position::new([x, y])
    }
}
