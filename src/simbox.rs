use rand::{rngs::ThreadRng, Rng};

use crate::params::MAX_PARTICLES_PER_CELL;
use crate::{
    params::DIMENSION,
    particle::{IsParticle, Particle, ParticleId},
    position::{DimVec, Position},
};

type CellId = usize;
type Cell = [ParticleId; MAX_PARTICLES_PER_CELL];
// cells_per_axis * cells_per_axis * MAX_PARTICLES_PER_CELL
// type CellGrid = Array<ParticleId, Ix3>;

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

// fn in_cyclic_range(p: f64, lower: f64, upper: f64) -> bool {
//     if lower < upper {
//         // no wrap around
//         p >= lower && p < upper
//     } else {
//         // wrap around
//         p >= lower || p < upper
//     }
// }

// struct CellGrid {
//     inner: Vec<ParticleId>
// }

// impl CellGrid {
//     fn new(cells_per_axis: [usize; DIMENSION], perticles: &[Particle]) -> Self {
//         // let tenants = Array::<ParticleId, Ix3>::zeros((cells_per_axis[0], cells_per_axis[1], MAX_PARTICLES_PER_CELL));
//         Self{inner: [] }
//     }
// }

// assumed to be periodic in every dimension
// dimension[i] = cell_dimensions[i] * cells_per_axis[i]
// Note: cell_dimensions * cells_per_axis must = dimensions
pub struct SimBox {
    dimensions: DimVec,
    cells_per_axis: [usize; DIMENSION],
    cell_dimensions: DimVec,
    cells: CellGrid,
}

// x_idx * cells_per_axis[1] * MAX_PARTICLES_PER_CELL + y_idx * MAX_PARTICLES_PER_CELL
// ParticleIdx::MAX = empty
impl SimBox {
    pub fn new(
        dimensions: [f64; DIMENSION],
        cells_per_axis: [usize; DIMENSION],
        cell_dimensions: [f64; DIMENSION],
        particles: &[Particle],
    ) -> Self {
        // assert_eq!(dimensions = cells_per_axis * cell_dimensions)

        // let tenants = CellGrid::new(cells_per_axis, particles);
        let cells = vec![[ParticleId::MAX; 4]; cells_per_axis[0] * cells_per_axis[1]];
        let mut r = Self {
            dimensions: DimVec::new(dimensions),
            cells_per_axis,
            cell_dimensions: DimVec::new(cell_dimensions),
            cells,
        };

        for p in particles.iter() {
            let cell_id = r.get_cell_id(p.pos());
            r.insert_p_into_cell(p.id(), cell_id);
        }
        r
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    // TODO: does negative wrong
    pub fn get_cell_id(&self, pos: Position) -> CellId {
        // adjust so all indexes are positive
        // indexes start at bottom left and scan bottom to top
        let x_idx =
            ((pos.x() + (self.dimensions.x() / 2.0)) / self.cell_dimensions.x()).floor() as usize;
        let y_idx =
            ((pos.y() + (self.dimensions.y() / 2.0)) / self.cell_dimensions.y()).floor() as usize;
        // println!("get_cell_id({:?}) ({:?},{:?}) = {:?} cell = [{:?}:{:?}[{:?}{:?}]", pos, x_idx, y_idx, x_idx * self.cells_per_axis[1] + y_idx, x_idx as f64 * self.cell_dimensions.x() - (self.dimensions.x() / 2.0),  (x_idx+1) as f64 * self.cell_dimensions.x() - (self.dimensions.x() / 2.0), y_idx as f64 * self.cell_dimensions.y() - (self.dimensions.y() / 2.0),  (y_idx+1) as f64 * self.cell_dimensions.y() - (self.dimensions.y() / 2.0) );
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

    // TODO: improve by keeping counts
    // TODO: improve by directly using an iterator and removing cloned
    // TODO: improve by using a neighbors data structure?
    pub fn get_neighbors(&self, p: &Particle) -> Vec<ParticleId> {
        let mut r: Vec<ParticleId> = Vec::new();

        let x_dim = self.cell_dimensions.x();
        let y_dim = self.cell_dimensions.y();

        r.extend(self.get_neighbor(p, -x_dim, y_dim)); // up left
        r.extend(self.get_neighbor(p, 0.0, y_dim)); // up
        r.extend(self.get_neighbor(p, x_dim, y_dim)); // up right
        r.extend(self.get_neighbor(p, -x_dim, 0.0)); // left
        r.extend(self.get_cell(p.pos())); // centre
        r.extend(self.get_neighbor(p, x_dim, 0.0)); // right
        r.extend(self.get_neighbor(p, -x_dim, -y_dim)); // bottom left
        r.extend(self.get_neighbor(p, 0.0, -y_dim)); // bottom center
        r.extend(self.get_neighbor(p, x_dim, -y_dim)); // bottom right
        r.iter()
            .filter(|&&p_id| p_id != ParticleId::MAX)
            .cloned()
            .collect()

        // let neighborhood = self.get_neighborhood(particle);
        // particles
        //     .iter()
        //     .filter(|p| self.in_neighborhood(neighborhood, p))
        //     .map(|p| p.id())
        //     .collect()
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

    // pub fn get_cell(&self, particle: &Particle) -> Position {
    //     let pos = particle.pos();
    //     let cell_dims = self.cell_dimensions;
    //     let x = (pos.x() / cell_dims.x()).floor() * cell_dims.x();
    //     let y = (pos.y() / cell_dims.y()).floor() * cell_dims.y();
    //     self.map_pos_into_box(Position::new([x, y]))
    // }

    // // get bottom left of neighborhood
    // pub fn get_neighborhood(&self, particle: &Particle) -> Position {
    //     let cell = self.get_cell(particle);
    //     self.map_pos_into_box(cell - self.cell_dimensions)
    // }

    // // neighborhood is bottom left corner of neighborhood
    // pub fn in_neighborhood(&self, neighborhood: Position, p: &Particle) -> bool {
    //     let left = neighborhood.x();
    //     let right = left + self.cell_dimensions.x() * 3.0;
    //     let right = map_into_range(right, self.min_x(), self.max_x()); // wrap around
    //     let down = neighborhood.y();
    //     let up = neighborhood.y() + self.cell_dimensions.y() * 3.0;
    //     let up = map_into_range(up, self.min_y(), self.max_y()); // wrap around
    //                                                              // right or up isn't in cyclic box?
    //     in_cyclic_range(p.pos().x(), left, right) && in_cyclic_range(p.pos().y(), down, up)
    // }

    // pub fn get_neighbors(&self, particles: &[Particle], particle: &Particle) -> Vec<ParticleId> {
    //     let neighborhood = self.get_neighborhood(particle);
    //     particles
    //         .iter()
    //         .filter(|p| self.in_neighborhood(neighborhood, p))
    //         .map(|p| p.id())
    //         .collect()
    // }

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
        // let x = map_into_range(sep.x(), self.min_x() *0.5, self.max_x() * 0.5);
        // let y = map_into_range(sep.y(), self.min_y() *0.5, self.max_y() * 0.5);
        // Position::new([x, y])
    }

    pub fn random_pos(&self, rng: &mut ThreadRng) -> Position {
        let x = rng.gen_range(self.min_x()..self.max_x());
        let y = rng.gen_range(self.min_y()..self.max_y());
        Position::new([x, y])
    }
}
