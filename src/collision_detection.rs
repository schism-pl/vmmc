use crate::{
    consts::NC_HALF_DIAG_LEN,
    particle::IsParticle,
    position::{PosDifference, Position},
    simbox::SimBox,
};
use std::f64::consts::FRAC_PI_2;

// TODO: needs lots of testing

impl SimBox {
    /// Calculates the vertices of the square.
    /// Postcondition: all returned positions are in box
    fn nc_vertices<P: IsParticle>(&self, p: &P) -> Vec<Position> {
        // TODO: assert is nanocube
        let diag_vec = p.or().scale_by(NC_HALF_DIAG_LEN);
        let diag_vec_rotated_90 = diag_vec.rotated_by(FRAC_PI_2);

        let c0 = self.map_pos_into_box(p.pos() + diag_vec);
        let c1 = self.map_pos_into_box(p.pos() + diag_vec_rotated_90);
        let c2 = self.map_pos_into_box(p.pos() - diag_vec);
        let c3 = self.map_pos_into_box(p.pos() - diag_vec_rotated_90);

        vec![c0, c1, c2, c3]
    }

    /// Checks if two squares overlap using the Separating Axis Theorem.
    pub fn nc_overlaps_with<P1: IsParticle, P2: IsParticle>(&self, a: &P1, b: &P2) -> bool {
        let a_vertices = self.nc_vertices(a);
        let b_vertices = self.nc_vertices(b);
        let a_normals = self.square_normals(&a_vertices);
        let b_normals = self.square_normals(&b_vertices);

        a_normals
            .iter()
            .all(|axis| overlap_on_axis(&a_vertices, &b_vertices, axis))
            && b_normals
                .iter()
                .all(|axis| overlap_on_axis(&a_vertices, &b_vertices, axis))
    }

    fn square_normals(&self, vertices: &[Position]) -> [PosDifference; 4] {
        let mut normals = [Position::new([0.0, 0.0]); 4];

        for i in 0..4 {
            let next_index = (i + 1) % 4; // wrap around to the first vertex after the last
            let side_vector = self.sep_in_box(vertices[next_index], vertices[i]);
            normals[i] = Position::new([side_vector.y(), -side_vector.x()]);
        }

        normals
    }
}

// TODO: no way this works
/// Checks if the projections of the vertices onto the axis overlap.
fn overlap_on_axis(vertices1: &[Position], vertices2: &[Position], axis: &Position) -> bool {
    let pa = vertices1[0].scalar_projected_on(*axis);
    let mut min_a = pa;
    let mut max_a = pa;

    let pb = vertices2[0].scalar_projected_on(*axis);
    let mut min_b = pb;
    let mut max_b = pb;

    for &vertex in vertices1.iter().skip(1) {
        let p = vertex.scalar_projected_on(*axis);
        min_a = min_a.min(p);
        max_a = max_a.max(p);
    }
    for &vertex in vertices2.iter().skip(1) {
        let p = vertex.scalar_projected_on(*axis);
        min_b = min_b.min(p);
        max_b = max_b.max(p);
    }

    max_a >= min_b && max_b >= min_a
}
