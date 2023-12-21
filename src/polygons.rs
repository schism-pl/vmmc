use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::{
    particle::{IsParticle, Particle, ParticleId},
    vmmc::Vmmc,
};

pub type PolygonId = usize;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Polygon {
    // id: PolygonId,
    vertices: Vec<ParticleId>,
}

impl Polygon {
    fn new(vertices: Vec<ParticleId>) -> Self {
        Self {
            // id,
            vertices,
        }
    }

    pub fn vertices(&self) -> &[ParticleId] {
        &self.vertices
    }

    pub fn edge_iter(&self) -> PolygonEdgeIterator {
        PolygonEdgeIterator::new(&self.vertices)
    }
}

pub struct PolygonEdgeIterator<'a> {
    vertices: &'a [ParticleId],
    index: usize,
}

impl<'a> PolygonEdgeIterator<'a> {
    fn new(vertices: &'a [ParticleId]) -> Self {
        Self { vertices, index: 0 }
    }
}

impl<'a> Iterator for PolygonEdgeIterator<'a> {
    type Item = (ParticleId, ParticleId);
    fn next(&mut self) -> Option<Self::Item> {
        let idx = self.index;
        if idx == self.vertices.len() {
            return None;
        }

        let edge = if idx == self.vertices.len() - 1 {
            Some((self.vertices[idx], self.vertices[0]))
        } else {
            Some((self.vertices[idx], self.vertices[idx + 1]))
        };
        self.index += 1;
        edge
    }
}

// Note: counter-clockwise-edges uniquely map to polygons
// TODO: figure out the trig here
pub fn tightest_neighbor(vmmc: &Vmmc, p0: &Particle, p1: &Particle) -> Option<ParticleId> {
    let mut lowest_cos = 1.0; // we only want tightest_angle that goes clockwise
    let mut tightest = ParticleId::MAX;
    // vector form p0 to p1
    let dist0 = vmmc.simbox().sep_in_box(p1.pos(), p0.pos());
    let nd0 = dist0.scalar_div(dist0.norm());
    for &p2_id in vmmc.determine_interactions(p1).iter() {
        let p2 = vmmc.particle(p2_id);
        if p0.id() != p2.id() {
            // vector from p1 to p2
            let dist1 = vmmc.simbox().sep_in_box(p2.pos(), p1.pos());
            let nd1 = dist1.scalar_div(dist1.norm());

            // dot product of 2 unit vectors = sin(theta)
            let sin_theta = nd0.cross_prod(nd1);
            // dot product of 2 unit vectors = cos(theta)
            let cos_theta = nd0.dot_prod(nd1);
            if cos_theta < lowest_cos && sin_theta < 0.0 {
                lowest_cos = cos_theta;
                tightest = p2_id;
            }
        }
    }
    if tightest == ParticleId::MAX {
        None
    } else {
        Some(tightest)
    }
}

pub fn calc_polygons(vmmc: &Vmmc, max_vertices: usize) -> Vec<Polygon> {
    let mut next_p = HashMap::new(); // { (p_id, p_id) -> p_id}
    for p in vmmc.particles().iter() {
        for &neighbor_id in vmmc.determine_interactions(p).iter() {
            let neighbor = vmmc.particle(neighbor_id);
            // let next = tightest_neighbor(vmmc, p, neighbor);
            if let Some(next) = tightest_neighbor(vmmc, p, neighbor) {
                next_p.insert((p.id(), neighbor.id()), next);
            }
        }
    }

    let mut polygons = Vec::new();
    let mut edge_to_polygon = HashMap::new(); // { (p_id, p_id) -> PolygonId }
    for &(src, dst) in next_p.keys() {
        if edge_to_polygon.contains_key(&(src, dst)) {
            continue; // we already know this edge!
        }

        let mut visited = vec![src, dst];
        for idx in visited.len()..max_vertices + 1 {
            let curr_src = visited[idx - 2];
            let curr_dst = visited[idx - 1];
            if curr_dst == src {
                // polygon found!
                visited.pop(); // remove duplicate vertex
                let id = polygons.len();
                let polygon = Polygon::new(visited.clone());
                // propagate polygon info
                for edge in polygon.edge_iter() {
                    edge_to_polygon.insert(edge, id);
                }
                assert!(polygon.vertices.len() >= 3);
                polygons.push(polygon);
                break;
            }

            // if edge has
            let edge = (curr_src, curr_dst);
            if next_p.contains_key(&edge) {
                let next_vertex = next_p[&(curr_src, curr_dst)];
                visited.push(next_vertex);
            } else {
                break;
            }
        }
    }
    polygons
}

pub fn calc_polygon_count(vmmc: &Vmmc, max_vertices: usize) -> usize {
    calc_polygons(vmmc, max_vertices).len()
}

// maps shapes to bond distribution
pub fn calc_bond_distribution(vmmc: &Vmmc) -> Vec<Vec<usize>> {
    let mut bond_counts_per_shape = Vec::new();
    for shape in vmmc.simbox().shapes() {
        let mut bond_counts = vec![0; shape.patches().len() + 1];
        for p in vmmc.particles().iter() {
            let bond_count = vmmc.determine_interactions(p).len();
            bond_counts[bond_count] += 1;
        }
        bond_counts_per_shape.push(bond_counts);
    }
    bond_counts_per_shape
}
