use petgraph::algo::all_simple_paths;
use petgraph::prelude::*;
use petgraph::{graph::UnGraph, Graph};

use crate::{particle::ParticleId, vmmc::Vmmc};

pub struct TilingGraph {
    // Undirected unweighted graph where vertices are ParticleId
    g: UnGraph<ParticleId, ()>,
}

impl TilingGraph {
    // TODO: can I make the graph Ungraph<(), ()> and have node index match particle index?
    pub fn from_vmmc(vmmc: &Vmmc) -> Self {
        let mut g = UnGraph::<ParticleId, ()>::new_undirected();

        // add nodes
        let mut nodes = Vec::new();
        for p in vmmc.particles() {
            let node = g.add_node(p.id());
            nodes.push(node);
        }

        // add edges for all interacting particles
        for p in vmmc.particles() {
            let interactions =
                vmmc.potential()
                    .determine_interactions(&vmmc.simbox(), &vmmc.particles(), p);
            for &neighbor_id in interactions.iter() {
                g.update_edge(nodes[p.id() as usize], nodes[neighbor_id as usize], ());
            }
        }

        Self { g }
    }

    // // returns polygons size 1-12
    // pub fn get_polygon_counts(&self) -> [usize; 12] {
    //     let result = [0; 12];
    //     for node_idx in self.g.node_indices(){
    //         for cycle in all_simple_paths(self.g, node_idx, node_idx, 2, 11){
    //             result[cycle.len()] += 1;
    //         }
    //     }
    //     result
    // }
}
