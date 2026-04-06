use petgraph::graph::UnGraph;

use super::UnitCell;

/// Tiling #01
/// Returns a petgraph UnGraph describing the unit cell of a 2-vertex tiling (4_4; 3_3.4_2). 
///
/// This is a 2-vertex tiling where:
/// - Vertex 0 has configuration 4_4 (4-sided polygons, 4 meeting at vertex)
/// - Vertex 1 has configuration 3_3.4_2 (3-sided and 4-sided polygons)
///
/// The graph structure is:
///     0 --- 1
///
/// Returns: UnitCell where vertices represent different vertex configurations
pub fn tiling_4_4_3_3_4_2() -> UnitCell {
    super::util::create_line(vec![4, 3]) // Unitcell is a square and a triangle
}

/// Tiling #07
/// Primitive unitcell defined
/// Returns a petgraph UnGraph describing the unit cell of a 2-vertex tiling (3.6.3.6; 3_2.6_2).
///
/// This is a 2-vertex tiling where:
/// - Vertex 0 has configuration 3.6.3.6 (trihexagonal pattern)
/// - Vertex 1 has configuration 3_2.6_2 (triangles and hexagons)
///
/// The graph structure is:
///     0 --- 1 --- 2 (Vertex 0 only connects to Vertex 1)
///                   (Vertex 1 connects to 2 difference vertices)
///
/// Returns: UnitCell where vertices represent different vertex configurations
pub fn tiling_3_6_3_6_3_2_6_2() -> UnitCell {
    super::util::create_line(vec![3, 6, 3]) // Triangle-hexagon-triangle pattern
}

/// Tiling #11
/// Returns a petgraph UnGraph describing the unit cell of a 2-vertex tiling (3_6; 3_2.6_2).
///
/// This is a 2-vertex tiling where:
/// - Vertex 0 has configuration 3_6 (triangles, 6 meeting at vertex)
/// - Vertex 1 has configuration 3_2.6_2 (triangles and hexagons)
///
/// The graph structure is:
///     0 (3) --- 1 (6) (Vertex 0 connects to 6 Vertex 1's)
///     2 (3) --- 1 (6) (Vertex 1 connects to 3 Vertex 1's and 1 Vertex 0) 
///     3 (3) --- 1 (6)
///
/// Returns: UnitCell where vertices represent different vertex configurations
pub fn tiling_3_6_3_2_6_2() -> UnitCell {
    let mut graph = UnGraph::new_undirected();

    // Create vertices: 3 triangles and 1 hexagon
    let triangle1 = graph.add_node(3);
    let hexagon = graph.add_node(6);
    let triangle2 = graph.add_node(3);
    let triangle3 = graph.add_node(3);

    // Connect all triangles to the shared hexagon
    graph.add_edge(triangle1, hexagon, ());
    graph.add_edge(triangle2, hexagon, ());
    graph.add_edge(triangle3, hexagon, ());

    graph
}

/// Tiling #14
/// Returns a petgraph UnGraph describing the unit cell of a 2-vertex tiling (3_2.6_2; 3_4.6).
///
/// This is a 2-vertex tiling where:
/// - Vertex 0 has configuration 3_2.6_2 (2 triangles and 2 hexagons)
/// - Vertex 1 has configuration 3_4.6 (4 triangles and 1 hexagon)
///
/// The graph structure is:
///     0 (3) --- 4 (6) (Vertex 0 has 3 triangles and connects to 4 other Vertex 1s)
///     1 (3) --- 4 (6)
///     2 (3) --- 4 (6)
///     3 (3) --- 4 (6)
///
/// Returns: UnitCell where vertices represent different vertex configurations
pub fn tiling_3_2_6_2_3_4_6() -> UnitCell {
    let mut graph = UnGraph::new_undirected();

    // Create vertices: 4 triangles and 1 hexagon
    let triangle1 = graph.add_node(3);
    let triangle2 = graph.add_node(3);
    let triangle3 = graph.add_node(3);
    let triangle4 = graph.add_node(3);
    let hexagon = graph.add_node(6);

    // Connect all triangles to the shared hexagon
    graph.add_edge(triangle1, hexagon, ());
    graph.add_edge(triangle2, hexagon, ());
    graph.add_edge(triangle3, hexagon, ());
    graph.add_edge(triangle4, hexagon, ());

    graph
}


/// Tiling #15
/// Returns a petgraph UnGraph describing the unit cell of a 2-vertex tiling (3_6; 3_2.4.3.4)
///
/// This is a 2-vertex tiling where:
/// - Vertex 0 has configuration 3_6 (6 triangles)
/// - Vertex 1 has configuration 3_2.4.3.4 (2 triangles, 1 square, 1 triangle, and 1 square)
///
/// The graph structure is:
///     0 

/// Unitcells undergoing testing

/// Tiling #07
/// Full unitcell defined
/// Under development
pub fn tiling_07() -> UnitCell {
    super::util::create_line(vec![6, 3, 6, 3]) // Hexagon-triangle-hexagon-triangle pattern
}

#[cfg(test)]
mod tests {
    use super::*;

/// Node count = number of polygons
/// Edge count = number of interfacing edges

    #[test]
    fn test_tiling_4_4_3_3_4_2() {
        let tiling = tiling_4_4_3_3_4_2();
        assert_eq!(tiling.node_count(), 2);
        assert_eq!(tiling.edge_count(), 1);
    }

    #[test]
    fn test_tiling_3_6_3_6_3_2_6_2() {
        let tiling = tiling_3_6_3_6_3_2_6_2();
        assert_eq!(tiling.node_count(), 3);
        assert_eq!(tiling.edge_count(), 2);
    }


    #[test]
    fn test_tiling_3_6_3_2_6_2() {
        let tiling = tiling_3_6_3_2_6_2();
        assert_eq!(tiling.node_count(), 4);
        assert_eq!(tiling.edge_count(), 3);
    }

    #[test]
    fn test_tiling_3_2_6_2_3_4_6() {
        let tiling = tiling_3_2_6_2_3_4_6();
        assert_eq!(tiling.node_count(), 5);
        assert_eq!(tiling.edge_count(), 4);
    }

    //Under development
    #[test]
    fn test_tiling_07() {
        let tiling = tiling_07();
        assert_eq!(tiling.node_count(), 4);
        assert_eq!(tiling.edge_count(), 4);
    }


}
