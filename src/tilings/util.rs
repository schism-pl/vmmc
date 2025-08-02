use petgraph::graph::UnGraph;

use super::UnitCell;

/// Creates a cycle graph with the given polygon sizes.
///
/// # Arguments
/// * `polygon_sizes` - Vector of polygon sizes (number of sides)
///
/// # Returns
/// UnitCell representing a cycle of polygons
pub fn create_cycle(polygon_sizes: Vec<usize>) -> UnitCell {
    let mut graph = UnGraph::new_undirected();

    if polygon_sizes.is_empty() {
        return graph;
    }

    // Add vertices
    let vertices: Vec<_> = polygon_sizes
        .iter()
        .map(|&size| graph.add_node(size))
        .collect();

    // Connect vertices in sequence
    for i in 0..vertices.len() {
        let current = vertices[i];
        let next = vertices[(i + 1) % vertices.len()];
        graph.add_edge(current, next, ());
    }

    graph
}

/// Creates a line graph with the given polygon sizes.
///
/// # Arguments
/// * `polygon_sizes` - Vector of polygon sizes (number of sides)
///
/// # Returns
/// UnitCell representing a line of polygons
pub fn create_line(polygon_sizes: Vec<usize>) -> UnitCell {
    let mut graph = UnGraph::new_undirected();

    if polygon_sizes.is_empty() {
        return graph;
    }

    // Add vertices
    let vertices: Vec<_> = polygon_sizes
        .iter()
        .map(|&size| graph.add_node(size))
        .collect();

    // Connect vertices in sequence (no closing edge)
    for i in 0..vertices.len() - 1 {
        graph.add_edge(vertices[i], vertices[i + 1], ());
    }

    graph
}

/// Creates a triangle graph with the given polygon sizes.
///
/// # Arguments
/// * `polygon_sizes` - Vector of exactly 3 polygon sizes
///
/// # Returns
/// UnitCell representing a triangle of polygons
pub fn create_triangle(polygon_sizes: Vec<usize>) -> UnitCell {
    assert_eq!(
        polygon_sizes.len(),
        3,
        "Triangle must have exactly 3 vertices"
    );
    create_cycle(polygon_sizes)
}

/// Creates a square graph with the given polygon sizes.
///
/// # Arguments
/// * `polygon_sizes` - Vector of exactly 4 polygon sizes
///
/// # Returns
/// UnitCell representing a square of polygons
pub fn create_square(polygon_sizes: Vec<usize>) -> UnitCell {
    assert_eq!(
        polygon_sizes.len(),
        4,
        "Square must have exactly 4 vertices"
    );
    create_cycle(polygon_sizes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_helper_functions() {
        // Test cycle creation
        let cycle = create_cycle(vec![3, 4, 6, 4]);
        assert_eq!(cycle.node_count(), 4);
        assert_eq!(cycle.edge_count(), 4);

        // Test line creation
        let line = create_line(vec![3, 6]);
        assert_eq!(line.node_count(), 2);
        assert_eq!(line.edge_count(), 1);

        // Test triangle creation
        let triangle = create_triangle(vec![6, 6, 6]);
        assert_eq!(triangle.node_count(), 3);
        assert_eq!(triangle.edge_count(), 3);

        // Test square creation
        let square = create_square(vec![4, 4, 4, 4]);
        assert_eq!(square.node_count(), 4);
        assert_eq!(square.edge_count(), 4);
    }
}
