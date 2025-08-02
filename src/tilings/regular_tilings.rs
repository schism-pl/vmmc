use super::UnitCell;

/// Returns a petgraph UnGraph describing the unit cell of a regular square tiling (4.4.4.4).
///
/// The square tiling unit cell consists of 4 vertices (representing 4-sided polygons)
/// connected in a square pattern. Each vertex represents a square tile.
///
/// The graph structure is:
///     0 --- 1
///     |     |
///     |     |
///     3 --- 2
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_4_4() -> UnitCell {
    super::util::create_square(vec![4, 4, 4, 4])
}

/// Returns a petgraph UnGraph describing the unit cell of a regular triangle tiling (3.3.3.3.3.3).
///
/// The triangle tiling unit cell consists of 6 vertices (representing 3-sided polygons)
/// connected in a triangular pattern. Each vertex represents a triangular tile.
///
/// The graph structure is:
///       0
///      / \
///     5   1
///    /     \
///   4       2
///    \     /
///     3 - 3
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_6() -> UnitCell {
    super::util::create_cycle(vec![3, 3, 3, 3, 3, 3])
}

/// Returns a petgraph UnGraph describing the unit cell of a regular hexagon tiling (6.6.6).
///
/// The hexagon tiling unit cell consists of 3 vertices (representing 6-sided polygons)
/// connected in a triangular pattern. Each vertex represents a hexagonal tile.
///
/// The graph structure is:
///     0
///    / \
///   2 - 1
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_6_3() -> UnitCell {
    super::util::create_triangle(vec![6, 6, 6])
}
