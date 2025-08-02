use super::UnitCell;

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.6.3.6).
///
/// This is the trihexagonal tiling - alternating triangles and hexagons.
/// The unit cell consists of 2 vertices (3-sided and 6-sided polygons).
///
/// The graph structure is:
///     0 --- 1
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_6_3_6() -> UnitCell {
    super::util::create_line(vec![3, 6])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.4.6.4).
///
/// This is the snub hexagonal tiling.
/// The unit cell consists of 4 vertices (3-sided, 4-sided, 6-sided, 4-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2 --- 3
///     |                     |
///     +---------------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_4_6_4() -> UnitCell {
    super::util::create_cycle(vec![3, 4, 6, 4])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.3.4.3.4).
///
/// This is the elongated triangular tiling.
/// The unit cell consists of 5 vertices (3-sided, 3-sided, 4-sided, 3-sided, 4-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2 --- 3 --- 4
///     |                           |
///     +---------------------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_3_4_3_4() -> UnitCell {
    super::util::create_cycle(vec![3, 3, 4, 3, 4])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.3.3.4.4).
///
/// This is the snub square tiling.
/// The unit cell consists of 5 vertices (3-sided, 3-sided, 3-sided, 4-sided, 4-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2 --- 3 --- 4
///     |                           |
///     +---------------------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_3_3_4_4() -> UnitCell {
    super::util::create_cycle(vec![3, 3, 3, 4, 4])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.3.3.3.6).
///
/// This is the snub trihexagonal tiling.
/// The unit cell consists of 5 vertices (3-sided, 3-sided, 3-sided, 3-sided, 6-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2 --- 3 --- 4
///     |                           |
///     +---------------------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_3_3_3_6() -> UnitCell {
    super::util::create_cycle(vec![3, 3, 3, 3, 6])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (4.8.8).
///
/// This is the truncated square tiling.
/// The unit cell consists of 3 vertices (4-sided, 8-sided, 8-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2
///     |             |
///     +-------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_4_8_8() -> UnitCell {
    super::util::create_triangle(vec![4, 8, 8])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (3.12.12).
///
/// This is the truncated trihexagonal tiling.
/// The unit cell consists of 3 vertices (3-sided, 12-sided, 12-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2
///     |             |
///     +-------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_3_12_12() -> UnitCell {
    super::util::create_triangle(vec![3, 12, 12])
}

/// Returns a petgraph UnGraph describing the unit cell of semiregular tiling (4.6.12).
///
/// This is the truncated hexagonal tiling.
/// The unit cell consists of 3 vertices (4-sided, 6-sided, 12-sided polygons).
///
/// The graph structure is:
///     0 --- 1 --- 2
///     |             |
///     +-------------+
///
/// Returns: UnitCell where vertices represent polygon sizes
pub fn tiling_4_6_12() -> UnitCell {
    super::util::create_triangle(vec![4, 6, 12])
}
