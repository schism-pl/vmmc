// Tilings module
// This module contains functionality related to tiling patterns and unit cell analysis

use petgraph::graph::UnGraph;

/// Type alias for unit cell graphs
pub type UnitCell = UnGraph<usize, ()>;

// Submodules
pub mod regular_tilings;
pub mod semiregular_tilings;
pub mod two_vertex_tilings;
pub mod util;

// Re-export commonly used items
pub use regular_tilings::{tiling_3_6, tiling_4_4, tiling_6_3};
pub use semiregular_tilings::{
    tiling_3_12_12, tiling_3_3_3_3_6, tiling_3_3_3_4_4, tiling_3_3_4_3_4, tiling_3_4_6_4,
    tiling_3_6_3_6, tiling_4_6_12, tiling_4_8_8,
};
pub use two_vertex_tilings::{
    tiling_3_2_6_2_3_4_6, tiling_3_6_3_2_6_2, tiling_3_6_3_6_3_2_6_2, tiling_4_4_3_3_4_2,
};

/// Returns a tiling based on the canonical name string.
///
/// Supported canonical names:
/// Regular tilings:
/// - "4.4.4.4" or "4_4" -> Square tiling
/// - "3.3.3.3.3.3" or "3_6" -> Triangle tiling  
/// - "6.6.6" or "6_3" -> Hexagon tiling
///
/// Semiregular tilings:
/// - "3.6.3.6" -> Trihexagonal tiling
/// - "3.4.6.4" -> Snub hexagonal tiling
/// - "3.3.4.3.4" -> Elongated triangular tiling
/// - "3.3.3.4.4" -> Snub square tiling
/// - "3.3.3.3.6" -> Snub trihexagonal tiling
/// - "4.8.8" -> Truncated square tiling
/// - "3.12.12" -> Truncated trihexagonal tiling
/// - "4.6.12" -> Truncated hexagonal tiling
///
/// 2-vertex tilings:
/// - "4_4; 3_3.4_2" -> 2-vertex tiling with different vertex configurations
/// - "3.6.3.6; 3_2.6_2" -> 2-vertex tiling with trihexagonal and mixed configurations
/// - "3_6; 3_2.6_2" -> 2-vertex tiling with triangular and mixed configurations
/// - "3_2.6_2; 3_4.6" -> 2-vertex tiling with mixed triangle-hexagon configurations
///
/// Returns: UnitCell corresponding to the tiling, or None if the string is not recognized
pub fn tiling_from_str(name: &str) -> Option<UnitCell> {
    match name {
        // Regular tilings
        "4.4.4.4" | "4_4" => Some(tiling_4_4()),
        "3.3.3.3.3.3" | "3_6" => Some(tiling_3_6()),
        "6.6.6" | "6_3" => Some(tiling_6_3()),

        // Semiregular tilings
        "3.6.3.6" => Some(tiling_3_6_3_6()),
        "3.4.6.4" => Some(tiling_3_4_6_4()),
        "3.3.4.3.4" => Some(tiling_3_3_4_3_4()),
        "3.3.3.4.4" => Some(tiling_3_3_3_4_4()),
        "3.3.3.3.6" => Some(tiling_3_3_3_3_6()),
        "4.8.8" => Some(tiling_4_8_8()),
        "3.12.12" => Some(tiling_3_12_12()),
        "4.6.12" => Some(tiling_4_6_12()),

        // 2-vertex tilings
        "4_4; 3_3.4_2" => Some(tiling_4_4_3_3_4_2()),
        "3.6.3.6; 3_2.6_2" => Some(tiling_3_6_3_6_3_2_6_2()),
        "3_6; 3_2.6_2" => Some(tiling_3_6_3_2_6_2()),
        "3_2.6_2; 3_4.6" => Some(tiling_3_2_6_2_3_4_6()),

        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tiling_from_str() {
        // Test regular tilings
        assert_eq!(tiling_from_str("4.4.4.4"), Some(tiling_4_4()));
        assert_eq!(tiling_from_str("3.3.3.3.3.3"), Some(tiling_3_6()));
        assert_eq!(tiling_from_str("6.6.6"), Some(tiling_6_3()));

        // Test semiregular tilings
        assert_eq!(tiling_from_str("3.6.3.6"), Some(tiling_3_6_3_6()));
        assert_eq!(tiling_from_str("3.4.6.4"), Some(tiling_3_4_6_4()));
        assert_eq!(tiling_from_str("3.3.4.3.4"), Some(tiling_3_3_4_3_4()));
        assert_eq!(tiling_from_str("3.3.3.4.4"), Some(tiling_3_3_3_4_4()));
        assert_eq!(tiling_from_str("3.3.3.3.6"), Some(tiling_3_3_3_3_6()));
        assert_eq!(tiling_from_str("4.8.8"), Some(tiling_4_8_8()));
        assert_eq!(tiling_from_str("3.12.12"), Some(tiling_3_12_12()));
        assert_eq!(tiling_from_str("4.6.12"), Some(tiling_4_6_12()));

        // Test 2-vertex tilings
        assert_eq!(tiling_from_str("4_4; 3_3.4_2"), Some(tiling_4_4_3_3_4_2()));
        assert_eq!(
            tiling_from_str("3.6.3.6; 3_2.6_2"),
            Some(tiling_3_6_3_6_3_2_6_2())
        );
        assert_eq!(tiling_from_str("3_6; 3_2.6_2"), Some(tiling_3_6_3_2_6_2()));
        assert_eq!(
            tiling_from_str("3_2.6_2; 3_4.6"),
            Some(tiling_3_2_6_2_3_4_6())
        );

        // Test invalid input
        assert_eq!(tiling_from_str("invalid"), None);
    }
}
