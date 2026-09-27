//! Low-level geometric primitives and robust predicates.

pub mod clip;
pub mod inflate;
pub mod point;
pub mod polygon;
pub mod predicates;

pub use point::Point2D;
pub use polygon::{BoundingBox, Polygon};
pub use predicates::{Orientation, incircle, orient2d};

/// Tolerance used only for non-topological "is this effectively zero" tests.
///
/// Never used for orientation, containment, or ordering decisions — those go through
/// [`predicates::orient2d`], which is exact.
pub const EPS: f64 = 1e-12;
