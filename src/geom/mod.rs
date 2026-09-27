//! Low-level geometric primitives and robust predicates.

use alloc::vec::Vec;

use crate::error::PathPlanError;

pub mod clearance;
pub mod point;
pub mod polygon;
pub mod predicates;

pub use point::Point2D;
pub use polygon::{BoundingBox, Polygon};
pub use predicates::{Orientation, incircle, orient2d};

/// The free space a route is planned in: the workspace, less the obstacles.
///
/// `margin` is **not** baked into this type. It is applied one step later, by eroding each
/// decomposition cell, because that is exact and the alternative — offsetting the obstacle
/// polygons — is not. See `docs/SPEC.md` section 3.3.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSpace {
    /// The workspace bounding box. A route stays inside it.
    pub workspace: BoundingBox,
    /// The obstacles. Their interiors are forbidden.
    pub obstacles: Vec<Polygon>,
}

impl FreeSpace {
    /// Constructs a free space, validating the workspace and the margin-independent inputs.
    ///
    /// # Errors
    /// Returns [`PathPlanError::DegenerateWorkspace`] if the bounding box has no positive extent.
    pub fn new(workspace: BoundingBox, obstacles: Vec<Polygon>) -> Result<Self, PathPlanError> {
        if workspace.is_degenerate() {
            return Err(PathPlanError::DegenerateWorkspace);
        }
        Ok(Self { workspace, obstacles })
    }

    /// Whether `p` is in the free space, ignoring `margin`.
    pub fn contains(&self, p: Point2D) -> bool {
        self.workspace.contains(p) && !self.obstacles.iter().any(|o| o.contains(p))
    }
}

/// Tolerance used only for non-topological "is this effectively zero" tests.
///
/// Never used for orientation, containment, or ordering decisions — those go through
/// [`predicates::orient2d`], which is exact.
pub const EPS: f64 = 1e-12;
