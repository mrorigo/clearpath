//! The rectilinear router. Milestone M7; a placeholder that reports the honest reason.

use crate::error::PathPlanError;
use crate::geom::Point2D;
use crate::spline::CubicBezierSegment;

use alloc::vec::Vec;

/// A rectilinear polyline: consecutive points share an abscissa or an ordinate exactly.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrthogonalPolyline {
    /// The vertices, from start to goal.
    pub points: Vec<Point2D>,
}

impl OrthogonalPolyline {
    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the polyline has no vertices.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Whether every consecutive pair differs on exactly one axis.
    pub fn is_axis_aligned(&self) -> bool {
        self.points.windows(2).all(|w| {
            (w[0].x == w[1].x) ^ (w[0].y == w[1].y)
        })
    }
}

/// Not yet implemented.
pub fn route_orthogonal(
    _scratch: &mut crate::router::scratch::Scratch,
    _req: &crate::router::RouteRequest,
) -> Result<OrthogonalPolyline, PathPlanError> {
    Err(PathPlanError::NoPathFound)
}

const _: Option<CubicBezierSegment> = None;
