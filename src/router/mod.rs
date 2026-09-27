//! The public routing API.

pub mod orthogonal;
pub mod scratch;
pub mod smooth;

use alloc::vec::Vec;

pub use orthogonal::OrthogonalPolyline;
pub use smooth::route_smooth;

use crate::error::PathPlanError;
use crate::geom::clearance::Clearance;
use crate::geom::{BoundingBox, FreeSpace, Point2D, Polygon};
use crate::spline::CubicBezierSegment;

/// A route endpoint, with an optional forced departure or arrival direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortConstraint {
    /// The point the route starts or ends at.
    pub point: Point2D,
    /// A forced tangent *direction*; the magnitude is derived from the neighbouring knot spacing.
    ///
    /// Honoured only if the resulting tangent is admissible at that knot (section 4.6);
    /// otherwise the tangent is projected and the direction is not met exactly. `None` leaves the
    /// direction to the solver.
    pub direction: Option<Point2D>,
}

impl PortConstraint {
    /// A port with no forced direction.
    pub fn free(point: Point2D) -> Self {
        Self { point, direction: None }
    }

    /// A port with a forced direction. The direction must be finite and non-zero; that is checked
    /// when the request is validated, not here, so a zero direction is reported as one
    /// `InvalidConfig` rather than a second error kind.
    pub fn directed(point: Point2D, direction: Point2D) -> Self {
        Self { point, direction: Some(direction) }
    }
}

/// Tunables. Defaults are the ones a caller gets from [`Config::default`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Config {
    /// Clearance required around every obstacle and inside the workspace. Non-negative and finite.
    pub margin: f64,
    /// The tangent solver's bias towards the seed directions. Non-negative and finite.
    ///
    /// Zero is *singular*: any constant tangent has identically zero curvature, so the
    /// minimum-curvature objective is flat along that whole direction and the solve has no unique
    /// answer. The default is 1.0; zero is accepted and falls back to the seeds.
    pub tangent_bias: f64,
    /// How many damping steps the containment repair may take. Greater than zero.
    pub max_repair_iters: u32,
    /// The per-step damping factor, in `(0, 1)`.
    pub repair_dampening: f64,
    /// The cap on arc segments per offset corner. Greater than zero.
    pub max_arc_segments: u32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            margin: 0.0,
            tangent_bias: 1.0,
            max_repair_iters: 16,
            repair_dampening: 0.5,
            max_arc_segments: 64,
        }
    }
}

impl Config {
    /// Checks the ranges, reporting which one is wrong.
    pub fn validate(&self) -> Result<(), PathPlanError> {
        if !self.margin.is_finite() || self.margin < 0.0 {
            return Err(PathPlanError::InvalidConfig("margin must be finite and >= 0"));
        }
        if !self.tangent_bias.is_finite() || self.tangent_bias < 0.0 {
            return Err(PathPlanError::InvalidConfig("tangent_bias must be finite and >= 0"));
        }
        if self.max_repair_iters == 0 {
            return Err(PathPlanError::InvalidConfig("max_repair_iters must be > 0"));
        }
        if !(self.repair_dampening > 0.0 && self.repair_dampening < 1.0) {
            return Err(PathPlanError::InvalidConfig("repair_dampening must be in (0, 1)"));
        }
        if self.max_arc_segments == 0 {
            return Err(PathPlanError::InvalidConfig("max_arc_segments must be > 0"));
        }
        Ok(())
    }
}

/// A routing request.
#[derive(Clone, Debug)]
pub struct RouteRequest {
    /// The region the route stays inside, and which the margin erodes.
    pub workspace: BoundingBox,
    /// Where the route starts.
    pub start: PortConstraint,
    /// Where the route ends.
    pub goal: PortConstraint,
    /// The forbidden polygons. Rectilinear for a nonzero `margin` (section 3.3).
    pub obstacles: Vec<Polygon>,
    /// The tunables.
    pub config: Config,
}

impl RouteRequest {
    /// A request with no obstacles and default tunables.
    pub fn new(workspace: BoundingBox, start: Point2D, goal: Point2D) -> Self {
        Self {
            workspace,
            start: PortConstraint::free(start),
            goal: PortConstraint::free(goal),
            obstacles: Vec::new(),
            config: Config::default(),
        }
    }

    /// The free space of the request, before the margin is applied.
    pub fn free_space(&self) -> Result<FreeSpace, PathPlanError> {
        FreeSpace::new(self.workspace, self.obstacles.clone())
    }

    /// The clearance predicate for the inflated domain this request routes in.
    ///
    /// The *un-eroded* workspace goes in: [`Clearance::new`] applies the margin itself. Eroding
    /// here as well shrank the workspace by twice the margin, which is not a safety problem — the
    /// obstacle test is what carries the margin — but it rejects routes that do fit, by reporting
    /// endpoints as outside the workspace.
    pub fn clearance(&self) -> Clearance {
        Clearance::new(&self.obstacles, self.workspace, self.config.margin)
    }

    /// Checks the request's invariants, before any geometry is built.
    pub fn validate(&self) -> Result<(), PathPlanError> {
        self.config.validate()?;
        if self.workspace.is_degenerate() {
            return Err(PathPlanError::DegenerateWorkspace);
        }
        for port in [&self.start, &self.goal] {
            if !port.point.is_finite() {
                return Err(PathPlanError::InvalidConfig("endpoint must be finite"));
            }
            if let Some(d) = port.direction {
                if !d.is_finite() || d.is_zero() {
                    return Err(PathPlanError::InvalidConfig("port direction must be finite and non-zero"));
                }
            }
        }
        if self.start.point == self.goal.point {
            return Err(PathPlanError::DegenerateEndpoints);
        }
        for port in [&self.start, &self.goal] {
            if !self.workspace.contains(port.point) {
                return Err(PathPlanError::EndpointOutsideWorkspace { point: port.point });
            }
            for obstacle in &self.obstacles {
                // A point exactly on a boundary counts as inside: an endpoint on an obstacle is
                // not a request worth answering.
                if obstacle.contains(port.point) {
                    return Err(PathPlanError::EndpointInObstacle { point: port.point });
                }
            }
        }
        Ok(())
    }
}

/// Which route family to produce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RouteKind {
    /// A `C^1` chain of cubic Bezier segments.
    Smooth,
    /// A rectilinear polyline.
    Orthogonal,
}

/// Either kind of route.
#[derive(Clone, Debug, PartialEq)]
pub enum Route {
    /// Cubic Bezier segments.
    Smooth(Vec<CubicBezierSegment>),
    /// A rectilinear polyline.
    Orthogonal(OrthogonalPolyline),
}

/// A reusable planner. Holds scratch buffers so a warm query does not allocate.
#[derive(Debug, Default)]
pub struct PathPlanner {
    pub(crate) scratch: scratch::Scratch,
}

impl PathPlanner {
    /// A planner with empty scratch.
    pub fn new() -> Self {
        Self { scratch: scratch::Scratch::default() }
    }

    /// A `C^1` chain of cubic Bezier segments from `start` to `goal`.
    pub fn route_smooth(
        &mut self,
        req: &RouteRequest,
    ) -> Result<Vec<CubicBezierSegment>, PathPlanError> {
        route_smooth(&mut self.scratch, req)
    }

    /// A rectilinear polyline from `start` to `goal`.
    pub fn route_orthogonal(
        &mut self,
        req: &RouteRequest,
    ) -> Result<OrthogonalPolyline, PathPlanError> {
        orthogonal::route_orthogonal(&mut self.scratch, req)
    }

    /// Either kind of route.
    pub fn route(&mut self, req: &RouteRequest, kind: RouteKind) -> Result<Route, PathPlanError> {
        match kind {
            RouteKind::Smooth => self.route_smooth(req).map(Route::Smooth),
            RouteKind::Orthogonal => self.route_orthogonal(req).map(Route::Orthogonal),
        }
    }
}
