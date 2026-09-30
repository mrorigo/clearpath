//! Error type for the crate.

use alloc::fmt;

use crate::geom::Point2D;

/// Why an obstacle polygon was rejected at ingest.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidObstacleReason {
    /// Fewer than 3 distinct vertices.
    TooFewVertices,
    /// A vertex appears more than once in the ring.
    RepeatedVertex,
    /// Two consecutive vertices are bitwise equal.
    ZeroLengthEdge,
    /// Some coordinate is NaN or infinite.
    NaNOrInfinite,
    /// A pair of non-adjacent edges crosses.
    SelfIntersecting,
    /// The ring encloses zero area.
    DegenerateArea,
    /// The shoelace area overflowed, so the ring's orientation is not a number.
    ///
    /// The sum is a running total of cross terms, so a coordinate beyond roughly
    /// `sqrt(f64::MAX)` overflows it. The result is `NaN`, and a ring accepted with a
    /// `NaN` area has no usable orientation.
    AreaOverflow,
}

/// Everything that can go wrong during a route query.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum PathPlanError {
    /// The start or goal point is inside an obstacle, or exactly on its boundary.
    EndpointInObstacle {
        /// The offending point.
        point: Point2D,
    },
    /// The start or goal point is not inside the workspace bounding box.
    EndpointOutsideWorkspace {
        /// The offending point.
        point: Point2D,
    },
    /// Start and goal are the same point, so no curve with a well-defined tangent exists.
    DegenerateEndpoints,
    /// No feasible path exists between start and goal.
    NoPathFound,
    /// The workspace bounding box is empty or inverted.
    DegenerateWorkspace,
    /// A nonzero `margin` was requested with an obstacle that is not rectilinear.
    ///
    /// The `margin` model is exact for rectilinear geometry and is deliberately not approximated
    /// for anything else: realising a margin on a general polygon needs the obstacle's offset
    /// boundary trimmed against itself, and an untrimmed approximation produces geometry that
    /// passes within `margin` of the very obstacle it came from. A nonzero margin with a
    /// non-rectilinear obstacle is reported instead. See `docs/SPEC.md` section 3.3.
    MarginUnsupportedGeometry {
        /// The margin that was requested.
        margin: f64,
    },

    /// An obstacle polygon was malformed.
    InvalidObstacle {
        /// Why it was rejected.
        reason: InvalidObstacleReason,
    },
    /// A `Config` field was out of range.
    InvalidConfig(&'static str),
}

impl fmt::Display for PathPlanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EndpointInObstacle { point } => {
                write!(f, "point ({}, {}) is inside or on an obstacle", point.x, point.y)
            }
            Self::EndpointOutsideWorkspace { point } => {
                write!(f, "point ({}, {}) is outside the workspace", point.x, point.y)
            }
            Self::DegenerateEndpoints => write!(f, "start and goal are the same point"),
            Self::NoPathFound => write!(f, "no feasible path exists between start and goal"),
            Self::DegenerateWorkspace => write!(f, "workspace bounds are empty or inverted"),
            Self::MarginUnsupportedGeometry { margin } => write!(
                f,
                "margin {margin} requires rectilinear obstacles (axis-aligned edges only)"
            ),
            Self::InvalidObstacle { reason } => write!(f, "invalid obstacle: {reason:?}"),
            Self::InvalidConfig(why) => write!(f, "invalid configuration: {why}"),
        }
    }
}

impl core::error::Error for PathPlanError {}
