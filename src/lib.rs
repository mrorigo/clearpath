//! `pathplan` — deterministic 2D obstacle-avoidance path planning and smooth spline fitting.
//!
//! See `docs/SPEC.md` for the normative specification this crate implements.
//!
//! The crate is `#![no_std]` and needs only `alloc`. Enable the default `std` feature for the
//! ergonomic `std`-backed build; build with `--no-default-features` for bare-metal targets.
//!
//! # Example
//!
//! ```
//! use pathplan::{BoundingBox, Point2D, Polygon};
//!
//! // A 2x2 square obstacle.
//! let obstacle = Polygon::new(vec![
//!     Point2D::new(4.0, 4.0),
//!     Point2D::new(6.0, 4.0),
//!     Point2D::new(6.0, 6.0),
//!     Point2D::new(4.0, 6.0),
//! ]).unwrap();
//!
//! assert!(obstacle.contains(Point2D::new(5.0, 5.0)));
//! assert!(!obstacle.contains(Point2D::new(7.0, 5.0)));
//! assert_eq!(obstacle.bounds(), BoundingBox::new(Point2D::new(4.0, 4.0), Point2D::new(6.0, 6.0)));
//! ```

#![no_std]
#![deny(unsafe_code)]
#![deny(missing_docs)]

extern crate alloc;

#[cfg(test)]
extern crate std;

pub mod decomp;
pub mod error;
pub mod corridor;
pub mod funnel;
pub mod spline;
pub mod geom;

pub use error::{InvalidObstacleReason, PathPlanError};
pub use geom::{BoundingBox, Orientation, Point2D, Polygon};

/// `sqrt` without `std`.
///
/// `f64::sqrt` exists in `core`, but routing it through `libm` unconditionally keeps one numeric
/// path rather than two that could differ between builds. The angular functions the orthogonal
/// router will need (`sin`, `cos`, `atan2`) are added with it, for the same reason.
#[inline]
pub(crate) fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}

/// `ceil` without `std`.
#[inline]
pub(crate) fn ceil(x: f64) -> f64 {
    libm::ceil(x)
}
