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

pub mod error;
pub mod geom;

pub use error::{InvalidObstacleReason, PathPlanError};
pub use geom::{BoundingBox, Orientation, Point2D, Polygon};

/// Transcendentals go through `libm` unconditionally rather than behind `cfg(feature = "std")`.
/// `f64::sqrt`, `sin` and `cos` exist in `core`, but `acos`, `atan2` and `ceil` do not, and a
/// single code path is worth more here than the marginal dependency cost.
macro_rules! libm_fn {
    ($name:ident, $path:path, $($arg:ident),*) => {
        #[doc = concat!("`", stringify!($name), "`, via `libm`.")]
        #[inline]
        pub(crate) fn $name($($arg: f64),*) -> f64 {
            $path($($arg),*)
        }
    };
}

libm_fn!(sqrt, libm::sqrt, x);
libm_fn!(sin, libm::sin, x);
libm_fn!(cos, libm::cos, x);
libm_fn!(acos, libm::acos, x);
libm_fn!(atan2, libm::atan2, y, x);
libm_fn!(ceil, libm::ceil, x);
