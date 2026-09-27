//! Fixtures shared by the integration tests.
//!
//! These were five identical copies of a box constructor, and two copies of a distance function —
//! one of which re-implemented [`clearpath::geom::clearance::distance_to_ring`], a public function
//! the library already exports.

#![allow(dead_code)] // each test binary uses a different subset

use clearpath::geom::Point2D;
use clearpath::{BoundingBox, Polygon};

/// An axis-aligned box from `(min_x, min_y, max_x, max_y)`.
pub fn box_poly(c: [f64; 4]) -> Polygon {
    Polygon::new(vec![
        Point2D::new(c[0], c[1]),
        Point2D::new(c[2], c[1]),
        Point2D::new(c[2], c[3]),
        Point2D::new(c[0], c[3]),
    ])
    .unwrap()
}

/// Shorthand for a point.
pub fn p(x: f64, y: f64) -> Point2D {
    Point2D::new(x, y)
}

/// The 100x100 workspace most fixtures use.
pub const WORKSPACE: BoundingBox = BoundingBox {
    min: Point2D { x: 0.0, y: 0.0 },
    max: Point2D { x: 100.0, y: 100.0 },
};
