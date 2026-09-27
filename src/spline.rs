//! Cubic Bezier segments and the spline solver.

use crate::geom::point::Point2D;
use crate::geom::BoundingBox;

/// One cubic Bezier segment, with its four control points in counter-clockwise-free order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicBezierSegment {
    /// Start of the curve.
    pub p0: Point2D,
    /// First control point.
    pub p1: Point2D,
    /// Second control point.
    pub p2: Point2D,
    /// End of the curve.
    pub p3: Point2D,
}

impl CubicBezierSegment {
    /// The point at parameter `t`.
    pub fn evaluate(&self, t: f64) -> Point2D {
        let u = 1.0 - t;
        let (a, b) = (u * u * u, 3.0 * u * u * t);
        let (c, d) = (3.0 * u * t * t, t * t * t);
        self.p0 * a + self.p1 * b + self.p2 * c + self.p3 * d
    }

    /// The derivative at parameter `t`.
    pub fn tangent(&self, t: f64) -> Point2D {
        let u = 1.0 - t;
        self.p0 * (-3.0 * u * u)
            + self.p1 * (6.0 * u * t - 3.0 * u * u)
            + self.p2 * (3.0 * u * u - 6.0 * u * t)
            + self.p3 * (3.0 * t * t)
    }

    /// The axis-aligned bounding box of the control hull, which contains the curve.
    pub fn bounding_box(&self) -> BoundingBox {
        BoundingBox::new(self.p0, self.p3)
            .with(self.p1)
            .with(self.p2)
    }
}
