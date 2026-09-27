//! Polygons and bounding boxes.

use alloc::vec::Vec;

use super::point::Point2D;
use super::predicates::{Orientation, orient2d};
use crate::error::{InvalidObstacleReason, PathPlanError};

/// An axis-aligned bounding box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingBox {
    /// Lower corner, inclusive.
    pub min: Point2D,
    /// Upper corner, inclusive.
    pub max: Point2D,
}

impl BoundingBox {
    /// Constructs a box from two corners, normalising them so `min <= max`.
    pub fn new(a: Point2D, b: Point2D) -> Self {
        Self { min: Point2D::new(a.x.min(b.x), a.y.min(b.y)), max: Point2D::new(a.x.max(b.x), a.y.max(b.y)) }
    }

    /// Width of the box.
    #[inline]
    pub fn width(&self) -> f64 {
        self.max.x - self.min.x
    }

    /// Height of the box.
    #[inline]
    pub fn height(&self) -> f64 {
        self.max.y - self.min.y
    }

    /// Area of the box.
    #[inline]
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    /// Whether both extents are strictly positive.
    #[inline]
    pub fn is_degenerate(&self) -> bool {
        !(self.width() > 0.0 && self.height() > 0.0)
    }

    /// Whether `p` is inside the closed box.
    #[inline]
    pub fn contains(&self, p: Point2D) -> bool {
        p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y
    }

    /// The four corners in counter-clockwise order, starting at `min`.
    pub fn corners(&self) -> [Point2D; 4] {
        [
            self.min,
            Point2D::new(self.max.x, self.min.y),
            self.max,
            Point2D::new(self.min.x, self.max.y),
        ]
    }

    /// The box shrunk inward by `d` on all four sides.
    pub fn eroded(&self, d: f64) -> Self {
        Self {
            min: Point2D::new(self.min.x + d, self.min.y + d),
            max: Point2D::new(self.max.x - d, self.max.y - d),
        }
    }
}

/// A closed, simple polygon, stored counter-clockwise.
///
/// "Simple" means the boundary does not cross itself and does not touch itself other than at
/// the shared endpoints of adjacent edges. Self-intersecting rings are rejected at construction
/// rather than being given an arbitrary winding-number interpretation, because a self-intersecting
/// ring has no well-defined inside.
///
/// A point that lies exactly on the boundary is reported as **inside** by
/// [`Polygon::contains`]. That is deliberate: it is the conservative answer for a path planner,
/// since it means an endpoint can never sit on an obstacle.
#[derive(Clone, Debug, PartialEq)]
pub struct Polygon {
    vertices: Vec<Point2D>,
}

impl Polygon {
    /// Validates a ring and stores it counter-clockwise.
    ///
    /// Rejects rings with non-finite coordinates, fewer than three distinct vertices, a
    /// zero-length edge, any repeated vertex, self-intersection, or zero enclosed area. A
    /// clockwise ring is accepted and reversed, so callers need not care about winding.
    pub fn new(mut vertices: Vec<Point2D>) -> Result<Self, PathPlanError> {
        if vertices.iter().any(|p| !p.is_finite()) {
            return Err(Self::err(InvalidObstacleReason::NaNOrInfinite));
        }
        for i in 0..vertices.len() {
            if vertices[i] == vertices[(i + 1) % vertices.len()] {
                return Err(Self::err(InvalidObstacleReason::ZeroLengthEdge));
            }
        }
        let distinct = deduped_len(&vertices);
        if distinct < 3 {
            return Err(Self::err(InvalidObstacleReason::TooFewVertices));
        }
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                if vertices[i] == vertices[j] {
                    return Err(Self::err(InvalidObstacleReason::RepeatedVertex));
                }
            }
        }
        if is_self_intersecting(&vertices) {
            return Err(Self::err(InvalidObstacleReason::SelfIntersecting));
        }
        if signed_area(&vertices) == 0.0 {
            return Err(Self::err(InvalidObstacleReason::DegenerateArea));
        }
        if signed_area(&vertices) < 0.0 {
            vertices.reverse();
        }
        Ok(Self { vertices })
    }

    /// Wraps an already-validated counter-clockwise ring without re-checking.
    ///
    /// # Panics
    /// Panics if the ring is degenerate, i.e. has fewer than three distinct vertices. For
    /// caller-supplied data prefer [`Polygon::new`], which reports rather than panics.
    pub fn new_unchecked(vertices: Vec<Point2D>) -> Self {
        assert!(deduped_len(&vertices) >= 3, "degenerate ring");
        Self { vertices }
    }

    /// The ring, counter-clockwise.
    #[inline]
    pub fn vertices(&self) -> &[Point2D] {
        &self.vertices
    }

    /// Number of vertices.
    #[inline]
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// Whether the ring has no vertices. Always `false` for a validated `Polygon`.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// The axis-aligned bounding box of the ring.
    pub fn bounds(&self) -> BoundingBox {
        let mut b = BoundingBox { min: Point2D::new(f64::INFINITY, f64::INFINITY), max: Point2D::new(f64::NEG_INFINITY, f64::NEG_INFINITY) };
        for p in &self.vertices {
            b.min.x = b.min.x.min(p.x);
            b.min.y = b.min.y.min(p.y);
            b.max.x = b.max.x.max(p.x);
            b.max.y = b.max.y.max(p.y);
        }
        b
    }

    /// Whether `p` is inside the polygon or exactly on its boundary.
    pub fn contains(&self, p: Point2D) -> bool {
        if !self.bounds().contains(p) {
            return false;
        }
        for i in 0..self.vertices.len() {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % self.vertices.len()];
            if on_segment(a, b, p) {
                return true;
            }
        }
        // Crossing count against the ray `p -> +x`. The half-open rule on `y` makes a ray
        // through a vertex count once. An edge that crosses the scanline contributes iff the
        // crossing is strictly right of `p`, which for an upward edge means `p` is left of
        // `a -> b` and for a downward edge means `p` is right of it.
        let mut inside = false;
        for i in 0..self.vertices.len() {
            let a = self.vertices[i];
            let b = self.vertices[(i + 1) % self.vertices.len()];
            let upward = a.y <= p.y && p.y < b.y;
            let downward = b.y <= p.y && p.y < a.y;
            if upward && is_left(a, b, p) || downward && is_right(a, b, p) {
                inside = !inside;
            }
        }
        inside
    }

    /// The area of the polygon; positive because the ring is counter-clockwise.
    pub fn area(&self) -> f64 {
        signed_area(&self.vertices)
    }

    fn err(reason: InvalidObstacleReason) -> PathPlanError {
        PathPlanError::InvalidObstacle { reason }
    }
}

#[inline]
fn is_left(a: Point2D, b: Point2D, c: Point2D) -> bool {
    orient2d(a, b, c) == Orientation::CounterClockwise
}

#[inline]
fn is_right(a: Point2D, b: Point2D, c: Point2D) -> bool {
    orient2d(a, b, c) == Orientation::Clockwise
}

/// Whether `p` lies on the closed segment `a`..`b`.
pub fn on_segment(a: Point2D, b: Point2D, p: Point2D) -> bool {
    if orient2d(a, b, p) != Orientation::Collinear {
        return false;
    }
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

/// The signed area of a ring; positive for counter-clockwise.
fn signed_area(vertices: &[Point2D]) -> f64 {
    let mut acc = 0.0;
    for i in 0..vertices.len() {
        let a = vertices[i];
        let b = vertices[(i + 1) % vertices.len()];
        acc += a.cross(b);
    }
    acc / 2.0
}

fn deduped_len(vertices: &[Point2D]) -> usize {
    let mut n = 0usize;
    for i in 0..vertices.len() {
        let mut seen = false;
        for j in 0..i {
            if vertices[i] == vertices[j] {
                seen = true;
                break;
            }
        }
        if !seen {
            n += 1;
        }
    }
    n
}

/// Whether any two non-adjacent edges of the ring cross, or any non-adjacent pair of edges shares
/// a point.
///
/// Quadratic in the ring length. Obstacle rings are small (a box is 4), so this is not a hot
/// path; if a caller ever passes thousand-vertex rings this becomes the ingest bottleneck and
/// should be replaced with a sweep.
fn is_self_intersecting(vertices: &[Point2D]) -> bool {
    let n = vertices.len();
    if n < 4 {
        return false;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let adjacent = j == i + 1 || (i == 0 && j == n - 1);
            let a0 = vertices[i];
            let a1 = vertices[(i + 1) % n];
            let b0 = vertices[j];
            let b1 = vertices[(j + 1) % n];
            if adjacent {
                // Adjacent edges legitimately share exactly one endpoint; anything else is a
                // self-touch. The shared endpoint is `a1 == b0` for consecutive edges and
                // `a0 == b1` for the wrap-around pair.
                let shared = if j == i + 1 { a1 == b0 } else { a0 == b1 };
                if !shared {
                    return true;
                }
            } else if segments_intersect(a0, a1, b0, b1) {
                return true;
            }
        }
    }
    false
}

/// Whether the closed segments `a0a1` and `b0b1` share at least one point, including touching.
pub fn segments_intersect(a0: Point2D, a1: Point2D, b0: Point2D, b1: Point2D) -> bool {
    let d1 = orient2d(a0, a1, b0);
    let d2 = orient2d(a0, a1, b1);
    let d3 = orient2d(b0, b1, a0);
    let d4 = orient2d(b0, b1, a1);

    if d1 != d2 && d3 != d4 {
        return true;
    }
    // At least one endpoint lies on the other segment (covers all collinear and touching cases).
    (d1.is_collinear() && on_segment(a0, a1, b0))
        || (d2.is_collinear() && on_segment(a0, a1, b1))
        || (d3.is_collinear() && on_segment(b0, b1, a0))
        || (d4.is_collinear() && on_segment(b0, b1, a1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn square() -> Vec<Point2D> {
        vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(0.0, 2.0),
        ]
    }

    #[test]
    fn accepts_a_square_and_is_ccw() {
        let poly = Polygon::new(square()).unwrap();
        assert_eq!(poly.area(), 4.0);
        assert_eq!(poly.len(), 4);
    }

    #[test]
    fn reverses_clockwise_input() {
        let mut v = square();
        v.reverse();
        let poly = Polygon::new(v).unwrap();
        assert_eq!(poly.area(), 4.0);
        assert_eq!(poly.vertices()[0], Point2D::new(0.0, 0.0));
    }

    #[test]
    fn rejects_non_finite() {
        let mut v = square();
        v[1].x = f64::NAN;
        assert_eq!(
            Polygon::new(v).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::NaNOrInfinite }
        );
    }

    #[test]
    fn rejects_zero_length_edge() {
        let v = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(0.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(0.0, 2.0),
        ];
        assert_eq!(
            Polygon::new(v).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::ZeroLengthEdge }
        );
    }

    #[test]
    fn rejects_too_few_vertices() {
        let v = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0), Point2D::new(0.0, 1.0), Point2D::new(0.0, 0.0)];
        assert_eq!(
            Polygon::new(v).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::ZeroLengthEdge }
        );
        assert_eq!(
            Polygon::new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::TooFewVertices }
        );
    }

    #[test]
    fn rejects_repeated_vertex() {
        // (2, 0) appears at indices 1 and 4, which are not adjacent.
        let v = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(0.0, 2.0),
            Point2D::new(2.0, 0.0),
        ];
        assert_eq!(
            Polygon::new(v).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::RepeatedVertex }
        );
    }

    #[test]
    fn rejects_self_intersection() {
        // A bowtie.
        let v = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(0.0, 2.0),
        ];
        assert_eq!(
            Polygon::new(v).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::SelfIntersecting }
        );
    }

    #[test]
    fn rejects_zero_area() {
        let v = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(0.0, 1.0),
        ];
        // area = 1.0, not zero: use an exactly-collinear ring instead.
        let collinear = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
        ];
        assert!(Polygon::new(v).is_ok());
        assert_eq!(
            Polygon::new(collinear).unwrap_err(),
            PathPlanError::InvalidObstacle { reason: InvalidObstacleReason::DegenerateArea }
        );
    }

    #[test]
    fn contains_interior_exterior_and_boundary() {
        let poly = Polygon::new(square()).unwrap();
        assert!(poly.contains(Point2D::new(1.0, 1.0)));
        assert!(!poly.contains(Point2D::new(3.0, 1.0)));
        // Boundary points count as inside, so an endpoint can never sit on an obstacle.
        assert!(poly.contains(Point2D::new(0.0, 1.0)));
        assert!(poly.contains(Point2D::new(1.0, 0.0)));
        assert!(poly.contains(Point2D::new(2.0, 2.0)));
    }

    #[test]
    fn bounds_and_area() {
        let poly = Polygon::new(square()).unwrap();
        assert_eq!(poly.bounds(), BoundingBox::new(Point2D::ZERO, Point2D::new(2.0, 2.0)));
        assert!((poly.area() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn segment_intersection() {
        let a0 = Point2D::new(0.0, 0.0);
        let a1 = Point2D::new(4.0, 0.0);
        assert!(segments_intersect(a0, a1, Point2D::new(1.0, -1.0), Point2D::new(1.0, 1.0)));
        assert!(!segments_intersect(a0, a1, Point2D::new(0.0, 1.0), Point2D::new(4.0, 1.0)));
        // touching counts
        assert!(segments_intersect(a0, a1, Point2D::new(2.0, 0.0), Point2D::new(2.0, 1.0)));
    }

    #[test]
    fn bounding_box_helpers() {
        let b = BoundingBox::new(Point2D::new(2.0, 4.0), Point2D::new(0.0, 0.0));
        assert_eq!(b.min, Point2D::ZERO);
        assert_eq!(b.max, Point2D::new(2.0, 4.0));
        assert_eq!(b.width(), 2.0);
        assert_eq!(b.height(), 4.0);
        assert_eq!(b.area(), 8.0);
        assert!(!b.is_degenerate());
        assert!(b.contains(Point2D::new(1.0, 1.0)));
        assert!(!b.contains(Point2D::new(3.0, 1.0)));
        assert!(b.eroded(0.5) == BoundingBox::new(Point2D::new(0.5, 0.5), Point2D::new(1.5, 3.5)));
        assert_eq!(b.corners()[0], Point2D::ZERO);
        assert_eq!(b.corners()[2], b.max);
    }

    #[test]
    fn degenerate_box_detected() {
        let b = BoundingBox { min: Point2D::new(1.0, 1.0), max: Point2D::new(1.0, 5.0) };
        assert!(b.is_degenerate());
    }
}
