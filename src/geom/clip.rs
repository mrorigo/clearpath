//! Clipping a polygon against a rectangle.

use alloc::vec::Vec;

use super::point::Point2D;
use super::polygon::BoundingBox;

/// Clips `subject` to `clip` with Sutherland–Hodgman, against the rectangle's four half-planes.
///
/// The output ring is *weakly simple*: for a non-convex subject, Sutherland–Hodgman traverses
/// "bridge" edges twice. The enclosed region is still correct, which is all that is needed
/// because the output is only ever used as a forbidden region for containment tests, never as
/// geometry to be walked along.
///
/// Consecutive duplicate vertices are removed, so a caller can validate the result with
/// [`crate::Polygon::new`] without tripping its zero-length-edge check.
pub fn clip_to_rect(subject: &[Point2D], clip: BoundingBox) -> Vec<Point2D> {
    let mut ring = subject.to_vec();
    ring = clip_half_plane(ring, |p| p.x >= clip.min.x, |a, b| intersect_vertical(a, b, clip.min.x));
    if ring.is_empty() {
        return ring;
    }
    ring = clip_half_plane(ring, |p| p.x <= clip.max.x, |a, b| intersect_vertical(a, b, clip.max.x));
    if ring.is_empty() {
        return ring;
    }
    ring = clip_half_plane(ring, |p| p.y >= clip.min.y, |a, b| intersect_horizontal(a, b, clip.min.y));
    if ring.is_empty() {
        return ring;
    }
    ring = clip_half_plane(ring, |p| p.y <= clip.max.y, |a, b| intersect_horizontal(a, b, clip.max.y));
    drop_consecutive_duplicates(&ring)
}

/// Clips `ring` to the half-plane for which `inside` holds, inserting crossings found by
/// `intersect`.
fn clip_half_plane<F, G>(ring: Vec<Point2D>, inside: F, intersect: G) -> Vec<Point2D>
where
    F: Fn(Point2D) -> bool,
    G: Fn(Point2D, Point2D) -> Point2D,
{
    let mut out: Vec<Point2D> = Vec::with_capacity(ring.len() + 4);
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let a_in = inside(a);
        let b_in = inside(b);
        if a_in {
            out.push(a);
        }
        if a_in != b_in {
            out.push(intersect(a, b));
        }
    }
    out
}

fn intersect_vertical(a: Point2D, b: Point2D, x: f64) -> Point2D {
    if a.x == b.x {
        return Point2D::new(x, a.y);
    }
    let t = ((x - a.x) / (a.x - b.x)).clamp(0.0, 1.0);
    let y = (a.y + t * (b.y - a.y)).clamp(a.y.min(b.y), a.y.max(b.y));
    Point2D::new(x, y)
}

fn intersect_horizontal(a: Point2D, b: Point2D, y: f64) -> Point2D {
    if a.y == b.y {
        return Point2D::new(a.x, y);
    }
    let t = ((y - a.y) / (a.y - b.y)).clamp(0.0, 1.0);
    let x = (a.x + t * (b.x - a.x)).clamp(a.x.min(b.x), a.x.max(b.x));
    Point2D::new(x, y)
}

/// Removes runs of bitwise-equal consecutive vertices, including across the wrap-around.
fn drop_consecutive_duplicates(ring: &[Point2D]) -> Vec<Point2D> {
    let mut out: Vec<Point2D> = Vec::with_capacity(ring.len());
    for &p in ring {
        if out.last() != Some(&p) {
            out.push(p);
        }
    }
    while out.len() > 1 && out[0] == out[out.len() - 1] {
        out.pop();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn box_poly(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2D> {
        vec![
            Point2D::new(x0, y0),
            Point2D::new(x1, y0),
            Point2D::new(x1, y1),
            Point2D::new(x0, y1),
        ]
    }

    #[test]
    fn fully_inside_is_unchanged() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let out = clip_to_rect(&box_poly(1.0, 1.0, 2.0, 2.0), rect);
        assert_eq!(out.len(), 4);
    }

    #[test]
    fn fully_outside_is_empty() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        assert!(clip_to_rect(&box_poly(20.0, 20.0, 21.0, 21.0), rect).is_empty());
    }

    #[test]
    fn straddling_the_right_edge_is_clipped() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let out = clip_to_rect(&box_poly(5.0, 1.0, 15.0, 2.0), rect);
        assert!(out.iter().all(|p| p.x <= 10.0));
        assert!(out.iter().any(|p| p.x == 10.0));
    }

    #[test]
    fn touching_only_yields_nothing() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        assert!(clip_to_rect(&box_poly(10.0, 0.0, 20.0, 10.0), rect).len() < 3);
    }

    #[test]
    fn no_duplicate_consecutive_vertices() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        for (x0, y0, x1, y1) in [(0.0, 0.0, 5.0, 5.0), (3.0, 3.0, 30.0, 30.0), (-5.0, -5.0, 15.0, 15.0)] {
            let out = clip_to_rect(&box_poly(x0, y0, x1, y1), rect);
            for i in 0..out.len() {
                assert_ne!(out[i], out[(i + 1) % out.len()], "duplicate at {i}");
            }
        }
    }

    #[test]
    fn l_shape_stays_inside() {
        let rect = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let l = vec![
            Point2D::new(1.0, 1.0),
            Point2D::new(9.0, 1.0),
            Point2D::new(9.0, 3.0),
            Point2D::new(3.0, 3.0),
            Point2D::new(3.0, 9.0),
            Point2D::new(1.0, 9.0),
        ];
        let out = clip_to_rect(&l, rect);
        assert!(out.iter().all(|p| rect.contains(*p)));
        assert!(out.len() >= 3);
    }
}
