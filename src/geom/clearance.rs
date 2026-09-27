//! The clearance predicate: is a point or a polygon in the free space of the inflated domain?
//!
//! This is where the section 4.1 guarantee is actually decided, so it is written as an exact test
//! against the obstacles rather than as a sample of a region. An earlier version tested only the
// midpoints of a control hull's four edges, and a hull that clipped an obstacle's corner passed —
//! which is a collision, and the kind that only shows up in the field.

use alloc::vec::Vec;

use super::point::Point2D;
use super::polygon::{BoundingBox, Polygon, segments_intersect};
use super::predicates::{Orientation, orient2d};

/// The free space of the inflated domain, as a predicate.
///
/// `margin` is the same value passed to the decomposition, so this tests the domain the planner
/// actually routed in rather than a separately-maintained copy of the rules.
#[derive(Clone, Debug)]
pub struct Clearance {
    obstacles: Vec<Polygon>,
    /// The eroded workspace: a point must be inside this.
    workspace: BoundingBox,
}

impl Clearance {
    /// Builds the predicate for the domain the planner used.
    pub fn new(obstacles: &[Polygon], workspace: BoundingBox, margin: f64) -> Self {
        Self { obstacles: obstacles.to_vec(), workspace: workspace.eroded(margin) }
    }

    /// The eroded workspace.
    pub fn workspace(&self) -> BoundingBox {
        self.workspace
    }

    /// The obstacles, as the planner saw them.
    pub fn obstacles(&self) -> &[Polygon] {
        &self.obstacles
    }

    /// Whether `p` is in the free space.
    ///
    /// A distance test, not a containment test: at `margin == 0` a taut path's knots lie exactly on
    /// obstacle boundaries, and the guarantee is `>= margin`, so the boundary is in.
    pub fn is_free(&self, p: Point2D) -> bool {
        self.workspace.contains(p)
            && self.obstacles.iter().all(|o| {
                let ring = o.vertices();
                let d = distance_to_ring(p, ring);
                // `Polygon::contains` reports the boundary as inside, and at `margin == 0` a taut
                // path's knots lie exactly on boundaries. So the interior test has to be
                // "strictly inside", which is `contains` *and* a positive distance — a point on the
                // boundary has distance zero and an interior point does not.
                let strictly_inside = o.contains(p) && d > 0.0;
                !strictly_inside && d >= 0.0
            })
    }

    /// Whether the convex hull of `points` is in the free space.
    ///
    /// Exact: the hull's edges are tested against every obstacle edge, and every obstacle vertex is
    /// tested for being inside the hull. Sampling the hull instead is what let a corner clip
    /// through.
    pub fn hull_is_free(&self, points: &[Point2D]) -> bool {
        let hull = convex_hull(points);
        if hull.len() < 3 {
            // A degenerate hull is its points; the caller has already checked them.
            return hull.iter().all(|p| self.is_free(*p));
        }
        for p in &hull {
            if !self.is_free(*p) {
                return false;
            }
        }
        for obstacle in &self.obstacles {
            let ring = obstacle.vertices();
            for v in ring {
                if point_in_convex(&hull, *v) {
                    return false;
                }
            }
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                for k in 0..hull.len() {
                    let (c, d) = (hull[k], hull[(k + 1) % hull.len()]);
                    if segments_intersect(a, b, c, d) {
                        return false;
                    }
                }
            }
        }
        true
    }
}

/// Minimum distance from `p` to a ring's boundary.
pub fn distance_to_ring(p: Point2D, ring: &[Point2D]) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let ab = b - a;
        let len2 = ab.norm_squared();
        let t = if len2 == 0.0 { 0.0 } else { ((p - a).dot(ab) / len2).clamp(0.0, 1.0) };
        best = best.min(p.distance(a + ab * t));
    }
    best
}

/// The convex hull of a small point set, counter-clockwise, without repeating the first point.
///
/// Andrew's monotone chain. Used on a handful of control points, so the sort is by orientation
/// against the chosen extreme and ties fall back to the coordinates, which keeps it a
/// deterministic function of the input.
pub fn convex_hull(points: &[Point2D]) -> Vec<Point2D> {
    let mut pts: Vec<Point2D> = points.to_vec();
    pts.sort_by(|a, b| {
        a.x.total_cmp(&b.x).then_with(|| a.y.total_cmp(&b.y))
    });
    pts.dedup();
    if pts.len() < 3 {
        return pts;
    }
    let lower = half(&pts);
    let upper = half(&pts.iter().rev().copied().collect::<Vec<_>>());
    // `lower` runs from the first point to the last, `upper` from the last back to the first, so
    // each contributes both duplicated endpoints. Dropping the last of each closes the ring.
    let mut out: Vec<Point2D> = lower[..lower.len() - 1].to_vec();
    out.extend_from_slice(&upper[..upper.len() - 1]);
    out
}

fn half(sorted: &[Point2D]) -> Vec<Point2D> {
    let mut out: Vec<Point2D> = Vec::with_capacity(sorted.len());
    for &p in sorted {
        while out.len() >= 2
            && orient2d(out[out.len() - 2], out[out.len() - 1], p) != Orientation::CounterClockwise
        {
            out.pop();
        }
        out.push(p);
    }
    out
}

/// Whether `p` is inside the convex polygon `poly`, given that it may lie on the boundary.
pub fn point_in_convex(poly: &[Point2D], p: Point2D) -> bool {
    if poly.len() < 3 {
        return poly.contains(&p);
    }
    for i in 0..poly.len() {
        if orient2d(poly[i], poly[(i + 1) % poly.len()], p) == Orientation::Clockwise {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn square() -> Polygon {
        Polygon::new(vec![
            Point2D::new(40.0, 40.0),
            Point2D::new(60.0, 40.0),
            Point2D::new(60.0, 60.0),
            Point2D::new(40.0, 60.0),
        ])
        .unwrap()
    }

    fn domain(margin: f64) -> Clearance {
        Clearance::new(
            &[square()],
            BoundingBox::new(Point2D::ZERO, Point2D::new(100.0, 100.0)),
            margin,
        )
    }

    #[test]
    fn free_and_blocked_points() {
        let c = domain(0.0);
        assert!(c.is_free(Point2D::new(10.0, 10.0)));
        assert!(!c.is_free(Point2D::new(50.0, 50.0)));
        // A knot on the boundary is in the free space at margin 0.
        assert!(c.is_free(Point2D::new(40.0, 50.0)));
    }

    #[test]
    fn the_margin_moves_the_workspace_bound() {
        let c = domain(10.0);
        assert!(!c.is_free(Point2D::new(5.0, 50.0)));
        assert!(c.is_free(Point2D::new(15.0, 50.0)));
    }

    #[test]
    fn a_hull_that_clips_a_corner_is_rejected() {
        // The case a midpoint sample misses: a triangle whose only contact with the square is a
        // sliver of one edge's interior, with both endpoints clear.
        let c = domain(0.0);
        let clipping = vec![
            Point2D::new(38.0, 20.0),
            Point2D::new(38.0, 80.0),
            Point2D::new(70.0, 50.0),
        ];
        assert!(!c.hull_is_free(&clipping), "a hull crossing the square must be rejected");
    }

    #[test]
    fn a_clear_hull_is_accepted() {
        let c = domain(0.0);
        let clear = vec![
            Point2D::new(5.0, 20.0),
            Point2D::new(20.0, 20.0),
            Point2D::new(80.0, 20.0),
            Point2D::new(95.0, 20.0),
        ];
        assert!(c.hull_is_free(&clear));
    }

    #[test]
    fn a_hull_containing_an_obstacle_vertex_is_rejected() {
        // The whole obstacle inside the hull, so no edge crosses it.
        let c = domain(0.0);
        let surrounding = vec![
            Point2D::new(10.0, 10.0),
            Point2D::new(90.0, 10.0),
            Point2D::new(90.0, 90.0),
            Point2D::new(10.0, 90.0),
        ];
        assert!(!c.hull_is_free(&surrounding));
    }

    #[test]
    fn convex_hull_of_a_bowtie_is_the_hull() {
        let hull = convex_hull(&[
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
        ]);
        assert_eq!(hull.len(), 4);
        for w in hull.windows(2) {
            assert_eq!(orient2d(w[0], w[1], Point2D::new(0.5, 0.5)), Orientation::CounterClockwise);
        }
    }

    #[test]
    fn convex_hull_of_collinear_points_is_a_segment() {
        let hull = convex_hull(&[Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0), Point2D::new(2.0, 2.0)]);
        assert!(hull.len() <= 2);
    }

    #[test]
    fn distance_to_ring() {
        let ring: Vec<Point2D> = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(10.0, 0.0),
            Point2D::new(10.0, 10.0),
            Point2D::new(0.0, 10.0),
        ];
        assert_eq!(super::distance_to_ring(Point2D::new(5.0, 5.0), &ring), 5.0);
        assert_eq!(super::distance_to_ring(Point2D::new(-1.0, 5.0), &ring), 1.0);
    }
}
