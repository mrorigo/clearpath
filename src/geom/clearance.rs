//! The clearance predicate: is a point or a polygon in the free space of the inflated domain?
//!
//! This is where the section 4.1 guarantee is actually decided, so it is written as an exact test
//! against the obstacles rather than as a sample of a region. An earlier version tested only the
// midpoints of a control hull's four edges, and a hull that clipped an obstacle's corner passed —
//! which is a collision, and the kind that only shows up in the field.

use alloc::vec::Vec;

use super::point::Point2D;
use super::polygon::{BoundingBox, Polygon, segments_cross_properly};
use super::predicates::{Orientation, orient2d};

/// The free space of the inflated domain, as a predicate.
///
/// `margin` is the same value passed to the decomposition, so this tests the domain the planner
/// actually routed in rather than a separately-maintained copy of the rules.
#[derive(Clone, Copy, Debug)]
pub struct Clearance<'a> {
    /// Borrowed, not owned. Cloning the obstacles here cost one allocation per obstacle on every
    /// query — two hundred of them on a two-hundred-box route, for a value the caller already
    /// holds and never mutates.
    obstacles: &'a [Polygon],
    /// The eroded workspace: a point must be inside this.
    workspace: BoundingBox,
    /// The margin, which the obstacles are held off by.
    margin: f64,
}

impl<'a> Clearance<'a> {
    /// Builds the predicate for the domain the planner used.
    pub fn new(obstacles: &'a [Polygon], workspace: BoundingBox, margin: f64) -> Self {
        Self { obstacles, workspace: workspace.eroded(margin), margin }
    }

    /// The eroded workspace.
    pub fn workspace(&self) -> BoundingBox {
        self.workspace
    }

    /// The obstacles, as the planner saw them.
    pub fn obstacles(&self) -> &[Polygon] {
        self.obstacles
    }

    /// The tolerance used when comparing against a bounding box, in distance units.
    fn slack(&self) -> f64 {
        1e-9 * self.margin.abs().max(1.0)
    }

    /// Whether `p` is at least `margin` away (less `slack`) from the segment `c..d`.
    ///
    /// Rejects on the segment's own box first, so the common case — a vertex nowhere near this edge
    /// of the hull — costs four comparisons rather than a projection and a square root.
    fn point_near_segment(&self, p: Point2D, c: Point2D, d: Point2D, slack: f64) -> bool {
        let (dx, dy) = (d.x - c.x, d.y - c.y);
        let lo_x = c.x.min(d.x) - self.margin;
        let hi_x = c.x.max(d.x) + self.margin;
        let lo_y = c.y.min(d.y) - self.margin;
        let hi_y = c.y.max(d.y) + self.margin;
        if p.x < lo_x || p.x > hi_x || p.y < lo_y || p.y > hi_y {
            return true;
        }
        let len2 = dx * dx + dy * dy;
        let q = if len2 == 0.0 {
            c
        } else {
            let t = ((p.x - c.x) * dx + (p.y - c.y) * dy) / len2;
            let t = t.clamp(0.0, 1.0);
            Point2D::new(c.x + t * dx, c.y + t * dy)
        };
        q.distance(p) >= self.margin - slack
    }

    /// Whether `p` is in the free space.
    ///
    /// A distance test, not a containment test: at `margin == 0` a taut path's knots lie exactly on
    /// obstacle boundaries, and the guarantee is `>= margin`, so the boundary is in.
    pub fn is_free(&self, p: Point2D) -> bool {
        // A tolerance on the workspace bound, for the same reason the obstacle test has one: a
        // curve sampled at a knot that lies exactly on the eroded boundary can evaluate a fraction
        // of an ulp outside it, and reporting that as a collision is a rounding artefact dressed up
        // as a geometric fact.
        let slack = 1e-9 * self.workspace.width().abs().max(self.workspace.height().abs()).max(1.0);
        let ws = self.workspace;
        if p.x < ws.min.x - slack
            || p.x > ws.max.x + slack
            || p.y < ws.min.y - slack
            || p.y > ws.max.y + slack
        {
            return false;
        }
        // The margin is the *distance* an obstacle is held off by, so it is a distance test and not
        // an outset geometry: the obstacles here are the caller's originals, and eroding the
        // workspace alone would leave them un-grown. A small relative slack absorbs the rounding in
        // `distance_to_ring` at exactly `margin`.
        let slack = 1e-9 * self.margin.abs().max(1.0);
        for o in self.obstacles {
            // A point further than `margin` from an obstacle's bounding box is further than
            // `margin` from the obstacle, so the per-edge loop below can be skipped. On a route
            // with many boxes most obstacles are nowhere near any given hull, and this is four
            // comparisons against a per-edge loop of a square root each.
            let b = o.bounds();
            if p.x < b.min.x - self.margin - slack
                || p.x > b.max.x + self.margin + slack
                || p.y < b.min.y - self.margin - slack
                || p.y > b.max.y + self.margin + slack
            {
                continue;
            }
            let d = distance_to_ring(p, o.vertices());
            // `Polygon::contains` reports the boundary as inside, and at `margin == 0` a taut path's
            // knots lie exactly on boundaries. So the interior test has to be "strictly inside",
            // which is `contains` *and* a positive distance.
            if (o.contains(p) && d > 0.0) || d < self.margin - slack {
                return false;
            }
        }
        true
    }

    /// Whether the convex hull of `points` is in the free space.
    ///
    /// Exact: the hull's edges are tested against every obstacle edge, and every obstacle vertex is
    /// tested for being inside the hull. Sampling the hull instead is what let a corner clip
    /// through.
    pub fn hull_is_free(&self, points: &[Point2D]) -> bool {
        let mut buffer = [Point2D::ZERO; 8];
        let count = hull_into(points, &mut buffer);
        let hull = &buffer[..count];
        for p in hull {
            if !self.is_free(*p) {
                return false;
            }
        }
        if hull.len() < 2 {
            return true;
        }
        // The hull's own box, to skip whole obstacles.
        let mut hbox = BoundingBox::new(hull[0], hull[0]);
        for p in hull {
            hbox = hbox.with(*p);
        }
        for obstacle in self.obstacles {
            let b = obstacle.bounds();
            if hbox.max.x < b.min.x - self.margin - self.slack()
                || hbox.min.x > b.max.x + self.margin + self.slack()
                || hbox.max.y < b.min.y - self.margin - self.slack()
                || hbox.min.y > b.max.y + self.margin + self.slack()
            {
                continue;
            }
            let ring = obstacle.vertices();
            for k in 0..hull.len() {
                // The hull's closing polyline, which for a two-point hull is the segment itself.
                // That case matters: a zero-tangent segment is exactly two points, and checking only
                // its endpoints would accept a segment that passes through an obstacle.
                let (c, d) = (hull[k], hull[(k + 1) % hull.len()]);
                // The local rejection. The hull's box passed this obstacle, which says very
                // little: the hull is long and thin, so its box is most of the workspace. This says
                // whether the edge *itself* can come within `margin` of the obstacle's box — grown by
                // the margin, because the check below is about the margin and not about contact. An
                // un-grown box here skips exactly the edges that graze: a curve can be clear of the
                // obstacle and still inside its margin, and that is the case the margin exists for.
                let grown = BoundingBox {
                    min: Point2D::new(b.min.x - self.margin, b.min.y - self.margin),
                    max: Point2D::new(b.max.x + self.margin, b.max.y + self.margin),
                };
                if !segment_meets_box(c, d, grown) {
                    continue;
                }
                for i in 0..ring.len() {
                    if segments_cross_properly(ring[i], ring[(i + 1) % ring.len()], c, d) {
                        return false;
                    }
                }
                // The margin is a *metric* property and the crossing test is topological, so the
                // edges need their own clearance check. Without it a hull can have four clear
                // vertices and an edge that passes within `margin` of an obstacle, and the curve
                // between those vertices is then too close — which is the whole thing the check
                // exists to prevent.
                //
                // For a segment that does not cross an obstacle edge, the closest approach to that
                // edge is attained at one of its endpoints, so testing every obstacle *vertex*
                // against the segment is sufficient.
                let slack = 1e-9 * self.margin.abs().max(1.0);
                for v in ring {
                    if !self.point_near_segment(*v, c, d, slack) {
                        return false;
                    }
                }
            }
            // A vertex *strictly* inside the hull means the obstacle is enclosed, or partly
            // covered. On the boundary is contact, which the guarantee permits. Meaningless for a
            // degenerate hull, so it is skipped there.
            if hull.len() >= 3 {
                for v in ring {
                    if point_strictly_in_convex(hull, *v) {
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

/// The convex hull of up to four points, counter-clockwise, without repeating the first point.
///
/// Every caller passes a Bezier segment's four control points, and this was one allocation and one
/// sort per hull test — and the repair runs a hull test per segment per iteration. Four points fit
/// in a fixed array, so the sort is an insertion sort over that and nothing is allocated.
pub fn convex_hull(points: &[Point2D]) -> Vec<Point2D> {
    let mut buffer = [Point2D::ZERO; 8];
    let count = hull_into(points, &mut buffer);
    buffer[..count].to_vec()
}

/// The convex hull of up to four points, into a caller-supplied buffer, returning the count.
///
/// `hull_is_free` runs once per segment per repair iteration, and it allocated three `Vec`s each
/// time — two chains and the output. Four input points can produce at most four hull vertices, so
/// a fixed buffer removes the allocations entirely and the repair stops being allocation-bound.
pub fn hull_into(points: &[Point2D], out: &mut [Point2D]) -> usize {
    debug_assert!(out.len() >= 8);
    let mut pts: [Point2D; 4] = [Point2D::ZERO; 4];
    let n = points.len().min(4);
    pts[..n].copy_from_slice(&points[..n]);
    for i in 1..n {
        let key = pts[i];
        let mut j = i;
        while j > 0 && less(&pts[j - 1], &key) {
            pts[j] = pts[j - 1];
            j -= 1;
        }
        pts[j] = key;
    }
    // Deduplicate: the control points of a zero-tangent segment repeat.
    let mut len = 0usize;
    for i in 0..n {
        let p = pts[i];
        if len == 0 || pts[len - 1] != p {
            pts[len] = p;
            len += 1;
        }
    }
    let sorted = &pts[..len];
    if len < 3 {
        for (i, p) in sorted.iter().enumerate() {
            out[i] = *p;
        }
        return len;
    }
    // Two monotone chains over a fixed buffer: `lower` from the first point to the last, `upper`
    // back again, so each carries both duplicated endpoints and the ring closes by dropping the
    // last of each.
    let mut lower: [Point2D; 4] = [Point2D::ZERO; 4];
    let mut upper: [Point2D; 4] = [Point2D::ZERO; 4];
    let (mut ln, mut un) = (0usize, 0usize);
    for p in sorted.iter() {
        while ln >= 2 && orient2d(lower[ln - 2], lower[ln - 1], *p) != Orientation::CounterClockwise {
            ln -= 1;
        }
        lower[ln] = *p;
        ln += 1;
    }
    for p in sorted.iter().rev() {
        while un >= 2 && orient2d(upper[un - 2], upper[un - 1], *p) != Orientation::CounterClockwise {
            un -= 1;
        }
        upper[un] = *p;
        un += 1;
    }
    let mut count = 0usize;
    for p in lower.iter().take(ln.saturating_sub(1)) {
        out[count] = *p;
        count += 1;
    }
    for p in upper.iter().take(un.saturating_sub(1)) {
        out[count] = *p;
        count += 1;
    }
    count
}

/// Lexicographic order on the bit patterns, so the sort is a total order and the hull is a
/// deterministic function of the input.
fn less(a: &Point2D, b: &Point2D) -> bool {
    a.x.total_cmp(&b.x) == core::cmp::Ordering::Less
        || (a.x == b.x && a.y.total_cmp(&b.y) == core::cmp::Ordering::Less)
}

/// Whether the closed segment `a..b` meets the axis-aligned box, by the slab method.
///
/// This is the cheap local test that the hull's own bounding box cannot be. A control hull on a
/// long route is long and thin, so its box covers most of the workspace and every obstacle passes
/// it — but the hull's individual *edges* are mostly short, and this rejects the obstacles each one
/// cannot reach before any segment-versus-segment work happens. Measured on a real route, that was
/// the difference between six thousand segment tests per query and a few hundred.
pub fn segment_meets_box(a: Point2D, b: Point2D, extent: BoundingBox) -> bool {
    let mut t0 = 0.0f64;
    let mut t1 = 1.0f64;
    for (p, d, lo, hi) in [
        (a.x, b.x - a.x, extent.min.x, extent.max.x),
        (a.y, b.y - a.y, extent.min.y, extent.max.y),
    ] {
        if d == 0.0 {
            // Parallel to this slab: outside it means no intersection at all.
            if p < lo || p > hi {
                return false;
            }
        } else {
            let (mut near, mut far) = ((lo - p) / d, (hi - p) / d);
            if near > far {
                core::mem::swap(&mut near, &mut far);
            }
            t0 = t0.max(near);
            t1 = t1.min(far);
            if t0 > t1 {
                return false;
            }
        }
    }
    true
}

/// Whether `p` is inside the convex polygon `poly`, given that it may lie on the boundary.
pub fn point_in_convex(poly: &[Point2D], p: Point2D) -> bool {
    if poly.len() < 3 {
        return poly.contains(&p);
    }
    (0..poly.len())
        .all(|i| orient2d(poly[i], poly[(i + 1) % poly.len()], p) != Orientation::Clockwise)
}

/// Whether `p` is strictly inside the convex polygon `poly`; the boundary does not count.
pub fn point_strictly_in_convex(poly: &[Point2D], p: Point2D) -> bool {
    poly.len() >= 3
        && (0..poly.len())
            .all(|i| orient2d(poly[i], poly[(i + 1) % poly.len()], p) == Orientation::CounterClockwise)
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

    fn domain(obstacles: &[Polygon], margin: f64) -> Clearance<'_> {
        Clearance::new(obstacles, BoundingBox::new(Point2D::ZERO, Point2D::new(100.0, 100.0)), margin)
    }

    #[test]
    fn free_and_blocked_points() {
        let obstacles = [square()];
        let c = domain(&obstacles, 0.0);
        assert!(c.is_free(Point2D::new(10.0, 10.0)));
        assert!(!c.is_free(Point2D::new(50.0, 50.0)));
        // A knot on the boundary is in the free space at margin 0.
        assert!(c.is_free(Point2D::new(40.0, 50.0)));
    }

    #[test]
    fn the_margin_moves_the_workspace_bound() {
        let obstacles = [square()];
        let c = domain(&obstacles, 10.0);
        assert!(!c.is_free(Point2D::new(5.0, 50.0)));
        assert!(c.is_free(Point2D::new(15.0, 50.0)));
    }

    #[test]
    fn a_hull_that_clips_a_corner_is_rejected() {
        // The case a midpoint sample misses: a triangle whose only contact with the square is a
        // sliver of one edge's interior, with both endpoints clear.
        let obstacles = [square()];
        let c = domain(&obstacles, 0.0);
        let clipping = vec![
            Point2D::new(38.0, 20.0),
            Point2D::new(38.0, 80.0),
            Point2D::new(70.0, 50.0),
        ];
        assert!(!c.hull_is_free(&clipping), "a hull crossing the square must be rejected");
    }

    #[test]
    fn a_clear_hull_is_accepted() {
        let obstacles = [square()];
        let c = domain(&obstacles, 0.0);
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
        let obstacles = [square()];
        let c = domain(&obstacles, 0.0);
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
