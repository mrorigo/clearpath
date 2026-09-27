//! The margin model: offsetting obstacles outward and eroding the workspace inward.
//!
//! This is the only place in the crate where `margin` is applied. See `docs/SPEC.md` section
//! 3.3.

use alloc::vec::Vec;

use super::clip::clip_to_rect;
use super::point::Point2D;
use super::polygon::{BoundingBox, Polygon};
use crate::error::PathPlanError;
use crate::{acos, atan2, ceil, cos, sin};

/// Relative radius by which the offset radius exceeds `margin`.
///
/// Corner arcs are approximated with chords inscribed in the true arc, which lie *inside* it: a
/// chord subtending `step` passes at `R * cos(step / 2)` from the corner vertex. To keep the
/// clearance at `margin` the arc radius is `R = margin * (1 + SAGITTA_RATIO)` and the step is
/// bounded by `2 * acos(1 / (1 + SAGITTA_RATIO))`, so that
/// `R * cos(step / 2) >= margin` exactly.
///
/// The cost is `theta / (2 * acos(1 / (1 + SAGITTA_RATIO)))` arc segments per convex corner, i.e.
/// six for a right angle at this ratio.
pub const SAGITTA_RATIO: f64 = 1e-2;

/// The free space a route is planned in, after `margin` has been applied.
#[derive(Clone, Debug, PartialEq)]
pub struct OffsetDomain {
    /// The workspace, eroded inward by `margin` on all four sides.
    pub workspace: BoundingBox,
    /// The obstacles, offset outward by `margin` and clipped to `workspace`.
    pub obstacles: Vec<Polygon>,
}

/// Applies the margin model: erodes `workspace` and offsets and clips every obstacle.
///
/// `margin` must be finite and `>= 0`. An obstacle that does not reach the eroded workspace
/// disappears from the domain, which is correct: it no longer constrains anything. Note that
/// the erosion does not usually *swallow* an obstacle that sits inside the workspace — an
/// obstacle at distance `d` from a wall survives iff `d >= 2 * margin`, so a large margin makes
/// the workspace smaller, not the obstacles.
///
/// An obstacle whose offset ring self-intersects — which happens when `margin` is large relative
/// to a narrow notch in that obstacle — is reported as
/// [`PathPlanError::InvalidObstacle`] with
/// [`InvalidObstacleReason::MarginTooLarge`](crate::InvalidObstacleReason::MarginTooLarge)
/// rather than being dropped. Silently dropping it would let a route pass through the obstacle.
pub fn build_domain(
    workspace: BoundingBox,
    obstacles: &[Polygon],
    margin: f64,
    max_arc_segments: u32,
) -> Result<OffsetDomain, PathPlanError> {
    if !margin.is_finite() || margin < 0.0 {
        return Err(PathPlanError::InvalidConfig("margin must be finite and >= 0"));
    }
    if workspace.is_degenerate() {
        return Err(PathPlanError::DegenerateWorkspace);
    }
    let eroded = workspace.eroded(margin);
    if eroded.is_degenerate() {
        return Err(PathPlanError::WorkspaceEroded { margin });
    }

    let mut out: Vec<Polygon> = Vec::with_capacity(obstacles.len());
    for (index, obstacle) in obstacles.iter().enumerate() {
        let grown = offset_polygon(obstacle.vertices(), margin, max_arc_segments);
        let clipped = clip_to_rect(&grown, eroded);
        if has_fewer_than_three_distinct_vertices(&clipped) {
            // The offset obstacle does not reach the eroded workspace.
            continue;
        }
        let poly = Polygon::new(clipped).map_err(|_| {
            let _ = index;
            PathPlanError::InvalidObstacle { reason: crate::InvalidObstacleReason::MarginTooLarge }
        })?;
        debug_assert_eq!(poly.area().signum(), 1.0);
        out.push(poly);
    }
    Ok(OffsetDomain { workspace: eroded, obstacles: out })
}

fn has_fewer_than_three_distinct_vertices(ring: &[Point2D]) -> bool {
    if ring.len() < 3 {
        return true;
    }
    let mut distinct = 0usize;
    for (i, p) in ring.iter().enumerate() {
        if !ring[..i].contains(p) {
            distinct += 1;
        }
    }
    distinct < 3
}

/// The outward offset of a counter-clockwise ring by `radius`.
///
/// The result is the boundary of the `radius`-neighbourhood of the ring's enclosed region:
/// each edge is translated along its right-hand normal, each convex corner is filled with a
/// circular arc centred on the corner vertex, and each reflex corner is cut back to where the two
/// neighbouring offset lines cross. That last point matters: a reflex corner's tangent points sit
/// at distance `radius` from the vertex, and the chord between them passes closer than `radius`,
/// so emitting the chord would break the clearance guarantee.
///
/// `radius == 0` returns the input unchanged.
pub fn offset_polygon(vertices: &[Point2D], radius: f64, max_arc_segments: u32) -> Vec<Point2D> {
    if radius == 0.0 || vertices.len() < 3 {
        return vertices.to_vec();
    }
    let n = vertices.len();
    // Over-inflate so the inscribed chords still clear the original ring by `radius`.
    let r = radius * (1.0 + SAGITTA_RATIO);
    // Right-hand normal of a direction; the ring is counter-clockwise so outward is to the right.
    let right_normal = |d: Point2D| Point2D::new(d.y, -d.x);

    let mut out: Vec<Point2D> = Vec::with_capacity(n * 8);
    for i in 0..n {
        let prev = vertices[(i + n - 1) % n];
        let cur = vertices[i];
        let next = vertices[(i + 1) % n];

        let (Some(d_in), Some(d_out)) =
            ((cur - prev).normalize(), (next - cur).normalize())
        else {
            continue;
        };

        let turn = atan2(d_in.cross(d_out), d_in.dot(d_out));
        let start = cur + right_normal(d_in) * r;
        let end = cur + right_normal(d_out) * r;

        if turn <= 0.0 {
            // Reflex or collinear corner. The offset boundary is where the two neighbouring
            // offset lines cross; for a collinear corner that is `start` itself.
            let crossing = if turn == 0.0 {
                start
            } else {
                line_intersection(start, d_in, end, d_out).unwrap_or(start)
            };
            out.push(crossing);
            continue;
        }

        let max_step = 2.0 * acos(1.0 / (1.0 + SAGITTA_RATIO));
        let segments = if max_step > 0.0 {
            (ceil(turn / max_step) as u32).clamp(1, max_arc_segments.max(1))
        } else {
            1
        };
        let a0 = atan2(start.y - cur.y, start.x - cur.x);
        // `segments + 1` samples: the two tangent points, with sub-arc midpoints between them.
        for k in 0..=segments {
            let ang = a0 + turn * (k as f64 / segments as f64);
            out.push(cur + Point2D::new(cos(ang), sin(ang)) * r);
        }
    }
    out
}

/// Intersection of the two lines `p + t*u` and `q + s*v`, or `None` if they are parallel.
fn line_intersection(p: Point2D, u: Point2D, q: Point2D, v: Point2D) -> Option<Point2D> {
    let denom = u.cross(v);
    if denom == 0.0 {
        return None;
    }
    let t = (q - p).cross(v) / denom;
    Some(p + u * t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
        Polygon::new(vec![
            Point2D::new(x0, y0),
            Point2D::new(x1, y0),
            Point2D::new(x1, y1),
            Point2D::new(x0, y1),
        ])
        .unwrap()
    }

    fn l_shape() -> Polygon {
        Polygon::new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(3.0, 0.0),
            Point2D::new(3.0, 1.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(1.0, 3.0),
            Point2D::new(0.0, 3.0),
        ])
        .unwrap()
    }

    /// Minimum distance from `p` to the boundary of a ring.
    fn dist_to_ring(p: Point2D, ring: &[Point2D]) -> f64 {
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

    #[test]
    fn zero_radius_is_identity() {
        let s = square(0.0, 0.0, 2.0, 2.0);
        assert_eq!(offset_polygon(s.vertices(), 0.0, 64), s.vertices().to_vec());
    }

    #[test]
    fn straight_edges_and_arcs_clear_the_radius() {
        // The load-bearing property: every point of the offset ring is at distance >= radius
        // from the original ring. This is the guarantee the whole margin model rests on.
        for radius in [0.01, 0.5, 1.0, 3.0] {
            let s = square(0.0, 0.0, 2.0, 2.0);
            let grown = offset_polygon(s.vertices(), radius, 64);
            for p in &grown {
                let d = dist_to_ring(*p, s.vertices());
                assert!(d >= radius * (1.0 - 1e-12), "radius {radius}: {p:?} is only {d} out");
            }
        }
    }

    #[test]
    fn chords_midway_between_arc_samples_clear_the_radius() {
        // The critical case the sagitta bound exists for: a point on the interior of a chord,
        // which is closer to the corner vertex than either endpoint.
        for radius in [0.1, 1.0, 5.0] {
            let s = square(0.0, 0.0, 2.0, 2.0);
            let grown = offset_polygon(s.vertices(), radius, 64);
            for i in 0..grown.len() {
                let a = grown[i];
                let b = grown[(i + 1) % grown.len()];
                for f in [0.25f64, 0.5, 0.75] {
                    let d = dist_to_ring(a.lerp(b, f), s.vertices());
                    assert!(
                        d >= radius * (1.0 - 1e-12),
                        "radius {radius}: chord interior {a:?}..{b:?} only {d} out"
                    );
                }
            }
        }
    }

    #[test]
    fn offset_contains_the_true_offset() {
        // Points on the true offset boundary (a cardinal direction from an edge midpoint) must
        // be inside the offset ring, or the offset is too small.
        let s = square(0.0, 0.0, 2.0, 2.0);
        let radius = 0.5;
        let grown = Polygon::new(offset_polygon(s.vertices(), radius, 64)).unwrap();
        for p in [
            Point2D::new(-radius, 1.0),
            Point2D::new(2.0 + radius, 1.0),
            Point2D::new(1.0, -radius),
            Point2D::new(1.0, 2.0 + radius),
        ] {
            assert!(grown.contains(p), "true offset point {p:?} escaped the offset ring");
        }
    }

    #[test]
    fn reflex_corner_is_cut_to_the_offset_lines() {
        // The L shape's inner corner is reflex. The offset boundary there must pass at distance
        // >= radius from that corner, which a chord between the tangent points would not.
        let radius = 0.3;
        let l = l_shape();
        let inner = Point2D::new(1.0, 1.0);
        let grown = offset_polygon(l.vertices(), radius, 64);
        for i in 0..grown.len() {
            let a = grown[i];
            let b = grown[(i + 1) % grown.len()];
            for f in [0.0f64, 0.5, 1.0] {
                let p = a.lerp(b, f);
                assert!(
                    p.distance(inner) >= radius - 1e-12,
                    "offset boundary passes {p:?} within {} of the reflex corner",
                    p.distance(inner)
                );
            }
        }
    }

    #[test]
    fn offset_ring_stays_simple_for_a_modest_margin() {
        // `Polygon::new` rejects a self-intersecting ring, so accepting it is the assertion.
        for radius in [0.01, 0.2, 0.5, 1.0] {
            let square_ring = offset_polygon(square(0.0, 0.0, 2.0, 2.0).vertices(), radius, 64);
            assert!(Polygon::new(square_ring).is_ok(), "square offset at {radius}");
            let l_ring = offset_polygon(l_shape().vertices(), radius, 64);
            assert!(Polygon::new(l_ring).is_ok(), "L offset at {radius}");
        }
    }

    #[test]
    fn right_angle_costs_six_arc_segments() {
        let s = square(0.0, 0.0, 2.0, 2.0);
        let grown = offset_polygon(s.vertices(), 1.0, 64);
        // Seven samples per corner: two tangent points plus five sub-arc midpoints.
        assert_eq!(grown.len(), 4 * 7);
    }

    #[test]
    fn arc_segment_cap_is_respected() {
        let s = square(0.0, 0.0, 2.0, 2.0);
        let grown = offset_polygon(s.vertices(), 1.0, 2);
        assert_eq!(grown.len(), 4 * 3);
    }

    #[test]
    fn domain_erodes_the_workspace() {
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let domain = build_domain(ws, &[], 1.0, 64).unwrap();
        assert_eq!(domain.workspace, BoundingBox::new(Point2D::new(1.0, 1.0), Point2D::new(9.0, 9.0)));
        assert!(domain.obstacles.is_empty());
    }

    #[test]
    fn domain_offsets_and_clips_obstacles() {
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(20.0, 20.0));
        let obstacles = vec![square(4.0, 4.0, 6.0, 6.0)];
        let domain = build_domain(ws, &obstacles, 0.5, 64).unwrap();
        assert_eq!(domain.obstacles.len(), 1);
        let b = domain.obstacles[0].bounds();
        assert!(b.min.x < 4.0 && b.max.x > 6.0, "obstacle not grown: {b:?}");
    }

    #[test]
    fn a_large_margin_shrinks_the_workspace_rather_than_the_obstacles() {
        // An obstacle hugging a wall is *not* swallowed: its offset reaches deeper into the
        // eroded workspace than the wall does. Pinning this down, because the intuitive
        // expectation is the opposite and getting it wrong would silently unblock a route.
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let obstacles = vec![square(0.0, 0.0, 1.0, 1.0)];
        let domain = build_domain(ws, &obstacles, 2.0, 64).unwrap();
        assert_eq!(domain.obstacles.len(), 1);
        let b = domain.obstacles[0].bounds();
        assert!(b.max.x > 2.0, "offset should reach past the eroded wall at x=2: {b:?}");
    }

    #[test]
    fn obstacle_outside_the_workspace_is_dropped() {
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        let obstacles = vec![square(20.0, 20.0, 30.0, 30.0)];
        let domain = build_domain(ws, &obstacles, 0.0, 64).unwrap();
        assert!(domain.obstacles.is_empty());
    }

    #[test]
    fn margin_larger_than_the_workspace_is_an_error() {
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        assert_eq!(build_domain(ws, &[], 5.0, 64).unwrap_err(), PathPlanError::WorkspaceEroded { margin: 5.0 });
    }

    #[test]
    fn negative_or_nan_margin_is_rejected() {
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(10.0, 10.0));
        assert!(matches!(build_domain(ws, &[], -1.0, 64), Err(PathPlanError::InvalidConfig(_))));
        assert!(matches!(build_domain(ws, &[], f64::NAN, 64), Err(PathPlanError::InvalidConfig(_))));
    }

    #[test]
    fn degenerate_workspace_is_rejected() {
        let ws = BoundingBox { min: Point2D::ZERO, max: Point2D::ZERO };
        assert_eq!(build_domain(ws, &[], 0.0, 64).unwrap_err(), PathPlanError::DegenerateWorkspace);
    }

    #[test]
    fn a_margin_that_swallows_a_notch_is_reported_not_dropped() {
        // A narrow slot: a large enough offset ring self-intersects. The obstacle must not
        // silently vanish from the domain.
        let slot = Polygon::new(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(10.0, 0.0),
            Point2D::new(10.0, 1.0),
            Point2D::new(2.0, 1.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(10.0, 2.0),
            Point2D::new(10.0, 3.0),
            Point2D::new(0.0, 3.0),
        ])
        .unwrap();
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(20.0, 20.0));
        // margin 0.6 fills the 1.0-wide slot exactly; 1.0 would over-shoot and self-intersect.
        assert!(build_domain(ws, core::slice::from_ref(&slot), 0.4, 64).is_ok());
        let err = build_domain(ws, &[slot], 1.5, 64).unwrap_err();
        assert_eq!(
            err,
            PathPlanError::InvalidObstacle { reason: crate::InvalidObstacleReason::MarginTooLarge }
        );
    }
}
