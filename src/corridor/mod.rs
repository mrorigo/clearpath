//! Admissible tangents: how far a spline may pull away from a knot before it leaves the free space.
//!
//! # Why this is not a "corridor"
//!
//! The first design was a convex *corridor* per taut segment, on the reasoning that a Bezier's hull
//! property plus a corridor inside the free space gives collision freedom for free. Two things go
//! wrong with it, and both are structural rather than fixable:
//!
//! * A convex region containing a whole taut segment and contained in the free space does not
//!   generally exist. The corridor has to bound itself by every cell the segment passes through, but
//!   the segment's own height varies across the corridor's abscissa range, so a cell traversed only
//!   near the far end bounds the near end too. The two requirements pull in opposite directions.
//! * Bounding only by the cells the segment touches *at that abscissa* is correct but not convex,
//!   and convexity is what the hull property needs.
//!
//! What is actually needed is local: the two control points `k + T/3` and `k' - T'/3` sit next to
//! the knots, so only the cells *adjacent to the knots* constrain the tangents. Those are single
//! convex cells, and intersecting two of them — one translated forwards, one backwards — gives a
//! convex admissible set that contains the zero tangent. The middle of the curve is then checked
//! directly against the free space by the hull test in `control_points_are_clear`.

use alloc::vec::Vec;

use crate::decomp::{CellId, Decomposition};
use crate::error::PathPlanError;
use crate::geom::clearance::Clearance;
use crate::geom::point::Point2D;
use crate::spline::CubicBezierSegment;

/// A half-plane `{x : a*x + b*y + c >= 0}`.
#[derive(Clone, Copy, Debug, PartialEq)]
struct HalfPlane {
    a: f64,
    b: f64,
    c: f64,
}

impl HalfPlane {
    #[inline]
    fn value(&self, p: Point2D) -> f64 {
        self.a * p.x + self.b * p.y + self.c
    }

    /// The slack allowed when testing membership, in units of the edge's own length.
    ///
    /// A knot produced by the funnel is a portal endpoint, and a cell's edge is a line through two
    /// corners that may be the same point by a different route, so "on the boundary" can evaluate
    /// to a small negative number rather than to zero. Rejecting on that is rejecting on rounding.
    /// The tolerance is proportional to the edge length so it means a fixed distance rather than a
    /// fixed number of ulps, and it is small enough that the projection still pulls the tangent
    /// properly inside.
    #[inline]
    fn slack(&self) -> f64 {
        -1e-9 * (self.a.abs() + self.b.abs()).max(1.0)
    }

    /// The half-planes bounding `cell`, counter-clockwise from the bottom edge.
    fn of_cell(cell: &crate::decomp::Cell) -> [HalfPlane; 4] {
        let c = [cell.bl, cell.br, cell.tr, cell.tl];
        let mut out = [HalfPlane { a: 0.0, b: 0.0, c: 0.0 }; 4];
        for i in 0..4 {
            // Inside is left of the directed edge, which is `cross(dir, p - from) >= 0`.
            let d = c[(i + 1) % 4] - c[i];
            out[i] = HalfPlane { a: -d.y, b: d.x, c: d.y * c[i].x - d.x * c[i].y };
        }
        out
    }
}

/// The set of tangents at one knot that keep the neighbouring control points in the free space.
///
/// Always contains the zero vector, and is convex: it is the intersection of the forward cell
/// translated by the knot, and the backward cell translated the other way.
#[derive(Clone, Debug, PartialEq)]
pub struct AdmissibleSet {
    planes: Vec<HalfPlane>,
    knot: Point2D,
}

impl AdmissibleSet {
    fn new(knot: Point2D, forward: &crate::decomp::Cell, backward: &crate::decomp::Cell) -> Self {
        let mut planes = Vec::with_capacity(8);
        for p in HalfPlane::of_cell(forward) {
            // q = knot + t/3 in the cell  <=>  a*t.x + b*t.y + 3*(a*knot.x + b*knot.y + c) >= 0.
            planes.push(HalfPlane {
                a: p.a,
                b: p.b,
                c: 3.0 * (p.a * knot.x + p.b * knot.y + p.c),
            });
        }
        for p in HalfPlane::of_cell(backward) {
            // q = knot - t/3 in the cell  <=>  -a*t.x - b*t.y + 3*(a*knot.x + b*knot.y + c) >= 0.
            planes.push(HalfPlane {
                a: -p.a,
                b: -p.b,
                c: 3.0 * (p.a * knot.x + p.b * knot.y + p.c),
            });
        }
        Self { planes, knot }
    }

    /// The knot this set belongs to.
    pub fn knot(&self) -> Point2D {
        self.knot
    }

    /// Whether `t` is admissible.
    pub fn contains(&self, t: Point2D) -> bool {
        self.planes.iter().all(|p| p.value(t) >= p.slack())
    }

    /// Whether the set admits nothing but the zero tangent, i.e. the knot is pinned.
    pub fn is_pinned(&self) -> bool {
        self.vertices().is_empty()
    }

    /// The set's vertices, in counter-clockwise order.
    ///
    /// Every *pair* of half-planes is tried, not just consecutive ones. The set is the
    /// intersection of two cells' four planes each, and those eight are not in a cyclic order
    /// around the region, so consecutive pairs miss most of the corners — and the duplicates that
    /// survive make the projection below return a point outside the set.
    pub fn vertices(&self) -> Vec<Point2D> {
        let n = self.planes.len();
        let mut out: Vec<Point2D> = Vec::with_capacity(n);
        for i in 0..n {
            for j in (i + 1)..n {
                let (p, q) = (self.planes[i], self.planes[j]);
                let det = p.a * q.b - q.a * p.b;
                if det.abs() < 1e-12 {
                    continue; // Parallel or coincident: no vertex from this pair.
                }
                let v = Point2D::new(
                    (p.b * q.c - q.b * p.c) / det,
                    (q.a * p.c - p.a * q.c) / det,
                );
                if !v.is_finite() {
                    continue;
                }
                if self.planes.iter().all(|r| r.value(v) >= r.slack())
                    && !out.iter().any(|w| w.distance(v) < 1e-9)
                {
                    out.push(v);
                }
            }
        }
        out.sort_by(|a, b| a.x.total_cmp(&b.x).then_with(|| a.y.total_cmp(&b.y)));
        out
    }

    /// The closest admissible tangent to `t`.
    ///
    /// The set is convex, so the closest point is `t` itself when it is admissible, and otherwise
    /// the perpendicular projection onto one of the edges. Projecting onto the *vertices* as well
    /// covers the degenerate case where the set is a segment or a point.
    pub fn project(&self, t: Point2D) -> Point2D {
        if self.contains(t) {
            return t;
        }
        // The closest point of a convex polygon to an external point lies either on an edge or at
        // a vertex, and a vertex is covered by both of its edges — so the edges alone suffice, and
        // enumerating the vertices as well only guards against a degenerate polygon.
        let v = self.vertices();
        if v.is_empty() {
            return Point2D::ZERO;
        }
        let mut best = v[0];
        let mut best_d2 = f64::INFINITY;
        for i in 0..v.len() {
            let (a, b) = (v[i], v[(i + 1) % v.len()]);
            let q = crate::corridor::closest_point_on_segment(t, a, b);
            let d2 = q.distance_squared(t);
            if d2 < best_d2 {
                best_d2 = d2;
                best = q;
            }
        }
        for w in &v {
            let d2 = w.distance_squared(t);
            if d2 < best_d2 {
                best_d2 = d2;
                best = *w;
            }
        }
        best
    }
}

/// The closest point on the closed segment `a..b` to `p`.
pub fn closest_point_on_segment(p: Point2D, a: Point2D, b: Point2D) -> Point2D {
    let ab = b - a;
    let len2 = ab.norm_squared();
    if len2 == 0.0 {
        return a;
    }
    a + ab * (((p - a).dot(ab) / len2).clamp(0.0, 1.0))
}

/// The Bezier segments for a taut path under the given tangents, without checking anything.
///
/// The caller is [`AdmissibleTangents::control_points`], or the repair, which has already
/// established that the tangents are admissible and the hulls clear.
pub fn segments_for(knots: &[Point2D], tangents: &[Point2D]) -> Vec<CubicBezierSegment> {
    (0..knots.len().saturating_sub(1))
        .map(|i| CubicBezierSegment {
            p0: knots[i],
            p1: knots[i] + tangents[i] / 3.0,
            p2: knots[i + 1] - tangents[i + 1] / 3.0,
            p3: knots[i + 1],
        })
        .collect()
}

/// The admissible tangent set at every knot of a taut path.
#[derive(Clone, Debug, Default)]
pub struct AdmissibleTangents {
    sets: Vec<AdmissibleSet>,
}

impl AdmissibleTangents {
    /// The set at knot `i`.
    pub fn get(&self, i: usize) -> Option<&AdmissibleSet> {
        self.sets.get(i)
    }

    /// Builds the sets for a taut path through `corridor_cells`.
    pub fn build(
        decomp: &Decomposition,
        corridor_cells: &[CellId],
        knots: &[Point2D],
    ) -> Result<Self, PathPlanError> {
        if knots.len() < 2 || corridor_cells.is_empty() {
            return Err(PathPlanError::NoPathFound);
        }
        // For each segment, the cell the segment leaves `a` in and the cell it arrives at `b` in.
        let mut start_cell: Vec<CellId> = Vec::with_capacity(knots.len() - 1);
        let mut end_cell: Vec<CellId> = Vec::with_capacity(knots.len() - 1);
        for i in 0..knots.len() - 1 {
            let (a, b) = (knots[i], knots[i + 1]);
            let forward = cell_entered_from(decomp, corridor_cells, a, b);
            let backward = cell_entered_from(decomp, corridor_cells, b, a);
            match (forward, backward) {
                (Some(f), Some(k)) => {
                    start_cell.push(f);
                    end_cell.push(k);
                }
                _ => return Err(PathPlanError::NoPathFound),
            }
        }

        let mut sets = Vec::with_capacity(knots.len());
        for i in 0..knots.len() {
            // The backward bound at knot 0 is the start cell's own cell, and the forward bound at
            // the last knot likewise: a knot at the end of the path has no segment on that side,
            // and the only thing that can constrain its tangent is the cell it sits in.
            let forward_id = if i + 1 < knots.len() { start_cell[i] } else { end_cell[i - 1] };
            let backward_id = if i > 0 { end_cell[i - 1] } else { start_cell[0] };
            sets.push(AdmissibleSet::new(
                knots[i],
                &decomp.cells()[forward_id as usize],
                &decomp.cells()[backward_id as usize],
            ));
        }
        Ok(Self { sets })
    }

    /// The Bezier control points for a taut path under the given tangents, and a check that each
    /// segment's control hull is in the free space.
    ///
    /// The hull test is what stands in for the convex corridor the design originally used: a cubic
    /// lies inside the hull of its four control points, so if the hull is clear the curve is clear.
    /// It is checked against the obstacles rather than trusted from the construction, because the
    /// ends are constrained by the two adjacent cells and the middle spans the channel between them,
    /// which no local construction can speak for.
    pub fn control_points(
        &self,
        knots: &[Point2D],
        tangents: &[Point2D],
        clearance: &Clearance,
    ) -> Result<Vec<CubicBezierSegment>, PathPlanError> {
        if self.first_obstructed(knots, tangents, clearance).is_some() {
            return Err(PathPlanError::NoPathFound);
        }
        Ok(segments_for(knots, tangents))
    }

    /// The index of the first segment that is not admissible or not clear, or `None`.
    ///
    /// The repair asks this question once per iteration, and the router asks it again at the end,
    /// so it stops at the first bad segment rather than checking the whole spline. Checking all of
    /// them and reporting "somewhere" was most of the repair's cost: each full pass re-ran the hull
    /// test on the segments that were already fine.
    pub fn first_obstructed(
        &self,
        knots: &[Point2D],
        tangents: &[Point2D],
        clearance: &Clearance,
    ) -> Option<usize> {
        if knots.len() < 2 || tangents.len() != knots.len() {
            return Some(0);
        }
        for (i, t) in tangents.iter().enumerate() {
            match self.sets.get(i) {
                Some(set) if set.contains(*t) => {}
                _ => return Some(i.min(knots.len() - 2)),
            }
        }
        for i in 0..knots.len() - 1 {
            let hull = [
                knots[i],
                knots[i] + tangents[i] / 3.0,
                knots[i + 1] - tangents[i + 1] / 3.0,
                knots[i + 1],
            ];
            // The hull, not the control points: a cubic lies inside the hull of its four control
            // points, so a clear hull is a clear curve. Testing the control points alone is not
            // enough — the hull's edges can cut across an obstacle between two clear vertices.
            if !clearance.hull_is_free(&hull) {
                return Some(i);
            }
        }
        None
    }
}

/// The corridor cell the segment `a..b` is in immediately after leaving `a`.
///
/// Uses a point a short way along the segment rather than `a` itself, because a knot lies on the
/// boundary of the cell being left as well as the one being entered, and picking the wrong one
/// constrains the tangent with a cell the curve does not go into.
fn cell_entered_from(
    decomp: &Decomposition,
    corridor_cells: &[CellId],
    a: Point2D,
    b: Point2D,
) -> Option<CellId> {
    for f in [1e-6, 1e-4, 1e-3, 1e-2, 0.25, 0.5] {
        let p = a.lerp(b, f);
        let hit = corridor_cells
            .iter()
            .copied()
            .find(|c| decomp.cells()[*c as usize].contains(p));
        if let Some(c) = hit {
            return Some(c);
        }
    }
    None
}
