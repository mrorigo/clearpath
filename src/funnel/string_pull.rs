//! Lee & Preparata's funnel algorithm.
//!
//! Given a corridor of convex cells, the shortest path from `start` to `goal` inside the corridor
//! is the taut string through its portals. The funnel algorithm finds it in `$O(N)$` in the number
//! of portals, maintaining a wedge with its apex at the last committed knot and tightening it
//! portal by portal; when one side of the wedge crosses the other, the knot is committed and the
//! apex jumps there.

use alloc::vec;
use alloc::vec::Vec;

use super::Corridor;
use crate::decomp::Decomposition;
use crate::error::PathPlanError;
use crate::geom::point::Point2D;
use crate::ceil;

/// The taut string through a corridor: the ordered knots, including `start` and `goal`.
///
/// Invariant L2 (section 4.7): every knot lies on the boundary of one of the corridor's cells, and
/// every point of every segment lies inside the union of the corridor's cells. The M4 gate tests
/// that directly.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TautPath {
    /// The knots, from `start` to `goal`.
    pub knots: Vec<Point2D>,
}

impl TautPath {
    /// Number of knots.
    #[inline]
    pub fn len(&self) -> usize {
        self.knots.len()
    }

    /// Whether the path has no knots at all. A taut path always has at least two (start and goal),
    /// so this exists only to satisfy the usual container convention.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.knots.is_empty()
    }
}

/// Twice the signed area of `a`, `b`, `c`, scaled by `s`.
///
/// The funnel's inequalities are written against a *negated* cross product, and the negation has to
/// be undone for a corridor that runs right to left: mirroring `x` flips the orientation of every
/// triangle, so travelling in `-x` needs the opposite sign. `s` is `-1.0` while travelling in `+x`
/// and `+1.0` while travelling in `-x`.
///
/// This is a sign convention, not a predicate improvement, and getting it wrong is quiet: the
/// funnel still returns a plausible polyline, it just leaves the corridor. The L2 gate is what
/// catches that, which is why it samples the emitted path against the corridor's own cells rather
/// than trusting the funnel.
///
/// The area itself is deliberately not computed exactly. The funnel's decisions are about the sign
/// of a turn in a sequence of points lying exactly on cell boundaries, where a rounding error would
/// move a knot by a measurable amount.
#[inline]
fn triarea2(s: f64, a: Point2D, b: Point2D, c: Point2D) -> f64 {
    s * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y))
}

/// Pulls the taut string through `corridor` from `start` to `goal`.
///
/// The first and last portals are the degenerate single points `{start}` and `{goal}`, per Lee &
/// Preparata; the portal sequence is otherwise the corridor's own portals, in order.
pub fn string_pull(
    decomp: &Decomposition,
    corridor: &Corridor,
    start: Point2D,
    goal: Point2D,
) -> Result<TautPath, PathPlanError> {
    if corridor.cells.is_empty() {
        return Err(PathPlanError::NoPathFound);
    }
    if corridor.cells.len() == 1 {
        // Start and goal share a convex cell, so the straight segment is optimal and stays inside.
        return Ok(TautPath { knots: vec![start, goal] });
    }

    // Portals in travel order, with `left`/`right` resolved against the direction of travel
    // *through that portal*, not globally.
    //
    // A corridor is not necessarily monotone in x. The cell graph fans out — a cell's right side
    // can overlap several cells of the next slab — so A* can legitimately produce a corridor that
    // steps right and then steps back left, which is how it serves two endpoints in the same slab.
    // A single global orientation is then wrong for half the portals, and the funnel commits knots
    // from the wrong chain: the path still looks like a plausible path, but it leaves the corridor.
    //
    // `left` is always the upper endpoint and `right` the lower one, and the *sign* passed to
    // `triarea2` carries the direction of travel instead. Splitting it that way — rather than
    // swapping which endpoint is called `left` — is what keeps the two conventions from drifting
    // apart, which they did while both were carrying part of the meaning.
    let mut xs: Vec<f64> = Vec::with_capacity(corridor.portals.len() + 2);
    xs.push(start.x);
    for id in &corridor.portals {
        xs.push(decomp.portals()[*id as usize].x);
    }
    xs.push(goal.x);

    // The orientation sign for each portal, precomputed so that the loop does not need an
    // out-of-range lookahead at the final (degenerate goal) portal, which inherits the direction of
    // the portal before it.
    let mut signs: Vec<f64> = Vec::with_capacity(xs.len());
    for k in 0..xs.len() {
        signs.push(if k + 1 < xs.len() && xs[k + 1] >= xs[k] { -1.0 } else { 1.0 });
    }
    if let Some(last) = signs.len().checked_sub(2).and_then(|i| signs.get(i).copied()) {
        if let Some(slot) = signs.last_mut() {
            *slot = last;
        }
    }

    let mut portals: Vec<(Point2D, Point2D)> = Vec::with_capacity(xs.len());
    for k in 0..xs.len() {
        if k == 0 || k + 1 == xs.len() {
            let point = if k == 0 { start } else { goal };
            portals.push((point, point));
            continue;
        }
        let p = decomp.portals()[corridor.portals[k - 1] as usize];
        portals.push((p.top(), p.bottom()));
    }

    let mut path: Vec<Point2D> = Vec::new();
    // `apex` is the last committed knot; `left`/`right` are the funnel's current sides, and the
    // `*_index` values are where in the portal list each side came from, so the scan can be
    // restarted at the apex. The sides are scalars rather than a stack of points: the funnel
    // only ever needs the current pair, and keeping them as indices is a use-after-free waiting to
    // happen on a restart.
    let mut apex = start;
    // Which portal produced each side's current point. The scan restarts at `from + 1` after a
    // commit, so this has to be the *producing portal*, not the length of the chain: the two
    // differ by however many times each side has been retightened, and using the length restarts
    // one portal too far — which silently drops a knot at a slot's end and lets a segment leave the
    // corridor through the obstacle it was supposed to go around.
    let mut left_from = 0usize;
    let mut right_from = 0usize;
    let mut left = start;
    let mut right = start;

    let mut i = 0usize;
    while i < portals.len() {
        let (new_left, new_right) = portals[i];
        // The sign for the portal being processed: `-x` travel inverts every orientation.
        let s = signs[i];


        // Tighten the right side. A right point that is not strictly right of the apex narrows the
        // wedge. If it also crosses past the left side, the left chain is the shorter way round, so
        // the left point becomes a knot and the wedge restarts from it.
        if triarea2(s, apex, right, new_right) <= 0.0 {
            if apex == right || triarea2(s, apex, left, new_right) > 0.0 {
                right = new_right;
                right_from = i;
            } else {
                path.push(left);
                apex = left;
                i = left_from + 1;
                left = apex;
                right = apex;
                left_from = i;
                right_from = i;
                continue;
            }
        }

        // Mirror image for the left side.
        if triarea2(s, apex, left, new_left) >= 0.0 {
            if apex == left || triarea2(s, apex, right, new_left) < 0.0 {
                left = new_left;
                left_from = i;
            } else {
                path.push(right);
                apex = right;
                i = right_from + 1;
                left = apex;
                right = apex;
                left_from = i;
                right_from = i;
                continue;
            }
        }

        i += 1;
    }

    // The funnel emits only the *committed corners*; the start and goal are prepended and appended
    // here. Omitting the start is a silent bug that looks like a path that begins at a corner.
    path.push(goal);
    path.insert(0, start);
    let path = drop_redundant_knots(path)?;

    // Safety net. The funnel is proven for a sleeve whose portals are ordered along the path, and
    // the search is arranged to hand it one (monotone moves, with the endpoints' abscissae forced
    // into the event set). A corridor can still fail that precondition — the endpoints sharing an
    // abscissa, which no event line can separate — and then the funnel returns a path that *looks*
    // plausible and leaves the corridor. Refusing is the only safe answer: the alternative is a
    // collision.
    if !inside_corridor(decomp, &corridor.cells, &path.knots) {
        return Err(PathPlanError::NoPathFound);
    }
    Ok(path)
}

/// Whether every point of the polyline lies inside one of `cells`.
fn inside_corridor(
    decomp: &Decomposition,
    cells: &[crate::decomp::CellId],
    knots: &[Point2D],
) -> bool {
    // Two samples per unit length: the failure this guards against is a segment cutting across an
    // obstacle edge, so the density only has to be fine enough to hit the edge.
    for seg in knots.windows(2) {
        let steps = (ceil(seg[0].distance(seg[1]) * 2.0) as usize).max(1);
        for k in 0..=steps {
            // The endpoint has to be taken verbatim rather than lerped to: `a + (b - a) * 1.0` is
            // not always exactly `b`, and a knot that lands a few ULP outside a cell's abscissa
            // fails the exact containment test. The knots are portal endpoints, which are exact.
            let p = if k == steps { seg[1] } else { seg[0].lerp(seg[1], k as f64 / steps as f64) };
            if !cells.iter().any(|c| decomp.cells()[*c as usize].contains(p)) {
                return false;
            }
        }
    }
    true
}

/// Removes knots that add nothing: duplicates and points collinear with their neighbours.
///
/// A collinear interior knot is redundant, and leaving it in would make Stage 4 emit a spline
/// segment that is a straight line, and would make the section 9.3 straight-line test ambiguous.
fn drop_redundant_knots(mut knots: Vec<Point2D>) -> Result<TautPath, PathPlanError> {
    if knots.len() < 2 {
        return Err(PathPlanError::NoPathFound);
    }
    // A portal's two endpoints coincide when the corridor has a pinch, so exact duplicates are
    // expected rather than exceptional.
    knots.dedup();
    if knots.len() < 2 {
        // Start and goal coincide after dedup: the caller should have rejected this earlier, but
        // reporting it here keeps the funnel total.
        return Ok(TautPath { knots });
    }
    let mut out: Vec<Point2D> = Vec::with_capacity(knots.len());
    out.push(knots[0]);
    for w in knots.windows(3) {
        let (a, b, c) = (w[0], w[1], w[2]);
        let degenerate = (b - a).is_zero() || (c - b).is_zero() || (b - a).cross(c - b).abs() <= 0.0;
        if !degenerate {
            out.push(b);
        }
    }
    out.push(knots[knots.len() - 1]);
    out.dedup();
    Ok(TautPath { knots: out })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decomp::cell::{Cell, Portal, PortalId};
    use crate::decomp::Decomposition;
    use alloc::vec;

    /// A hand-built corridor reproducing the shape that broke the funnel: a tall cell, a tall
    /// cell, a horizontal *slot*, a tall cell, a tall cell. The slot is the point of the test — the
    /// taut path has to touch both of the slot's ends, and a funnel that misses either produces a
    /// segment that leaves the corridor through the obstacle between them.
    fn slot_corridor() -> (Decomposition, Corridor) {
        let cells = vec![
            Cell { bl: Point2D::new(25.0, 47.0), br: Point2D::new(31.0, 47.0), tr: Point2D::new(31.0, 99.0), tl: Point2D::new(25.0, 99.0) },
            Cell { bl: Point2D::new(23.0, 1.0), br: Point2D::new(25.0, 1.0), tr: Point2D::new(25.0, 99.0), tl: Point2D::new(23.0, 99.0) },
            Cell { bl: Point2D::new(17.0, 39.0), br: Point2D::new(23.0, 39.0), tr: Point2D::new(23.0, 41.0), tl: Point2D::new(17.0, 41.0) },
            Cell { bl: Point2D::new(15.0, 1.0), br: Point2D::new(17.0, 1.0), tr: Point2D::new(17.0, 99.0), tl: Point2D::new(15.0, 99.0) },
            Cell { bl: Point2D::new(9.0, 47.0), br: Point2D::new(15.0, 47.0), tr: Point2D::new(15.0, 99.0), tl: Point2D::new(9.0, 99.0) },
        ];
        let portals = vec![
            Portal { x: 25.0, lo: 47.0, hi: 99.0, left: 1, right: 0 },
            Portal { x: 23.0, lo: 39.0, hi: 41.0, left: 2, right: 1 },
            Portal { x: 17.0, lo: 39.0, hi: 41.0, left: 3, right: 2 },
            Portal { x: 15.0, lo: 47.0, hi: 99.0, left: 4, right: 3 },
        ];
        let d = Decomposition {
            xs: vec![9.0, 15.0, 17.0, 23.0, 25.0, 31.0],
            slab_cells: Vec::new(),
            cells,
            portals,
            right_index: crate::decomp::cell::SideIndex::build(5, |_| &[]),
            left_index: crate::decomp::cell::SideIndex::build(5, |_| &[]),
        };
        let corridor = Corridor {
            cells: vec![0, 1, 2, 3, 4],
            portals: vec![0 as PortalId, 1, 2, 3],
        };
        (d, corridor)
    }

    #[test]
    fn the_funnel_touches_both_ends_of_a_slot() {
        let (d, corridor) = slot_corridor();
        let start = Point2D::new(30.0, 50.0);
        let goal = Point2D::new(10.0, 50.0);
        let path = string_pull(&d, &corridor, start, goal).unwrap();
        let want = vec![
            start,
            Point2D::new(25.0, 47.0),
            Point2D::new(23.0, 41.0),
            Point2D::new(17.0, 41.0),
            Point2D::new(15.0, 47.0),
            goal,
        ];
        assert_eq!(path.knots, want);
    }

    #[test]
    fn the_funnel_stays_inside_the_slot_corridor() {
        let (d, corridor) = slot_corridor();
        let start = Point2D::new(30.0, 50.0);
        let goal = Point2D::new(10.0, 50.0);
        let path = string_pull(&d, &corridor, start, goal).unwrap();
        for seg in path.knots.windows(2) {
            let steps = ((seg[0].distance(seg[1]) * 20.0).ceil() as usize).max(1);
            for k in 0..=steps {
                let p = seg[0].lerp(seg[1], k as f64 / steps as f64);
                assert!(
                    corridor.cells.iter().any(|c| d.cells()[*c as usize].contains(p)),
                    "sample {p:?} on {:?}..{:?} leaves the corridor",
                    seg[0],
                    seg[1]
                );
            }
        }
    }
}
