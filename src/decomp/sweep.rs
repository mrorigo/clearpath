//! The sweep that builds the vertical decomposition.

use alloc::vec::Vec;

use super::cell::{Cell, CellBuilder, CellId, Portal};
use super::Decomposition;
use crate::error::PathPlanError;
use crate::geom::point::Point2D;
use crate::geom::polygon::{BoundingBox, Polygon, rectilinear_rectangles};
use crate::geom::predicates::{Orientation, orient2d};
use crate::geom::FreeSpace;

/// A boundary edge: one edge of an offset obstacle, or one of the workspace walls.
#[derive(Clone, Copy, Debug)]
struct Edge {
    a: Point2D,
    b: Point2D,
    /// Which obstacle this edge belongs to, or `None` for a workspace wall.
    obstacle: Option<u32>,
    /// Index of the edge within its obstacle, used to keep the sort stable and deterministic.
    index: u32,
}

impl Edge {
    /// The x range the edge covers. A vertical edge has an empty range and is active in no slab.
    #[inline]
    fn x_span(&self) -> (f64, f64) {
        let lo = self.a.x.min(self.b.x);
        let hi = self.a.x.max(self.b.x);
        (lo, hi)
    }

    /// The y of the edge at `x`, by exact linear interpolation.
    #[inline]
    fn y_at(&self, x: f64) -> f64 {
        let dx = self.b.x - self.a.x;
        if dx == 0.0 {
            return self.a.y;
        }
        self.a.y + (self.b.y - self.a.y) * ((x - self.a.x) / dx)
    }

    /// Whether the edge covers the whole open slab `(lo, hi)`.
    #[inline]
    fn covers(&self, lo: f64, hi: f64) -> bool {
        let (a, b) = self.x_span();
        a <= lo && hi <= b
    }
}

/// Builds the vertical decomposition of the free space, with `margin` applied.
///
/// `margin` is realised *before* the sweep, by insetting the forbidden geometry: the workspace is
/// eroded and every rectilinear obstacle is partitioned into rectangles and each rectangle is
/// inset. That is exact — see [`rectilinear_rectangles`] for why the inset of a union of
/// rectangles is the erosion of the union — and it leaves the sweep working on plain
/// axis-aligned boxes, which is both the cheapest and the most robust thing to decompose.
///
/// Eroding the *cells* instead, or offsetting the obstacle rings, is what this replaces. Both are
/// wrong: a cell's own faces do not account for obstacle features just beyond a portal, and
/// offsetting a ring's boundary needs trimming between offset features, without which a `margin`
/// large relative to a feature produces geometry that cuts through the obstacle it came from.
/// A `margin` against a non-rectilinear obstacle is therefore rejected rather than approximated;
/// see [`PathPlanError::MarginUnsupportedGeometry`].
pub fn decompose(space: &FreeSpace, margin: f64) -> Result<Decomposition, PathPlanError> {
    if !margin.is_finite() || margin < 0.0 {
        return Err(PathPlanError::InvalidConfig("margin must be finite and >= 0"));
    }
    let ws = space.workspace;
    if ws.is_degenerate() {
        return Err(PathPlanError::DegenerateWorkspace);
    }
    let (workspace, forbidden) = apply_margin(space, margin)?;
    if workspace.is_degenerate() {
        // The margin is larger than the workspace on at least one axis, so there is no free space
        // left. An empty decomposition routes as `NoPathFound`, which is the honest answer.
        return Ok(Decomposition::default());
    }
    let edges = collect_edges(&workspace, &forbidden);
    let xs = event_coordinates(&edges, ws);

    // Degenerate free space: no cells at all.
    if xs.len() < 2 {
        return Ok(Decomposition::default());
    }


    // Probe points and the active edge list for each slab.
    let mut probes: Vec<Vec<(usize, Point2D)>> = Vec::with_capacity(xs.len() - 1);
    let mut active: Vec<Vec<Edge>> = Vec::with_capacity(xs.len() - 1);
    for k in 0..xs.len() - 1 {
        let (lo, hi) = (xs[k], xs[k + 1]);
        let mid = 0.5 * (lo + hi);
        let mut list: Vec<Edge> =
            edges.iter().copied().filter(|e| e.covers(lo, hi)).collect();
        sort_by_y_at(&mut list, mid);
        probes.push(gap_probes(&list, mid));
        active.push(list);
    }

    // Cells, per slab.
    let mut builder = CellBuilder::default();
    let mut slab_cells: Vec<Vec<CellId>> = Vec::with_capacity(xs.len() - 1);
    for k in 0..xs.len() - 1 {
        let mut ids: Vec<CellId> = Vec::with_capacity(probes[k].len());
        for (gap, probe) in &probes[k] {
            if is_free(&forbidden, *probe) && workspace.contains(*probe) {
                ids.push(builder.push_cell(cell_from_gap(&active[k], *gap, xs[k], xs[k + 1])));
            }
        }
        slab_cells.push(ids);
    }

    // Portals: at every event line, join the cells of the two adjacent slabs over their vertical
    // overlap. Both lists are sorted and disjoint in y, so this is a linear merge.
    for k in 1..xs.len() - 1 {
        let x = xs[k];
        connect_slabs(&mut builder, &slab_cells[k - 1], &slab_cells[k], x, true);
    }

    let mut cells = builder.cells;
    // Post-pass: a portal whose overlap collapsed to zero length carries no information and would
    // make the funnel see a zero-width gate. Dropping it also drops the adjacency, which is what
    // we want: two cells that touch only at a point are not connected.
    let mut remap: Vec<u32> = alloc::vec![u32::MAX; builder.portals.len()];
    let mut kept: Vec<Portal> = Vec::with_capacity(builder.portals.len());
    for (old, portal) in builder.portals.iter().enumerate() {
        if portal.hi > portal.lo {
            remap[old] = kept.len() as u32;
            kept.push(*portal);
        }
    }
    for cell in &mut cells {
        cell.left = cell.left.and_then(|p| {
            let new = remap[p as usize];
            (new != u32::MAX).then_some(new)
        });
        cell.right = cell.right.and_then(|p| {
            let new = remap[p as usize];
            (new != u32::MAX).then_some(new)
        });
    }

    Ok(Decomposition { xs, slab_cells, cells, portals: kept })
}

/// Every edge of every offset obstacle, plus the four workspace walls.
///
/// The walls are what keep the sweep inside the workspace: the bottom and top walls are always
/// active, so the gaps beyond them are never free, and the side walls close the x range.
fn collect_edges(ws: &BoundingBox, obstacles: &[Polygon]) -> Vec<Edge> {
    let mut edges: Vec<Edge> = Vec::new();
    for (i, p) in obstacles.iter().enumerate() {
        let v = p.vertices();
        for k in 0..v.len() {
            edges.push(Edge { a: v[k], b: v[(k + 1) % v.len()], obstacle: Some(i as u32), index: k as u32 });
        }
    }
    let corners = ws.corners();
    for k in 0..4 {
        edges.push(Edge {
            a: corners[k],
            b: corners[(k + 1) % 4],
            obstacle: None,
            index: 4 + k as u32,
        });
    }
    edges
}

/// The sorted, distinct event x coordinates: every vertex abscissa, plus the workspace's.
fn event_coordinates(edges: &[Edge], ws: BoundingBox) -> Vec<f64> {
    let mut xs: Vec<f64> = Vec::with_capacity(edges.len() + 2);
    for e in edges {
        xs.push(e.a.x);
        xs.push(e.b.x);
    }
    xs.push(ws.min.x);
    xs.push(ws.max.x);
    sort_f64_ascending(&mut xs);
    xs.dedup();
    // Drop event lines that the workspace walls make unreachable; a zero-width remainder would
    // produce a degenerate slab.
    xs.retain(|x| x.is_finite());
    xs
}

/// Ascending sort of `f64` as a total order on the bit patterns.
///
/// `f64::partial_cmp().unwrap()` panics on NaN, and comparing the numbers themselves is not a
/// total order on the bit patterns (`-0.0 == 0.0` but their bits differ), so neither can satisfy
/// the D2 determinism requirement of section 4.4.
///
/// The key is the standard IEEE-754 flip — negatives have every bit inverted, non-negatives get
/// the sign bit set — which orders *unsigned*, monotonically in the numeric value, with `-0.0`
/// before `+0.0` and NaN at the extremes. Comparing it as a signed integer would put every
/// non-negative value before every negative one, because the transform moves positives into the
/// upper half of the unsigned range.
pub(crate) fn sort_f64_ascending(values: &mut [f64]) {
    values.sort_by_key(|x| f64_sort_key(*x));
}

/// The total-order key for a `f64`. See [`sort_f64_ascending`].
#[inline]
pub(crate) fn f64_sort_key(x: f64) -> u64 {
    let bits = x.to_bits();
    if bits >> 63 == 0 { bits | (1u64 << 63) } else { !bits }
}

/// Sorts active edges bottom to top at abscissa `x`, with a deterministic tie-break.
///
/// The primary key is the exact orientation of the two edges at `x` — not a float comparison of
/// two separately-computed intercepts, which would disagree with the true order when the
/// intercepts are close.
fn sort_by_y_at(edges: &mut [Edge], x: f64) {
    edges.sort_by(|e1, e2| {
        // Normalise `e1` to run in the +x direction first. The sign of `orient2d` depends on the
        // direction of the directed edge, so comparing two edges whose stored endpoints happen to
        // run in opposite directions would report both as "less than" each other — a comparator
        // that is not antisymmetric, which `sort_by` is entitled to turn into an arbitrary order.
        let (p, q) = if e1.a.x <= e1.b.x { (e1.a, e1.b) } else { (e1.b, e1.a) };
        // `e2`'s intercept at `x` lies strictly above `e1` exactly when it is left of `p -> q`.
        let o = orient2d(p, q, Point2D::new(x, e2.y_at(x)));
        match o {
            Orientation::Clockwise => core::cmp::Ordering::Greater,
            Orientation::CounterClockwise => core::cmp::Ordering::Less,
            Orientation::Collinear => {
                // Collinear at this x. Break by the obstacles' identity so overlapping obstacles
                // keep a fixed relative order, then by edge index.
                (e1.obstacle, e1.index).cmp(&(e2.obstacle, e2.index))
            }
        }
    });
}

/// A probe point strictly inside each gap between consecutive active edges, plus one below the
/// lowest and one above the highest, each tagged with the gap index it belongs to.
///
/// The lowest active edge is always the workspace's bottom wall, so the below-gap probe is
/// outside the workspace and the above-gap probe is outside its top wall; both are classified by
/// the same free-space test and rejected. Carrying the gap index avoids recovering it by
/// re-scanning the active list, which cannot then disagree with the list it was built from.
fn gap_probes(active: &[Edge], x: f64) -> Vec<(usize, Point2D)> {
    let n = active.len();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let below = i.checked_sub(1).map_or(f64::NEG_INFINITY, |j| active[j].y_at(x));
        let above = active.get(i).map_or(f64::INFINITY, |e| e.y_at(x));
        let y = match (below.is_finite(), above.is_finite()) {
            (true, true) => 0.5 * (below + above),
            // Below the lowest edge: one unit under it. Above the highest: one unit over it.
            (false, true) => above - 1.0,
            (true, false) => below + 1.0,
            (false, false) => 0.0,
        };
        out.push((i, Point2D::new(x, y)));
    }
    out
}

/// The workspace eroded by `margin` and the forbidden geometry inset by `margin`.
///
/// Exact for rectilinear obstacles. See [`decompose`] for why, and for the rejection of a
/// non-rectilinear obstacle under a nonzero `margin`.
fn apply_margin(space: &FreeSpace, margin: f64) -> Result<(BoundingBox, Vec<Polygon>), PathPlanError> {
    let workspace = space.workspace.eroded(margin);
    let mut forbidden: Vec<Polygon> = Vec::with_capacity(space.obstacles.len());
    for obstacle in &space.obstacles {
        if margin == 0.0 {
            forbidden.push(obstacle.clone());
            continue;
        }
        let parts = rectilinear_rectangles(obstacle.vertices()).ok_or(
            PathPlanError::MarginUnsupportedGeometry { margin },
        )?;
        for part in parts {
            let grown = part.outset(margin);
            forbidden.push(Polygon::new(grown.corners().to_vec()).map_err(|_| {
                PathPlanError::MarginUnsupportedGeometry { margin }
            })?);
        }
    }
    Ok((workspace, forbidden))
}

/// Whether `p` is outside every forbidden polygon.
fn is_free(forbidden: &[Polygon], p: Point2D) -> bool {
    !forbidden.iter().any(|o| o.contains(p))
}

/// The cell for gap `i` of `active`: the region between `active[i-1]` and `active[i]` over
/// `[lo_x, hi_x]`. The gap above the last edge has a virtual lower boundary at `y = 0` and the
/// gap below the first has a virtual upper one; both are classified as outside the workspace and
/// rejected before reaching here, so the virtual values never surface.
fn cell_from_gap(active: &[Edge], i: usize, lo_x: f64, hi_x: f64) -> Cell {
    let lo_edge = i.checked_sub(1).and_then(|j| active.get(j)).copied();
    let hi_edge = active.get(i).copied();
    Cell {
        bl: Point2D::new(lo_x, lo_edge.map_or(0.0, |e| e.y_at(lo_x))),
        br: Point2D::new(hi_x, lo_edge.map_or(0.0, |e| e.y_at(hi_x))),
        tr: Point2D::new(hi_x, hi_edge.map_or(0.0, |e| e.y_at(hi_x))),
        tl: Point2D::new(lo_x, hi_edge.map_or(0.0, |e| e.y_at(lo_x))),
        left: None,
        right: None,
    }
}

/// Joins the cells of two adjacent slabs across the event line at `x`, creating a portal for every
/// pair with a positive-length vertical overlap. `attach_right` records the portal on the left
/// cell's right side and the right cell's left side.
fn connect_slabs(
    builder: &mut CellBuilder,
    left: &[CellId],
    right: &[CellId],
    x: f64,
    attach_right: bool,
) {
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() && j < right.len() {
        let l = builder.cells[left[i] as usize].y_range_at(x);
        let r = builder.cells[right[j] as usize].y_range_at(x);
        let lo = l.0.max(r.0);
        let hi = l.1.min(r.1);
        if hi > lo {
            let id = builder.push_portal(Portal { x, lo, hi, left: left[i], right: right[j] });
            if attach_right {
                builder.cells[left[i] as usize].right = Some(id);
                builder.cells[right[j] as usize].left = Some(id);
            }
        }
        // Advance whichever cell ends first; equal ends advance both.
        if l.1 < r.1 {
            i += 1;
        } else if r.1 < l.1 {
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::FreeSpace;
    use alloc::vec;
    use alloc::vec::Vec;

    fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
        Polygon::new(vec![
            Point2D::new(x0, y0),
            Point2D::new(x1, y0),
            Point2D::new(x1, y1),
            Point2D::new(x0, y1),
        ])
        .unwrap()
    }

    fn space(obstacles: Vec<Polygon>) -> FreeSpace {
        FreeSpace::new(
            BoundingBox::new(Point2D::ZERO, Point2D::new(100.0, 100.0)),
            obstacles,
        )
        .unwrap()
    }

    #[test]
    fn a_margin_erodes_the_workspace_walls() {
        let d = decompose(&space(vec![]), 5.0).unwrap();
        assert_eq!(d.locate(Point2D::new(1.0, 50.0)), None, "too close to the west wall");
        assert_eq!(d.locate(Point2D::new(99.0, 50.0)), None, "too close to the east wall");
        assert!(d.locate(Point2D::new(50.0, 50.0)).is_some());
    }

    #[test]
    fn a_margin_erodes_away_from_obstacles() {
        let d = decompose(&space(vec![square(40.0, 40.0, 60.0, 60.0)]), 5.0).unwrap();
        assert!(d.locate(Point2D::new(50.0, 38.0)).is_none(), "within margin of the box");
        assert!(d.locate(Point2D::new(50.0, 33.0)).is_some(), "clear of the box by the margin");
    }

    #[test]
    fn a_gap_narrower_than_twice_the_margin_closes() {
        // A full-height vertical barrier with a 10-unit gap in it.
        let barrier = || {
            space(vec![square(40.0, 0.0, 44.0, 45.0), square(40.0, 55.0, 44.0, 100.0)])
        };
        let open = decompose(&barrier(), 2.0).unwrap();
        assert!(open.locate(Point2D::new(42.0, 50.0)).is_some(), "a 6-unit passage should survive");
        let shut = decompose(&barrier(), 6.0).unwrap();
        assert!(shut.locate(Point2D::new(42.0, 50.0)).is_none(), "a 10-unit gap cannot pass");
        // The right-hand region still exists as free space; it is just unreachable.
        assert!(shut.locate(Point2D::new(80.0, 50.0)).is_some());
    }

    #[test]
    fn a_margin_larger_than_the_workspace_yields_no_cells() {
        let d = decompose(&space(vec![]), 60.0).unwrap();
        assert_eq!(d.cell_count(), 0);
        assert!(d.locate(Point2D::new(50.0, 50.0)).is_none());
    }

    #[test]
    fn a_margin_on_a_sloped_obstacle_is_rejected_not_approximated() {
        let triangle = Polygon::new(vec![
            Point2D::new(10.0, 10.0),
            Point2D::new(40.0, 10.0),
            Point2D::new(10.0, 40.0),
        ])
        .unwrap();
        // Exact at zero margin, refused above it.
        assert!(decompose(&space(vec![triangle.clone()]), 0.0).is_ok());
        assert_eq!(
            decompose(&space(vec![triangle]), 1.0).unwrap_err(),
            PathPlanError::MarginUnsupportedGeometry { margin: 1.0 }
        );
    }

    #[test]
    fn a_margin_on_an_l_shaped_obstacle_is_exact() {
        // An L is rectilinear, so a margin against it is supported, not refused.
        let l = Polygon::new(vec![
            Point2D::new(10.0, 10.0),
            Point2D::new(60.0, 10.0),
            Point2D::new(60.0, 30.0),
            Point2D::new(30.0, 30.0),
            Point2D::new(30.0, 60.0),
            Point2D::new(10.0, 60.0),
        ])
        .unwrap();
        let d = decompose(&space(vec![l]), 2.0).unwrap();
        // The arm's right face is x = 30 and the base's top face is y = 30, so a 2-unit margin
        // pushes the notch's free region out to x > 32 and y > 32.
        assert!(d.locate(Point2D::new(35.0, 35.0)).is_some(), "in the notch, clear of both faces");
        assert!(d.locate(Point2D::new(31.0, 35.0)).is_none(), "1 unit from the arm");
        assert!(d.locate(Point2D::new(35.0, 31.0)).is_none(), "1 unit from the base");
        assert!(d.locate(Point2D::new(33.0, 35.0)).is_some(), "3 units from the arm");
        // The base's interior is still solid, and its top face keeps its standoff.
        assert!(d.locate(Point2D::new(50.0, 20.0)).is_none(), "inside the base");
        assert!(d.locate(Point2D::new(50.0, 31.0)).is_none(), "1 unit above the base");
    }

    #[test]
    fn an_empty_domain_is_a_single_strip() {
        let d = decompose(&space(vec![]), 0.0).unwrap();
        assert!(d.cell_count() >= 1);
        assert!(d.locate(Point2D::new(50.0, 50.0)).is_some());
    }

    #[test]
    fn a_box_splits_free_space() {
        let d = decompose(&space(vec![square(40.0, 40.0, 60.0, 60.0)]), 0.0).unwrap();
        assert_eq!(d.locate(Point2D::new(50.0, 50.0)), None);
        assert!(d.locate(Point2D::new(50.0, 20.0)).is_some());
        assert!(d.locate(Point2D::new(20.0, 50.0)).is_some());
        assert!(d.locate(Point2D::new(80.0, 50.0)).is_some());
    }



    /// A full-height vertical barrier with a gap in it. Returns the decomposition and whether a
    /// cell in the gap exists, and whether the left half reaches the right half by walking
    /// right-portals.

    #[test]
    fn cells_are_convex_and_non_negative_in_area() {
        let d = decompose(&space(vec![square(40.0, 40.0, 60.0, 60.0), square(10.0, 10.0, 20.0, 20.0)]), 3.0).unwrap();
        for cell in d.cells() {
            let c = cell.corners();
            for i in 0..4 {
                let o = orient2d(c[i], c[(i + 1) % 4], c[(i + 2) % 4]);
                assert!(
                    o == Orientation::CounterClockwise || o == Orientation::Collinear,
                    "cell {c:?} is not convex"
                );
            }
            assert!(cell.area() >= 0.0);
        }
    }

    #[test]
    fn events_are_sorted_and_distinct() {
        let d = decompose(&space(vec![square(40.0, 40.0, 60.0, 60.0), square(10.0, 10.0, 20.0, 20.0)]), 0.0).unwrap();
        for w in d.xs.windows(2) {
            assert!(w[0] < w[1], "event xs not strictly increasing: {w:?}");
        }
    }

    #[test]
    fn sort_f64_is_a_total_order() {
        let mut v = vec![1.0, -0.0, 0.0, f64::NAN, -1.0, 1.0, 2.0, -2.0];
        sort_f64_ascending(&mut v);
        assert_eq!(v[0], -2.0);
        assert_eq!(v[1], -1.0);
        assert_eq!(v[2].to_bits(), (-0.0f64).to_bits());
        assert_eq!(v[3].to_bits(), 0.0f64.to_bits());
        assert_eq!(v[4], 1.0);
        assert_eq!(v[5], 1.0);
        assert_eq!(v[6], 2.0);
        assert!(v[7].is_nan());
    }

    #[test]
    fn sort_f64_keys_are_monotonic_in_the_numeric_value() {
        // The specific failure this guards: the IEEE flip orders *unsigned*, so comparing the
        // keys as signed integers would place every non-negative value before every negative one.
        let ascending = [-3.0f64, -0.5, 0.0, 0.5, 3.0];
        let keys: Vec<u64> = ascending.iter().map(|x| f64_sort_key(*x)).collect();
        assert!(keys.windows(2).all(|w| w[0] < w[1]), "{keys:?}");
        let mut probe = vec![3.0, -0.5, 0.5, 0.0, -3.0];
        sort_f64_ascending(&mut probe);
        assert_eq!(probe, ascending.to_vec());
    }

    #[test]
    fn degenerate_workspace_is_rejected() {
        let bad = FreeSpace {
            workspace: BoundingBox::new(Point2D::ZERO, Point2D::ZERO),
            obstacles: vec![],
        };
        assert_eq!(decompose(&bad, 0.0).unwrap_err(), PathPlanError::DegenerateWorkspace);
    }

    #[test]
    fn a_negative_or_nan_margin_is_rejected() {
        let s = space(vec![]);
        assert!(matches!(decompose(&s, -1.0), Err(PathPlanError::InvalidConfig(_))));
        assert!(matches!(decompose(&s, f64::NAN), Err(PathPlanError::InvalidConfig(_))));
    }
}
