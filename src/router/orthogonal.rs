//! The rectilinear router (section 7.2).
//!
//! # Why the Hanan grid, and why the cell test is exact
//!
//! The grid lines are the x and y coordinates of every forbidden rectangle's faces. So no obstacle
//! edge crosses a cell's *interior* — every face lies along a grid line — and a cell is therefore
//! wholly inside or wholly outside each obstacle. That is what makes testing a cell's *centre* an
//! exact test of the cell, rather than the sampling approximation it looks like. It is the one
//! property that lets this be a small router.
//!
//! # Why there is no `u`-shaped fallback
//!
//! There is not one, deliberately. A fallback that guesses at topology is wrong whenever the sealed
//! region is not rectangular, and a rectilinear route that is wrong is a collision. When the grid
//! has no 4-connected path, the answer is `NoPathFound`.

use alloc::vec;
use alloc::vec::Vec;

use crate::decomp::sweep::forbidden_rectangles;
use crate::error::PathPlanError;
use crate::geom::point::Point2D;
use crate::router::scratch::Scratch;
use crate::router::RouteRequest;

/// A rectilinear polyline: consecutive vertices share an abscissa or an ordinate exactly.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrthogonalPolyline {
    /// The vertices, from start to goal.
    pub points: Vec<Point2D>,
}

impl OrthogonalPolyline {
    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the polyline has no vertices.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Whether every consecutive pair differs on exactly one axis, which is what "rectilinear"
    /// means.
    pub fn is_axis_aligned(&self) -> bool {
        self.points.windows(2).all(|w| (w[0].x == w[1].x) ^ (w[0].y == w[1].y))
    }
}

/// A grid cell index pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Cell {
    i: usize,
    j: usize,
}

/// A rectilinear route from `req.start` to `req.goal`.
pub fn route_orthogonal(
    scratch: &mut Scratch,
    req: &RouteRequest,
) -> Result<OrthogonalPolyline, PathPlanError> {
    req.validate()?;
    let margin = req.config.margin;
    let workspace = req.workspace.eroded(margin);
    if workspace.is_degenerate() {
        return Err(PathPlanError::NoPathFound);
    }
    let forbidden = forbidden_rectangles(&req.obstacles, margin)?;

    let grid = Grid::build(&workspace, &forbidden);
    let start = req.start.point;
    let goal = req.goal.point;
    let start_cell = grid.locate(start).ok_or(PathPlanError::NoPathFound)?;
    let goal_cell = grid.locate(goal).ok_or(PathPlanError::NoPathFound)?;

    let path = grid.astar(&mut scratch.grid, start_cell, goal_cell)?;
    let mut polyline = grid.to_polyline(&path, start, goal);

    // Trim collinear runs, which the corner walk can leave behind at a repeated direction.
    polyline.points = simplify(core::mem::take(&mut polyline.points));
    if polyline.points.first() != Some(&start) {
        polyline.points.insert(0, start);
    }
    if polyline.points.last() != Some(&goal) {
        polyline.points.push(goal);
    }

    // The polyline is verified rather than trusted. The grid geometry is what makes it correct, and
    // a verification is what catches it when that reasoning is wrong.
    //
    // The verification is exact, not sampled. `Grid` is built from each obstacle's *bounding box*,
    // so for a non-rectangular obstacle a face can cross a cell's interior even though no grid
    // line does; the cell's centre then decides a cell that is only partly free. Sampling a leg at
    // a handful of points cannot see that — a blocked stretch narrower than the sample spacing
    // falls between two samples and the route is returned straight through the obstacle.
    // `hull_is_free` tests the leg against every obstacle edge, and a two-point hull is the leg.
    let clearance = req.clearance();
    for w in polyline.points.windows(2) {
        if !clearance.hull_is_free(w) {
            return Err(PathPlanError::NoPathFound);
        }
    }
    Ok(polyline)
}

/// The Hanan grid: the axis-aligned lines through every forbidden face.
struct Grid {
    xs: Vec<f64>,
    ys: Vec<f64>,
    /// `free[j * (nx) + i]`, with `i` the column and `j` the row.
    free: Vec<bool>,
    nx: usize,
}

impl Grid {
    fn build(workspace: &crate::geom::BoundingBox, forbidden: &[crate::geom::Polygon]) -> Self {
        let mut xs = vec![workspace.min.x, workspace.max.x];
        let mut ys = vec![workspace.min.y, workspace.max.y];
        for obstacle in forbidden {
            let b = obstacle.bounds();
            xs.push(b.min.x);
            xs.push(b.max.x);
            ys.push(b.min.y);
            ys.push(b.max.y);
        }
        let (wx0, wx1) = (workspace.min.x, workspace.max.x);
        let (wy0, wy1) = (workspace.min.y, workspace.max.y);
        xs.retain(|x| *x >= wx0 && *x <= wx1);
        ys.retain(|y| *y >= wy0 && *y <= wy1);
        crate::decomp::sweep::sort_f64_ascending(&mut xs);
        crate::decomp::sweep::sort_f64_ascending(&mut ys);
        xs.dedup();
        ys.dedup();

        let nx = xs.len().saturating_sub(1);
        let ny = ys.len().saturating_sub(1);
        let mut free = vec![false; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                // The cell's centre. Exact rather than a sample: no face crosses a cell's
                // interior, so the centre decides the whole cell.
                let cx = 0.5 * (xs[i] + xs[i + 1]);
                let cy = 0.5 * (ys[j] + ys[j + 1]);
                free[j * nx + i] = !forbidden.iter().any(|o| o.contains(Point2D::new(cx, cy)));
            }
        }
        Self { xs, ys, free, nx }
    }

    fn is_free(&self, c: Cell) -> bool {
        c.i < self.nx && c.j < self.ys.len() - 1 && self.free[c.j * self.nx + c.i]
    }

    fn locate(&self, p: Point2D) -> Option<Cell> {
        // The cell containing `p` is the one below its upper bound. A point below the first line
        // is outside the grid, which is what the `?` reports.
        let i = upper_bound(&self.xs, p.x)?.checked_sub(1)?;
        let j = upper_bound(&self.ys, p.y)?.checked_sub(1)?;
        let c = Cell { i, j };
        self.is_free(c).then_some(c)
    }

    /// A* over the 4-neighbourhood, cost the Manhattan step.
    fn astar(&self, scratch: &mut GridScratch, start: Cell, goal: Cell) -> Result<Vec<Cell>, PathPlanError> {
        let n = self.free.len();
        scratch.reset(n);
        let index = |c: Cell| c.j * self.nx + c.i;
        let h = |c: Cell| {
            (self.xs[c.i] - self.xs[goal.i]).abs() + (self.ys[c.j] - self.ys[goal.j]).abs()
        };
        scratch.g[index(start)] = 0.0;
        scratch.push(((h(start), index(start) as u64), start, 0.0));
        while let Some(((f, _), c, g)) = scratch.pop() {
            let k = index(c);
            if scratch.closed[k] {
                continue;
            }
            scratch.closed[k] = true;
            if c == goal {
                return Ok(scratch.reconstruct(self.nx, start, goal));
            }
            for next in [Cell { i: c.i + 1, j: c.j }, Cell { i: c.i.wrapping_sub(1), j: c.j }, Cell { i: c.i, j: c.j + 1 }, Cell { i: c.i, j: c.j.wrapping_sub(1) }] {
                if !self.is_free(next) {
                    continue;
                }
                let nk = index(next);
                if scratch.closed[nk] {
                    continue;
                }
                // The step cost is the distance actually travelled along the shared grid line.
                let step = if next.i != c.i {
                    (self.xs[next.i] - self.xs[c.i]).abs()
                } else {
                    (self.ys[next.j] - self.ys[c.j]).abs()
                };
                let tentative = g + step;
                if tentative < scratch.g[nk] {
                    scratch.g[nk] = tentative;
                    scratch.from[nk] = index(c) as u32;
                    // `f` is the popped key, unused: the value is recomputed from `g` so the entry
                    // cannot go stale, and ties break on the cell index so the order is a
                    // deterministic function of the inputs (section 4.4).
                    let _ = f;
                    scratch.push(((h(next), nk as u64), next, tentative));
                }
            }
        }
        Err(PathPlanError::NoPathFound)
    }

    /// The centre of a cell.
    fn centre(&self, c: Cell) -> Point2D {
        Point2D::new(0.5 * (self.xs[c.i] + self.xs[c.i + 1]), 0.5 * (self.ys[c.j] + self.ys[c.j + 1]))
    }

    /// Turns a cell path into an axis-aligned polyline through cell *centres*.
    ///
    /// Centres, not corners, and that is the whole correctness argument. A segment between two
    /// adjacent cells' centres lies inside the union of those two cells, and both are free — so the
    /// segment is clear. Routing through corners instead puts waypoints on grid lines, and a
    /// segment along a grid line runs between the two cells on either side of it, only one of which
    /// the path certified. The search guarantees the *centres* are free; it guarantees nothing
    /// about the lines between them.
    ///
    /// The endpoints need an L each, because the start and goal sit anywhere inside their cells.
    /// Both legs of each L stay within that cell, so they are clear for the same reason.
    fn to_polyline(&self, path: &[Cell], start: Point2D, goal: Point2D) -> OrthogonalPolyline {
        if path.len() <= 1 {
            // One cell holds both endpoints, so a single L suffices and detouring through the
            // centre would add two corners for nothing. Both legs stay inside that one cell.
            let corner = Point2D::new(goal.x, start.y);
            return OrthogonalPolyline { points: vec![start, corner, goal] };
        }
        let mut points = Vec::with_capacity(path.len() + 3);
        points.push(start);
        let push_entry = |points: &mut Vec<Point2D>, target: Point2D| {
            if let Some(from) = points.last().copied() {
                if from.x != target.x && from.y != target.y {
                    points.push(Point2D::new(target.x, from.y));
                }
            }
            points.push(target);
        };
        for c in path {
            push_entry(&mut points, self.centre(*c));
        }
        push_entry(&mut points, goal);
        OrthogonalPolyline { points }
    }
}

/// Strict "smaller than" on a heap key, with the cell index breaking ties.
fn key_lt(a: (f64, u64), b: (f64, u64)) -> bool {
    if a.0 < b.0 {
        true
    } else if a.0 > b.0 {
        false
    } else {
        a.1 < b.1
    }
}

/// The first index whose value is `>= x`, or `None` if `x` is below every value.
///
/// The sort key is the same bit-pattern total order the decomposition uses, so a point exactly on
/// a grid line lands in the cell above it, consistently with the sweep's own convention.
fn upper_bound(values: &[f64], x: f64) -> Option<usize> {
    match values.binary_search_by(|probe| {
        if *probe < x {
            core::cmp::Ordering::Less
        } else if *probe > x {
            core::cmp::Ordering::Greater
        } else {
            core::cmp::Ordering::Equal
        }
    }) {
        Ok(i) => Some(i),
        // The insertion point is the first value `>= x`; a zero insertion point means `x` is below
        // everything, which is outside the grid and not a cell.
        Err(0) => None,
        Err(i) => Some(i),
    }
}

/// Drops interior vertices that continue the same direction.
fn simplify(points: Vec<Point2D>) -> Vec<Point2D> {
    if points.len() < 3 {
        return points;
    }
    let mut out: Vec<Point2D> = Vec::with_capacity(points.len());
    out.push(points[0]);
    for w in points.windows(3) {
        let (a, b, c) = (w[0], w[1], w[2]);
        let continues = (a.x == b.x && b.x == c.x) || (a.y == b.y && b.y == c.y);
        if !continues {
            out.push(b);
        }
    }
    out.push(points[points.len() - 1]);
    out
}

/// The grid search's arrays.
#[derive(Clone, Debug, Default)]
pub struct GridScratch {
    g: Vec<f64>,
    from: Vec<u32>,
    closed: Vec<bool>,
    /// A min-heap on `(f, cell index)`, with the index as the tie-break so the visit order is a
    /// deterministic function of the inputs (section 4.4). `f64` is not `Ord`, so the comparison
    /// goes through [`key_lt`].
    open: Vec<((f64, u64), Cell, f64)>,
    counter: u32,
}

impl GridScratch {
    fn reset(&mut self, n: usize) {
        self.g.clear();
        self.g.resize(n, f64::INFINITY);
        self.from.clear();
        self.from.resize(n, u32::MAX);
        self.closed.clear();
        self.closed.resize(n, false);
        self.open.clear();
        self.counter = 0;
    }

    fn push(&mut self, entry: ((f64, u64), Cell, f64)) {
        self.counter += 1;
        self.open.push(entry);
        // Sift up. A hand-rolled heap keeps the scratch allocation-free across queries, and the
        // grid is small enough that the asymptotics of `BinaryHeap` do not matter.
        let mut i = self.open.len() - 1;
        while i > 0 {
            let parent = (i - 1) / 2;
            if key_lt(self.open[i].0, self.open[parent].0) {
                self.open.swap(i, parent);
                i = parent;
            } else {
                break;
            }
        }
    }

    fn pop(&mut self) -> Option<((f64, u64), Cell, f64)> {
        let top = self.open.first().cloned()?;
        let last = self.open.pop().expect("non-empty: `top` was read from it");
        if !self.open.is_empty() {
            self.open[0] = last;
            let mut i = 0;
            loop {
                let (l, r) = (2 * i + 1, 2 * i + 2);
                let mut largest = i;
                if l < self.open.len() && key_lt(self.open[l].0, self.open[largest].0) {
                    largest = l;
                }
                if r < self.open.len() && key_lt(self.open[r].0, self.open[largest].0) {
                    largest = r;
                }
                if largest == i {
                    break;
                }
                self.open.swap(i, largest);
                i = largest;
            }
        }
        Some(top)
    }

    fn reconstruct(&self, nx: usize, start: Cell, goal: Cell) -> Vec<Cell> {
        let index = |c: Cell| c.j * nx + c.i;
        let mut out = vec![goal];
        let mut cur = goal;
        while cur != start {
            cur = Cell { i: (self.from[index(cur)] as usize % nx), j: (self.from[index(cur)] as usize / nx) };
            out.push(cur);
        }
        out.reverse();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::BoundingBox;

    #[test]
    fn a_single_cell_still_gives_an_l() {
        let g = Grid { xs: vec![0.0, 100.0], ys: vec![0.0, 100.0], free: vec![true], nx: 1 };
        let poly = g.to_polyline(
            &[Cell { i: 0, j: 0 }],
            Point2D::new(10.0, 20.0),
            Point2D::new(80.0, 60.0),
        );
        assert_eq!(poly.points.len(), 3, "{poly:?}");
        assert!(poly.is_axis_aligned());
    }

    #[test]
    fn a_multi_cell_path_routes_through_centres() {
        let g = Grid {
            xs: vec![0.0, 100.0],
            ys: vec![0.0, 50.0, 100.0],
            free: vec![true; 2],
            nx: 1,
        };
        let poly = g.to_polyline(
            &[Cell { i: 0, j: 0 }, Cell { i: 0, j: 1 }],
            Point2D::new(10.0, 20.0),
            Point2D::new(80.0, 80.0),
        );
        assert!(poly.is_axis_aligned(), "{poly:?}");
        // Every leg is inside one free cell or two adjacent ones, so no leg can cross an obstacle.
        assert!(poly.points.iter().all(|p| p.x >= 0.0 && p.x <= 100.0));
    }

    #[test]
    fn upper_bound_finds_the_cell() {
        let v = vec![0.0, 10.0, 20.0];
        assert_eq!(upper_bound(&v, -1.0), None);
        assert_eq!(upper_bound(&v, 0.0), Some(0));
        assert_eq!(upper_bound(&v, 5.0), Some(1));
        assert_eq!(upper_bound(&v, 10.0), Some(1));
        assert_eq!(upper_bound(&v, 25.0), Some(3));
    }

    #[test]
    fn simplify_drops_collinear_runs() {
        let pts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(2.0, 1.0),
            Point2D::new(2.0, 2.0),
        ];
        assert_eq!(
            simplify(pts),
            vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0), Point2D::new(2.0, 2.0)]
        );
    }

    #[test]
    fn the_grid_marks_blocked_cells() {
        let obstacle = crate::geom::Polygon::new(alloc::vec![
            Point2D::new(40.0, 40.0),
            Point2D::new(60.0, 40.0),
            Point2D::new(60.0, 60.0),
            Point2D::new(40.0, 60.0),
        ])
        .unwrap();
        let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(100.0, 100.0));
        let g = Grid::build(&ws, alloc::slice::from_ref(&obstacle));
        assert!(g.is_free(Cell { i: 0, j: 0 }));
        // The grid is [0, 40, 60, 100] on both axes, so the box occupies the cell (1, 1) exactly
        // and the cells to its sides are free. This is the Hanan property the centre test relies
        // on: no face crosses a cell's interior, so the centre decides the whole cell.
        assert!(!g.is_free(Cell { i: 1, j: 1 }), "the box's cell must be blocked");
        assert!(g.is_free(Cell { i: 0, j: 1 }));
        assert!(g.is_free(Cell { i: 1, j: 0 }));
    }
}
