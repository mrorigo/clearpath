//! A* over the cell adjacency graph.

use alloc::boxed::Box;
use alloc::collections::BinaryHeap;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

use super::Corridor;
use crate::decomp::{CellId, Decomposition, PortalId};
use crate::error::PathPlanError;
use crate::geom::point::Point2D;

/// An open-list entry: `f`, then insertion order, then the cell.
///
/// The `f` comparison uses `f64::total_cmp`, and the tie-break is the insertion counter, so the
/// order the search visits cells in is a deterministic function of the inputs (section 4.4, D2).
/// Nothing here depends on a pointer or a hash iteration order.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Open {
    f: f64,
    counter: u32,
    cell: CellId,
}

impl Eq for Open {}

impl Ord for Open {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed: `BinaryHeap` is a max-heap and A* wants the smallest `f` first.
        other
            .f
            .total_cmp(&self.f)
            .then_with(|| other.counter.cmp(&self.counter))
    }
}

impl PartialOrd for Open {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Reusable scratch for the search, so a `PathPlanner` can avoid reallocating per query
/// (section 8.1).
#[derive(Clone, Debug, Default)]
pub struct SearchScratch {
    g_score: Vec<f64>,
    came_from: Vec<u32>,
    closed: Vec<bool>,
    open: BinaryHeap<Open>,
    counter: u32,
}

impl SearchScratch {
    /// Prepares the scratch for a decomposition with `cell_count` cells. The `Vec`s are cleared
    /// and resized rather than reallocated, so a warm query allocates nothing.
    pub(crate) fn reset(&mut self, cell_count: usize) {
        self.g_score.clear();
        self.g_score.resize(cell_count, f64::INFINITY);
        self.came_from.clear();
        self.came_from.resize(cell_count, u32::MAX);
        self.closed.clear();
        self.closed.resize(cell_count, false);
        self.open.clear();
        self.counter = 0;
    }
}

/// Finds a cell corridor from the cell containing `start` to the cell containing `goal`.
///
/// The search is restricted to moves in the direction of travel, so the corridor is **monotone in
/// x**. The cell graph fans out — a cell's right side overlaps several cells of the next slab — so
/// an unrestricted search can return a corridor that steps right and then back left, and the funnel
/// is only proven for a sleeve whose portals are ordered along the path. Monotonicity removes that
/// case by construction rather than by hoping the tie-breaks avoid it.
///
/// When the two points are in the same slab no monotone corridor exists, so the search is retried
/// without the restriction. That case is only reachable when the two abscissae are *equal* — the
/// caller is expected to have forced them into the sweep's event set otherwise, which puts them in
/// different slabs by construction. It is still covered, by the containment check in
/// [`string_pull`](super::string_pull::string_pull), which refuses to return a path that leaves the
/// corridor rather than returning a wrong one.
///
/// Returns [`PathPlanError::NoPathFound`] if either point is in no cell, or if the two cells are
/// in different connected components. Connectivity is not precomputed; it falls out of the
/// exhausted open list.
///
/// The cost is `$O(E \log V)$` over the cells and portals, where `$E$` counts portals — distinct
/// from the funnel's `$O(N)` in section 6.3, and tested separately.
pub fn search(
    decomp: &Decomposition,
    start: Point2D,
    goal: Point2D,
    scratch: &mut SearchScratch,
) -> Result<Corridor, PathPlanError> {
    let start_cell = decomp.locate(start).ok_or(PathPlanError::NoPathFound)?;
    let goal_cell = decomp.locate(goal).ok_or(PathPlanError::NoPathFound)?;

    if start_cell == goal_cell {
        return Ok(Corridor { cells: vec![start_cell], portals: Vec::new() });
    }
    if let Some(found) = run(decomp, start_cell, goal_cell, goal, true, scratch) {
        return Ok(found);
    }
    run(decomp, start_cell, goal_cell, goal, false, scratch).ok_or(PathPlanError::NoPathFound)
}

/// The search itself. `monotone` restricts the expansion to the direction of travel.
fn run(
    decomp: &Decomposition,
    start_cell: CellId,
    goal_cell: CellId,
    goal_point: Point2D,
    monotone: bool,
    scratch: &mut SearchScratch,
) -> Option<Corridor> {
    scratch.reset(decomp.cell_count());
    // `+1` if the goal lies to the right of the whole start cell, `-1` otherwise. Comparing against
    // the cell rather than the start *point* matters: a start point near the cell's right edge with
    // the goal just left of it would otherwise be sent the wrong way.
    let travel = if goal_point.x > decomp.cells()[start_cell as usize].br.x { 1.0 } else { -1.0 };
    scratch.g_score[start_cell as usize] = 0.0;
    let start_f = decomp.cells()[start_cell as usize].centroid().distance(goal_point);
    let counter = scratch.next_counter();
    scratch.open.push(Open { f: start_f, counter, cell: start_cell });

    while let Some(Open { cell, .. }) = scratch.open.pop() {
        if scratch.closed[cell as usize] {
            continue;
        }
        scratch.closed[cell as usize] = true;
        if cell == goal_cell {
            return Some(reconstruct(decomp, scratch, start_cell, goal_cell));
        }
        // The edge cost is the portal length: the only way out of a cell is across a portal, and a
        // longer portal is a wider gate that is more expensive to have routed through.
        let expand: Box<dyn Iterator<Item = (PortalId, CellId)>> = if !monotone {
            Box::new(decomp.right_neighbours(cell).chain(decomp.left_neighbours(cell)))
        } else if travel > 0.0 {
            Box::new(decomp.right_neighbours(cell))
        } else {
            Box::new(decomp.left_neighbours(cell))
        };
        for (portal, next) in expand {
            if scratch.closed[next as usize] {
                continue;
            }
            let step = decomp.portals()[portal as usize].length();
            let tentative = scratch.g_score[cell as usize] + step;
            if tentative < scratch.g_score[next as usize] {
                scratch.g_score[next as usize] = tentative;
                scratch.came_from[next as usize] = cell;
                // Admissible: no route can beat the straight-line distance from this centroid to
                // the goal, so the true cost is never less than the remaining distance.
                let h = decomp.cells()[next as usize].centroid().distance(goal_point);
                let counter = scratch.next_counter();
                scratch.open.push(Open { f: tentative + h, counter, cell: next });
            }
        }
    }
    None
}

impl SearchScratch {
    /// A strictly increasing counter, so equal-`f` entries break ties by insertion order.
    fn next_counter(&mut self) -> u32 {
        self.counter += 1;
        self.counter
    }
}

/// Walks `came_from` back from the goal, emitting the cell and portal sequence in forward order.
fn reconstruct(
    decomp: &Decomposition,
    scratch: &SearchScratch,
    start: CellId,
    goal: CellId,
) -> Corridor {
    let mut cells: Vec<CellId> = vec![goal];
    let mut current = goal;
    while current != start {
        current = scratch.came_from[current as usize] as CellId;
        cells.push(current);
    }
    cells.reverse();
    // Recover each portal from the cell's own adjacency rather than storing it during the search,
    // so the scratch holds no per-cell portal list.
    let mut portals: Vec<PortalId> = Vec::with_capacity(cells.len().saturating_sub(1));
    for pair in cells.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let found = decomp
            .right_portals(a)
            .iter()
            .chain(decomp.left_portals(a).iter())
            .find(|p| {
                let portal = decomp.portals()[**p as usize];
                (portal.left == b && portal.right == a) || (portal.right == b && portal.left == a)
            })
            .copied()
            .expect("came_from is only ever set across a portal");
        portals.push(found);
    }
    Corridor { cells, portals }
}
