//! The vertical (trapezoidal) decomposition of the offset free space.
//!
//! A sweep in `+x` partitions the free space into convex cells, one per gap between consecutive
//! active boundary edges in each slab. See `docs/SPEC.md` section 6.2 and lemma L1.

pub mod cell;
pub mod sweep;

use alloc::vec::Vec;

use crate::geom::Point2D;

pub use cell::{Cell, CellId, Portal, PortalId};

/// A partition of the free space into convex cells.
///
/// Invariant L1: the cells are convex, have pairwise disjoint interiors, and their union is the
/// free space. `tests::decomposition_is_a_valid_partition` checks this exhaustively on small
/// inputs against a direct free-space predicate.
#[derive(Clone, Debug, Default)]
pub struct Decomposition {
    /// Sorted, distinct event x coordinates. There are `xs.len() - 1` slabs.
    pub(crate) xs: Vec<f64>,
    /// Cells grouped by slab; `slab_cells[k]` holds the ids of the cells in slab `k`, ordered
    /// bottom to top.
    pub(crate) slab_cells: Vec<Vec<CellId>>,
    pub(crate) cells: Vec<Cell>,
    pub(crate) portals: Vec<Portal>,
}

impl Decomposition {
    /// The cells.
    #[inline]
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// The portals.
    #[inline]
    pub fn portals(&self) -> &[Portal] {
        &self.portals
    }

    /// Number of cells.
    #[inline]
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Number of portals.
    #[inline]
    pub fn portal_count(&self) -> usize {
        self.portals.len()
    }

    /// The cell to the right of `cell` across its right portal, if it has one.
    #[inline]
    pub fn right_neighbour(&self, cell: CellId) -> Option<(PortalId, CellId)> {
        let p = self.cells[cell as usize].right?;
        Some((p, self.portals[p as usize].right))
    }

    /// The cell to the left of `cell` across its left portal, if it has one.
    #[inline]
    pub fn left_neighbour(&self, cell: CellId) -> Option<(PortalId, CellId)> {
        let p = self.cells[cell as usize].left?;
        Some((p, self.portals[p as usize].left))
    }

    /// The cell containing `p`, or `None` if `p` is outside the free space.
    ///
    /// `O(log n)`: binary search for the slab, then a linear scan of that slab's cells, which are
    /// sorted and disjoint in y.
    pub fn locate(&self, p: Point2D) -> Option<CellId> {
        if !self.xs.is_empty() {
            let k = slab_index(&self.xs, p.x)?;
            for &id in &self.slab_cells[k] {
                let cell = &self.cells[id as usize];
                let (lo, hi) = cell.y_range_at(p.x);
                if p.y >= lo && p.y <= hi && cell.contains(p) {
                    return Some(id);
                }
            }
        }
        None
    }
}

/// The index of the slab whose closed x range contains `x`, or `None` if `x` is outside the
/// decomposition entirely.
pub(crate) fn slab_index(xs: &[f64], x: f64) -> Option<usize> {
    if xs.len() < 2 || x < xs[0] || x > xs[xs.len() - 1] {
        return None;
    }
    // The largest k with xs[k] <= x, then clamp: x on the last event line belongs to the last
    // slab; x elsewhere belongs to the slab that starts at or before it.
    let mut lo = 0usize;
    let mut hi = xs.len() - 1;
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if xs[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(lo.min(xs.len() - 2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Point2D;
    use alloc::vec;

    #[test]
    fn slab_index_finds_the_right_slab() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        assert_eq!(slab_index(&xs, 0.0), Some(0));
        assert_eq!(slab_index(&xs, 0.5), Some(0));
        assert_eq!(slab_index(&xs, 1.0), Some(1));
        assert_eq!(slab_index(&xs, 1.5), Some(1));
        assert_eq!(slab_index(&xs, 2.5), Some(2));
        assert_eq!(slab_index(&xs, 3.0), Some(2));
        assert_eq!(slab_index(&xs, 3.5), None);
        assert_eq!(slab_index(&xs, -0.5), None);
    }

    #[test]
    fn slab_index_of_a_degenerate_decomposition() {
        assert_eq!(slab_index(&[], 0.0), None);
        assert_eq!(slab_index(&[0.0], 0.0), None);
    }

    #[test]
    fn locate_on_an_empty_decomposition_is_none() {
        let d = Decomposition::default();
        assert_eq!(d.locate(Point2D::new(1.0, 1.0)), None);
        assert_eq!(d.cell_count(), 0);
        assert_eq!(d.portal_count(), 0);
        assert!(d.cells().is_empty());
    }
}
