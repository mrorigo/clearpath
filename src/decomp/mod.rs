//! The vertical (trapezoidal) decomposition of the offset free space.
//!
//! A sweep in `+x` partitions the free space into convex cells, one per gap between consecutive
//! active boundary edges in each slab. See `docs/SPEC.md` section 6.2 and lemma L1.

pub mod cell;
pub mod sweep;

use alloc::vec::Vec;

use crate::geom::Point2D;

use cell::SideIndex;

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
    /// Cells grouped by slab in one flat buffer, ordered bottom to top within a slab.
    pub(crate) slab_cells: Vec<CellId>,
    /// Where each slab's run in `slab_cells` ends; slab `k` is `[offsets[k-1], offsets[k])`, and
    /// `offsets` has one more entry than there are slabs.
    pub(crate) slab_cell_offsets: Vec<u32>,
    pub(crate) cells: Vec<Cell>,
    pub(crate) portals: Vec<Portal>,
    /// Portals to the right of each cell, as a compressed list.
    pub(crate) right_index: SideIndex,
    /// Portals to the left of each cell, as a compressed list.
    pub(crate) left_index: SideIndex,
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

    /// The ids of the cells in slab `k`, as a range into [`Decomposition::slab_cells`].
    ///
    /// Flat rather than a `Vec` per slab: a decomposition has one slab per distinct event
    /// abscissa, so the nested form was one allocation per slab on every query.
    pub fn slab_run(&self, k: usize) -> (usize, usize) {
        let lo = if k == 0 { 0 } else { self.slab_cell_offsets[k - 1] as usize };
        let hi = if k < self.slab_cell_offsets.len() {
            self.slab_cell_offsets[k] as usize
        } else {
            self.slab_cells.len()
        };
        (lo, hi)
    }

    /// The portals on the right side of `cell`.
    ///
    /// Usually zero or one, but can be several: a wide cell's right side can overlap several
    /// narrower cells in the next slab.
    #[inline]
    pub fn right_portals(&self, cell: CellId) -> &[PortalId] {
        self.right_index.get(cell)
    }

    /// The portals on the left side of `cell`.
    #[inline]
    pub fn left_portals(&self, cell: CellId) -> &[PortalId] {
        self.left_index.get(cell)
    }

    /// The portals joining `cell` to the cell on the other side of each.
    ///
    /// A convenience for the common case; a caller that needs to walk the whole fan uses
    /// [`Decomposition::right_portals`].
    pub fn right_neighbours(&self, cell: CellId) -> impl Iterator<Item = (PortalId, CellId)> + '_ {
        self.right_portals(cell)
            .iter()
            .map(move |p| (*p, self.portals[*p as usize].right))
    }

    /// The cells on the other side of each of `cell`'s left portals.
    pub fn left_neighbours(&self, cell: CellId) -> impl Iterator<Item = (PortalId, CellId)> + '_ {
        self.left_portals(cell)
            .iter()
            .map(move |p| (*p, self.portals[*p as usize].left))
    }

    /// The cell containing `p`, or `None` if `p` is outside the free space.
    ///
    /// `O(log n)`: binary search for the slab, then a linear scan of that slab's cells, which are
    /// sorted and disjoint in y.
    pub fn locate(&self, p: Point2D) -> Option<CellId> {
        if !self.xs.is_empty() {
            let k = slab_index(&self.xs, p.x)?;
            let (lo, hi) = self.slab_run(k);
            for &id in &self.slab_cells[lo..hi] {
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
