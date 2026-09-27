//! Finding a corridor and pulling a taut string through it.

pub mod cell_search;
pub mod string_pull;

use alloc::vec::Vec;

use crate::decomp::{CellId, PortalId};

/// A sequence of cells from the start cell to the goal cell, and the portals joining them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Corridor {
    /// The cells, in order from the start to the goal. Never empty.
    pub cells: Vec<CellId>,
    /// The portal between `cells[i]` and `cells[i + 1]`, so `portals.len() == cells.len() - 1`.
    pub portals: Vec<PortalId>,
}

impl Corridor {
    /// Number of cells in the corridor.
    #[inline]
    pub fn len(&self) -> usize {
        self.cells.len()
    }

    /// Whether the corridor holds no cells. A successful search never produces one.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Whether the corridor is a single cell, i.e. start and goal are in the same cell.
    #[inline]
    pub fn is_trivial(&self) -> bool {
        self.cells.len() <= 1
    }
}
