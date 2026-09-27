//! Reusable scratch buffers, so a warm query allocates less (section 8.1).
//!
//! The decomposition is *not* yet reused: `decompose` builds a fresh one, and the struct is rebuilt
//! per query because the sweep's cell count depends on the obstacles. That is the largest
//! remaining allocation in a warm query, and moving it here needs `Decomposition` to be
//! fillable in place. Recorded as a known gap rather than pretended away.

use alloc::vec::Vec;

use crate::corridor::AdmissibleSet;
use crate::decomp::PortalId;
use crate::funnel::cell_search::SearchScratch;
use crate::geom::Point2D;

/// Everything a query allocates after the decomposition, kept alive between queries.
#[derive(Clone, Debug, Default)]
pub struct Scratch {
    /// The cell search's arrays: `g_score`, `came_from`, `closed` and the open list.
    pub search: SearchScratch,
    /// The corridor's cells.
    pub corridor_cells: Vec<u32>,
    /// The corridor's portals.
    pub corridor_portals: Vec<PortalId>,
    /// The taut path's knots.
    pub knots: Vec<Point2D>,
    /// The solved tangents.
    pub tangents: Vec<Point2D>,
    /// The admissible set at each knot.
    pub admissible: Vec<AdmissibleSet>,
}

impl Scratch {
    /// A scratch with empty buffers.
    pub fn new() -> Self {
        Self::default()
    }
}
