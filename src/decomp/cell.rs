//! Cells: the convex trapezoids of the vertical decomposition.

use alloc::vec::Vec;

use crate::geom::point::Point2D;

/// Index of a cell in a [`Decomposition`](super::Decomposition).
pub type CellId = u32;

/// Index of a portal in a [`Decomposition`](super::Decomposition).
pub type PortalId = u32;

/// A vertical portal: the shared boundary between two horizontally adjacent cells.
///
/// Portals are always vertical segments. That is what makes the cells convex and is what lets the
/// funnel (§6.3) treat the corridor as a portal sequence.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Portal {
    /// The x coordinate of the portal.
    pub x: f64,
    /// Lower endpoint's y.
    pub lo: f64,
    /// Upper endpoint's y.
    pub hi: f64,
    /// The cell to the left of the portal.
    pub left: CellId,
    /// The cell to the right of the portal.
    pub right: CellId,
}

impl Portal {
    /// The length of the portal along y.
    #[inline]
    pub fn length(&self) -> f64 {
        self.hi - self.lo
    }

    /// Lower endpoint.
    #[inline]
    pub fn bottom(&self) -> Point2D {
        Point2D::new(self.x, self.lo)
    }

    /// Upper endpoint.
    #[inline]
    pub fn top(&self) -> Point2D {
        Point2D::new(self.x, self.hi)
    }
}

/// A convex cell: the region between two straight boundary edges, over one slab of the sweep.
///
/// The four corners are counter-clockwise starting bottom-left. A cell may be a triangle, when the
/// two boundary edges meet inside the slab.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell {
    /// Bottom-left corner.
    pub bl: Point2D,
    /// Bottom-right corner.
    pub br: Point2D,
    /// Top-right corner.
    pub tr: Point2D,
    /// Top-left corner.
    pub tl: Point2D,
    /// The portal on the left edge, if any.
    pub left: Option<PortalId>,
    /// The portal on the right edge, if any.
    pub right: Option<PortalId>,
}

impl Cell {
    /// The centroid, used as the A* node cost (§6.3a).
    #[inline]
    pub fn centroid(&self) -> Point2D {
        Point2D::new(
            (self.bl.x + self.br.x + self.tr.x + self.tl.x) / 4.0,
            (self.bl.y + self.br.y + self.tr.y + self.tl.y) / 4.0,
        )
    }

    /// The x range spanned by the cell.
    #[inline]
    pub fn x_range(&self) -> (f64, f64) {
        (self.bl.x, self.br.x)
    }

    /// The y of the lower boundary at `x`, linearly interpolated along `bl -> br`.
    #[inline]
    pub fn bottom_at(&self, x: f64) -> f64 {
        interpolate(self.bl, self.br, x)
    }

    /// The y of the upper boundary at `x`, linearly interpolated along `tl -> tr`.
    #[inline]
    pub fn top_at(&self, x: f64) -> f64 {
        interpolate(self.tl, self.tr, x)
    }

    /// The vertical extent of the cell at `x`.
    #[inline]
    pub fn y_range_at(&self, x: f64) -> (f64, f64) {
        (self.bottom_at(x), self.top_at(x))
    }

    /// Whether `p` is in the closed cell. Uses exact orientation, never an epsilon.
    pub fn contains(&self, p: Point2D) -> bool {
        use crate::geom::predicates::{Orientation, orient2d};
        if p.x < self.bl.x || p.x > self.br.x {
            return false;
        }
        let (lo, hi) = self.y_range_at(p.x);
        if p.y < lo || p.y > hi {
            return false;
        }
        // A point inside the bounding trapezoid is inside the cell; the orientation checks keep
        // this honest for the degenerate (triangular) cases where the bounding test is not tight.
        let o = [
            orient2d(self.bl, self.br, p),
            orient2d(self.br, self.tr, p),
            orient2d(self.tr, self.tl, p),
            orient2d(self.tl, self.bl, p),
        ];
        o.iter().all(|d| *d == Orientation::CounterClockwise || *d == Orientation::Collinear)
    }

    /// The cell's area, by the shoelace formula.
    pub fn area(&self) -> f64 {
        let pts = [self.bl, self.br, self.tr, self.tl];
        let mut acc = 0.0;
        for i in 0..4 {
            acc += pts[i].cross(pts[(i + 1) % 4]);
        }
        acc / 2.0
    }

    /// The cell's four corners, counter-clockwise from bottom-left.
    #[inline]
    pub fn corners(&self) -> [Point2D; 4] {
        [self.bl, self.br, self.tr, self.tl]
    }
}

/// The y coordinate of the segment `a -> b` at abscissa `x`.
#[inline]
fn interpolate(a: Point2D, b: Point2D, x: f64) -> f64 {
    let dx = b.x - a.x;
    if dx == 0.0 {
        return a.y;
    }
    a.y + (b.y - a.y) * ((x - a.x) / dx)
}

/// A one-off working set for the sweep, so a `Decomposition` is built without reallocating per
/// slab.
#[derive(Debug, Default)]
pub(crate) struct CellBuilder {
    pub cells: Vec<Cell>,
    pub portals: Vec<Portal>,
}

impl CellBuilder {
    pub(crate) fn push_cell(&mut self, cell: Cell) -> CellId {
        self.cells.push(cell);
        (self.cells.len() - 1) as CellId
    }

    pub(crate) fn push_portal(&mut self, portal: Portal) -> PortalId {
        self.portals.push(portal);
        (self.portals.len() - 1) as PortalId
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_trapezoid() -> Cell {
        Cell {
            bl: Point2D::new(0.0, 0.0),
            br: Point2D::new(1.0, 0.0),
            tr: Point2D::new(1.0, 1.0),
            tl: Point2D::new(0.0, 1.0),
            left: None,
            right: None,
        }
    }

    #[test]
    fn centroid_and_ranges() {
        let c = unit_trapezoid();
        assert_eq!(c.centroid(), Point2D::new(0.5, 0.5));
        assert_eq!(c.x_range(), (0.0, 1.0));
        assert_eq!(c.y_range_at(0.5), (0.0, 1.0));
        assert_eq!(c.area(), 1.0);
        assert_eq!(c.corners()[0], Point2D::new(0.0, 0.0));
    }

    #[test]
    fn slanted_boundaries_interpolate() {
        let c = Cell {
            bl: Point2D::new(0.0, 0.0),
            br: Point2D::new(2.0, 0.0),
            tr: Point2D::new(2.0, 4.0),
            tl: Point2D::new(0.0, 2.0),
            left: None,
            right: None,
        };
        assert_eq!(c.bottom_at(1.0), 0.0);
        assert_eq!(c.top_at(0.0), 2.0);
        assert_eq!(c.top_at(1.0), 3.0);
        assert_eq!(c.top_at(2.0), 4.0);
        assert_eq!(c.area(), 6.0);
    }

    #[test]
    fn contains_interior_exterior_boundary() {
        let c = unit_trapezoid();
        assert!(c.contains(Point2D::new(0.5, 0.5)));
        assert!(c.contains(Point2D::new(0.0, 0.5)));
        assert!(c.contains(Point2D::new(1.0, 0.5)));
        assert!(c.contains(Point2D::new(0.5, 0.0)));
        assert!(!c.contains(Point2D::new(0.5, 1.5)));
        assert!(!c.contains(Point2D::new(1.5, 0.5)));
    }

    #[test]
    fn a_vertical_side_reports_its_endpoint_height() {
        // A cell whose x range collapses to a line: interpolation has no x to work with, so it
        // must return the stored height rather than dividing by zero.
        let c = Cell {
            bl: Point2D::new(1.0, 0.0),
            br: Point2D::new(1.0, 5.0),
            tr: Point2D::new(1.0, 7.0),
            tl: Point2D::new(1.0, 5.0),
            left: None,
            right: None,
        };
        assert_eq!(c.bottom_at(1.0), 0.0);
        assert_eq!(c.top_at(1.0), 5.0);
    }

    #[test]
    fn portal_endpoints_and_length() {
        let p = Portal { x: 3.0, lo: 1.0, hi: 4.0, left: 0, right: 1 };
        assert_eq!(p.length(), 3.0);
        assert_eq!(p.bottom(), Point2D::new(3.0, 1.0));
        assert_eq!(p.top(), Point2D::new(3.0, 4.0));
    }

    #[test]
    fn builder_returns_consecutive_ids() {
        let mut b = CellBuilder::default();
        assert_eq!(b.push_cell(unit_trapezoid()), 0);
        assert_eq!(b.push_cell(unit_trapezoid()), 1);
        assert_eq!(b.cells.len(), 2);
        b.push_portal(Portal { x: 0.0, lo: 0.0, hi: 1.0, left: 0, right: 1 });
        assert_eq!(b.portals.len(), 1);
    }
}
