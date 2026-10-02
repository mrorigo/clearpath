//! Low-level geometric primitives and robust predicates.

use alloc::vec::Vec;

use crate::error::PathPlanError;

pub mod clearance;
pub mod point;
pub mod polygon;
pub mod predicates;

pub use point::Point2D;
pub use polygon::{BoundingBox, Polygon};
pub use predicates::{Orientation, incircle, orient2d};

/// The free space a route is planned in: the workspace, less the obstacles.
///
/// `margin` is **not** baked into this type. It is applied one step later, by eroding each
/// decomposition cell, because that is exact and the alternative — offsetting the obstacle
/// polygons — is not. See `docs/SPEC.md` section 3.3.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSpace {
    /// The workspace bounding box. A route stays inside it.
    pub workspace: BoundingBox,
    /// The obstacles. Their interiors are forbidden.
    pub obstacles: Vec<Polygon>,
}

impl FreeSpace {
    /// Constructs a free space, validating the workspace and the margin-independent inputs.
    ///
    /// # Errors
    /// Returns [`PathPlanError::DegenerateWorkspace`] if the bounding box has no positive extent.
    pub fn new(workspace: BoundingBox, obstacles: Vec<Polygon>) -> Result<Self, PathPlanError> {
        if workspace.is_degenerate() {
            return Err(PathPlanError::DegenerateWorkspace);
        }
        Ok(Self { workspace, obstacles })
    }

    /// Whether `p` is in the free space, ignoring `margin`.
    pub fn contains(&self, p: Point2D) -> bool {
        self.workspace.contains(p) && !self.obstacles.iter().any(|o| o.contains(p))
    }
}

/// A line as `(slope, intercept)`.
pub type Line = (f64, f64);

/// Whether the segment `a..b` lies inside the band between `lower` and `upper` over
/// `[x0, x1]`.
///
/// Everything here is affine in `x`, so the question is a pair of half-interval intersections
/// rather than a search: the segment is in the band where `y_segment - lower >= 0` and
/// `upper - y_segment >= 0`, and each of those is affine.
///
/// This exists because sampling a segment to decide whether it is inside a region is a *sample*,
/// not a proof, and the sample density that makes it convincing over a 1000-unit segment is
/// expensive enough to dominate the query. A cell of the decomposition is a trapezoid, so its
/// vertical bounds are lines, and the containment question is exact in closed form.
pub fn segment_within_band(
    a: Point2D,
    b: Point2D,
    lower: Line,
    upper: Line,
    x0: f64,
    x1: f64,
) -> bool {
    let (u, v) = (x0.min(x1), x0.max(x1));
    if u > v {
        return false;
    }
    // A vertical segment is not a function of `x`, so it has no slope-intercept form. Representing
    // it as one at slope 0 substitutes the *horizontal* line `y = a.y` for it, which silently
    // discards `b.y` and reduces the test to the single ordinate of `a` -- so a segment that
    // leaves the band at its far end is accepted. Its abscissa range is the single value `a.x`,
    // and the band over that abscissa is a closed ordinate interval, so both endpoints are tested
    // against both bounds there.
    if b.x == a.x {
        // A vertical segment over a single abscissa reduces to its ordinate range, so it is inside
        // the band exactly when the low end is at or above the lower bound and the high end is at
        // or below the upper bound. Both comparisons reuse `above` at a degenerate interval
        // [a.x, a.x], which keeps the boundary-tie tolerance identical to the sloped path.
        let (lo, hi) = (a.y.min(b.y), a.y.max(b.y));
        return above((0.0, lo), lower, a.x, a.x)
            && above(upper, (0.0, hi), a.x, a.x);
    }
    let slope = (b.y - a.y) / (b.x - a.x);
    let segment = (slope, a.y - slope * a.x);
    above(segment, lower, u, v) && above(upper, segment, u, v)
}

/// Whether `probe >= bound` throughout `[lo, hi]`, both affine.
///
/// The difference of the two lines is affine, so its minimum over the interval is at an endpoint
/// and the test is two evaluations. Dividing for a root and comparing against it does not work: a
/// segment that *starts* exactly on a cell edge makes the root land on the endpoint to within a
/// rounding error in the wrong direction, and the containment test then fails for a segment that is
/// exactly on the boundary — which is the case the whole margin model is built to allow.
fn above(probe: Line, bound: Line, lo: f64, hi: f64) -> bool {
    let m = probe.0 - bound.0;
    let c = probe.1 - bound.1;
    let at = |x: f64| {
        // A tolerance proportional to the magnitudes involved, so it is a distance rather than a
        // count of ulps. It is needed because both lines are *reconstructions*: a segment endpoint
        // that lies exactly on a cell edge is on it bit for bit, but re-evaluating the segment as
        // `slope * x + intercept` at that same abscissa can come out a fraction of an ulp the other
        // side. Admitting a boundary tie is the safe direction — the margin model already allows a
        // curve to touch what it was pushed off.
        let v = m * x + c;
        let scale = m.abs() * x.abs() + c.abs() + 1.0;
        v >= -1e-9 * scale
    };
    if m >= 0.0 {
        at(lo)
    } else {
        at(hi)
    }
}

/// Tolerance used only for non-topological "is this effectively zero" tests.
///
/// Never used for orientation, containment, or ordering decisions — those go through
/// [`predicates::orient2d`], which is exact.
pub const EPS: f64 = 1e-12;

#[cfg(test)]
mod band_tests {
    use super::*;

    fn line(m: f64, b: f64) -> Line {
        (m, b)
    }

    #[test]
    fn a_segment_inside_a_band_is_accepted() {
        assert!(segment_within_band(
            Point2D::new(0.0, 5.0),
            Point2D::new(10.0, 7.0),
            line(0.0, 0.0),
            line(0.0, 10.0),
            0.0,
            10.0
        ));
    }

    #[test]
    fn a_segment_leaving_the_band_is_rejected() {
        assert!(!segment_within_band(
            Point2D::new(0.0, 5.0),
            Point2D::new(10.0, 15.0),
            line(0.0, 0.0),
            line(0.0, 10.0),
            0.0,
            10.0
        ));
    }

    /// The case that was failing: a segment whose endpoint lies *exactly* on a cell edge, so the
    /// containment is a tie at that abscissa.
    #[test]
    fn a_segment_starting_exactly_on_the_lower_bound_is_accepted() {
        let a = Point2D::new(15.769401554340265, 77.23059844565974);
        let b = Point2D::new(90.0, 90.0);
        let lower = line(0.0, 77.23059844565974);
        let upper = line(0.0, 97.76940155434026);
        assert!(segment_within_band(a, b, lower, upper, 15.769401554340265, 16.230598445659734));
    }

    /// The same, with a sloped lower bound, which is where a root comparison loses a ULP.
    #[test]
    fn a_segment_starting_exactly_on_a_sloped_bound_is_accepted() {
        let a = Point2D::new(1.0, 1.0);
        let b = Point2D::new(9.0, 9.0);
        // The lower bound is the line y = x, which `a` lies on exactly.
        assert!(segment_within_band(a, b, line(1.0, 0.0), line(0.0, 20.0), 1.0, 9.0));
    }

    #[test]
    fn a_segment_touching_a_bound_is_accepted_either_side() {
        let a = Point2D::new(0.0, 5.0);
        let b = Point2D::new(10.0, 5.0);
        assert!(segment_within_band(a, b, line(0.0, 5.0), line(0.0, 9.0), 0.0, 10.0));
        assert!(segment_within_band(a, b, line(0.0, 1.0), line(0.0, 5.0), 0.0, 10.0));
    }

    #[test]
    fn a_vertical_segment_is_tested_by_its_ordinate_range() {
        // A vertical segment's abscissa range is a single value, but its *ordinate* range is the
        // whole segment, so both endpoints have to be inside the band.
        assert!(segment_within_band(
            Point2D::new(5.0, 2.0),
            Point2D::new(5.0, 8.0),
            line(0.0, 1.0),
            line(0.0, 9.0),
            5.0,
            5.0
        ));
        assert!(!segment_within_band(
            Point2D::new(5.0, 0.0),
            Point2D::new(5.0, 8.0),
            line(0.0, 1.0),
            line(0.0, 9.0),
            5.0,
            5.0
        ));
    }

    /// A vertical segment is not a function of `x`, and its abscissa range is the single value
    /// `a.x`. Writing it as a slope-intercept pair at slope 0 therefore substitutes the *horizontal*
    /// line `y = a.y` for it and discards `b.y` altogether, so the test comes down to the one
    /// ordinate of `a` and a segment that leaves the band at its far end is accepted. The
    /// ordinate range of the segment is the whole segment, so a vertical knot is inside the band
    /// at `a.x` only when *both* ends are inside it, whichever end is which.
    ///
    /// All of the negative cases live in one function on purpose. The broken code reads a single
    /// ordinate, so a lone case is satisfied by whatever that one ordinate happens to be, and a
    /// second case with a different one could pass in a way that says nothing about the first. Held
    /// together they cannot: each one puts the violation at a different end and in a different
    /// direction, so no choice of a single ordinate clears them all. The last group needs the
    /// bounds themselves to be sloped, because a horizontal bound has the same value everywhere,
    /// and so cannot tell a bound read at the segment's own abscissa from one read at zero, which
    /// is the mistake that survives when the ordinate range is fixed but the abscissa is not.
    #[test]
    fn a_vertical_segment_is_rejected_when_either_end_leaves_the_band() {
        let flat_lower = line(0.0, 3.0);
        let flat_upper = line(0.0, 9.0);
        let inside = Point2D::new(5.0, 4.0);
        let high = Point2D::new(5.0, 20.0);
        // Out through the upper bound, and the same segment with its ends swapped: `a..b` and
        // `b..a` are one segment, so neither order may be accepted.
        assert!(!segment_within_band(inside, high, flat_lower, flat_upper, 5.0, 5.0));
        assert!(!segment_within_band(high, inside, flat_lower, flat_upper, 5.0, 5.0));
        // Down through the lower bound, from an ordinate that is inside the band, so the far end
        // is the only one that can betray it.
        let top = Point2D::new(5.0, 8.0);
        let deep = Point2D::new(5.0, -50.0);
        assert!(!segment_within_band(top, deep, flat_lower, flat_upper, 5.0, 5.0));
        // Wholly below the band, with no endpoint near it. A band checked at the wrong abscissa is
        // still a band, so this is what catches a bound read at the wrong place.
        let far = Point2D::new(5.0, -100.0);
        assert!(!segment_within_band(far, deep, flat_lower, flat_upper, 5.0, 5.0));
        // A sloped band at the same abscissa: y = x below and y = 2x + 5 above, so at x = 5 the
        // band is y in [5, 15].
        let sloped_lower = line(1.0, 0.0);
        let sloped_upper = line(2.0, 5.0);
        let in_band = Point2D::new(5.0, 6.0);
        let under = Point2D::new(5.0, 4.0);
        let over = Point2D::new(5.0, 16.0);
        let within = Point2D::new(5.0, 14.0);
        assert!(segment_within_band(in_band, within, sloped_lower, sloped_upper, 5.0, 5.0));
        assert!(!segment_within_band(in_band, over, sloped_lower, sloped_upper, 5.0, 5.0));
        assert!(!segment_within_band(in_band, under, sloped_lower, sloped_upper, 5.0, 5.0));
    }

    /// The positive control for the test above: a vertical segment whose both ends are inside the
    /// band is inside it, so widening the test to the whole ordinate range does not start
    /// rejecting the legitimate vertical knots the decomposition produces. Without it, a fix that
    /// rejected every vertical segment would satisfy every negative case.
    #[test]
    fn a_vertical_segment_fully_inside_the_band_is_accepted() {
        let lower = line(0.0, 3.0);
        let upper = line(0.0, 9.0);
        let a = Point2D::new(5.0, 4.0);
        let b = Point2D::new(5.0, 8.0);
        assert!(segment_within_band(a, b, lower, upper, 5.0, 5.0));
    }

    #[test]
    fn a_partial_range_is_tested_only_where_it_is_asked() {
        // Over [0, 1] the segment is inside; over the whole range it is not.
        let a = Point2D::new(0.0, 0.5);
        let b = Point2D::new(10.0, 5.0);
        assert!(segment_within_band(a, b, line(0.0, 0.0), line(0.0, 10.0), 0.0, 1.0));
        assert!(!segment_within_band(a, b, line(0.0, 0.0), line(0.0, 1.0), 0.0, 10.0));
    }
}
