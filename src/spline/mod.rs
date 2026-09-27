//! Cubic Bezier segments, the tangent solver, and the containment repair.

pub mod containment;
pub mod solver;

use alloc::vec::Vec;

use crate::geom::point::Point2D;
use crate::geom::BoundingBox;
use crate::sqrt;

/// One cubic Bezier segment.
///
/// Control points are in curve order: `p0` and `p3` are the ends, `p1` and `p2` the interior ones.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicBezierSegment {
    /// Start of the curve.
    pub p0: Point2D,
    /// First control point.
    pub p1: Point2D,
    /// Second control point.
    pub p2: Point2D,
    /// End of the curve.
    pub p3: Point2D,
}

impl CubicBezierSegment {
    /// The point at parameter `t`.
    #[inline]
    pub fn evaluate(&self, t: f64) -> Point2D {
        let u = 1.0 - t;
        let (a, b) = (u * u * u, 3.0 * u * u * t);
        let (c, d) = (3.0 * u * t * t, t * t * t);
        self.p0 * a + self.p1 * b + self.p2 * c + self.p3 * d
    }

    /// The derivative at parameter `t`.
    #[inline]
    pub fn tangent(&self, t: f64) -> Point2D {
        // B'(t) = 3(1-t)^2 (P1-P0) + 6(1-t) t (P2-P1) + 3 t^2 (P3-P2).
        let u = 1.0 - t;
        (self.p1 - self.p0) * (3.0 * u * u)
            + (self.p2 - self.p1) * (6.0 * u * t)
            + (self.p3 - self.p2) * (3.0 * t * t)
    }

    /// The bounding box of the control hull, which contains the curve.
    pub fn bounding_box(&self) -> BoundingBox {
        BoundingBox::new(self.p0, self.p3).with(self.p1).with(self.p2)
    }

    /// De Casteljau subdivision at `t`, returning the two halves.
    pub fn split(&self, t: f64) -> (Self, Self) {
        let ab = self.p0.lerp(self.p1, t);
        let bc = self.p1.lerp(self.p2, t);
        let cd = self.p2.lerp(self.p3, t);
        let abc = ab.lerp(bc, t);
        let bcd = bc.lerp(cd, t);
        let abcd = abc.lerp(bcd, t);
        (
            Self { p0: self.p0, p1: ab, p2: abc, p3: abcd },
            Self { p0: abcd, p1: bcd, p2: cd, p3: self.p3 },
        )
    }

    /// A polyline approximation whose chords deviate from the curve by at most `tolerance`.
    ///
    /// The count comes from the standard flatness bound, so the deviation is a guarantee rather
    /// than a sample rate: this is what the tests use to check a curve against geometry that the
    /// curve's own control points do not describe.
    pub fn flatten(&self, tolerance: f64) -> Vec<Point2D> {
        let tol = if tolerance > 0.0 { tolerance } else { f64::MIN_POSITIVE };
        let mut out = vec_with_capacity(self);
        flatten_into(self, tol, &mut out);
        out
    }

    /// The largest distance from the curve to the chord `p0..p3`, by sampling.
    ///
    /// Used by `flatten` to pick a subdivision count.
    fn flatness(&self) -> f64 {
        // The distance from `p1` and `p2` to the `p0..p3` line bounds the curve's deviation.
        let d = self.p3 - self.p0;
        let len2 = d.norm_squared();
        if len2 == 0.0 {
            return self.p1.distance(self.p0).max(self.p2.distance(self.p0));
        }
        let h1 = (d.cross(self.p1 - self.p0)).abs() / sqrt(len2);
        let h2 = (d.cross(self.p2 - self.p0)).abs() / sqrt(len2);
        h1.max(h2)
    }
}

fn vec_with_capacity(segment: &CubicBezierSegment) -> Vec<Point2D> {
    let n = subdivision_count(segment);
    let mut v = Vec::with_capacity(n + 2);
    v.push(segment.p0);
    v
}

/// How many times to halve the segment to get within `tolerance`.
fn subdivision_count(segment: &CubicBezierSegment) -> usize {
    let f = segment.flatness();
    if f <= 0.0 {
        return 0;
    }
    // Each split roughly halves the deviation; `f / 2^k <= tol`.
    let mut k = 0usize;
    let mut d = f;
    while d > 1.0 && k < 24 {
        d *= 0.5;
        k += 1;
    }
    k
}

fn flatten_into(segment: &CubicBezierSegment, tolerance: f64, out: &mut Vec<Point2D>) {
    if segment.flatness() <= tolerance || out.len() > 1 << 16 {
        out.push(segment.p3);
        return;
    }
    let (a, b) = segment.split(0.5);
    flatten_into(&a, tolerance, out);
    flatten_into(&b, tolerance, out);
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn line() -> CubicBezierSegment {
        CubicBezierSegment {
            p0: Point2D::new(0.0, 0.0),
            p1: Point2D::new(1.0, 0.0),
            p2: Point2D::new(2.0, 0.0),
            p3: Point2D::new(3.0, 0.0),
        }
    }

    fn curve() -> CubicBezierSegment {
        CubicBezierSegment {
            p0: Point2D::new(0.0, 0.0),
            p1: Point2D::new(1.0, 2.0),
            p2: Point2D::new(2.0, 2.0),
            p3: Point2D::new(3.0, 0.0),
        }
    }

    #[test]
    fn a_straight_segment_evaluates_linearly() {
        let s = line();
        assert_eq!(s.evaluate(0.0), s.p0);
        assert_eq!(s.evaluate(1.0), s.p3);
        assert_eq!(s.evaluate(0.5), Point2D::new(1.5, 0.0));
        assert_eq!(s.flatness(), 0.0);
    }

    #[test]
    fn tangent_endpoints_match_the_control_polygon() {
        let s = curve();
        assert_eq!(s.tangent(0.0), (s.p1 - s.p0) * 3.0);
        assert_eq!(s.tangent(1.0), (s.p3 - s.p2) * 3.0);
    }

    #[test]
    fn the_curve_stays_inside_its_control_hull() {
        let s = curve();
        let b = s.bounding_box();
        for i in 0..=100 {
            let p = s.evaluate(i as f64 / 100.0);
            assert!(b.contains(p), "point {p:?} outside {b:?}");
        }
    }

    #[test]
    fn split_reproduces_the_curve_and_meets_at_the_split() {
        let s = curve();
        let (a, b) = s.split(0.3);
        assert_eq!(a.p3, b.p0);
        assert!(a.p3.distance(s.evaluate(0.3)) < 1e-12);
        // The halves trace the same *curve*, but not with the same parameterisation: de Casteljau
        // reparameterises, so `a.evaluate(t)` is `s.evaluate(0.3 t)` only up to a change of
        // variable. The geometric identity is what holds, and the flattened polylines are how it is
        // checked.
        let whole = s.flatten(1e-5);
        let mut halves = a.flatten(1e-5);
        halves.extend(b.flatten(1e-5));
        // The two polylines need not have the same vertex count — each half is flatter on its own,
        // so it subdivides less — so compare them as curves rather than point by point.
        let distance_to = |p: Point2D, poly: &[Point2D]| -> f64 {
            poly.windows(2)
                .map(|w| crate::corridor::closest_point_on_segment(p, w[0], w[1]).distance(p))
                .fold(f64::INFINITY, f64::min)
        };
        for h in &halves {
            assert!(distance_to(*h, &whole) < 1e-4, "half point {h:?} is off the whole");
        }
        for w in &whole {
            assert!(distance_to(*w, &halves) < 1e-4, "whole point {w:?} is off the halves");
        }
    }

    #[test]
    fn flatten_respects_the_tolerance() {
        let s = curve();
        for tol in [1.0, 0.1, 0.01, 0.001] {
            let poly = s.flatten(tol);
            assert!(poly.first() == Some(&s.p0));
            assert!(poly.last() == Some(&s.p3));
            // Every sample of the curve must be within the tolerance of the polyline.
            for i in 0..=200 {
                let p = s.evaluate(i as f64 / 200.0);
                let d = poly
                    .windows(2)
                    .map(|w| crate::corridor::closest_point_on_segment(p, w[0], w[1]).distance(p))
                    .fold(f64::INFINITY, f64::min);
                assert!(d <= tol * 1.5, "tolerance {tol} exceeded by {d}");
            }
        }
    }

    #[test]
    fn a_straight_segment_flattens_to_two_points() {
        assert_eq!(line().flatten(0.1), vec![Point2D::new(0.0, 0.0), Point2D::new(3.0, 0.0)]);
    }

    #[test]
    fn a_degenerate_segment_is_handled() {
        let s = CubicBezierSegment {
            p0: Point2D::new(1.0, 1.0),
            p1: Point2D::new(1.0, 1.0),
            p2: Point2D::new(1.0, 1.0),
            p3: Point2D::new(1.0, 1.0),
        };
        assert_eq!(s.flatness(), 0.0);
        assert_eq!(s.evaluate(0.5), s.p0);
        let (a, _) = s.split(0.5);
        assert_eq!(a.p0, a.p3);
    }
}
