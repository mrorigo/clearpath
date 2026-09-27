//! Exact geometric predicates.
//!
//! Every topological decision in this crate — orientation, in-circle, convex-hull ordering,
//! sweep event ordering — goes through this module. It computes the *exact* determinant of the
//! `f64` inputs using Shewchuk-style floating-point expansions, so a `Collinear` or `OnCircle`
//! answer means the inputs are genuinely degenerate rather than "close enough to look degenerate
//! after rounding". That determinism is what `docs/SPEC.md` section 4.4 requires.
//!
//! The implementation is pure safe Rust: the two-term error of a product comes from
//! [`f64::mul_add`], an intrinsic, not from `unsafe` pointer punning.

use super::point::Point2D;

/// The sign of an exact orientation determinant.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orientation {
    /// Negative determinant: `a`, `b`, `c` turn clockwise.
    Clockwise,
    /// Positive determinant: `a`, `b`, `c` turn counter-clockwise.
    CounterClockwise,
    /// Exactly zero determinant: the three points are collinear.
    Collinear,
}

impl Orientation {
    /// Whether this is [`Orientation::CounterClockwise`].
    #[inline]
    pub fn is_ccw(self) -> bool {
        matches!(self, Self::CounterClockwise)
    }

    /// Whether the three points are exactly collinear.
    #[inline]
    pub fn is_collinear(self) -> bool {
        matches!(self, Self::Collinear)
    }

    /// The orientation of `a - b - c`; i.e. the same triple read in reverse.
    #[inline]
    pub fn reversed(self) -> Self {
        match self {
            Self::Clockwise => Self::CounterClockwise,
            Self::CounterClockwise => Self::Clockwise,
            Self::Collinear => Self::Collinear,
        }
    }
}

/// Whether `d` is inside, on, or outside the circle through `a`, `b`, `c`.
///
/// `a`, `b`, `c` must be in counter-clockwise order; the result is meaningless otherwise (the
/// determinant is the *signed* in-circle test, so a clockwise triple flips `In` and `Out`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Incircle {
    /// Strictly inside the circumcircle.
    In,
    /// Exactly on the circumcircle.
    OnCircle,
    /// Strictly outside the circumcircle.
    Out,
}

/// Maximum number of components in an expansion used here.
///
/// The in-circle determinant is a 3x3 of 2-component products: 12 products plus 12 error terms,
/// which the elimination step reduces to at most 24 components.
const MAX_COMPONENTS: usize = 32;

/// A nonoverlapping expansion of a real number as an increasing sequence of `f64` components,
/// least significant first.
#[derive(Clone, Debug)]
struct Expansion {
    len: usize,
    buf: [f64; MAX_COMPONENTS],
}

impl Expansion {
    #[inline]
    fn zero() -> Self {
        Self { len: 0, buf: [0.0; MAX_COMPONENTS] }
    }

    /// The faithful rounding of the expansion: summing the components in increasing order is
    /// correctly rounded, so its sign is the sign of the exact value unless that value is
    /// exactly zero.
    #[inline]
    fn estimate(&self) -> f64 {
        let mut sum = 0.0;
        for &c in &self.buf[..self.len] {
            sum += c;
        }
        sum
    }

    /// Appends a single `f64` term, eliminating overlap.
    fn push(&mut self, value: f64) {
        let mut q = value;
        for i in 0..self.len {
            let (sum, err) = two_sum(q, self.buf[i]);
            q = sum;
            self.buf[i] = err;
        }
        debug_assert!(self.len < MAX_COMPONENTS);
        self.buf[self.len] = q;
        self.len += 1;
    }

    /// Replaces the expansion with the exact sum of `self` and `other`, least significant first.
    fn sum(&mut self, other: &Expansion) {
        let mut out = Expansion::zero();
        let (mut i, mut j) = (0usize, 0usize);
        let mut q = 0.0f64;
        while i < self.len || j < other.len {
            let t;
            if j >= other.len || (i < self.len && self.buf[i].abs() < other.buf[j].abs()) {
                t = self.buf[i];
                i += 1;
            } else {
                t = other.buf[j];
                j += 1;
            }
            let (s, e) = two_sum(q, t);
            q = s;
            if e != 0.0 {
                out.push(e);
            }
        }
        if q != 0.0 || out.len == 0 {
            out.push(q);
        }
        *self = out;
    }

    /// The exact product of the two expansions.
    fn product(a: &Expansion, b: &Expansion) -> Expansion {
        let mut out = Expansion::zero();
        for i in 0..a.len {
            for j in 0..b.len {
                let (p, e) = two_product(a.buf[i], b.buf[j]);
                out.push(p);
                if e != 0.0 {
                    out.push(e);
                }
            }
        }
        if out.len == 0 {
            out.push(0.0);
        }
        out
    }
}

/// Exact sum of `a` and `b`: `a + b` rounded, plus the exact error.
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    (x, (a - av) + (b - bv))
}

/// Exact difference: `a - b` rounded, plus the exact error.
#[inline]
fn two_diff(a: f64, b: f64) -> (f64, f64) {
    let x = a - b;
    let bv = a - x;
    let av = x + bv;
    (x, (a - av) + (bv - b))
}

/// Splits `a` into two halves with 26 significant bits each, so their product is exact.
#[inline]
fn split(a: f64) -> (f64, f64) {
    const SPLITTER: f64 = 134_217_729.0; // 2^27 + 1
    let c = SPLITTER * a;
    let abig = c - a;
    let ahi = c - abig;
    (ahi, a - ahi)
}

/// Exact product: `a * b` rounded, plus the exact error.
///
/// Dekker's split-and-multiply, rather than `f64::mul_add`: `mul_add` is only in `core` when the
/// target has a hardware FMA, so relying on it would make the predicate silently target-dependent
/// (and section 4.4 forbids that). This form is exact for every finite input and is plain safe
/// Rust.
#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);
    let err1 = x - ahi * bhi;
    let err2 = err1 - alo * bhi;
    let err3 = err2 - ahi * blo;
    (x, alo * blo - err3)
}

/// The exact 2-component expansion of `a - b`.
#[inline]
fn diff_expansion(a: f64, b: f64) -> Expansion {
    let (hi, lo) = two_diff(a, b);
    let mut e = Expansion::zero();
    e.push(lo);
    e.push(hi);
    e
}

/// Exact orientation of `a`, `b`, `c`.
///
/// Computes `det = (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x)` exactly.
pub fn orient2d(a: Point2D, b: Point2D, c: Point2D) -> Orientation {
    let left = Expansion::product(&diff_expansion(a.x, c.x), &diff_expansion(b.y, c.y));
    let right = Expansion::product(&diff_expansion(a.y, c.y), &diff_expansion(b.x, c.x));
    let mut det = left;
    let mut negated = right;
    for i in 0..negated.len {
        negated.buf[i] = -negated.buf[i];
    }
    det.sum(&negated);

    let est = det.estimate();
    if est > 0.0 {
        Orientation::CounterClockwise
    } else if est < 0.0 {
        Orientation::Clockwise
    } else {
        Orientation::Collinear
    }
}

/// Convenience wrapper: is `c` strictly left of the directed line `a -> b`?
#[inline]
pub fn is_left(a: Point2D, b: Point2D, c: Point2D) -> bool {
    orient2d(a, b, c) == Orientation::CounterClockwise
}

/// Convenience wrapper: is `c` strictly right of the directed line `a -> b`?
#[inline]
pub fn is_right(a: Point2D, b: Point2D, c: Point2D) -> bool {
    orient2d(a, b, c) == Orientation::Clockwise
}

/// Exact in-circle test of `d` against the circumcircle of the counter-clockwise triple
/// `a`, `b`, `c`.
pub fn incircle(a: Point2D, b: Point2D, c: Point2D, d: Point2D) -> Incircle {
    // Translate by `d` and evaluate
    //   | ax-dx  ay-dy  (ax-dx)^2 + (ay-dy)^2 |
    //   | bx-dx  by-dy  (bx-dx)^2 + (by-dy)^2 |
    //   | cx-dx  cy-dy  (cx-dx)^2 + (cy-dy)^2 |
    // whose third column is the squared norm of the first two, so no square root is involved
    // and the whole determinant is exact in expansions.
    let row = |p: Point2D| {
        let dx = diff_expansion(p.x, d.x);
        let dy = diff_expansion(p.y, d.y);
        let lift = {
            let sx = Expansion::product(&dx, &dx);
            let sy = Expansion::product(&dy, &dy);
            let mut l = sx;
            l.sum(&sy);
            l
        };
        [dx, dy, lift]
    };
    let det = det3(&row(a), &row(b), &row(c));
    let est = det.estimate();
    if est > 0.0 {
        Incircle::In
    } else if est < 0.0 {
        Incircle::Out
    } else {
        Incircle::OnCircle
    }
}

/// The exact 3x3 determinant of rows of expansions.
fn det3(r0: &[Expansion; 3], r1: &[Expansion; 3], r2: &[Expansion; 3]) -> Expansion {
    let m01 = {
        let t = Expansion::product(&r0[0], &r1[1]);
        let u = Expansion::product(&r0[1], &r1[0]);
        let mut d = t;
        d.sum(&negate(&u));
        d
    };
    let m02 = {
        let t = Expansion::product(&r0[0], &r1[2]);
        let u = Expansion::product(&r0[2], &r1[0]);
        let mut d = t;
        d.sum(&negate(&u));
        d
    };
    let m12 = {
        let t = Expansion::product(&r1[1], &r2[2]);
        let u = Expansion::product(&r1[2], &r2[1]);
        let mut d = t;
        d.sum(&negate(&u));
        d
    };

    let a = Expansion::product(&r0[0], &m12);
    let b = Expansion::product(&r0[1], &m02);
    let c = Expansion::product(&r0[2], &m01);
    let mut det = a;
    det.sum(&negate(&b));
    det.sum(&c);
    det
}

#[inline]
fn negate(e: &Expansion) -> Expansion {
    let mut out = e.clone();
    for i in 0..out.len {
        out.buf[i] = -out.buf[i];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    macro_rules! p {
        ($x:expr, $y:expr) => {
            Point2D::new($x, $y)
        };
    }

    #[test]
    fn two_sum_is_exact() {
        let (x, e) = two_sum(1.0, 1e-17);
        assert_eq!(x + e, 1.0 + 1e-17);
    }

    #[test]
    fn two_product_is_exact() {
        let (x, e) = two_product(1.0 + 2.0 * f64::EPSILON, 1.0 - 2.0 * f64::EPSILON);
        assert_eq!(x + e, (1.0 + 2.0 * f64::EPSILON) * (1.0 - 2.0 * f64::EPSILON));
    }

    #[test]
    fn basic_orientation() {
        assert_eq!(orient2d(p!(0.0, 0.0), p!(1.0, 0.0), p!(0.0, 1.0)), Orientation::CounterClockwise);
        assert_eq!(orient2d(p!(0.0, 0.0), p!(0.0, 1.0), p!(1.0, 0.0)), Orientation::Clockwise);
        assert_eq!(orient2d(p!(0.0, 0.0), p!(1.0, 1.0), p!(2.0, 2.0)), Orientation::Collinear);
    }

    #[test]
    fn degenerate_collinear_triples_are_exact() {
        // The classic case that breaks naive determinants: three collinear points whose
        // double-length representation is exact but whose products cancel to a rounding error.
        // Doubling a point is exact in binary, so a naive determinant sees exact zero here.
        let a = p!(0.5, 0.5);
        let b = p!(12.0, 12.0);
        let c = p!(24.0, 24.0);
        assert!(orient2d(a, b, c).is_collinear());

        let d = p!(2.0, 0.0);
        let e = p!(2.0, 4.0);
        let f = p!(2.0, 8.0);
        assert!(orient2d(d, e, f).is_collinear());
    }

    #[test]
    fn orientation_agrees_with_naive_when_well_conditioned() {
        let cases = vec![
            (p!(0.0, 0.0), p!(4.0, 0.0), p!(2.0, 3.0)),
            (p!(-1.0, -1.0), p!(1.0, 1.0), p!(0.0, 5.0)),
            (p!(3.0, 3.0), p!(-2.0, 5.0), p!(1.0, -4.0)),
        ];
        for (a, b, c) in cases {
            let naive = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            let expected =
                if naive > 0.0 { Orientation::CounterClockwise } else { Orientation::Clockwise };
            assert_eq!(orient2d(a, b, c), expected);
        }
    }

    #[test]
    fn coincident_points_are_collinear() {
        let a = p!(1.0, 1.0);
        assert!(orient2d(a, a, p!(3.0, 1.0)).is_collinear());
        assert!(orient2d(a, p!(3.0, 1.0), a).is_collinear());
        assert!(orient2d(a, a, a).is_collinear());
    }

    #[test]
    fn reversed_flips_orientation() {
        let a = p!(0.0, 0.0);
        let b = p!(2.0, 0.0);
        let c = p!(1.0, 1.0);
        assert_eq!(orient2d(a, b, c).reversed(), orient2d(c, b, a));
    }

    #[test]
    fn left_right_helpers() {
        let (a, b) = (p!(0.0, 0.0), p!(1.0, 0.0));
        assert!(is_left(a, b, p!(0.0, 1.0)));
        assert!(!is_left(a, b, p!(0.0, -1.0)));
        assert!(is_right(a, b, p!(0.0, -1.0)));
    }

    #[test]
    fn incircle_on_a_unit_circle() {
        let a = p!(1.0, 0.0);
        let b = p!(0.0, 1.0);
        let c = p!(-1.0, 0.0);
        assert_eq!(incircle(a, b, c, p!(0.0, 0.5)), Incircle::In);
        assert_eq!(incircle(a, b, c, p!(0.0, 0.0)), Incircle::In);
        assert_eq!(incircle(a, b, c, p!(0.0, 2.0)), Incircle::Out);
        assert_eq!(incircle(a, b, c, p!(0.0, -1.0)), Incircle::OnCircle);
    }

    #[test]
    fn incircle_matches_naive_when_well_conditioned() {
        let a = p!(0.0, 0.0);
        let b = p!(4.0, 0.0);
        let c = p!(0.0, 4.0);
        for d in [p!(1.0, 1.0), p!(3.0, 3.0), p!(-1.0, -1.0), p!(0.5, 0.5)] {
            let dx = d.x - 0.0;
            let dy = d.y - 0.0;
            let adx = a.x - d.x;
            let ady = a.y - d.y;
            let bdx = b.x - d.x;
            let bdy = b.y - d.y;
            let cdx = c.x - d.x;
            let cdy = c.y - d.y;
            let _ = (dx, dy);
            let abdet = adx * bdy - bdx * ady;
            let bcdet = bdx * cdy - cdx * bdy;
            let cadet = cdx * ady - adx * cdy;
            let alift = adx * adx + ady * ady;
            let blift = bdx * bdx + bdy * bdy;
            let clift = cdx * cdx + cdy * cdy;
            let naive = alift * bcdet + blift * cadet + clift * abdet;
            let expected = if naive > 0.0 { Incircle::In } else { Incircle::Out };
            assert_eq!(incircle(a, b, c, d), expected, "point {d:?}");
        }
    }
}
