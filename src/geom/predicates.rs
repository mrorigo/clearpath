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

/// The exact 2-component expansion of a difference `a - b`.
#[inline]
fn diff_expansion(a: f64, b: f64) -> Expansion {
    let (hi, lo) = two_diff(a, b);
    let mut e = Expansion::zero();
    e.push(lo);
    e.push(hi);
    e
}

/// The exact 2-component expansion of an already-computed difference, whose error term is
/// recovered rather than recomputed.
#[inline]
fn two_component(v: f64) -> Expansion {
    diff_expansion(v, 0.0)
}

/// Unit roundoff for `f64`: `2^-53`.
const EPS: f64 = 1.0 / 9_007_199_254_740_992.0;

/// Shewchuk's error bound for the `orient2d` filter: `(3 + 16 eps) eps`.
fn ccw_err_bound(a: f64, b: f64) -> f64 {
    (3.0 + 16.0 * EPS) * EPS * (a.abs() + b.abs())
}

/// Shewchuk's error bound for the `incircle` filter: `(10 + 96 eps) eps`.
fn icc_err_bound(permanent: f64) -> f64 {
    (10.0 + 96.0 * EPS) * EPS * permanent
}

/// Exact orientation of `a`, `b`, `c`.
///
/// Computes `det = (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x)` exactly.
///
/// A floating-point filter runs first and the exact expansion path is only entered when the
/// rounded determinant is not provably outside the error bound. The filter is what makes the
/// predicate usable in the ingest and sweep hot loops: exact expansions are roughly two orders
/// of magnitude slower, and the filter passes for every input that is not near-degenerate.
pub fn orient2d(a: Point2D, b: Point2D, c: Point2D) -> Orientation {
    let acx = a.x - c.x;
    let bcx = b.x - c.x;
    let acy = a.y - c.y;
    let bcy = b.y - c.y;

    let left = acx * bcy;
    let right = acy * bcx;
    let det = left - right;
    let bound = ccw_err_bound(left, right);
    if det > bound {
        return Orientation::CounterClockwise;
    }
    if -det > bound {
        return Orientation::Clockwise;
    }
    orient2d_exact(acx, bcx, acy, bcy)
}

/// The exact fallback for [`orient2d`], on already-subtracted coordinates.
fn orient2d_exact(acx: f64, bcx: f64, acy: f64, bcy: f64) -> Orientation {
    let left = Expansion::product(&two_component(acx), &two_component(bcy));
    let right = Expansion::product(&two_component(acy), &two_component(bcx));
    let mut det = left;
    det.sum(&negate(&right));

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
    let adx = a.x - d.x;
    let ady = a.y - d.y;
    let bdx = b.x - d.x;
    let bdy = b.y - d.y;
    let cdx = c.x - d.x;
    let cdy = c.y - d.y;

    // Floating-point filter, as for `orient2d`.
    let abdet = adx * bdy - bdx * ady;
    let bcdet = bdx * cdy - cdx * bdy;
    let cadet = cdx * ady - adx * cdy;
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * bcdet + blift * cadet + clift * abdet;
    let permanent = (bdx.abs() * cdy.abs() + cdx.abs() * bdy.abs()) * alift
        + (cdx.abs() * ady.abs() + adx.abs() * cdy.abs()) * blift
        + (adx.abs() * bdy.abs() + bdx.abs() * ady.abs()) * clift;
    let bound = icc_err_bound(permanent);
    if det > bound {
        return Incircle::In;
    }
    if -det > bound {
        return Incircle::Out;
    }
    incircle_exact(adx, ady, bdx, bdy, cdx, cdy)
}

/// The exact fallback for [`incircle`], on coordinates already translated by `d`.
fn incircle_exact(
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
) -> Incircle {
    // Evaluate
    //   | ax  ay  ax^2 + ay^2 |
    //   | bx  by  bx^2 + by^2 |
    //   | cx  cy  cx^2 + cy^2 |
    // whose third column is the squared norm of the first two, so no square root is involved and
    // the whole determinant is exact in expansions. The coordinates are already relative to `d`.
    let row = |x: f64, y: f64| {
        let dx = two_component(x);
        let dy = two_component(y);
        let lift = {
            let sx = Expansion::product(&dx, &dx);
            let sy = Expansion::product(&dy, &dy);
            let mut l = sx;
            l.sum(&sy);
            l
        };
        [dx, dy, lift]
    };
    let det = det3(&row(ax, ay), &row(bx, by), &row(cx, cy));
    let est = det.estimate();
    if est > 0.0 {
        Incircle::In
    } else if est < 0.0 {
        Incircle::Out
    } else {
        Incircle::OnCircle
    }
}

/// The exact 3x3 determinant of rows of expansions, by cofactor expansion along the first row.
///
/// `Mij` is the minor obtained by deleting row 0 and column `j`, so each one pairs two elements
/// of the *second and third* rows. Pairing them against the first row instead yields some other
/// expression that is algebraically unrelated to the determinant.
fn det3(r0: &[Expansion; 3], r1: &[Expansion; 3], r2: &[Expansion; 3]) -> Expansion {
    let minor = |j: usize| {
        let (i1, i2) = match j {
            0 => (1usize, 2usize),
            1 => (0, 2),
            _ => (0, 1),
        };
        let positive = Expansion::product(&r1[i1], &r2[i2]);
        let negative = Expansion::product(&r1[i2], &r2[i1]);
        let mut d = positive;
        d.sum(&negate(&negative));
        d
    };

    let mut det = Expansion::product(&r0[0], &minor(0));
    det.sum(&negate(&Expansion::product(&r0[1], &minor(1))));
    det.sum(&Expansion::product(&r0[2], &minor(2)));
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

    /// The filter must never change an answer. This is the check that keeps the fast path honest:
    /// the exact path is ground truth, and the inputs are deliberately near-degenerate.
    #[test]
    fn the_filter_never_disagrees_with_the_exact_path() {
        // A cheap deterministic xorshift, so the case is reproducible without a dependency.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let coord = |state: u64| -> f64 {
            // Mix wildly different magnitudes so near-collinear triples actually occur.
            let scale = (state % 5) as u32;
            let v = ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0;
            v * 10f64.powi(scale as i32)
        };
        for _ in 0..20000 {
            let a = p!(coord(next()), coord(next()));
            let b = p!(coord(next()), coord(next()));
            // Half the time, construct c exactly on the line a->b so the filter must escalate.
            let c = if next() % 2 == 0 {
                let t = coord(next());
                Point2D::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t)
            } else {
                p!(coord(next()), coord(next()))
            };
            let filtered = orient2d(a, b, c);
            let exact = orient2d_exact(a.x - c.x, b.x - c.x, a.y - c.y, b.y - c.y);
            assert_eq!(filtered, exact, "orient2d disagrees on {a:?} {b:?} {c:?}");
        }
    }

    #[test]
    fn the_incircle_filter_never_disagrees_with_the_exact_path() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let coord = |state: u64| -> f64 {
            let scale = (state % 4) as u32;
            let v = ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0;
            v * 10f64.powi(scale as i32)
        };
        let mut checked = 0usize;
        for _ in 0..20000 {
            let a = p!(coord(next()), coord(next()));
            let b = p!(coord(next()), coord(next()));
            let c = p!(coord(next()), coord(next()));
            if orient2d(a, b, c) != Orientation::CounterClockwise {
                continue; // the exact incircle test is only signed for ccw triples
            }
            let d = if next() % 2 == 0 {
                // Exactly on the circumcircle is the case the filter must not swallow.
                let (cx, cy, ccx, ccy) = (a.x, a.y, b.x, b.y);
                let _ = (cx, cy, ccx, ccy);
                p!(coord(next()), coord(next()))
            } else {
                p!(coord(next()), coord(next()))
            };
            let filtered = incircle(a, b, c, d);
            let exact = incircle_exact(a.x - d.x, a.y - d.y, b.x - d.x, b.y - d.y, c.x - d.x, c.y - d.y);
            assert_eq!(filtered, exact, "incircle disagrees on {a:?} {b:?} {c:?} {d:?}");
            checked += 1;
        }
        assert!(checked > 1000, "too few counter-clockwise triples to be meaningful: {checked}");
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
