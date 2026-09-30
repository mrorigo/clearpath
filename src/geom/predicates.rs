//! Exact geometric predicates.
//!
//! Every topological decision in this crate — orientation, in-circle, convex-hull ordering,
//! sweep event ordering — goes through this module. It computes the *exact* determinant of the
//! `f64` inputs using Shewchuk-style floating-point expansions, so a `Collinear` or `OnCircle`
//! answer means the inputs are genuinely degenerate rather than "close enough to look degenerate
//! after rounding". That determinism is what `docs/SPEC.md` section 4.4 requires.
//!
//! The implementation is pure safe Rust: the two-term error of a product comes from Dekker's
//! split-and-multiply over `f64`, not from `mul_add` (which is only in `core` when the target has a
//! hardware FMA, and would make the predicate target-dependent) and not from `unsafe` pointer
//! punning.

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
    ///
    /// A term that is not finite, and an expansion that has already filled its buffer, are both
    /// dropped rather than appended. This is a backstop against writing out of bounds, and it is
    /// what stops an overflowed expansion from aborting in a release build, where the
    /// `debug_assert!` below is compiled out.
    ///
    /// A term that is exactly zero is *not* dropped: zero is a finite component, and dropping it
    /// would truncate the expansion -- which is why `sum` guards its own pushes rather than relying
    /// on this method. Only `inf` and `NaN` are discarded, and discarding them leaves the finite
    /// terms in place, so an expansion that has overflowed degrades to a rounded answer rather than
    /// to a wrong topological decision.
    fn push(&mut self, value: f64) {
        if !value.is_finite() || self.len >= MAX_COMPONENTS {
            return;
        }
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

/// The degree of the `orient2d` determinant in the coordinates: a product of two differences.
const TWO_COMPONENT_DEGREE: u32 = 2;

/// The degree of the in-circle determinant in the coordinates.
///
/// Each row of the 3x3 is `(x, y, x^2 + y^2)`, whose entries are at most quadratic, so a cofactor
/// expansion multiplies three of them and the determinant is homogeneous of degree six. The
/// determinant is a *difference* of such products, so the terms cancel heavily and a smaller degree
/// would often suffice -- but "often" is not a correctness argument, and the exactness claim is
/// unconditional, so the rescale is budgeted for the full degree.
const SIX_COMPONENT_DEGREE: u32 = 6;

/// The exponent budget a scaled product is held to.
///
/// A degree-`d` determinant over coordinates of magnitude `m` is a sum of products of at most `d`
/// of them, so holding `m` at or below `2^(EXPANSION_PRODUCT_BUDGET / d)` keeps every product below
/// `2^1020` -- inside `f64`'s 1023-bit range, with the rest covering the error terms `two_product`
/// adds and the summation. The budget sits just under the point where the unscaled products would
/// overflow, so the rescale engages as late as possible and leaves ordinary geometry untouched.
const EXPANSION_PRODUCT_BUDGET: i32 = 1020;

/// A power-of-two scale that brings a magnitude-`largest` set of `degree` products back into range,
/// or `None` when the inputs are already in range.
///
/// A determinant of degree `d` over coordinates of magnitude at most `m` is a sum of products of at
/// most `d` of them, so holding `m` at or below `2^(EXPANSION_PRODUCT_BUDGET / d)` keeps every
/// product, and the error terms `two_product` adds, representable. The caller multiplies its inputs by
/// the returned factor; because the factor is a power of two, those multiplications are exact, so
/// the rescaling introduces no rounding at all. `None` keeps the common case free of extra work.
fn exact_rescale(largest: f64, degree: u32) -> Option<f64> {
    if !largest.is_finite() || largest == 0.0 {
        return None;
    }
    // The magnitude is bounded by 2 raised to the unbiased binary exponent, so bounding that
    // exponent bounds the products regardless of the mantissa. It is read from the bit pattern
    // because the integer log is not available in core for the bare-metal targets this crate builds
    // for, and a floating-point logarithm would be a rounding that must not decide an exactness
    // question. The exponent field is biased by 1023; a subnormal has a zero field and is far
    // below any limit, so it reports zero and correctly needs no rescaling.
    let biased = ((largest.to_bits() >> 52) & 0x7ff) as i32;
    let exponent = (biased - 1023).max(0);
    let limit = EXPANSION_PRODUCT_BUDGET / degree as i32;
    if exponent <= limit {
        return None;
    }
    // 2 raised to the negated shift, written straight into the exponent field so the factor is
    // exactly a power of two. The field is biased by 1023, and a negative biased exponent is
    // unreachable for a normal double, so the guard below covers the one case where the shift is
    // large enough to leave the normal range -- and returning `None` there is safe, because it
    // only happens for a subnormal input that never needed rescaling to begin with.
    let shift = exponent - limit;
    let scale_biased = 1023 - shift;
    if scale_biased <= 0 {
        return None;
    }
    debug_assert!(shift > 0);
    Some(f64::from_bits((scale_biased as u64) << 52))
}

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
    // The exact path represents the determinant as a sum of `f64` products, so it can only be
    // exact while the largest product is representable -- that is, while the coordinates are below
    // `sqrt(f64::MAX)`, about 1.34e154. Above that, `two_product` returns `x = inf` and `e = NaN`,
    // the expansion sums to `NaN`, and because both `est > 0.0` and `est < 0.0` are false for
    // `NaN` the predicate would report `Collinear` for a triple that is definitively not
    // collinear.
    //
    // Rescaling the inputs by a power of two fixes that without giving up exactness. The
    // determinant is homogeneous of degree two, so `det(s*p) = s^2 * det(p)`, and `s^2` is a
    // strictly positive multiple: it cannot change the sign, and it is represented exactly, so the
    // sign of the rescaled determinant is the sign of the true one. Dividing by a power of two is
    // itself exact for any input not driven to zero by underflow, and a power of two is used
    // precisely so no rounding is introduced in the division.
    //
    // The scale is taken from the already-subtracted differences rather than the raw coordinates:
    // those are what the products are formed from, and scaling them cannot overflow a subtraction
    // that has already happened.
    let largest = acx.abs().max(bcx.abs()).max(acy.abs()).max(bcy.abs());
    let (acx, bcx, acy, bcy) = match exact_rescale(largest, TWO_COMPONENT_DEGREE) {
        Some(scale) => (acx * scale, bcx * scale, acy * scale, bcy * scale),
        None => (acx, bcx, acy, bcy),
    };
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
    // The exact path squares its inputs into the third column, so the determinant is quartic in the
    // coordinates and overflows far earlier than `orient2d` does -- at `f64::MAX^(1/4)`, about
    // 1.16e77. Past that the expansion fills with `inf` and `NaN` components and, because
    // `MAX_COMPONENTS` is a fixed buffer whose `debug_assert!` guard is compiled out in release,
    // the overflow used to write past the end of the buffer and abort.
    //
    // The same rescaling argument applies: the in-circle determinant is homogeneous of degree
    // four, so scaling the inputs by a power of two multiplies the determinant by a strictly
    // positive exact power of two and cannot change its sign.
    let largest = ax
        .abs()
        .max(ay.abs())
        .max(bx.abs())
        .max(by.abs())
        .max(cx.abs())
        .max(cy.abs());
    // A 3x3 determinant of rows (x, y, x^2+y^2) with all entries of degree at most 2 is homogeneous
    // of degree 2 + 2 + 2 = 6, not 4: the squared-norm column contributes a factor of the whole
    // coordinate once per pairing, so each product in the cofactor expansion has degree up to 6.
    // Budgeting for degree 4 under-scales, and the residual overflow is what pushed the answer to
    // the wrong sign at the top of the range.
    let (ax, ay, bx, by, cx, cy) = match exact_rescale(largest, SIX_COMPONENT_DEGREE) {
        Some(scale) => (ax * scale, ay * scale, bx * scale, by * scale, cx * scale, cy * scale),
        None => (ax, ay, bx, by, cx, cy),
    };
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

    /// `orient2d`'s exact path represents the determinant as a sum of `f64` products, so it can
    /// only be exact while `magnitude^2` is representable -- that is, while the coordinates are
    /// below `sqrt(f64::MAX)` (about 1.34e154). Above that the products overflow to infinity, the
    /// expansion sums to `NaN`, and both `est > 0.0` and `est < 0.0` are false, so the predicate
    /// reported `Collinear` for triples that are definitively not collinear.
    #[test]
    fn orient2d_is_exact_when_the_determinant_overflows_f64() {
        // A right angle at the origin, scaled past the overflow threshold. `c` is strictly left of
        // `a -> b`, so the exact answer is counter-clockwise at every magnitude.
        for &magnitude in &[1e155, 1e160, 1e200, 1e300, 1.7e308] {
            let (a, b, c) = (p!(0.0, 0.0), p!(magnitude, 0.0), p!(0.0, magnitude));
            assert_eq!(
                orient2d(a, b, c),
                Orientation::CounterClockwise,
                "counter-clockwise right angle of size {magnitude:e} was not reported as such"
            );
            assert_eq!(
                orient2d(a, c, b),
                Orientation::Clockwise,
                "clockwise right angle of size {magnitude:e} was not reported as such"
            );
        }
    }

    /// The same overflow reached through a translation rather than through the origin, which is
    /// the shape the crate actually meets: the predicate subtracts `c` first, so the intermediate
    /// differences are what overflow, not the input coordinates.
    #[test]
    fn orient2d_is_exact_above_the_overflow_threshold_with_a_nonzero_origin() {
        for &magnitude in &[1e155, 1e200, 1e300] {
            let o = magnitude;
            let (a, b, c) = (p!(o, o), p!(2.0 * o, o), p!(o, 2.0 * o));
            assert_eq!(
                orient2d(a, b, c),
                Orientation::CounterClockwise,
                "translated counter-clockwise triple of size {magnitude:e} was not reported as such"
            );
            assert_eq!(
                orient2d(a, c, b),
                Orientation::Clockwise,
                "translated clockwise triple of size {magnitude:e} was not reported as such"
            );
        }
    }

    /// The prescale must not turn a genuinely degenerate triple into a non-degenerate one, and
    /// must not stop answering `Collinear` for collinear input at any magnitude.
    #[test]
    fn orient2d_still_reports_genuinely_collinear_triples_at_every_magnitude() {
        for &magnitude in &[1e0, 1e10, 1e150, 1e155, 1e200, 1e300] {
            let (a, b, c) = (p!(0.0, 0.0), p!(magnitude, magnitude), p!(3.0 * magnitude, 3.0 * magnitude));
            assert_eq!(
                orient2d(a, b, c),
                Orientation::Collinear,
                "collinear triple of size {magnitude:e} was not reported as collinear"
            );
            // A repeated point is degenerate in the same way.
            let (a, b, c) = (p!(magnitude, magnitude), p!(magnitude, magnitude), p!(0.0, 0.0));
            assert_eq!(
                orient2d(a, b, c),
                Orientation::Collinear,
                "repeated-point triple of size {magnitude:e} was not reported as collinear"
            );
        }
    }

    /// The prescale divides by a power of two, which is exact, so the predicate must agree with the
    /// unscaled answer for every magnitude where the unscaled answer is meaningful. This is the
    /// check that the fix does not quietly change behaviour below the threshold.
    #[test]
    fn orient2d_agrees_across_the_overflow_threshold() {
        let mut checked = 0usize;
        for exponent in 0..=308u32 {
            let magnitude = f64::from_bits((1023u64 + u64::from(exponent)) << 52);
            if !(magnitude.is_finite() && magnitude > 0.0) {
                continue;
            }
            // Exactly representable coordinates keep the ground truth exact: the determinant is
            // then a product of powers of two and its sign is not in doubt.
            let (a, b, c) = (p!(0.0, 0.0), p!(magnitude, 0.0), p!(0.0, magnitude));
            assert_eq!(
                orient2d(a, b, c),
                Orientation::CounterClockwise,
                "power-of-two magnitude 2^{exponent} lost its orientation"
            );
            checked += 1;
        }
        assert!(checked > 200, "too few magnitudes checked: {checked}");
    }

    /// `incircle`'s exact path squares its inputs into the third column, so it overflows at
    /// `f64::MAX^(1/4)` (about 1.16e77) rather than at `sqrt(f64::MAX)`. It did so before this
    /// fix by reporting `OnCircle` for points that are clearly inside or outside, and by
    /// exhausting the expansion buffer and panicking outright above about 5.6e102.
    #[test]
    fn incircle_is_exact_when_the_determinant_overflows_f64() {
        for &magnitude in &[1e78, 1e99, 1e150, 1e300] {
            let (a, b, c) = (p!(0.0, 0.0), p!(magnitude, 0.0), p!(0.0, magnitude));
            // The circumcircle of this right angle is centred at the midpoint of the hypotenuse,
            // (m/2, m/2), with radius m/sqrt(2), so that midpoint is strictly inside it. The
            // quarter point (m/4, m/4) is deliberately not used here: it sits exactly on the
            // circle, at half the radius, and that is the `OnCircle` case below.
            let inside = p!(magnitude / 2.0, magnitude / 2.0);
            assert_eq!(
                incircle(a, b, c, inside),
                Incircle::In,
                "interior point of size {magnitude:e} was not reported inside"
            );
            // And a point far beyond the hypotenuse is strictly outside it.
            let outside = p!(magnitude * 4.0, 0.0);
            assert_eq!(
                incircle(a, b, c, outside),
                Incircle::Out,
                "exterior point of size {magnitude:e} was not reported outside"
            );
        }
    }

    /// The four points `(0,0)`, `(r,0)`, `(0,r)`, `(r,r)` are exactly concyclic whenever `r` is a
    /// power of two, so the predicate must report `OnCircle` rather than collapsing to `In` or
    /// `Out` with an overflowed expansion. The magnitude below is a power of two for that reason:
    /// for an arbitrary `r` the fourth point is concyclic only up to rounding, so the exact answer
    /// is then legitimately `In` or `Out` rather than `OnCircle`.
    #[test]
    fn incircle_reports_on_circle_when_the_determinant_overflows_f64() {
        // Powers of two spanning the range where the quartic determinant overflows, from just
        // above the fourth root of the largest double up to the largest double itself.
        for exponent in [256u32, 300, 400, 500, 600, 700, 1020] {
            let magnitude = f64::from_bits(u64::from(1023 + exponent) << 52);
            let (a, b, c, d) = (
                p!(0.0, 0.0),
                p!(magnitude, 0.0),
                p!(0.0, magnitude),
                p!(magnitude, magnitude),
            );
            assert_eq!(
                incircle(a, b, c, d),
                Incircle::OnCircle,
                "concyclic point at 2^{exponent} was not reported on the circle"
            );
        }
    }

    /// The in-circle determinant is homogeneous of degree six, not four: each row of the 3x3 is
    /// `(x, y, x^2 + y^2)`, so a cofactor expansion multiplies three at-most-quadratic rows. A
    /// rescale budgeted for degree four under-scales, and the residual overflow at the very top of
    /// the `f64` range was enough to flip the reported sign. This case is the regression for that.
    #[test]
    fn incircle_is_exact_at_the_top_of_the_f64_range() {
        // The largest finite double, halved so that the point below is a usable probe. 8.9e307 is
        // near f64::MAX and its square would be about 1e616 if the rescale did not engage.
        let m = 8.9e307f64;
        let (a, b, c) = (p!(0.0, 0.0), p!(m, 0.0), p!(0.0, m));
        // The hypotenuse midpoint, strictly inside the circumcircle.
        assert_eq!(
            incircle(a, b, c, p!(m / 2.0, m / 2.0)),
            Incircle::In,
            "circumcentre of a right angle at the top of the range was not reported inside"
        );
        // A point three quarters of the way out along the x axis, strictly outside it. The factor
        // is applied as m/2 * 3 rather than m * 3 so the probe itself stays finite.
        assert_eq!(
            incircle(a, b, c, p!(m / 2.0 * 3.0, 0.0)),
            Incircle::Out,
            "exterior point at the top of the range was not reported outside"
        );
    }

    /// The expansion buffer is a fixed 32 slots and its `debug_assert!` guard is compiled out in
    /// release, so an overflowed expansion used to write past the end of the buffer and abort. This
    /// exercises the same magnitudes in a release build.
    #[test]
    fn incircle_does_not_overflow_the_expansion_buffer() {
        for &magnitude in &[1e103, 1e150, 1e300, 1.7e308] {
            let (a, b, c) = (p!(0.0, 0.0), p!(magnitude, 0.0), p!(0.0, magnitude));
            let d = p!(magnitude * 4.0, magnitude * 4.0);
            // The answer is what is under test; the point of the case is that this returns at all.
            let _ = incircle(a, b, c, d);
        }
    }
}
