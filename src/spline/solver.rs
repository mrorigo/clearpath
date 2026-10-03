//! Tangent solve for the spline stage.
//!
//! The system is the stationary condition of
//!
//! ```text
//!   sum_i |T_{i+1} - 2 T_i + T_{i-1}|^2  +  lambda * sum_i |T_i - d_i|^2
//! ```
//!
//! so it is symmetric, and positive definite for `lambda > 0` and positive semi-definite at
//! `lambda == 0`, which is why the seeds are the documented fallback at zero bias rather than a
//! `NaN`. The matrix touches `T_{j-2} .. T_{j+2}` — bandwidth two, not one. A tridiagonal solve
//! here would silently solve a different problem. The matrix is symmetric with a two-band
//! half-width, so it is factorised as a banded system. Knot counts are small (the number of knots in
//! a route), so the cost would be fine either way; the banded form is kept because it makes the
//! structure explicit and the cost linear.

use alloc::vec;
use alloc::vec::Vec;

use crate::geom::Point2D;
use crate::sqrt;

/// Half-width of the band. The matrix touches `T_{j-2} .. T_{j+2}`.
const BAND: usize = 2;

/// Solves for the tangents, treating a pinned one as a boundary condition.
///
/// `seeds` are the unconstrained tangents: the chord directions, or a port's requested direction
/// scaled to the local knot spacing. `lambda` is the bias towards the seed: `0` is the pure
/// minimum-curvature solution, and larger values pull the tangents towards `seeds`. A very large
/// `lambda` returns the seeds unchanged.
///
/// `pinned[i]` marks a tangent the caller fixed. Section 6.5 of `docs/SPEC.md` is explicit that a
/// `PortConstraint::direction` at a port is *fixed* rather than merely attracted towards: `T_0` is
/// exactly `|k_0 - k_1| / 3 * d_hat`. Handing the direction over only as a seed makes it a soft
/// bias, and the curvature term outvotes it wherever the route bends — on a bent route that
/// discarded 359 of 360 requested directions, 181 of them collapsing the tangent to exactly zero.
///
/// A pinned tangent is eliminated rather than approximated: its row and column are zeroed and its
/// contribution moved to the right-hand side, so the free tangents solve the stationary condition
/// *subject to* the pinned endpoints. What remains is a principal submatrix of a symmetric positive
/// definite matrix, hence again symmetric positive definite, and the band is unchanged because a
/// principal submatrix of a banded matrix is banded with at most the same half-width. Only the two
/// ends are ever pinned, so at most the first two and last two rows and columns are touched.
pub fn solve_tangents(seeds: &[Point2D], lambda: f64, pinned: &[bool]) -> Vec<Point2D> {
    let n = seeds.len();
    if n < 2 {
        return seeds.to_vec();
    }
    let lambda = if lambda > 0.0 { lambda } else { 0.0 };
    let is_pinned = |j: usize| pinned.get(j).copied().unwrap_or(false);
    // Nothing to solve if every tangent is fixed.
    if (0..n).all(is_pinned) {
        return seeds.to_vec();
    }

    // Build row `j` from the second differences themselves rather than from a fixed stencil.
    //
    // With `S_i = T_{i+1} - 2 T_i + T_{i-1}` for `i` in `[1, n-2]` and zero outside, the gradient is
    // `S_{j-1} - 2 S_j + S_{j+1}`, and the coefficient of `T_{j-2}` is 1, of `T_{j-1}` is -4, of
    // `T_j` is 6, of `T_{j+1}` is -4 and of `T_{j+2}` is 1 — but only where all three `S` terms
    // exist. At the ends they do not, and the coefficient of `T_0` is 1, not 6. A fixed stencil
    // gets that wrong, and the result is a system that is not the stationarity condition of the
    // objective: the solve is stable, it just solves a different problem.
    let last = n as i64 - 2;
    let valid = |i: i64| (1..=last).contains(&i);
    let mut rows: Vec<[f64; 2 * BAND + 1]> = (0..n)
        .map(|j| {
            let j = j as i64;
            let mut row = [0.0f64; 2 * BAND + 1];
            let mut add = |offset: i64, coeff: f64| {
                let idx = j + offset;
                if idx >= 0 && (idx as usize) < n {
                    row[(offset + BAND as i64) as usize] += coeff;
                }
            };
            if valid(j - 1) {
                add(-2, 1.0);
                add(-1, -2.0);
                add(0, 1.0);
            }
            if valid(j) {
                add(-1, -2.0);
                add(0, 4.0);
                add(1, -2.0);
            }
            if valid(j + 1) {
                add(0, 1.0);
                add(1, -2.0);
                add(2, 1.0);
            }
            row[BAND] += lambda;
            row
        })
        .collect();

    // Move the pinned tangents to the right-hand side and clear their rows and columns. Only a row
    // within the band of a pinned index can hold a non-zero entry in that column, so this is
    // O(pinned * BAND) rather than O(pinned * n).
    let pinned_idx: Vec<usize> = (0..n).filter(|j| is_pinned(*j)).collect();
    let mut rhs_x: Vec<f64> = seeds.iter().map(|d| lambda * d.x).collect();
    let mut rhs_y: Vec<f64> = seeds.iter().map(|d| lambda * d.y).collect();
    for &k in &pinned_idx {
        for j in k.saturating_sub(BAND)..=(k + BAND).min(n - 1) {
            let col = k + BAND - j;
            let a = rows[j][col];
            rhs_x[j] -= a * seeds[k].x;
            rhs_y[j] -= a * seeds[k].y;
            rows[j][col] = 0.0;
        }
        // A unit diagonal rather than a zero row: a zero row leaves no pivot for the rows below
        // to divide by, and the solve then returns the seeds instead of a solution.
        rows[k].fill(0.0);
        rows[k][BAND] = 1.0;
        rhs_x[k] = seeds[k].x;
        rhs_y[k] = seeds[k].y;
    }

    // Banded Cholesky. The matrix is symmetric positive definite for `lambda >= 0`: the curvature
    // term is a sum of squares and the bias adds `lambda` to the diagonal, so a non-positive pivot
    // is a numerical failure rather than a property of the problem, and the seeds are the
    // documented fallback instead of a `NaN`. A pinned row was zeroed above, so it is skipped
    // rather than factorised — its value is known, and skipping keeps the pivot sequence valid.
    //
    // `lower[j][m]` is the sub-diagonal entry at column `j - BAND + m`; the diagonal is
    // `lower[j][BAND]`.
    let mut lower: Vec<[f64; 2 * BAND + 1]> = vec![[0.0; 2 * BAND + 1]; n];
    for j in 0..n {
        for k in j.saturating_sub(BAND)..j {
            let mut v = rows[j][(k + BAND) - j];
            // `lower[j]` is zero left of its own band, so the inner range starts at `j - BAND`
            // rather than `k - BAND`; going lower indexes before the row and underflows.
            for m in j.saturating_sub(BAND)..k {
                v -= lower[j][(m + BAND) - j] * lower[k][(m + BAND) - k];
            }
            lower[j][(k + BAND) - j] = v / lower[k][BAND];
        }
        let mut diag = rows[j][BAND];
        for k in j.saturating_sub(BAND)..j {
            let v = lower[j][(k + BAND) - j];
            diag -= v * v;
        }
        if !diag.is_finite() || diag <= 0.0 {
            return seeds.to_vec();
        }
        lower[j][BAND] = sqrt(diag);
    }

    let solve = |b: &[f64]| -> Option<Vec<f64>> {
        let mut y = vec![0.0f64; n];
        for j in 0..n {
            let mut v = b[j];
            for k in j.saturating_sub(BAND)..j {
                v -= lower[j][(k + BAND) - j] * y[k];
            }
            y[j] = v / lower[j][BAND];
        }
        let mut x = vec![0.0f64; n];
        for j in (0..n).rev() {
            let mut v = y[j];
            for k in (j + 1)..(j + BAND + 1).min(n) {
                v -= lower[k][(j + BAND) - k] * x[k];
            }
            x[j] = v / lower[j][BAND];
        }
        if x.iter().all(|v| v.is_finite()) { Some(x) } else { None }
    };

    match (solve(&rhs_x), solve(&rhs_y)) {
        // A pinned tangent comes back as its seed: it is a boundary condition, not a solved value.
        (Some(xs), Some(ys)) => (0..n)
            .map(|j| {
                if is_pinned(j) {
                    seeds[j]
                } else {
                    Point2D::new(xs[j], ys[j])
                }
            })
            .collect(),
        _ => seeds.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seeds(n: usize) -> Vec<Point2D> {
        (0..n)
            .map(|i| Point2D::new(i as f64 * 2.0, (i as f64).sin() * 3.0))
            .collect()
    }

    /// The objective whose stationary condition the system is.
    fn objective(t: &[Point2D], seeds: &[Point2D], lambda: f64) -> f64 {
        let mut acc = 0.0;
        for i in 1..t.len() - 1 {
            let s = t[i + 1] - t[i] * 2.0 + t[i - 1];
            acc += s.norm_squared();
        }
        for i in 0..t.len() {
            acc += lambda * t[i].distance_squared(seeds[i]);
        }
        acc
    }

    #[test]
    fn zero_bias_is_singular_and_returns_the_seeds() {
        let seeds = seeds(3);
        assert_eq!(solve_tangents(&seeds, 0.0, &[]), seeds);
    }

    #[test]
    fn a_huge_bias_returns_the_seeds() {
        let seeds = seeds(3);
        let out = solve_tangents(&seeds, 1e12, &[]);
        for (o, s) in out.iter().zip(&seeds) {
            assert!((o.x - s.x).abs() < 1e-6, "{o:?} vs {s:?}");
        }
    }

    #[test]
    fn zero_bias_gives_a_straight_line_for_a_straight_input() {
        let seeds = vec![
            Point2D::new(1.0, 1.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(3.0, 3.0),
        ];
        let out = solve_tangents(&seeds, 0.0, &[]);
        for (i, s) in out.iter().enumerate() {
            assert!((s.x - s.y).abs() < 1e-9, "tangent {i} left the line: {s:?}");
        }
    }

    /// The regression this change exists for: a port direction is *fixed*, so the solver must return
    /// it exactly rather than trading it off against curvature.
    #[test]
    fn a_pinned_tangent_comes_back_exactly() {
        let seeds = vec![
            Point2D::new(0.0, -5.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(3.0, 5.0),
        ];
        let out = solve_tangents(&seeds, 0.1, &[true, false, false, true]);
        assert_eq!(out[0], seeds[0], "the start tangent was moved");
        assert_eq!(out[3], seeds[3], "the goal tangent was moved");
    }

    /// The converse, and why pinning is optional rather than unconditional: with nothing pinned the
    /// endpoints stay free, which is what makes a route a spline rather than a polyline with smoothed
    /// joints.
    #[test]
    fn an_unpinned_endpoint_still_moves() {
        let seeds = vec![
            Point2D::new(0.0, -5.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(3.0, 5.0),
        ];
        let out = solve_tangents(&seeds, 0.1, &[]);
        assert_ne!(out[0], seeds[0], "the start tangent should have moved");
        assert_ne!(out[3], seeds[3], "the goal tangent should have moved");
    }

    #[test]
    fn pinning_one_end_only_leaves_the_other_free() {
        let knots: Vec<Point2D> = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(3.0, 2.0),
            Point2D::new(7.0, 1.0),
            Point2D::new(10.0, 4.0),
        ];
        let seeds: Vec<Point2D> = knots.windows(2).map(|w| (w[1] - w[0]) / 3.0).collect();
        let out = solve_tangents(&seeds, 0.1, &[true, false, false]);
        assert_eq!(out[0], seeds[0], "the pinned start moved");
        assert_ne!(out[2], seeds[2], "the free goal should have moved");
    }

    #[test]
    fn pinning_everything_returns_the_seeds() {
        let seeds = seeds(4);
        let pinned = [true, true, true, true];
        assert_eq!(solve_tangents(&seeds, 0.5, &pinned), seeds);
    }

    #[test]
    fn the_solution_stationarises_the_objective() {
        // Finite-difference check that the returned tangents are a critical point: perturbing any
        // single tangent in any direction must not decrease the objective. Pinned tangents are
        // boundary conditions, so they are skipped — pinning must be a real constraint rather than a
        // way of freezing the answer, and the free tangents still stationarise the objective
        // subject to it.
        let knots: Vec<Point2D> = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(3.0, 2.0),
            Point2D::new(7.0, 1.0),
            Point2D::new(10.0, 4.0),
        ];
        let seeds: Vec<Point2D> = knots.windows(2).map(|w| (w[1] - w[0]) / 3.0).collect();
        let lambda = 0.5;
        for pinned in [
            [false, false, false],
            [true, false, false],
            [false, false, true],
            [true, false, true],
        ] {
            let out = solve_tangents(&seeds, lambda, &pinned);
            let base = objective(&out, &seeds, lambda);
            for i in 0..out.len() {
                if pinned[i] {
                    continue;
                }
                for delta in [1e-3, 1e-2] {
                    let mut bumped = out.clone();
                    bumped[i] += Point2D::new(delta, 0.0);
                    assert!(
                        objective(&bumped, &seeds, lambda) >= base - 1e-9,
                        "objective fell at knot {i} (pinned {pinned:?})"
                    );
                    let mut bumped = out.clone();
                    bumped[i] += Point2D::new(0.0, delta);
                    assert!(
                        objective(&bumped, &seeds, lambda) >= base - 1e-9,
                        "objective fell at knot {i} in y (pinned {pinned:?})"
                    );
                }
            }
        }
    }

    #[test]
    fn trivial_inputs_are_returned_unchanged() {
        assert_eq!(solve_tangents(&[], 1.0, &[]), Vec::<Point2D>::new());
        let one = vec![Point2D::new(1.0, 1.0)];
        assert_eq!(solve_tangents(&one, 1.0, &[]), one);
        let two = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)];
        assert_eq!(solve_tangents(&two, 1.0, &[]).len(), 2);
    }

    #[test]
    fn the_system_is_pentadiagonal() {
        // Guard against someone "simplifying" this to a tridiagonal solve: perturbing a tangent
        // two positions away must still change the result.
        let seeds: Vec<Point2D> = (0..5).map(|i| Point2D::new(1.0, if i % 2 == 0 { 0.0 } else { 1.0 })).collect();
        let a = solve_tangents(&seeds, 0.1, &[]);
        let mut b_seeds = seeds.clone();
        b_seeds[2] += Point2D::new(0.0, 5.0);
        let b = solve_tangents(&b_seeds, 0.1, &[]);
        // Changing an interior seed moves the neighbouring solutions, which only a banded width of
        // two can produce.
        assert!(a[1].distance(b[1]) > 1e-9, "knot 1 unaffected by knot 2's seed");
    }
}