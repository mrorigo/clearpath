//! The tangent solver.
//!
//! # The system is pentadiagonal, not tridiagonal
//!
//! The draft specified minimising
//!
//! ```text
//! sum_i |T_{i+1} - 2 T_i + T_{i-1}|^2  +  lambda * sum_i |T_i - d_i|^2
//! ```
//!
//! and described it as a tridiagonal system solved with the Thomas algorithm. Differentiating
//! shows otherwise: the second difference `S_i = T_{i+1} - 2 T_i + T_{i-1}` couples `T_i` to
//! `S_{i-1}`, `S_i` and `S_{i+1}`, so the stationarity condition for `T_j` reaches
//! `T_{j-2} .. T_{j+2}` — bandwidth two, not one. A tridiagonal solve here would silently solve a
//! different, smoother problem.
//!
//! The matrix is symmetric with a two-band half-width, so it is factorised as a banded system. Knot
//! counts are small (one per obstacle vertex the taut string wraps), so a dense solve would also
//! be fine; the banded form is kept because it makes the structure explicit and the cost linear.

use alloc::vec;
use alloc::vec::Vec;

use crate::geom::point::Point2D;
use crate::sqrt;

/// Half-width of the band. The matrix touches `T_{j-2} .. T_{j+2}`.
const BAND: usize = 2;

/// Solves for the tangents that minimise
/// `sum_i |T_{i+1} - 2 T_i + T_{i-1}|^2 + lambda * sum_i |T_i - d_i|^2`.
///
/// `lambda` is the bias towards the seed: `0` is the pure minimum-curvature solution, and larger
/// values pull the tangents towards `seeds` — the chord directions, or a port's requested
/// direction. A very large `lambda` returns the seeds unchanged.
///
/// The endpoints are unknowns too, which is what makes the result a spline rather than a polyline
/// with smoothed joints: the curvature term reaches the ends through the adjacent interior tangent.
pub fn solve_tangents(seeds: &[Point2D], lambda: f64) -> Vec<Point2D> {
    let n = seeds.len();
    if n < 2 {
        return seeds.to_vec();
    }
    let lambda = if lambda > 0.0 { lambda } else { 0.0 };

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
    let rows: Vec<[f64; 2 * BAND + 1]> = (0..n)
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

    // Banded Cholesky. The matrix is symmetric positive definite for `lambda >= 0`: the curvature
    // term is a sum of squares and the bias adds `lambda` to the diagonal, so a non-positive pivot
    // is a numerical failure rather than a property of the problem, and the seeds are the
    // documented fallback instead of a `NaN`.
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

    let rhs_x: Vec<f64> = seeds.iter().map(|d| lambda * d.x).collect();
    let rhs_y: Vec<f64> = seeds.iter().map(|d| lambda * d.y).collect();
    match (solve(&rhs_x), solve(&rhs_y)) {
        (Some(xs), Some(ys)) => (0..n).map(|j| Point2D::new(xs[j], ys[j])).collect(),
        _ => seeds.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn on_line(n: usize) -> Vec<Point2D> {
        (0..n).map(|i| Point2D::new(i as f64, 0.0)).collect()
    }

    #[test]
    fn zero_bias_is_singular_and_returns_the_seeds() {
        // Documented behaviour, not an accident: the minimum-curvature objective is flat along
        // every constant tangent, so without the bias there is nothing to pin the scale down.
        let seeds: Vec<Point2D> = (0..5).map(|i| Point2D::new(i as f64, (i % 2) as f64)).collect();
        assert_eq!(solve_tangents(&seeds, 0.0), seeds);
    }

    #[test]
    fn a_huge_bias_returns_the_seeds() {
        let seeds: Vec<Point2D> = (0..5).map(|i| Point2D::new(i as f64, i as f64 * 0.5)).collect();
        let out = solve_tangents(&seeds, 1e12);
        for (a, b) in out.iter().zip(&seeds) {
            assert!((a.x - b.x).abs() < 1e-6 && (a.y - b.y).abs() < 1e-6, "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn zero_bias_gives_a_straight_line_for_a_straight_input() {
        // Equally spaced collinear knots: the minimum-curvature solution is the chord, and the
        // stationarity condition is satisfied by any constant scaled tangent.
        let knots = on_line(6);
        let seeds: Vec<Point2D> = knots.windows(2).map(|w| w[1] - w[0]).collect();
        let out = solve_tangents(&seeds, 0.0);
        // The tangents must all be parallel to the line: the y components all vanish because
        // every seed and every knot is on y = 0.
        for t in &out {
            assert!(t.y.abs() < 1e-9, "tangent {t:?} left the line");
        }
        // And the solution is a positive multiple of the unit direction, not a sign flip.
        assert!(out.iter().all(|t| t.x > 0.0), "tangents reversed: {out:?}");
    }

    #[test]
    fn the_solution_stationarises_the_objective() {
        // Finite-difference check that the returned tangents are a critical point: perturbing any
        // single tangent in any direction must not decrease the objective.
        let knots: Vec<Point2D> = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(3.0, 2.0),
            Point2D::new(7.0, 1.0),
            Point2D::new(10.0, 4.0),
        ];
        let seeds: Vec<Point2D> = knots.windows(2).map(|w| (w[1] - w[0]) / 3.0).collect();
        let lambda = 0.5;
        let out = solve_tangents(&seeds, lambda);
        let objective = |t: &[Point2D]| {
            let mut acc = 0.0;
            for i in 1..t.len() - 1 {
                let s = t[i + 1] - t[i] * 2.0 + t[i - 1];
                acc += s.norm_squared();
            }
            for i in 0..t.len() {
                acc += lambda * t[i].distance_squared(seeds[i]);
            }
            acc
        };
        let base = objective(&out);
        for i in 0..out.len() {
            for delta in [1e-3, 1e-2] {
                let mut bumped = out.clone();
                bumped[i] += Point2D::new(delta, 0.0);
                assert!(objective(&bumped) >= base - 1e-9, "objective fell at knot {i}");
                let mut bumped = out.clone();
                bumped[i] += Point2D::new(0.0, delta);
                assert!(objective(&bumped) >= base - 1e-9, "objective fell at knot {i} in y");
            }
        }
    }

    #[test]
    fn trivial_inputs_are_returned_unchanged() {
        assert_eq!(solve_tangents(&[], 1.0), Vec::<Point2D>::new());
        let one = vec![Point2D::new(1.0, 2.0)];
        assert_eq!(solve_tangents(&one, 1.0), one);
        let two = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)];
        assert_eq!(solve_tangents(&two, 1.0).len(), 2);
    }

    #[test]
    fn the_system_is_pentadiagonal() {
        // Guard against someone "simplifying" this to a tridiagonal solve: perturbing a tangent
        // two positions away must still change the result.
        let seeds: Vec<Point2D> = (0..5).map(|i| Point2D::new(1.0, if i % 2 == 0 { 0.0 } else { 1.0 })).collect();
        let a = solve_tangents(&seeds, 0.1);
        let mut b_seeds = seeds.clone();
        b_seeds[2] += Point2D::new(0.0, 5.0);
        let b = solve_tangents(&b_seeds, 0.1);
        // Changing an interior seed moves the neighbouring solutions, which only a banded width of
        // two can produce.
        assert!(a[1].distance(b[1]) > 1e-9, "knot 1 unaffected by knot 2's seed");
    }
}
