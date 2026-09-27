//! Clamping the solved tangents until every control hull is in the free space.
//!
//! The solver produces tangents that minimise curvature, with no knowledge of the geometry. This
//! module brings them back inside it, and the order of operations is forced by two facts:
//!
//! * Projection is on the **tangent**, never on an individual control point. `$T_i$` is one value
//!   shared by two segments, so clamping a control point would break `$C^1$`; clamping the tangent
//!   cannot, because both segments read the same number afterwards.
//! * The loop therefore runs over the whole spline, not per segment, and it **cannot fail**. The
//!   zero tangents are always admissible — that is what "do nothing" means — and they make every
//!   segment the straight chord between its knots, which is inside the free space because the taut
//!   path is (L2). So the worst case is a piecewise-linear route, not a collision or an error.


use crate::corridor::AdmissibleTangents;
use crate::geom::clearance::Clearance;
use crate::geom::point::Point2D;

/// Clamps `tangents` until every segment's control hull is in the free space.
///
/// Returns `true` if the result is a valid spline, `false` if it fell back to the taut polyline.
///
/// # Why damping is not enough on its own
///
/// The obvious repair is to damp the tangents towards zero and re-check. It does not converge: a
/// tangent whose control point pokes a *positive* amount into an obstacle shrinks that amount
/// geometrically but never to zero, so the hull stays obstructed however many iterations run, and
/// the route ends up fully flattened when only one knot was at fault.
///
/// So the repair is targeted. When a segment's hull is not clear, the tangents at *its* two ends are
/// zeroed — that makes the segment the straight chord, which is in the free space because the taut
/// path is (L2) — and the rest of the spline is re-checked. Straight runs between two such corners
/// keep their smoothing, which is the whole point of the solve.
///
/// The loop is bounded, and the bound does not affect correctness: zeroing every remaining tangent
/// is always valid, and is the documented outcome when the bound is reached.
pub fn clamp_and_repair(
    admissible: &AdmissibleTangents,
    knots: &[Point2D],
    tangents: &mut [Point2D],
    clearance: &Clearance,
    max_iters: u32,
    dampening: f64,
) -> bool {
    project_all(admissible, tangents);
    if admissible.first_obstructed(knots, tangents, clearance).is_none() {
        return true;
    }
    // Damp once: for an overshoot of a few percent this is enough and it keeps the shape.
    for t in tangents.iter_mut() {
        *t = *t * dampening;
    }
    project_all(admissible, tangents);
    if admissible.first_obstructed(knots, tangents, clearance).is_none() {
        return true;
    }
    // Only the segments at or after `resume` can have changed, so the scan starts there.
    let mut resume = 0usize;
    for _ in 0..max_iters {
        // Find the first obstructed segment and flatten exactly that one. Stopping at the first,
        // and not re-testing segments that zeroing cannot have touched, is what keeps the repair
        // proportional to the number of *broken* segments rather than to the route's length.
        let Some(bad) = admissible.first_obstructed_from(knots, tangents, clearance, resume) else {
            return true;
        };
        tangents[bad] = Point2D::ZERO;
        tangents[bad + 1] = Point2D::ZERO;
        project_all(admissible, tangents);
        resume = bad.saturating_sub(1);
    }
    for t in tangents.iter_mut() {
        *t = Point2D::ZERO;
    }
    admissible.first_obstructed(knots, tangents, clearance).is_none()
}

/// Clamps every tangent into its knot's admissible set.
pub fn project_all(admissible: &AdmissibleTangents, tangents: &mut [Point2D]) {
    for (i, t) in tangents.iter_mut().enumerate() {
        if let Some(set) = admissible.get(i) {
            *t = set.project(*t);
        }
    }
}

