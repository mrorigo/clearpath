//! The smooth route: stages 1a through 4 of `docs/SPEC.md`.

use alloc::vec::Vec;

use crate::corridor::AdmissibleTangents;
use crate::decomp::sweep::decompose_with_guides;
use crate::error::PathPlanError;
use crate::funnel::cell_search::search;
use crate::funnel::string_pull::string_pull;
use crate::geom::point::Point2D;
use crate::router::scratch::Scratch;
use crate::router::RouteRequest;
use crate::spline::containment::clamp_and_repair;
use crate::spline::solver::solve_tangents;
use crate::spline::CubicBezierSegment;

/// Produces a `C^1` chain of cubic Bezier segments from `req.start` to `req.goal`.
///
/// The stages are the ones in the specification, and each one's output is the next one's input:
/// the margin model, the decomposition, the corridor search, the funnel, the admissible tangents,
/// the tangent solve, and the containment repair.
pub fn route_smooth(
    scratch: &mut Scratch,
    req: &RouteRequest,
) -> Result<Vec<CubicBezierSegment>, PathPlanError> {
    req.validate()?;
    let free_space = req.free_space()?;
    let clearance = req.clearance();

    // Stage 1: the inflated domain, with the endpoints' abscissae forced onto event lines so the
    // corridor can be monotone in x (section 6.3a).
    let decomp = decompose_with_guides(
        &free_space,
        req.config.margin,
        &[req.start.point.x, req.goal.point.x],
    )?;

    // Stage 2: a corridor, then the taut string through it.
    let corridor = search(&decomp, req.start.point, req.goal.point, &mut scratch.search)?;
    let path = string_pull(&decomp, &corridor, req.start.point, req.goal.point)?;

    // Stage 3: what each knot's tangent is allowed to be.
    let admissible = AdmissibleTangents::build(&decomp, &corridor.cells, &path.knots)?;

    // Stage 4: the tangents, then the repair.
    let seeds = seed_tangents(&path.knots, &req.start, &req.goal);
    let mut tangents = solve_tangents(&seeds, req.config.tangent_bias);
    clamp_and_repair(
        &admissible,
        &path.knots,
        &mut tangents,
        &clearance,
        req.config.max_repair_iters,
        req.config.repair_dampening,
    );

    admissible.control_points(&path.knots, &tangents, &clearance)
}

/// The unconstrained tangents: a third of the local chord, or a port's direction at that scale.
///
/// The factor of a third is not arbitrary — with it, a straight segment's control points sit at
/// `1/9` and `8/9` of the chord, which is the natural cubic's parameterisation and makes an
/// unrouted straight line reproduce itself.
fn seed_tangents(
    knots: &[Point2D],
    start: &crate::router::PortConstraint,
    goal: &crate::router::PortConstraint,
) -> Vec<Point2D> {
    let n = knots.len();
    let mut seeds: Vec<Point2D> = Vec::with_capacity(n);
    for i in 0..n {
        // The local chord: the longer of the two adjacent chords, so the tangent's magnitude
        // reflects the spacing it has to span.
        let incoming = if i > 0 { knots[i] - knots[i - 1] } else { Point2D::ZERO };
        let outgoing = if i + 1 < n { knots[i + 1] - knots[i] } else { Point2D::ZERO };
        let span = if incoming.length() >= outgoing.length() { incoming } else { outgoing };
        let requested = match i {
            0 => start.direction,
            j if j + 1 == n => goal.direction,
            _ => None,
        };
        seeds.push(match requested.and_then(|d| d.normalize()) {
            Some(unit) => unit * (span.length() / 3.0),
            None => span / 3.0,
        });
    }
    seeds
}
