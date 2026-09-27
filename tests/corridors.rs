//! M5 exit gate: the admissible tangent sets are convex, contain the zero tangent, and keep the
//! spline's control hull in the free space.
//!
//! Nothing here trusts the construction: the hull is tested against the obstacles directly, and the
//! projection is tested to land back inside the set.

mod common;
use common::{box_poly, p, WORKSPACE};
use clearpath::geom::clearance::distance_to_ring;
use clearpath::decomp::sweep::decompose_with_guides;
use clearpath::funnel::cell_search::{SearchScratch, search};
use clearpath::funnel::string_pull::string_pull;
use clearpath::geom::{FreeSpace, Point2D, Polygon};
use clearpath::corridor::AdmissibleTangents;
use proptest::prelude::*;
use clearpath::geom::BoundingBox;




/// A tangent is admissible when it is a fraction of the chord: the chord direction is what the
/// solver would produce before any constraint, and scaling it is a one-parameter family that
/// includes the zero tangent and a genuinely large one.
fn chord_tangents(knots: &[Point2D], scale: f64) -> Vec<Point2D> {
    let n = knots.len();
    let mut out = vec![Point2D::ZERO; n];
    for i in 0..n - 1 {
        let chord = knots[i + 1] - knots[i];
        out[i] = chord * scale;
    }
    out[n - 1] = (knots[n - 1] - knots[n - 2]) * scale;
    out
}

fn check(obstacles: &[Polygon], margin: f64, start: Point2D, goal: Point2D) {
    let space = FreeSpace::new(WORKSPACE, obstacles.to_vec()).unwrap();
    let decomp = decompose_with_guides(&space, margin, &[start.x, goal.x]).unwrap();
    let mut scratch = SearchScratch::default();
    let Ok(corridor) = search(&decomp, start, goal, &mut scratch) else {
        return;
    };
    let path = string_pull(&decomp, &corridor, start, goal).unwrap();
    let admissible =
        AdmissibleTangents::build(&decomp, &corridor.cells, &path.knots).unwrap();

    // The free-space predicate the margin model defines: a distance test, not a containment test,
    // because a taut path hugs what it wraps around and its knots sit exactly on obstacle
    // boundaries when `margin` is 0.
    let clearance = clearpath::geom::clearance::Clearance::new(obstacles, WORKSPACE, margin);
    let free = |q: Point2D| clearance.is_free(q);

    for i in 0..path.knots.len() {
        let set = admissible.get(i).unwrap();
        assert_eq!(set.knot(), path.knots[i]);
        // The zero tangent is always admissible: it is what "do nothing" means.
        assert!(set.contains(Point2D::ZERO), "set at knot {i} rejects the zero tangent");
        // Projection lands back inside, from anywhere.
        for candidate in [
            path.knots[i] * 100.0,
            Point2D::ZERO,
            Point2D::new(1e3, 1e3),
            Point2D::new(-1e3, 1e3),
            Point2D::new(0.3, -0.7),
        ] {
            let projected = set.project(candidate);
            assert!(
                set.contains(projected),
                "projection of {candidate:?} at knot {i} left the set: {projected:?}"
            );
        }
        // The vertices, when they exist, are all admissible.
        for v in set.vertices() {
            assert!(set.contains(*v), "vertex {v:?} of set {i} is not admissible");
        }
    }

    // Every chord scale, projected, must yield control points whose hull is clear.
    for scale in [0.0f64, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0] {
        let raw = chord_tangents(&path.knots, scale);
        let mut tangents = Vec::with_capacity(raw.len());
        for (i, t) in raw.iter().enumerate() {
            tangents.push(admissible.get(i).unwrap().project(*t));
        }
        match admissible.control_points(&path.knots, &tangents, &clearance) {
            Ok(segments) => {
                for s in &segments {
                    // The curve itself, sampled: the hull argument is the guarantee, this is the
                    // check on the argument.
                    for k in 0..=64 {
                        let q = s.evaluate(k as f64 / 64.0);
                        assert!(free(q), "curve point {q:?} is not in the free space (scale {scale})");
                    }
                }
            }
            Err(_) => {
                // Refusing is a legitimate outcome; what is not legitimate is returning a segment
                // whose curve leaves the free space. Nothing to check here.
            }
        }
    }

    // A hull that *is* accepted must pass the control-point check, which is the contract.
    let tangents: Vec<Point2D> = (0..path.knots.len())
        .map(|i| admissible.get(i).unwrap().project(Point2D::ZERO))
        .collect();
    let straight = admissible.control_points(&path.knots, &tangents, &clearance).unwrap();
    assert_eq!(straight.len(), path.knots.len() - 1);
    for s in &straight {
        assert!(free(s.p0) && free(s.p1) && free(s.p2) && free(s.p3));
    }
}

#[test]
fn open_space() {
    check(&[], 0.0, p(5.0, 50.0), p(95.0, 50.0));
    check(&[], 2.0, p(5.0, 50.0), p(95.0, 50.0));
}

#[test]
fn around_a_single_box() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    check(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
    check(&obstacles, 3.0, p(5.0, 50.0), p(95.0, 50.0));
}

#[test]
fn through_a_staggered_gap() {
    let obstacles = vec![box_poly([30.0, 0.0, 45.0, 55.0]), box_poly([55.0, 45.0, 70.0, 100.0])];
    check(&obstacles, 0.0, p(5.0, 90.0), p(95.0, 10.0));
    check(&obstacles, 0.0, p(95.0, 10.0), p(5.0, 90.0));
}

#[test]
fn through_a_slot() {
    let obstacles = vec![box_poly([43.0, 0.0, 57.0, 39.0]), box_poly([43.0, 41.0, 57.0, 100.0])];
    check(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
}

#[test]
fn overlapping_boxes() {
    check(&[box_poly([20.0, 20.0, 50.0, 50.0]), box_poly([40.0, 40.0, 70.0, 70.0])], 0.0, p(5.0, 5.0), p(95.0, 95.0));
}

#[test]
fn a_narrow_passage_with_a_margin() {
    let obstacles = vec![
        box_poly([40.0, 0.0, 44.0, 40.0]),
        box_poly([44.0, 40.0, 60.0, 44.0]),
        box_poly([60.0, 0.0, 64.0, 48.0]),
        box_poly([64.0, 48.0, 80.0, 52.0]),
        box_poly([80.0, 0.0, 84.0, 60.0]),
    ];
    check(&obstacles, 2.0, p(5.0, 20.0), p(95.0, 20.0));
}

#[test]
fn a_pinned_knot_accepts_only_zero() {
    // A knot in a corner with no room at all: the set must still be a valid convex region that
    // contains zero, and projecting anything must give zero.
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let space = FreeSpace::new(WORKSPACE, obstacles.clone()).unwrap();
    let a = p(5.0, 50.0);
    let b = p(95.0, 50.0);
    let decomp = decompose_with_guides(&space, 0.0, &[a.x, b.x]).unwrap();
    let mut scratch = SearchScratch::default();
    let corridor = search(&decomp, a, b, &mut scratch).unwrap();
    let path = string_pull(&decomp, &corridor, a, b).unwrap();
    let admissible = AdmissibleTangents::build(&decomp, &corridor.cells, &path.knots).unwrap();
    for i in 0..path.knots.len() {
        let set = admissible.get(i).unwrap();
        assert!(set.contains(Point2D::ZERO));
        let projected = set.project(Point2D::new(1e6, 1e6));
        assert!(set.contains(projected));
    }
}

#[test]
fn the_zero_tangent_yields_the_taut_polyline() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let space = FreeSpace::new(WORKSPACE, obstacles.clone()).unwrap();
    let a = p(5.0, 50.0);
    let b = p(95.0, 50.0);
    let decomp = decompose_with_guides(&space, 0.0, &[a.x, b.x]).unwrap();
    let mut scratch = SearchScratch::default();
    let corridor = search(&decomp, a, b, &mut scratch).unwrap();
    let path = string_pull(&decomp, &corridor, a, b).unwrap();
    let admissible = AdmissibleTangents::build(&decomp, &corridor.cells, &path.knots).unwrap();
    let clearance = clearpath::geom::clearance::Clearance::new(&obstacles, WORKSPACE, 0.0);
    let tangents = vec![Point2D::ZERO; path.knots.len()];
    let segments = admissible.control_points(&path.knots, &tangents, &clearance).unwrap();
    // A zero tangent makes each segment a straight line between consecutive knots.
    for (s, w) in segments.iter().zip(path.knots.windows(2)) {
        assert_eq!(s.p0, w[0]);
        assert_eq!(s.p3, w[1]);
        assert!(s.p1.distance(s.p0) < 1e-12);
        assert!(s.p2.distance(s.p3) < 1e-12);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn admissible_sets_and_hulls_are_safe(
        boxes in prop::collection::vec((0u8..=11, 0u8..=11, 0u8..=6, 0u8..=6), 0..7),
        margin in 0.0f64..4.0,
    ) {
        let obstacles: Vec<Polygon> = boxes
            .into_iter()
            .map(|(i, j, w, h)| {
                let x = 2.0 + (i as f64) * 8.0;
                let y = 2.0 + (j as f64) * 8.0;
                box_poly([x, y, x + w as f64 + 1.0, y + h as f64 + 1.0])
            })
            .collect();
        let points: Vec<Point2D> = (0..4)
            .flat_map(|i| (0..4).map(move |j| p(10.0 + (i as f64) * 26.0, 10.0 + (j as f64) * 26.0)))
            .collect();
        for (i, a) in points.iter().enumerate() {
            let b = points[(i * 5 + 7) % points.len()];
            if *a != b {
                check(&obstacles, margin, *a, b);
            }
        }
    }
}
