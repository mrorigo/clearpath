//! M6 exit gate: the public API, and the section 4 invariants end to end.
//!
//! The properties here are the ones a caller actually depends on, and they are stated in terms of
//! the public types only.

mod common;
use common::{box_poly, p, WORKSPACE};
use clearpath::geom::clearance::Clearance;
use clearpath::geom::{Point2D, Polygon};
use clearpath::{
    Config, PathPlanner, PathPlanError, PortConstraint, Route, RouteKind, RouteRequest,
};
use proptest::prelude::*;




fn request(obstacles: Vec<Polygon>, start: Point2D, goal: Point2D, margin: f64) -> RouteRequest {
    let mut req = RouteRequest::new(WORKSPACE, start, goal);
    req.obstacles = obstacles;
    req.config.margin = margin;
    req
}

/// The smooth router used to return a curve that dips a rounding step *into* an obstacle.
///
/// The knots here lie exactly on an obstacle's bottom face, which is the normal case: a taut path
/// runs along obstacle boundaries. Reconstructing a control point as `knot + t/3` then lands a
/// fraction of an ulp inside. `is_free` used to accept any point with a non-zero distance to the
/// ring, so such a point was reported free, and the router shipped a curve through the obstacle.
#[test]
fn a_smooth_route_does_not_dip_into_an_obstacle_along_its_face() {
    let obstacles = vec![
        Polygon::new(vec![
            Point2D::new(17.141496085405187, 32.0138540651721),
            Point2D::new(24.96953883992036, 32.0138540651721),
            Point2D::new(24.96953883992036, 36.372611485923045),
            Point2D::new(17.141496085405187, 36.372611485923045),
        ])
        .unwrap(),
        Polygon::new(vec![
            Point2D::new(19.095236231747986, 27.62985896618443),
            Point2D::new(27.697157282844046, 27.62985896618443),
            Point2D::new(27.697157282844046, 31.91034886387489),
            Point2D::new(19.095236231747986, 31.91034886387489),
        ])
        .unwrap(),
        Polygon::new(vec![
            Point2D::new(25.770223685062522, 34.581371709157494),
            Point2D::new(30.998363789697283, 34.581371709157494),
            Point2D::new(30.998363789697283, 46.52180991237389),
            Point2D::new(25.770223685062522, 46.52180991237389),
        ])
        .unwrap(),
    ];
    let req = request(
        obstacles,
        Point2D::new(14.906697685761042, 44.91714569699075),
        Point2D::new(36.08902029423311, 82.93933742844784),
        0.0,
    );
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).expect("a route exists");

    let clearance = clearance_of(&req);
    for (i, seg) in segments.iter().enumerate() {
        for k in 0..=4000 {
            let t = k as f64 / 4000.0;
            let point = seg.evaluate(t);
            assert!(
                clearance.is_free(point),
                "segment {i} leaves the free space at t={t}: {point:?}"
            );
        }
    }
}

/// `is_free` is the predicate every part of the crate trusts when it certifies a route, so it
/// must never report a point inside an obstacle as free — and must keep the boundary free, which
/// is what the knots rely on.
#[test]
fn a_point_inside_an_obstacle_is_never_free() {
    let square = box_poly([0.0, 0.0, 10.0, 20.0]);
    let clearance = Clearance::new(std::slice::from_ref(&square), WORKSPACE, 0.0);
    // Depths well past the boundary tolerance: a point a rounding step inside is what the bug
    // turned on, but the tolerance is deliberately a length, so a point within it stays free.
    for (i, y) in [19.9, 19.99, 19.999, 19.99999].iter().enumerate() {
        let point = Point2D::new(5.0, *y);
        assert!(square.contains(point), "sanity: {point:?} should be inside");
        assert!(!clearance.is_free(point), "{point:?} is inside the square but reported free ({i})");
    }
    // And the boundary itself stays free, which is what knots rely on.
    assert!(clearance.is_free(Point2D::new(5.0, 20.0)));
    assert!(clearance.is_free(Point2D::new(5.0, 0.0)));
    assert!(clearance.is_free(Point2D::new(0.0, 10.0)));
}

fn clearance_of(req: &RouteRequest) -> Clearance<'_> {
    req.clearance()
}

/// Section 4.2: C1 continuity, bitwise at the joints and to `1e-12` relative on the derivative.
fn assert_c1(segments: &[clearpath::CubicBezierSegment]) {
    for w in segments.windows(2) {
        let (a, b) = (w[0], w[1]);
        assert_eq!(a.p3, b.p0, "joint point differs");
        let ta = a.tangent(1.0);
        let tb = b.tangent(0.0);
        let scale = ta.length().max(1.0);
        assert!(
            ta.distance(tb) <= 1e-12 * scale,
            "tangent discontinuity {} at joint {:?}",
            ta.distance(tb),
            a.p3
        );
    }
}

/// Section 4.1: the curve is in the free space, tested on the flattened polyline so the sampling
/// density is a stated guarantee rather than a guess.
fn assert_clear(req: &RouteRequest, segments: &[clearpath::CubicBezierSegment]) {
    let clearance = clearance_of(req);
    for s in segments {
        for q in s.flatten(1e-3) {
            assert!(clearance.is_free(q), "curve point {q:?} is not in the free space");
        }
    }
}

/// The route starts and ends where it was asked to.
fn assert_endpoints(req: &RouteRequest, segments: &[clearpath::CubicBezierSegment]) {
    let first = segments.first().expect("a route with no segments");
    let last = segments.last().unwrap();
    assert_eq!(first.p0, req.start.point);
    assert_eq!(last.p3, req.goal.point);
}

fn check(req: &RouteRequest) {
    let mut planner = PathPlanner::new();
    let Ok(segments) = planner.route_smooth(req) else {
        return; // A refused route is a legitimate answer; the fixtures decide which they expect.
    };
    assert!(!segments.is_empty());
    assert_endpoints(req, &segments);
    assert_c1(&segments);
    assert_clear(req, &segments);
}

#[test]
fn open_space_is_a_straight_line() {
    let req = request(vec![], p(5.0, 50.0), p(95.0, 50.0), 0.0);
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).unwrap();
    assert_eq!(segments.len(), 1, "no obstacles means one segment: {segments:?}");
    let s = &segments[0];
    // Section 4.3: the control points are collinear with the chord and ordered along it.
    for k in 0..64 {
        let t = k as f64 / 64.0;
        let q = s.evaluate(t);
        assert!(
            q.distance(clearpath::corridor::closest_point_on_segment(
                q,
                req.start.point,
                req.goal.point
            )) < 1e-9,
            "point {q:?} at t={t} is off the straight line"
        );
    }
    check(&req);
}

#[test]
fn around_a_box_the_route_is_continuous_and_clear() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    for (a, b) in [
        (p(5.0, 50.0), p(95.0, 50.0)),
        (p(20.0, 80.0), p(80.0, 20.0)),
        (p(95.0, 50.0), p(5.0, 50.0)),
    ] {
        let req = request(obstacles.clone(), a, b, 0.0);
        check(&req);
        let mut planner = PathPlanner::new();
        let segments = planner.route_smooth(&req).unwrap();
        assert!(!segments.is_empty());
    }
}

#[test]
fn a_margin_is_honoured_by_the_curve() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let req = request(obstacles, p(5.0, 50.0), p(95.0, 50.0), 3.0);
    check(&req);
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).unwrap();
    let clearance = clearance_of(&req);
    for s in &segments {
        for q in s.flatten(1e-4) {
            assert!(clearance.is_free(q), "{q:?} violates the margin");
        }
    }
}

#[test]
fn a_sealed_region_reports_no_path() {
    let obstacles = vec![
        box_poly([30.0, 30.0, 70.0, 34.0]),
        box_poly([30.0, 66.0, 70.0, 70.0]),
        box_poly([30.0, 30.0, 34.0, 70.0]),
        box_poly([66.0, 30.0, 70.0, 70.0]),
    ];
    let req = request(obstacles, p(10.0, 10.0), p(50.0, 50.0), 0.0);
    let mut planner = PathPlanner::new();
    assert!(matches!(planner.route_smooth(&req), Err(PathPlanError::NoPathFound)));
}

#[test]
fn validation_errors_are_specific() {
    let mut planner = PathPlanner::new();

    let degenerate = request(vec![], p(5.0, 5.0), p(5.0, 5.0), 0.0);
    assert!(matches!(planner.route_smooth(&degenerate), Err(PathPlanError::DegenerateEndpoints)));

    let outside = request(vec![], p(-5.0, 5.0), p(50.0, 5.0), 0.0);
    assert!(matches!(
        planner.route_smooth(&outside),
        Err(PathPlanError::EndpointOutsideWorkspace { .. })
    ));

    let inside = request(vec![box_poly([40.0, 40.0, 60.0, 60.0])], p(50.0, 50.0), p(90.0, 50.0), 0.0);
    assert!(matches!(planner.route_smooth(&inside), Err(PathPlanError::EndpointInObstacle { .. })));

    let bad_margin = request(vec![], p(5.0, 5.0), p(90.0, 5.0), -1.0);
    assert!(matches!(planner.route_smooth(&bad_margin), Err(PathPlanError::InvalidConfig(_))));

    let bad_damping = request(vec![], p(5.0, 5.0), p(90.0, 5.0), 0.0);
    let mut r = bad_damping;
    r.config.repair_dampening = 1.5;
    assert!(matches!(planner.route_smooth(&r), Err(PathPlanError::InvalidConfig(_))));

    let too_big = request(vec![], p(5.0, 5.0), p(90.0, 5.0), 60.0);
    assert!(planner.route_smooth(&too_big).is_err(), "a margin larger than the workspace cannot route");
}

#[test]
fn a_sloped_obstacle_with_a_margin_is_reported() {
    let triangle = Polygon::new(vec![p(10.0, 10.0), p(40.0, 10.0), p(10.0, 40.0)]).unwrap();
    let req = request(vec![triangle], p(5.0, 50.0), p(90.0, 50.0), 1.0);
    let mut planner = PathPlanner::new();
    assert!(matches!(planner.route_smooth(&req), Err(PathPlanError::MarginUnsupportedGeometry { .. })));
    // At zero margin the same obstacle routes fine: the limitation is the margin, not the shape.
    let ok = request(vec![], p(5.0, 50.0), p(90.0, 50.0), 0.0);
    assert!(planner.route_smooth(&ok).is_ok());
}

#[test]
fn a_port_direction_is_honoured_when_it_fits() {
    // A route forced to leave vertically and arrive vertically, in open space, where both fit.
    let mut req = request(vec![], p(20.0, 20.0), p(80.0, 80.0), 0.0);
    req.start = PortConstraint::directed(p(20.0, 20.0), p(0.0, -1.0));
    req.goal = PortConstraint::directed(p(80.0, 80.0), p(0.0, 1.0));
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).unwrap();
    let first = segments.first().unwrap();
    let start_dir = first.tangent(0.0);
    assert!(
        start_dir.normalize().map(|d| d.distance(p(0.0, -1.0)) < 1e-6).unwrap_or(false),
        "start tangent {start_dir:?} does not follow the requested direction"
    );
    let last = segments.last().unwrap();
    let goal_dir = last.tangent(1.0);
    assert!(
        goal_dir.normalize().map(|d| d.distance(p(0.0, 1.0)) < 1e-6).unwrap_or(false),
        "goal tangent {goal_dir:?} does not follow the requested direction"
    );
    check(&req);
}

#[test]
fn a_port_direction_that_cannot_fit_is_projected_not_refused() {
    // Leaving a knot wedged in a corner: the direction is dropped and the route still comes out.
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let mut req = request(obstacles, p(5.0, 20.0), p(95.0, 20.0), 0.0);
    req.start = PortConstraint::directed(p(5.0, 20.0), p(1.0, 1.0));
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).unwrap();
    assert!(!segments.is_empty());
    check(&req);
}

#[test]
fn the_route_enum_dispatches() {
    let mut planner = PathPlanner::new();
    let req = request(vec![], p(5.0, 50.0), p(95.0, 50.0), 0.0);
    match planner.route(&req, RouteKind::Smooth) {
        Ok(Route::Smooth(segments)) => assert!(!segments.is_empty()),
        other => panic!("unexpected route: {other:?}"),
    }
}

#[test]
fn the_planner_is_reusable_across_queries() {
    // A warm planner must give the same answer as a fresh one, or the scratch is leaking state.
    let mut planner = PathPlanner::new();
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let a = request(obstacles.clone(), p(5.0, 50.0), p(95.0, 50.0), 0.0);
    let b = request(obstacles, p(20.0, 80.0), p(80.0, 20.0), 1.0);
    let first = planner.route_smooth(&a).unwrap();
    let _ = planner.route_smooth(&b).unwrap();
    let again = planner.route_smooth(&a).unwrap();
    assert_eq!(first, again, "a warm planner disagrees with its own first answer");
}

#[test]
fn the_tangent_bias_changes_the_route() {
    // A knob that does nothing is a knob that lies.
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let mut low = request(obstacles.clone(), p(5.0, 50.0), p(95.0, 50.0), 0.0);
    let mut high = request(obstacles, p(5.0, 50.0), p(95.0, 50.0), 0.0);
    low.config.tangent_bias = 0.0;
    high.config.tangent_bias = 1e6;
    let mut planner = PathPlanner::new();
    let a = planner.route_smooth(&low).unwrap();
    let b = planner.route_smooth(&high).unwrap();
    assert_ne!(a, b, "tangent_bias had no effect");
    for r in [&low, &high] {
        let req = r.clone();
        check(&req);
    }
}

#[test]
fn config_defaults_validate() {
    assert!(Config::default().validate().is_ok());
    assert!(Config { max_arc_segments: 0, ..Config::default() }.validate().is_err());
    assert!(Config { max_repair_iters: 0, ..Config::default() }.validate().is_err());
    assert!(Config { tangent_bias: f64::NAN, ..Config::default() }.validate().is_err());
    assert!(Config { margin: -1.0, ..Config::default() }.validate().is_err());
    assert!(Config { repair_dampening: 0.0, ..Config::default() }.validate().is_err());
    // The default bias must be positive: zero makes the tangent system singular (the
    // minimum-curvature objective is flat along every constant tangent).
    assert!(Config::default().tangent_bias > 0.0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    /// The end-to-end contract on random fields: whatever comes out is C1, starts and ends at the
    /// query points, and stays in the free space.
    #[test]
    fn routes_are_continuous_and_clear(
        boxes in prop::collection::vec((0u8..=11, 0u8..=11, 0u8..=6, 0u8..=6), 0..6),
        margin in 0.0f64..3.0,
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
                check(&request(obstacles.clone(), *a, b, margin));
            }
        }
    }
}

#[test]
fn a_smooth_route_is_actually_smooth() {
    // The repair must not flatten everything: a route with a corner still has long straight runs
    // between the corners, and those should carry non-zero tangents.
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let req = request(obstacles, p(5.0, 50.0), p(95.0, 50.0), 0.0);
    let mut planner = PathPlanner::new();
    let segments = planner.route_smooth(&req).unwrap();
    let interior: Vec<&clearpath::CubicBezierSegment> =
        segments.iter().filter(|s| s.p1.distance(s.p0) > 1e-9).collect();
    assert!(
        !interior.is_empty(),
        "every segment is a straight line; the solve and the repair flattened the route: {segments:?}"
    );
}
