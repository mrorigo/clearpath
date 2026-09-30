//! M7 exit gate: the rectilinear router.
//!
//! The contract is narrow on purpose: an axis-aligned polyline from the start to the goal, in the
//! free space. Everything else — optimality, aesthetics — is not claimed.

mod common;
use common::{box_poly, p, WORKSPACE};
use clearpath::geom::{Point2D, Polygon};
use clearpath::{Config, PathPlanError, PathPlanner, PortConstraint, Route, RouteKind, RouteRequest};
use proptest::prelude::*;




fn request(obstacles: Vec<Polygon>, start: Point2D, goal: Point2D, margin: f64) -> RouteRequest {
    let mut req = RouteRequest::new(WORKSPACE, start, goal);
    req.obstacles = obstacles;
    req.config.margin = margin;
    req
}

fn check(req: &RouteRequest) {
    let mut planner = PathPlanner::new();
    let Ok(poly) = planner.route_orthogonal(req) else {
        return;
    };
    let clearance = req.clearance();

    assert!(poly.points.len() >= 2, "a route needs at least two vertices: {poly:?}");
    assert_eq!(poly.points.first(), Some(&req.start.point), "does not start at the start");
    assert_eq!(poly.points.last(), Some(&req.goal.point), "does not end at the goal");
    assert!(poly.is_axis_aligned(), "not rectilinear: {:?}", poly.points);

    // No repeated vertices, and no zero-length segments.
    for w in poly.points.windows(2) {
        assert!(w[0] != w[1], "zero-length segment in {poly:?}");
    }
    for i in 0..poly.points.len() {
        assert!(
            !poly.points[..i].contains(&poly.points[i]),
            "repeated vertex in {poly:?}"
        );
    }

    // Every point of every segment is in the free space. Sampled, and the corner points are
    // checked exactly too.
    for w in poly.points.windows(2) {
        assert!(clearance.is_free(w[0]), "vertex {:?} is not free", w[0]);
        assert!(clearance.is_free(w[1]), "vertex {:?} is not free", w[1]);
        for k in 1..16 {
            let q = w[0].lerp(w[1], k as f64 / 16.0);
            assert!(clearance.is_free(q), "point {q:?} on {:?}..{:?} is not free", w[0], w[1]);
        }
    }
}

#[test]
fn open_space_is_an_l() {
    // A rectilinear route between two points in an empty box: one corner, no more.
    let req = request(vec![], p(10.0, 20.0), p(80.0, 60.0), 0.0);
    let mut planner = PathPlanner::new();
    let poly = planner.route_orthogonal(&req).unwrap();
    assert!(poly.points.len() <= 3, "expected at most one corner, got {:?}", poly.points);
    assert!(poly.is_axis_aligned());
    check(&req);
}

/// The router's verification used to sample each leg at nine points, so a blocked stretch
/// narrower than an eighth of the leg fell between two samples and the route was returned
/// straight through the obstacle. These are the inputs that produced it, found by fuzzing.
#[test]
fn a_route_never_crosses_a_non_rectangular_obstacle() {
    // A triangle and a quadrilateral: neither is its own bounding box, so a face of each
    // crosses the interior of a grid cell whose centre the router certifies as free.
    let obstacles = vec![
        Polygon::new(vec![
            Point2D::new(62.2985272601378, 76.96932122183367),
            Point2D::new(58.65098857647418, 80.12628183000274),
            Point2D::new(57.04354354209111, 76.24987005449628),
            Point2D::new(59.80890537115175, 71.66867667327912),
        ])
            .unwrap(),
        Polygon::new(vec![
            Point2D::new(37.939392606205324, 71.73262863860343),
            Point2D::new(38.26792768762064, 80.92843389733514),
            Point2D::new(31.673446815380757, 72.77132666893101),
        ])
        .unwrap(),
    ];
    for (a, b) in [
        (
            Point2D::new(78.32082755447614, 6.355276226875372),
            Point2D::new(39.83035360063738, 76.9996360975841),
        ),
        (
            Point2D::new(17.821504441688013, 33.06899301462444),
            Point2D::new(78.21979712547615, 91.5848254339149),
        ),
    ] {
        let req = request(obstacles.clone(), a, b, 0.0);
        let mut planner = PathPlanner::new();
        let Ok(poly) = planner.route_orthogonal(&req) else {
            continue;
        };
        // The contract: every returned route is in the free space. `check` samples, so this
        // asserts the property directly against the exact predicate for each leg.
        let clearance = req.clearance();
        for w in poly.points.windows(2) {
            assert!(
                clearance.hull_is_free(w),
                "leg {:?}..{:?} is not free, in route {poly:?}",
                w[0],
                w[1]
            );
        }
        check(&req);
    }
}

/// The same defect, reduced to a shape a reader can check by eye: a thin sliver of an
/// obstacle that pokes into a cell the grid certified on its centre.
#[test]
fn a_thin_protrusion_between_samples_is_still_refused() {
    // A narrow triangle whose apex is the only part inside the bounding box, so the grid
    // places no line at the apex and no sample of the long leg lands on it.
    let tri = Polygon::new(vec![p(49.0, 49.0), p(51.0, 49.0), p(50.0, 51.0)]).unwrap();
    let req = request(vec![tri], p(20.0, 20.0), p(80.0, 80.0), 0.0);
    let mut planner = PathPlanner::new();
    if let Ok(poly) = planner.route_orthogonal(&req) {
        let clearance = req.clearance();
        for w in poly.points.windows(2) {
            assert!(clearance.hull_is_free(w), "leg {:?}..{:?} is not free", w[0], w[1]);
        }
    }
}

#[test]
fn around_a_box() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    for (a, b) in [
        (p(5.0, 50.0), p(95.0, 50.0)),
        (p(5.0, 5.0), p(95.0, 95.0)),
        (p(20.0, 80.0), p(80.0, 20.0)),
    ] {
        let req = request(obstacles.clone(), a, b, 0.0);
        let mut planner = PathPlanner::new();
        assert!(planner.route_orthogonal(&req).is_ok(), "no route from {a:?} to {b:?}");
        check(&req);
    }
}

#[test]
fn through_a_staggered_gap() {
    let obstacles = vec![box_poly([30.0, 0.0, 45.0, 55.0]), box_poly([55.0, 45.0, 70.0, 100.0])];
    let req = request(obstacles, p(5.0, 90.0), p(95.0, 10.0), 0.0);
    let mut planner = PathPlanner::new();
    assert!(planner.route_orthogonal(&req).is_ok());
    check(&req);
}

#[test]
fn a_margin_is_honoured() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let req = request(obstacles, p(5.0, 50.0), p(95.0, 50.0), 3.0);
    check(&req);
}

#[test]
fn a_sealed_region_reports_no_path() {
    // A ring with a hole, with the hole's interior reachable only through the ring itself.
    let obstacles = vec![
        box_poly([30.0, 30.0, 70.0, 34.0]),
        box_poly([30.0, 66.0, 70.0, 70.0]),
        box_poly([30.0, 30.0, 34.0, 70.0]),
        box_poly([66.0, 30.0, 70.0, 70.0]),
    ];
    let req = request(obstacles, p(10.0, 10.0), p(50.0, 50.0), 0.0);
    let mut planner = PathPlanner::new();
    assert!(matches!(planner.route_orthogonal(&req), Err(PathPlanError::NoPathFound)));
}

#[test]
fn a_diagonal_only_passage_is_closed() {
    // Two boxes whose corners touch leave no 4-connected opening, so a rectilinear route must
    // fail rather than squeeze through the point of contact. A smooth route may also fail here;
    // what matters is that neither claims a route that does not exist.
    let obstacles = vec![box_poly([40.0, 0.0, 60.0, 50.0]), box_poly([50.0, 50.0, 70.0, 100.0])];
    let req = request(obstacles, p(5.0, 5.0), p(95.0, 95.0), 0.0);
    let mut planner = PathPlanner::new();
    match planner.route_orthogonal(&req) {
        Ok(poly) => {
            // If a route exists, it must go around, so it has to leave the diagonal band.
            assert!(poly.is_axis_aligned());
            check(&req);
        }
        Err(PathPlanError::NoPathFound) => {}
        Err(other) => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn a_sloped_obstacle_with_a_margin_is_reported() {
    let triangle = Polygon::new(vec![p(10.0, 10.0), p(40.0, 10.0), p(10.0, 40.0)]).unwrap();
    let req = request(vec![triangle], p(5.0, 50.0), p(90.0, 50.0), 1.0);
    let mut planner = PathPlanner::new();
    assert!(matches!(
        planner.route_orthogonal(&req),
        Err(PathPlanError::MarginUnsupportedGeometry { .. })
    ));
}

#[test]
fn validation_errors_match_the_smooth_router() {
    // Both routers share the request validation, so they must report the same thing.
    let mut planner = PathPlanner::new();
    let degenerate = request(vec![], p(5.0, 5.0), p(5.0, 5.0), 0.0);
    assert!(matches!(planner.route_orthogonal(&degenerate), Err(PathPlanError::DegenerateEndpoints)));
    let outside = request(vec![], p(-5.0, 5.0), p(50.0, 5.0), 0.0);
    assert!(matches!(
        planner.route_orthogonal(&outside),
        Err(PathPlanError::EndpointOutsideWorkspace { .. })
    ));
    let bad = RouteRequest { config: Config { margin: -1.0, ..Config::default() }, ..request(vec![], p(5.0, 5.0), p(50.0, 5.0), 0.0) };
    assert!(matches!(planner.route_orthogonal(&bad), Err(PathPlanError::InvalidConfig(_))));
}

#[test]
fn the_route_enum_dispatches_to_the_orthogonal_branch() {
    let mut planner = PathPlanner::new();
    let req = request(vec![], p(10.0, 20.0), p(80.0, 60.0), 0.0);
    match planner.route(&req, RouteKind::Orthogonal) {
        Ok(Route::Orthogonal(poly)) => assert!(poly.is_axis_aligned()),
        other => panic!("unexpected route: {other:?}"),
    }
}

#[test]
fn a_port_direction_is_ignored_by_the_rectilinear_router() {
    // A rectilinear segment is axis-aligned, so an arbitrary requested direction cannot be
    // honoured. It is documented as not applying, and the route still comes out.
    let mut req = request(vec![], p(20.0, 20.0), p(80.0, 80.0), 0.0);
    req.start = PortConstraint::directed(p(20.0, 20.0), p(1.0, 1.0));
    let mut planner = PathPlanner::new();
    let poly = planner.route_orthogonal(&req).unwrap();
    assert!(poly.is_axis_aligned());
    check(&req);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn rectilinear_routes_are_axis_aligned_and_clear(
        boxes in prop::collection::vec((0u8..=10, 0u8..=10, 0u8..=6, 0u8..=6), 0..6),
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
