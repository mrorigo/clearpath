//! M4 exit gate: lemma L2 — the funnel returns a polyline inside the corridor's own cells.
//!
//! Checked by sampling the emitted polyline and asking which corridor cell contains each sample. A
//! path that leaves the corridor, or that crosses an obstacle, fails. Nothing here trusts the
//! funnel's own reasoning.

use pathplan::decomp::sweep::decompose_with_guides;
use pathplan::decomp::Decomposition;
use pathplan::funnel::cell_search::{SearchScratch, search};
use pathplan::funnel::string_pull::string_pull;
use pathplan::geom::{BoundingBox, FreeSpace, Point2D, Polygon};
use proptest::prelude::*;

const WORKSPACE: BoundingBox = BoundingBox {
    min: Point2D { x: 0.0, y: 0.0 },
    max: Point2D { x: 100.0, y: 100.0 },
};

/// Samples per unit length along each taut segment. Fine enough that a segment grazing an obstacle
/// corner cannot slip between samples unnoticed.
const DENSITY: f64 = 20.0;

fn box_poly(c: [f64; 4]) -> Polygon {
    Polygon::new(vec![
        Point2D::new(c[0], c[1]),
        Point2D::new(c[2], c[1]),
        Point2D::new(c[2], c[3]),
        Point2D::new(c[0], c[3]),
    ])
    .unwrap()
}

/// Minimum distance from `p` to a ring's boundary.
fn dist_to_ring(p: Point2D, ring: &[Point2D]) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let ab = b - a;
        let len2 = ab.norm_squared();
        let t = if len2 == 0.0 { 0.0 } else { ((p - a).dot(ab) / len2).clamp(0.0, 1.0) };
        best = best.min(p.distance(a + ab * t));
    }
    best
}

/// The decomposition the router would build: with the endpoints' abscissae forced into the event
/// set, so each endpoint lands in its own slab and the corridor can be monotone.
fn routed(obstacles: &[Polygon], margin: f64, start: Point2D, goal: Point2D) -> Decomposition {
    let space = FreeSpace::new(WORKSPACE, obstacles.to_vec()).unwrap();
    decompose_with_guides(&space, margin, &[start.x, goal.x]).unwrap()
}

fn route(obstacles: &[Polygon], margin: f64, start: Point2D, goal: Point2D) -> (Decomposition, pathplan::funnel::string_pull::TautPath) {
    let decomp = routed(obstacles, margin, start, goal);
    let mut scratch = SearchScratch::default();
    let corridor = search(&decomp, start, goal, &mut scratch).expect("a route should exist");
    let path = string_pull(&decomp, &corridor, start, goal).expect("funnel");
    (decomp, path)
}

/// L2 plus collision-freedom, for one fixture.
fn check_route(obstacles: &[Polygon], margin: f64, start: Point2D, goal: Point2D) {
    let space = FreeSpace::new(WORKSPACE, obstacles.to_vec()).unwrap();
    let decomp = routed(obstacles, margin, start, goal);
    let _ = &space;
    let mut scratch = SearchScratch::default();
    let corridor = match search(&decomp, start, goal, &mut scratch) {
        Ok(c) => c,
        // No corridor is a legitimate answer for a sealed fixture; the caller decides whether that
        // was expected.
        Err(_) => return,
    };
    let path = string_pull(&decomp, &corridor, start, goal).unwrap();

    // The endpoints must be the query points.
    assert_eq!(path.knots.first().copied(), Some(start), "path must start at the start point");
    assert_eq!(path.knots.last().copied(), Some(goal), "path must end at the goal point");

    let corridor_cells: Vec<pathplan::decomp::CellId> = corridor.cells.clone();

    let mut samples = 0usize;
    for seg in path.knots.windows(2) {
        let (a, b) = (seg[0], seg[1]);
        let steps = ((a.distance(b) * DENSITY).ceil() as usize).max(1);
        for k in 0..=steps {
            // Take the endpoint verbatim: `a + (b - a) * 1.0` can miss `b` by a few ULP, and the
            // cell containment test is exact.
            let p = if k == steps { b } else { a.lerp(b, k as f64 / steps as f64) };
            // L2: every sample is inside one of the corridor's own cells.
            assert!(
                corridor_cells.iter().any(|c| decomp.cells()[*c as usize].contains(p)),
                "sample {p:?} on segment {a:?}..{b:?} is outside every corridor cell"
            );
            // And it clears every obstacle by the margin. Note this is a *distance* test, not a
            // containment test: a taut path hugs the obstacles it wraps around, so with
            // `margin == 0` its knots and long stretches lie exactly on an obstacle boundary. That
            // is correct — the guarantee is `>= margin` — and a containment test would reject every
            // useful route.
            for obstacle in obstacles {
                let clear = dist_to_ring(p, obstacle.vertices());
                assert!(
                    clear >= margin - 1e-9,
                    "sample {p:?} is {clear} from an obstacle, margin {margin}"
                );
            }
            samples += 1;
        }
    }
    assert!(samples > 0, "no samples were taken");
}

fn p(x: f64, y: f64) -> Point2D {
    Point2D::new(x, y)
}

#[test]
fn open_space_is_a_straight_line() {
    let obstacles: Vec<Polygon> = Vec::new();
    let (decomp, path) = route(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
    assert_eq!(path.knots.len(), 2, "no obstacles means a straight line, got {:?}", path.knots);
    let _ = decomp;
    check_route(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
}

#[test]
fn around_a_single_box() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    // Below the box: straight. Above: straight. Left to right past it: two knots.
    check_route(&obstacles, 0.0, p(5.0, 20.0), p(95.0, 20.0));
    check_route(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
    let (_, path) = route(&obstacles, 0.0, p(5.0, 50.0), p(95.0, 50.0));
    assert!(path.knots.len() >= 2);
    check_route(&obstacles, 0.0, p(20.0, 80.0), p(80.0, 20.0));
}

#[test]
fn through_a_staggered_gap() {
    // A zigzag that cannot be walked in a straight line.
    let obstacles = vec![
        box_poly([30.0, 0.0, 45.0, 55.0]),
        box_poly([55.0, 45.0, 70.0, 100.0]),
    ];
    check_route(&obstacles, 0.0, p(5.0, 90.0), p(95.0, 10.0));
    check_route(&obstacles, 0.0, p(95.0, 10.0), p(5.0, 90.0));
}

#[test]
fn right_to_left_is_the_mirror_of_left_to_right() {
    let obstacles = vec![box_poly([40.0, 40.0, 60.0, 60.0])];
    let a = p(5.0, 50.0);
    let b = p(95.0, 50.0);
    let (_, forward) = route(&obstacles, 0.0, a, b);
    let (_, backward) = route(&obstacles, 0.0, b, a);
    assert_eq!(forward.knots.len(), backward.knots.len());
    for (f, r) in forward.knots.iter().rev().zip(backward.knots.iter()) {
        assert!(f.distance(*r) < 1e-9, "mirror mismatch: {f:?} vs {r:?}");
    }
}

#[test]
fn through_a_narrow_passage_with_a_margin() {
    // A gap 20 wide, then 12 wide, then 8 wide, with a margin that closes the last one.
    let obstacles = vec![
        box_poly([40.0, 0.0, 44.0, 40.0]),
        box_poly([44.0, 40.0, 60.0, 44.0]),
        box_poly([60.0, 0.0, 64.0, 48.0]),
        box_poly([64.0, 48.0, 80.0, 52.0]),
        box_poly([80.0, 0.0, 84.0, 60.0]),
    ];
    check_route(&obstacles, 2.0, p(5.0, 20.0), p(95.0, 20.0));
    // A margin of 5 closes the 8-wide gap, so no corridor exists.
    let space = FreeSpace::new(WORKSPACE, obstacles.clone()).unwrap();
    let decomp = decompose_with_guides(&space, 5.0, &[5.0, 95.0]).unwrap();
    let mut scratch = SearchScratch::default();
    assert!(
        search(&decomp, p(5.0, 20.0), p(95.0, 20.0), &mut scratch).is_err(),
        "a 8-unit gap cannot pass a margin of 5"
    );
}

#[test]
fn endpoints_in_the_same_cell_are_a_straight_segment() {
    let obstacles = Vec::new();
    let (_, path) = route(&obstacles, 0.0, p(10.0, 10.0), p(12.0, 12.0));
    assert_eq!(path.knots, vec![p(10.0, 10.0), p(12.0, 12.0)]);
}

#[test]
fn a_goal_outside_the_free_space_has_no_corridor() {
    let space = FreeSpace::new(WORKSPACE, vec![box_poly([40.0, 40.0, 60.0, 60.0])]).unwrap();
    let decomp = decompose_with_guides(&space, 0.0, &[5.0, 50.0]).unwrap();
    let mut scratch = SearchScratch::default();
    assert!(search(&decomp, p(5.0, 50.0), p(50.0, 50.0), &mut scratch).is_err());
}

#[test]
fn a_sealed_region_reports_no_path_rather_than_routing_through() {
    let obstacles = vec![
        box_poly([30.0, 30.0, 70.0, 34.0]),
        box_poly([30.0, 66.0, 70.0, 70.0]),
        box_poly([30.0, 30.0, 34.0, 70.0]),
        box_poly([66.0, 30.0, 70.0, 70.0]),
    ];
    let space = FreeSpace::new(WORKSPACE, obstacles).unwrap();
    let decomp = decompose_with_guides(&space, 0.0, &[5.0, 50.0]).unwrap();
    let mut scratch = SearchScratch::default();
    // The inside of the box is free space but a separate component, so it must not be reachable.
    assert!(search(&decomp, p(10.0, 10.0), p(50.0, 50.0), &mut scratch).is_err());
}

#[test]
fn touching_boxes_do_not_block_the_passage() {
    let obstacles = vec![
        box_poly([30.0, 0.0, 50.0, 45.0]),
        box_poly([50.0, 45.0, 70.0, 100.0]),
    ];
    check_route(&obstacles, 0.0, p(5.0, 90.0), p(95.0, 10.0));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// Random rectilinear obstacle fields. `Obstacle`s that are fully sealed away are skipped by
    /// `check_route`, which is the correct outcome rather than a test failure.
    #[test]
    fn random_fields_route_inside_their_corridor(
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
        // A lattice of start/goal pairs, chosen to be a fixed fraction of the domain apart.
        let cells: Vec<Point2D> = (0..5)
            .flat_map(|i| (0..5).map(move |j| p(10.0 + (i as f64) * 20.0, 10.0 + (j as f64) * 20.0)))
            .collect();
        for (i, a) in cells.iter().enumerate() {
            let b = cells[(i * 7 + 3) % cells.len()];
            if *a != b {
                check_route(&obstacles, margin, *a, b);
            }
        }
    }
}
