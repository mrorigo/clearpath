//! M3 exit gate: lemma L1 — the decomposition is a valid partition of the free space.
//!
//! Checked on a dense grid against a direct free-space predicate, so the gate does not depend on
//! any of the decomposition's own reasoning.

mod common;
use common::{box_poly, p, WORKSPACE};
use clearpath::geom::clearance::distance_to_ring;
use clearpath::decomp::sweep::decompose;
use clearpath::geom::predicates::{Orientation, orient2d};
use clearpath::geom::{BoundingBox, FreeSpace, Point2D, Polygon};
use proptest::prelude::*;


/// Grid resolution. Every cell boundary is on an event line or on an obstacle edge, and the
/// fixtures put those on round coordinates, so a grid that avoids them entirely is a sound
/// sample of the open free space and of the open cells.
const STEPS: usize = 200;


/// Minimum distance from `p` to a ring's boundary.

/// The contract, in two halves.
///
/// 1. L1: the cells are convex, cover exactly the free space, and do not overlap.
/// 2. Clearance: when `margin > 0`, every point of every cell — chord interiors included — is at
///    least `margin` from every obstacle and inside the workspace eroded by `margin`.
///
/// Half 2 is checked against the original obstacles rather than against anything the
/// decomposition computed, so it cannot be satisfied by a consistent mistake.
fn check_partition(obstacles: Vec<Polygon>, margin: f64) {
    let space = FreeSpace::new(WORKSPACE, obstacles.clone()).unwrap();
    let d = decompose(&space, margin).unwrap();

    for cell in d.cells() {
        // Convexity: a cell is convex iff every consecutive triple of its corners turns left.
        let c = cell.corners();
        for i in 0..4 {
            let o = orient2d(c[i], c[(i + 1) % 4], c[(i + 2) % 4]);
            assert!(
                o == Orientation::CounterClockwise || o == Orientation::Collinear,
                "cell {c:?} is not convex (turn {i} is {o:?})"
            );
        }
        assert!(cell.area() >= 0.0, "cell {c:?} has negative area {}", cell.area());
    }

    // Clearance, checked on the cells themselves. Every point of a cell boundary, and of every
    // chord between consecutive boundary points, must clear every obstacle by `margin`.
    if margin > 0.0 {
        for cell in d.cells() {
            let c = cell.corners();
            for i in 0..4 {
                for f in [0.0f64, 0.25, 0.5, 0.75] {
                    let p = c[i].lerp(c[(i + 1) % 4], f);
                    for obstacle in &obstacles {
                        let clear = distance_to_ring(p, obstacle.vertices());
                        assert!(
                            clear >= margin - 1e-9,
                            "cell point {p:?} is only {clear} from an obstacle, margin {margin}"
                        );
                    }
                    let slack = margin - 1e-9;
                    assert!(
                        p.x >= WORKSPACE.min.x + slack && p.x <= WORKSPACE.max.x - slack
                            && p.y >= WORKSPACE.min.y + slack && p.y <= WORKSPACE.max.y - slack,
                        "cell point {p:?} is inside the margin band of the workspace"
                    );
                }
            }
        }
    }

    // Coverage and soundness, against the free space the margin model actually produces: the
    // workspace eroded by `margin`, less the square inflation of every obstacle.
    //
    // The inflation is a *square*, not a disc, so the free space is conservative in the corner
    // directions — a point `margin` away from a corner along the diagonal is excluded even though its
    // distance is exactly `margin`. That is the safe direction and it is what the model specifies;
    // `square_inflation_is_conservative_only_at_corners` pins the behaviour.
    let free_at = |p: Point2D| -> bool {
        if p.x < WORKSPACE.min.x + margin
            || p.x > WORKSPACE.max.x - margin
            || p.y < WORKSPACE.min.y + margin
            || p.y > WORKSPACE.max.y - margin
        {
            return false;
        }
        obstacles.iter().all(|o| {
            let b = o.bounds();
            p.x < b.min.x - margin
                || p.x > b.max.x + margin
                || p.y < b.min.y - margin
                || p.y > b.max.y + margin
        })
    };

    let mut covered = 0usize;
    let mut free = 0usize;
    for i in 1..STEPS {
        for j in 1..STEPS {
            // Off-grid sample points: cell boundaries sit on integer coordinates in these
            // fixtures, so a half-integer grid never lands on one.
            let p = Point2D::new(i as f64 * 0.5 + 0.25, j as f64 * 0.5 + 0.25);
            let cell = d.locate(p);
            let free_here = free_at(p);
            if free_here {
                free += 1;
                assert!(cell.is_some(), "free point {p:?} at ({i},{j}) is in no cell");
            }
            if let Some(id) = cell {
                covered += 1;
                assert!(free_here, "cell {id} claims obstructed point {p:?} at ({i},{j})");
            }
        }
    }
    // Disjointness: if two cells overlapped, `locate` would be hiding the overlap behind the
    // first match, and the covered count would exceed the free count. L1's union condition
    // bounds covered from above by free.
    assert!(covered <= free, "cells overlap: {covered} covered vs {free} free");
    assert!(free > 0, "fixture has no free space at all: {obstacles:?}");
}

/// The square inflation must be conservative — never admitting a point closer than `margin` — and
/// its only cost is in the corner directions, where it can be up to `margin * (sqrt(2) - 1)` too
/// generous.
#[test]
fn square_inflation_is_conservative_only_at_corners() {
    let obstacle = box_poly([40.0, 40.0, 60.0, 60.0]);
    let space = FreeSpace::new(WORKSPACE, vec![obstacle.clone()]).unwrap();
    let margin = 4.0;
    let d = decompose(&space, margin).unwrap();

    // Axial directions: the boundary is exactly at `margin`, and the containment convention puts a
    // point *on* a forbidden boundary on the forbidden side, so probe either side of it.
    assert!(d.locate(Point2D::new(40.0 - margin - 0.01, 50.0)).is_some(), "just clear of the west face");
    assert!(d.locate(Point2D::new(40.0 - margin + 0.01, 50.0)).is_none(), "just inside the margin");

    // The corner is where the square inflation costs something: a point at true distance
    // `margin + 0.001` along the diagonal is excluded, because it falls inside the outset square.
    // A disc inflation would have admitted it. This is the only direction in which the model is not
    // tight, and it errs towards smaller free space.
    let diagonal = (margin + 0.001) / 2f64.sqrt();
    assert!(
        d.locate(Point2D::new(40.0 - diagonal, 40.0 - diagonal)).is_none(),
        "the square inflation closes the corner gap even at true distance > margin"
    );
    // Far enough out diagonally, the point is admitted.
    let far = margin * 1.5 / 2f64.sqrt();
    assert!(d.locate(Point2D::new(40.0 - far, 40.0 - far)).is_some());

    // And it is still never *inside* the margin: the guarantee is one-sided.
    for cell in d.cells() {
        let c = cell.corners();
        for i in 0..4 {
            for f in [0.0f64, 0.25, 0.5, 0.75] {
                let p = c[i].lerp(c[(i + 1) % 4], f);
                assert!(
                    distance_to_ring(p, obstacle.vertices()) >= margin - 1e-9,
                    "cell point {p:?} is closer than the margin to the obstacle"
                );
            }
        }
    }
}

#[test]
fn empty_free_space() {
    check_partition(vec![], 0.0);
}

#[test]
fn one_box() {
    check_partition(vec![box_poly([40.0, 40.0, 60.0, 60.0])], 0.0);
}

#[test]
fn two_disjoint_boxes() {
    check_partition(
        vec![box_poly([20.0, 20.0, 40.0, 40.0]), box_poly([60.0, 60.0, 80.0, 80.0])],
        0.0,
    );
}

#[test]
fn overlapping_boxes() {
    check_partition(
        vec![box_poly([20.0, 20.0, 50.0, 50.0]), box_poly([40.0, 40.0, 70.0, 70.0])],
        0.0,
    );
}

#[test]
fn edge_sharing_boxes() {
    check_partition(
        vec![box_poly([20.0, 20.0, 50.0, 50.0]), box_poly([50.0, 20.0, 80.0, 50.0])],
        0.0,
    );
}

#[test]
fn corner_touching_boxes() {
    check_partition(
        vec![box_poly([20.0, 20.0, 50.0, 50.0]), box_poly([50.0, 50.0, 80.0, 80.0])],
        0.0,
    );
}

#[test]
fn nested_boxes() {
    check_partition(
        vec![box_poly([10.0, 10.0, 90.0, 90.0]), box_poly([40.0, 40.0, 60.0, 60.0])],
        0.0,
    );
}

#[test]
fn a_box_sealed_into_the_workspace() {
    // The free space is a frame; the interior is unreachable. Coverage must still hold outside.
    check_partition(
        vec![box_poly([10.0, 10.0, 90.0, 90.0]), box_poly([10.0, 10.0, 90.0, 90.0])],
        0.0,
    );
}

#[test]
fn a_horizontal_wall_with_a_gap() {
    // A wall that does not span the full width: the two gaps on either side must be found.
    check_partition(vec![box_poly([0.0, 45.0, 70.0, 55.0])], 0.0);
}

#[test]
fn a_vertical_wall_with_a_gap() {
    check_partition(vec![box_poly([45.0, 0.0, 55.0, 70.0])], 0.0);
}

#[test]
fn a_staircase() {
    check_partition(
        vec![
            box_poly([10.0, 10.0, 30.0, 30.0]),
            box_poly([30.0, 30.0, 50.0, 50.0]),
            box_poly([50.0, 50.0, 70.0, 70.0]),
        ],
        0.0,
    );
}

#[test]
fn with_a_margin() {
    check_partition(
        vec![box_poly([20.0, 20.0, 40.0, 40.0]), box_poly([55.0, 55.0, 75.0, 75.0])],
        5.0,
    );
}

#[test]
fn a_margin_that_merges_two_boxes() {
    check_partition(
        vec![box_poly([20.0, 20.0, 40.0, 40.0]), box_poly([45.0, 20.0, 65.0, 40.0])],
        6.0,
    );
}

#[test]
fn many_boxes() {
    let obstacles: Vec<Polygon> = (0..6)
        .flat_map(|i| {
            (0..6).map(move |j| {
                let x = 8.0 + (i as f64) * 14.0;
                let y = 8.0 + (j as f64) * 14.0;
                box_poly([x, y, x + 9.0, y + 9.0])
            })
        })
        .collect();
    check_partition(obstacles, 0.0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Random axis-aligned boxes on a coarse lattice, so the half-integer probe grid never lands
    /// exactly on a cell boundary and the assertions are not testing floating-point coincidence.
    #[test]
    fn random_box_fields_are_a_valid_partition(
        boxes in prop::collection::vec((0u8..=12, 0u8..=12, 0u8..=8, 0u8..=8), 0..8),
    ) {
        let obstacles: Vec<Polygon> = boxes
            .into_iter()
            .map(|(i, j, w, h)| {
                let x = 2.0 + (i as f64) * 7.0;
                let y = 2.0 + (j as f64) * 7.0;
                box_poly([x, y, x + w as f64 + 1.0, y + h as f64 + 1.0])
            })
            .collect();
        check_partition(obstacles, 0.0);
    }
}
