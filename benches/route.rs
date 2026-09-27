//! Benchmarks.
//!
//! Two things are measured that a single end-to-end number would hide: where the time actually
//! goes (the stage split), and what the *warm* query costs, which is the figure section 8.2
//! names. The corpus is a fixed set of axis-aligned boxes on a lattice, so the numbers move only
//! when the implementation does.

#![allow(missing_docs)] // the criterion macros generate items without docs

use criterion::{Criterion, criterion_group, criterion_main};
use pathplan::corridor::AdmissibleTangents;
use pathplan::decomp::sweep::decompose_with_guides;
use pathplan::funnel::cell_search::search;
use pathplan::funnel::string_pull::string_pull;
use pathplan::geom::{BoundingBox, Point2D, Polygon};
use pathplan::spline::containment::clamp_and_repair;
use pathplan::spline::solver::solve_tangents;
use pathplan::{Config, PathPlanner, RouteRequest};
use std::time::Duration;

const WORKSPACE: BoundingBox = BoundingBox {
    min: Point2D { x: 0.0, y: 0.0 },
    max: Point2D { x: 1000.0, y: 1000.0 },
};

/// `n` boxes on a lattice, spread over the workspace.
fn corpus(n: usize) -> Vec<Polygon> {
    let pitch = (1000.0 / (n as f64).sqrt().max(1.0) * 0.55).max(12.0);
    (0..n)
        .map(|i| {
            let cols = (1000.0 / pitch).floor() as usize;
            let x = 8.0 + (i % cols.max(1)) as f64 * pitch;
            let y = 8.0 + (i / cols.max(1)) as f64 * pitch;
            Polygon::new(vec![
                Point2D::new(x, y),
                Point2D::new(x + 6.0, y),
                Point2D::new(x + 6.0, y + 6.0),
                Point2D::new(x, y + 6.0),
            ])
            .unwrap()
        })
        .collect()
}

fn request(obstacles: Vec<Polygon>, margin: f64) -> RouteRequest {
    let mut req = RouteRequest::new(WORKSPACE, Point2D::new(20.0, 500.0), Point2D::new(980.0, 500.0));
    req.obstacles = obstacles;
    req.config = Config { margin, ..Config::default() };
    req
}

fn bench_end_to_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("route_smooth");
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));
    for n in [10usize, 50, 200] {
        for margin in [0.0f64, 1.0] {
            let req = request(corpus(n), margin);
            let mut planner = PathPlanner::new();
            // One untimed call so the corpus is proven routable before it is timed.
            assert!(planner.route_smooth(&req).is_ok(), "corpus of {n} at margin {margin} has no route");
            group.bench_function(format!("{n}_boxes_margin_{margin}"), |b| {
                b.iter(|| {
                    let mut planner = PathPlanner::new();
                    planner.route_smooth(&req).unwrap()
                })
            });
        }
    }
    group.finish();
}

fn bench_stages(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage");
    group.warm_up_time(Duration::from_millis(400));
    group.measurement_time(Duration::from_secs(2));
    for n in [10usize, 50, 200] {
        let req = request(corpus(n), 0.0);
        let space = req.free_space().unwrap();
        let guides = [req.start.point.x, req.goal.point.x];
        let decomp = decompose_with_guides(&space, req.config.margin, &guides).unwrap();
        let mut scratch = pathplan::router::scratch::Scratch::default();
        let corridor = search(&decomp, req.start.point, req.goal.point, &mut scratch.search).unwrap();
        let path = string_pull(&decomp, &corridor, req.start.point, req.goal.point).unwrap();
        let admissible = AdmissibleTangents::build(&decomp, &corridor.cells, &path.knots).unwrap();
        let clearance = req.clearance();
        let seeds: Vec<Point2D> = (0..path.knots.len())
            .map(|i| {
                let inc = if i > 0 { path.knots[i] - path.knots[i - 1] } else { Point2D::ZERO };
                let out = if i + 1 < path.knots.len() { path.knots[i + 1] - path.knots[i] } else { Point2D::ZERO };
                let span = if inc.length() >= out.length() { inc } else { out };
                span / 3.0
            })
            .collect();

        let label = format!("{n}_boxes");
        group.bench_function(format!("decompose_{label}"), |b| {
            b.iter(|| decompose_with_guides(&space, 0.0, &guides).unwrap())
        });
        group.bench_function(format!("search_{label}"), |b| {
            b.iter(|| search(&decomp, req.start.point, req.goal.point, &mut scratch.search).unwrap())
        });
        group.bench_function(format!("funnel_{label}"), |b| {
            b.iter(|| string_pull(&decomp, &corridor, req.start.point, req.goal.point).unwrap())
        });
        group.bench_function(format!("solve_{label}"), |b| {
            b.iter(|| solve_tangents(&seeds, req.config.tangent_bias))
        });
        group.bench_function(format!("repair_{label}"), |b| {
            let t = solve_tangents(&seeds, req.config.tangent_bias);
            b.iter(|| {
                let mut t = t.clone();
                clamp_and_repair(&admissible, &path.knots, &mut t, &clearance, 16, 0.5)
            });
        });
    }
    group.finish();
}

fn bench_warm(c: &mut Criterion) {
    // The figure section 8.2 names: a planner that has already run, against a request of the same
    // shape. The difference from `route_smooth` above is the retained scratch.
    let mut group = c.benchmark_group("warm");
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));
    for n in [10usize, 50] {
        let req = request(corpus(n), 0.0);
        let mut planner = PathPlanner::new();
        planner.route_smooth(&req).expect("the fixture routes");
        group.bench_function(format!("{n}_boxes"), |b| {
            b.iter(|| planner.route_smooth(&req).unwrap())
        });
    }
    group.finish();
}

fn bench_orthogonal(c: &mut Criterion) {
    let mut group = c.benchmark_group("route_orthogonal");
    group.warm_up_time(Duration::from_millis(400));
    group.measurement_time(Duration::from_secs(2));
    for n in [10usize, 50] {
        let req = request(corpus(n), 0.0);
        let mut planner = PathPlanner::new();
        assert!(planner.route_orthogonal(&req).is_ok());
        group.bench_function(format!("{n}_boxes"), |b| {
            b.iter(|| planner.route_orthogonal(&req).unwrap())
        });
    }
    group.finish();
}

criterion_group!(benches, bench_end_to_end, bench_stages, bench_warm, bench_orthogonal);
criterion_main!(benches);
