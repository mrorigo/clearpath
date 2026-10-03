//! How many allocations does a warm query make, and where? (section 8.1)
//!
//! The spec claims a warm query allocates nothing. This measures it rather than asserting it, so the
//! claim and the code cannot drift apart silently — and the per-stage breakdown says which
//! allocations are worth removing, so a future pass has somewhere to start.

use clearpath::funnel::cell_search::SearchScratch;
use clearpath::geom::{BoundingBox, Point2D, Polygon};
use clearpath::{Config, PathPlanner, RouteRequest};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

static ALLOCS: AtomicUsize = AtomicUsize::new(0);
static BYTES: AtomicUsize = AtomicUsize::new(0);
static COUNTING: AtomicUsize = AtomicUsize::new(0);

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) == 1 {
            ALLOCS.fetch_add(1, Ordering::Relaxed);
            BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        }
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) == 1 {
            ALLOCS.fetch_add(1, Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

struct Counted(usize, usize);

impl Counted {
    /// The *minimum* of three runs.
    ///
    /// A `Vec`'s allocation count depends on the capacity it happens to have, so a single run
    /// measures as much about what ran before as about the stage itself. The minimum is the count a
    /// warm caller would pay.
    fn measure(mut f: impl FnMut()) -> Self {
        let mut best: Option<Self> = None;
        for _ in 0..3 {
            let one = Counted::once(&mut f);
            best = Some(match best {
                Some(b) if b.0 <= one.0 => b,
                _ => one,
            });
        }
        best.expect("three runs")
    }

    fn once(f: &mut dyn FnMut()) -> Self {
        ALLOCS.store(0, Ordering::Relaxed);
        BYTES.store(0, Ordering::Relaxed);
        COUNTING.store(1, Ordering::Relaxed);
        f();
        COUNTING.store(0, Ordering::Relaxed);
        Self(ALLOCS.load(Ordering::Relaxed), BYTES.load(Ordering::Relaxed))
    }
}

const WORKSPACE: BoundingBox = BoundingBox {
    min: Point2D { x: 0.0, y: 0.0 },
    max: Point2D { x: 1000.0, y: 1000.0 },
};

fn request(n: usize) -> RouteRequest {
    let mut obstacles: Vec<Polygon> = (0..n)
        .map(|i| {
            let x = 8.0 + (i % 20) as f64 * 48.0;
            let y = 8.0 + (i / 20) as f64 * 48.0;
            Polygon::new(vec![
                Point2D::new(x, y),
                Point2D::new(x + 6.0, y),
                Point2D::new(x + 6.0, y + 6.0),
                Point2D::new(x, y + 6.0),
            ])
            .unwrap()
        })
        .collect();
    obstacles.truncate(n);
    let mut req = RouteRequest::new(WORKSPACE, Point2D::new(20.0, 500.0), Point2D::new(980.0, 500.0));
    req.obstacles = obstacles;
    req.config = Config { margin: 0.0, ..Config::default() };
    req
}

#[test]
fn a_warm_query_allocates_less_than_a_cold_one() {
    for n in [10usize, 50] {
        let req = request(n);
        let mut planner = PathPlanner::new();
        planner.route_smooth(&req).expect("the fixture routes");

        let mut warm_fn = || {
            planner.route_smooth(&req).expect("the fixture routes");
        };
        let warm = Counted::once(&mut warm_fn);
        let cold = Counted::once(&mut || {
            PathPlanner::new().route_smooth(&req).expect("the fixture routes");
        });
        println!(
            "n={n} boxes: warm {} allocs / {} bytes, cold {} allocs / {} bytes",
            warm.0, warm.1, cold.0, cold.1
        );
        assert!(
            warm.0 <= cold.0,
            "{n} boxes: a warm query allocated more ({}) than a cold one ({})",
            warm.0,
            cold.0
        );
    }

    per_stage_allocation_breakdown();
}

/// Where the allocations are, by stage.
///
/// Called from the test above rather than being its own, because the counters are global and
/// `cargo test` runs test functions in parallel: two tests counting at once produce nonsense.
///
/// Kept because "reduce allocations" is only actionable if it is known *which* ones: the
/// decomposition is roughly half of them and scales with the number of event abscissae, while the
/// rest are one small `Vec` per stage and barely move the latency.
fn per_stage_allocation_breakdown() {
    use clearpath::corridor::AdmissibleTangents;
    use clearpath::decomp::sweep::decompose_with_guides;
    use clearpath::funnel::cell_search::search;
    use clearpath::funnel::string_pull::string_pull;
    use clearpath::spline::containment::clamp_and_repair;
    use clearpath::spline::solver::solve_tangents;

    for n in [10usize, 50] {
        let req = request(n);
        let space = req.free_space().unwrap();
        let guides = [req.start.point.x, req.goal.point.x];
        let clearance = req.clearance();

        let Counted(decompose, bytes) = Counted::measure(|| {
            let d = decompose_with_guides(&space, 0.0, &guides).unwrap();
            std::hint::black_box(&d);
        });
        let d = decompose_with_guides(&space, 0.0, &guides).unwrap();
        let Counted(search_allocs, _) = Counted::measure(|| {
            let mut sc = SearchScratch::default();
            let c = search(&d, req.start.point, req.goal.point, &mut sc).unwrap();
            std::hint::black_box(&c);
        });
        let mut scratch = SearchScratch::default();
        let c = search(&d, req.start.point, req.goal.point, &mut scratch).unwrap();
        let Counted(funnel_allocs, _) = Counted::measure(|| {
            let p = string_pull(&d, &c, req.start.point, req.goal.point).unwrap();
            std::hint::black_box(&p);
        });
        let path = string_pull(&d, &c, req.start.point, req.goal.point).unwrap();
        let Counted(admissible, _) = Counted::measure(|| {
            let a = AdmissibleTangents::build(&d, &c.cells, &path.knots).unwrap();
            std::hint::black_box(&a);
        });
        let adm = AdmissibleTangents::build(&d, &c.cells, &path.knots).unwrap();
        let seeds: Vec<Point2D> = (0..path.knots.len())
            .map(|i| {
                let inc = if i > 0 { path.knots[i] - path.knots[i - 1] } else { Point2D::ZERO };
                let out =
                    if i + 1 < path.knots.len() { path.knots[i + 1] - path.knots[i] } else { Point2D::ZERO };
                let span = if inc.length() >= out.length() { inc } else { out };
                span / 3.0
            })
            .collect();
        let Counted(solve, _) = Counted::measure(|| {
            let t = solve_tangents(&seeds, 0.0, &[]);
            std::hint::black_box(&t);
        });
        let Counted(repair, _) = Counted::measure(|| {
            let mut t = solve_tangents(&seeds, 0.0, &[]);
            clamp_and_repair(&adm, &path.knots, &mut t, &clearance, 16, 0.5, &[]);
            std::hint::black_box(&t);
        });
        let mut t = solve_tangents(&seeds, 0.0, &[]);
        clamp_and_repair(&adm, &path.knots, &mut t, &clearance, 16, 0.5, &[]);
        let Counted(final_check, _) = Counted::measure(|| {
            let s = adm.control_points(&path.knots, &t, &clearance).unwrap();
            std::hint::black_box(&s);
        });
        let Counted(clearance_cost, _) = Counted::measure(|| {
            let c = req.clearance();
            std::hint::black_box(c.workspace());
        });

        println!(
            "n={n:3}  decompose {decompose:3} ({bytes:5}B)  search {search_allocs:2}  funnel {funnel_allocs:2}  \
             admissible {admissible:2}  solve {solve}  repair {repair:2}  final {final_check:2}  \
             clearance {clearance_cost}",
        );
    }
}
