//! How many allocations does a warm query make? (section 8.1)
//!
//! The spec claims a warm query allocates nothing. This measures it rather than asserting it, so
//! the claim and the code cannot drift apart silently.

use pathplan::geom::{BoundingBox, Point2D, Polygon};
use pathplan::{Config, PathPlanner, RouteRequest};
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
    fn measure(f: impl FnOnce()) -> Self {
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

        let warm = Counted::measure(|| {
            planner.route_smooth(&req).expect("the fixture routes");
        });
        let cold = Counted::measure(|| {
            PathPlanner::new().route_smooth(&req).expect("the fixture routes");
        });
        println!(
            "{n} boxes: warm {} allocs / {} bytes, cold {} allocs / {} bytes",
            warm.0, warm.1, cold.0, cold.1
        );
        assert!(
            warm.0 <= cold.0,
            "{n} boxes: a warm query allocated more ({}) than a cold one ({})",
            warm.0,
            cold.0
        );
    }
}
