//! The two fuzz harnesses, promoted to permanent property tests.
//!
//! `probes/fuzz_router.rs` and `probes/fuzz_smooth.rs` found two real bugs during development:
//! a route check that sampled each leg at nine points and so missed a blocked stretch narrower
//! than an eighth of the leg, and a `Clearance::is_free` that accepted points inside an obstacle.
//! Both are fixed, and both regressions are pinned in `tests/orthogonal.rs` and
//! `tests/routing.rs`. What is here is the other half: the search that found them, kept running.
//!
//! The shape is the same as the probes. Random queries go in; the router's own answer is compared
//! against ground truth by densely sampling every leg of the returned route and testing every
//! sample with `Clearance::is_free` — the same predicate the library routes with, because that is
//! the contract. There is no second, weaker point-in-polygon test here: a check that disagrees
//! with the library's own notion of "free" is a check of something else.
//!
//! # Seeds
//!
//! The iteration count is small enough to run in CI, so the seeds have to do the work that
//! iteration count did in the probes. `PROPERTIES_SEED` pins the seed for a repro; without it the
//! seed is mixed from the clock and the process id, so the six repeated runs in `scripts/ci.sh`
//! explore different seeds. Each seed is deterministic, so a failure names the seed that produced
//! it and replays exactly.

mod common;
use common::{box_poly, p, WORKSPACE};
use clearpath::geom::clearance::Clearance;
use clearpath::geom::{Point2D, Polygon};
use clearpath::{Config, CubicBezierSegment, PathPlanner, RouteRequest};
use std::time::{SystemTime, UNIX_EPOCH};

/// The environment variable that pins the seed, for replaying a failure.
const SEED_VAR: &str = "PROPERTIES_SEED";

/// Samples per leg. A leg of a rectilinear route can be up to the full 100 units of the
/// workspace, so 64 samples is one every 1.6 units — well under the width of the thin
/// protrusions the grid-based router used to place its samples across.
const LEG_SAMPLES: usize = 64;

/// Samples per spline segment. A cubic Bezier between knots is shorter than a workspace
/// diagonal, so 64 samples is well inside the resolution the probes used and still cheap.
const SEGMENT_SAMPLES: usize = 64;

/// The seed this run explores, and the clock value it was mixed with, so a failure report
/// carries both.
struct Seed {
    value: u64,
    source: &'static str,
}

/// A clock- and pid-mixed seed, or the pinned one from `PROPERTIES_SEED`.
///
/// Mixing rather than using the raw clock matters: `cargo test` runs the test binaries in the
/// same second often enough that two of the six CI runs would otherwise get the same seed and
/// explore the same cases.
fn seed() -> Seed {
    if let Ok(pinned) = std::env::var(SEED_VAR) {
        if let Ok(value) = u64::from_str_radix(pinned.trim_start_matches("0x"), 16)
            .or_else(|_| pinned.parse::<u64>())
        {
            return Seed { value, source: "PROPERTIES_SEED" };
        }
    }
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(0);
    // `std::process::id` is a separate clock source, so it breaks ties between processes
    // started in the same nanosecond as well as the same second.
    let mixed = nanos
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add((std::process::id() as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));
    Seed { value: mixed, source: "clock" }
}

/// A xorshift64 PRNG. Small, and the same generator the probes used, so a seed that reproduces
/// in a probe reproduces here.
///
/// Distinct from the proptest strategy generator the other test files use: those draw one
/// arbitrary query per case, whereas these draw thousands of queries per case, and a
/// `TestRunner`'s failure output shrinks to one case rather than to the first bad query in a
/// sequence of them.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        // xorshift64 has a zero fixed point, so a zero seed is replaced rather than allowed to
        // emit a constant stream.
        Rng(if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed })
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Uniform in `[0, 1)`, from the top 53 bits so every value is exactly representable.
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.unit() * (hi - lo)
    }

    /// Uniform in `0..n`. `next() % n` is biased towards the low values, which for `n <= 7` is
    /// a 2% skew at worst and does not matter for choosing an obstacle shape.
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// A random convex polygon inscribed in a circle of radius `r` about `(cx, cy)`, so its bounding
/// box is strictly larger than the shape itself.
///
/// This is the case `Grid::build` cannot see: it places grid lines only at the bounding box
/// faces, so a face of the shape crosses the interior of a cell whose centre it certifies as
/// free. That is the shape that the nine-sample route check used to sail through.
fn random_convex(rng: &mut Rng, cx: f64, cy: f64, r: f64) -> Polygon {
    let n = 3 + rng.below(5) as usize; // 3..=7 vertices
    let base = rng.range(0.0, std::f64::consts::TAU);
    let step = std::f64::consts::TAU / n as f64;
    let pts: Vec<Point2D> = (0..n)
        .map(|i| {
            let a = base + step * i as f64;
            let rad = r * rng.range(0.35, 1.0);
            Point2D::new(cx + rad * a.cos(), cy + rad * a.sin())
        })
        .collect();
    // A near-degenerate ring is rejected by `Polygon::new`, and the box is a fair substitute:
    // the router sees a legal obstacle either way.
    Polygon::new(pts).unwrap_or_else(|_| box_poly([cx - r, cy - r, cx + r, cy + r]))
}

/// A random axis-aligned box with a corner in `[15, 85]` and sides in `[3, 12]`.
///
/// Smaller and more central than the orthogonal probe's boxes, because the smooth router's
/// funnel sees the workspace as a whole: a box that reaches an edge or a corner of the workspace
/// changes which side of the obstacle the corridor wraps, and a wide range of those would leave
/// the funnel refusing most queries and the sampler untested.
fn random_wall(rng: &mut Rng) -> Polygon {
    let x0 = rng.range(15.0, 85.0);
    let y0 = rng.range(15.0, 85.0);
    let w = rng.range(3.0, 12.0);
    let h = rng.range(3.0, 12.0);
    box_poly([x0, y0, (x0 + w).min(99.0), (y0 + h).min(99.0)])
}

/// A random point in `[1, 99]` squared, redrawn while it lands inside an obstacle.
///
/// Redrawing is what keeps validation from rejecting the query: `RouteRequest::validate`
/// reports `EndpointInObstacle`, and a query that never reaches the router tests nothing.
fn random_endpoint(rng: &mut Rng, obstacles: &[Polygon]) -> Point2D {
    let mut guard = 0;
    loop {
        let q = p(rng.range(1.0, 99.0), rng.range(1.0, 99.0));
        if !obstacles.iter().any(|o| o.contains(q)) {
            return q;
        }
        guard += 1;
        if guard >= 50 {
            // The workspace is 98 units wide and the obstacles at most a third of it, so 50
            // redraws is not reachable. Reported rather than looped on, in case it ever is.
            return q;
        }
    }
}

/// The counterexample, at full precision.
///
/// Every coordinate goes out through `{:?}`, never `{:.4}`: a rounded coordinate is a different
/// point, and re-running the case with it produces a different answer, which is the one thing a
/// counterexample must not do.
fn describe(req: &RouteRequest) -> String {
    let w = req.workspace;
    let mut out = format!(
        "\n  workspace: {:?}..{:?}\n  start    : {:?}\n  goal     : {:?}\n  margin   : {:?}\n  obstacles:",
        w.min,
        w.max,
        req.start.point,
        req.goal.point,
        req.config.margin
    );
    for (i, o) in req.obstacles.iter().enumerate() {
        let verts: Vec<String> = o
            .vertices()
            .iter()
            .map(|v| format!("Point2D::new({:?}, {:?})", v.x, v.y))
            .collect();
        out.push_str(&format!("\n    {}: Polygon::new(vec![{}]).unwrap(),", i, verts.join(", ")));
    }
    out
}

/// The first sample of `leg` that `clearance` rejects, and its parameter along the leg.
fn first_blocked_sample(clearance: &Clearance<'_>, a: Point2D, b: Point2D) -> Option<(f64, Point2D)> {
    (0..=LEG_SAMPLES)
        .map(|k| k as f64 / LEG_SAMPLES as f64)
        .map(|t| (t, a.lerp(b, t)))
        .find(|(_, q)| !clearance.is_free(*q))
}

/// The first sample of `segment` that `clearance` rejects, and its parameter along the segment.
fn first_blocked_sample_on_segment(
    clearance: &Clearance<'_>,
    segment: &CubicBezierSegment,
) -> Option<(f64, Point2D)> {
    (0..=SEGMENT_SAMPLES)
        .map(|k| k as f64 / SEGMENT_SAMPLES as f64)
        .map(|t| (t, segment.evaluate(t)))
        .find(|(_, q)| !clearance.is_free(*q))
}

/// Every returned rectilinear route is entirely in the free space.
///
/// Ported from `probes/fuzz_router.rs`. The probes mixed axis-aligned boxes with non-rectangular
/// convex shapes, and included a random margin on half the queries, because a margin is what
/// turns a router's answer into a claim about clearance and the obstacles must then be
/// rectilinear for the request to be answerable at all.
#[test]
fn every_orthogonal_route_is_free() {
    const ITERATIONS: usize = 2_000;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for it in 0..ITERATIONS {
        // 1..=2 obstacles, each a box or a convex shape, centred in the workspace.
        let n_obs = 1 + rng.below(2) as usize;
        let mut obstacles = Vec::with_capacity(n_obs);
        for _ in 0..n_obs {
            let cx = rng.range(20.0, 80.0);
            let cy = rng.range(20.0, 80.0);
            let r = rng.range(3.0, 18.0);
            let poly = if rng.next() % 2 == 0 {
                box_poly([cx - r, cy - r, cx + r, cy + r])
            } else {
                random_convex(&mut rng, cx, cy, r)
            };
            if WORKSPACE.contains(p(cx, cy)) {
                obstacles.push(poly);
            }
        }
        if obstacles.is_empty() {
            continue;
        }

        let start = p(rng.range(1.0, 99.0), rng.range(1.0, 99.0));
        let goal = p(rng.range(1.0, 99.0), rng.range(1.0, 99.0));
        let margin = if rng.next() % 2 == 0 { 0.0 } else { rng.range(0.5, 4.0) };

        let req = RouteRequest {
            workspace: WORKSPACE,
            obstacles: obstacles.clone(),
            start: clearpath::PortConstraint::free(start),
            goal: clearpath::PortConstraint::free(goal),
            config: Config { margin, ..Config::default() },
        };
        if req.validate().is_err() {
            continue;
        }

        let Ok(polyline) = planner.route_orthogonal(&req) else {
            refused += 1;
            continue;
        };
        routed += 1;

        // Ground truth: dense sampling of every leg, against every obstacle, through the
        // library's own predicate.
        let clearance = req.clearance();
        for (leg, w) in polyline.points.windows(2).enumerate() {
            if let Some((t, q)) = first_blocked_sample(&clearance, w[0], w[1]) {
                panic!(
                    "the orthogonal router returned a route that is not free\n\
                     seed {:#x} ({}) iteration {it}\n\
                     leg {leg} of {}: {:?} -> {:?} (length {:?})\n\
                     first blocked sample: t={t:?} point {q:?}\n\
                     polyline: {polyline:?}{}",
                    seed.value,
                    seed.source,
                    polyline.points.len().saturating_sub(1),
                    w[0],
                    w[1],
                    w[0].distance(w[1]),
                    describe(&req),
                );
            }
        }
    }

    println!(
        "orthogonal: {routed} routed, {refused} refused, {ITERATIONS} queries, seed {:#x} ({})",
        seed.value, seed.source
    );
}

/// Every returned spline is entirely in the free space.
///
/// Ported from `probes/fuzz_smooth.rs`. The probe restricted itself to vertical-wall boxes and
/// said why: those are the shapes the orthogonal grid handles exactly, so a refusal there is a
/// funnel or string-pull problem rather than a decomposition one. That restriction is kept.
///
/// The probe also counted refusals, and so does this. A funnel that refuses a route that exists
/// is a real defect, but a property test cannot assert the negative of it — there is no
/// in-library oracle for "a path exists here", and a crossing-number check written in the test
/// would be a second geometry engine whose disagreements would be indistinguishable from
/// router bugs. So refusals are reported, not asserted on.
#[test]
fn every_smooth_spline_is_free() {
    const ITERATIONS: usize = 2_000;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for it in 0..ITERATIONS {
        // 1..=3 vertical walls.
        let n_obs = 1 + rng.below(3) as usize;
        let obstacles: Vec<Polygon> = (0..n_obs).map(|_| random_wall(&mut rng)).collect();
        let start = random_endpoint(&mut rng, &obstacles);
        let goal = random_endpoint(&mut rng, &obstacles);

        let req = RouteRequest {
            workspace: WORKSPACE,
            obstacles: obstacles.clone(),
            start: clearpath::PortConstraint::free(start),
            goal: clearpath::PortConstraint::free(goal),
            config: Config::default(),
        };
        if req.validate().is_err() {
            continue;
        }

        let Ok(segments) = planner.route_smooth(&req) else {
            refused += 1;
            continue;
        };
        routed += 1;

        let clearance = req.clearance();
        for (i, segment) in segments.iter().enumerate() {
            if let Some((t, q)) = first_blocked_sample_on_segment(&clearance, segment) {
                let ends: Vec<String> = segments
                    .iter()
                    .map(|s| format!("{:?} -> {:?}", s.evaluate(0.0), s.evaluate(1.0)))
                    .collect();
                panic!(
                    "the smooth router returned a spline that is not free\n\
                     seed {:#x} ({}) iteration {it}\n\
                     segment {i} of {}: {:?} -> {:?}\n\
                     first blocked sample: t={t:?} point {q:?}\n\
                     segment ends:\n    {}\n{}",
                    seed.value,
                    seed.source,
                    segments.len(),
                    segment.evaluate(0.0),
                    segment.evaluate(1.0),
                    ends.join("\n    "),
                    describe(&req),
                );
            }
        }
    }

    println!(
        "smooth: {routed} routed, {refused} refused, {ITERATIONS} queries, seed {:#x} ({})",
        seed.value, seed.source
    );
}
