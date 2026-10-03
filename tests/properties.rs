//! The two fuzz harnesses, promoted to permanent property tests.
//!
//! `probes/fuzz_router.rs` and `probes/fuzz_smooth.rs` found two real bugs during development:
//! a route check that sampled each leg at nine points and so missed a blocked stretch narrower
//! than an eighth of the leg, and a `Clearance::is_free` that accepted points inside an obstacle.
//! Both are fixed, and both regressions are pinned in `tests/orthogonal.rs` and
//! `tests/routing.rs`. What is here is the other half: the search that found them, kept running.

mod common;

use common::{box_poly, p, WORKSPACE};
use clearpath::geom::clearance::Clearance;
use clearpath::geom::{Point2D, Polygon};
use clearpath::router::OrthogonalPolyline;
use clearpath::{Config, CubicBezierSegment, PathPlanner, RouteRequest};
use std::time::{SystemTime, UNIX_EPOCH};

/// Samples per polyline leg. At the workspace scale a leg is tens of units long, so this is
/// one sample per fraction of a unit — far finer than the nine-per-leg check that caused the
/// original bug, which is why the density here is a real choice and not an arbitrary number.
const LEG_SAMPLES: usize = 64;

/// Samples per spline segment.
const SEGMENT_SAMPLES: usize = 64;

/// A seed, and where it came from, so a failure can be replayed.
struct Seed {
    value: u64,
    source: &'static str,
}

/// The seed for this run.
///
/// `PROPERTIES_SEED` pins it for replay. Otherwise it is mixed from the clock *and* the process
/// id: `scripts/ci.sh` runs the suite six times, and a raw clock alone hands two of those runs
/// close enough together to explore the same cases.
fn seed() -> Seed {
    if let Ok(raw) = std::env::var("PROPERTIES_SEED") {
        let trimmed = raw.trim();
        let value = match trimmed.strip_prefix("0x") {
            Some(hex) => u64::from_str_radix(hex, 16),
            None => trimmed.parse::<u64>(),
        };
        if let Ok(value) = value {
            return Seed { value, source: "PROPERTIES_SEED" };
        }
    }
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_nanos()).unwrap_or(0);
    let value = (nanos as u64) ^ ((std::process::id() as u64) << 32);
    Seed { value: value | 1, source: "clock" }
}

/// xorshift64. Zero is replaced by the caller, because zero is this generator's fixed point.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// A value in `[0, n)`.
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// A value in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (self.next() as f64 / u64::MAX as f64) * (hi - lo)
    }
}

/// A convex polygon inside the disc of radius `r` about `(cx, cy)`.
///
/// Deliberately not rectilinear: a box is its own bounding box, so the router's grid splits
/// exactly at its faces and a box proves much less than it appears to.
fn random_convex(rng: &mut Rng, cx: f64, cy: f64, r: f64) -> Polygon {
    let n = 3 + rng.below(5) as usize;
    let phase = rng.range(0.0, core::f64::consts::TAU);
    let mut pts = Vec::with_capacity(n);
    for k in 0..n {
        let angle = phase + core::f64::consts::TAU * k as f64 / n as f64;
        let radius = r * (0.6 + 0.4 * (rng.next() as f64 / u64::MAX as f64));
        pts.push(p(cx + radius * angle.cos(), cy + radius * angle.sin()));
    }
    Polygon::new(pts).expect("a convex ring is a valid polygon")
}

/// A tall axis-aligned wall, the shape the orthogonal grid handles exactly.
fn random_wall(rng: &mut Rng) -> Polygon {
    let x = rng.range(20.0, 75.0);
    let width = rng.range(3.0, 10.0);
    let lo = rng.range(0.0, 40.0);
    let hi = lo + rng.range(25.0, 60.0);
    box_poly([x, lo, x + width, hi])
}

/// A point in the workspace that is not inside any of `obstacles`.
fn random_endpoint(rng: &mut Rng, obstacles: &[Polygon]) -> Point2D {
    for _ in 0..64 {
        let q = p(rng.range(2.0, 98.0), rng.range(2.0, 98.0));
        if obstacles.iter().all(|o| !o.contains(q)) {
            return q;
        }
    }
    p(2.0, 2.0)
}

/// A mix of box and non-rectangular convex obstacles, centred inside the workspace.
///
/// The mix is deliberate. A box is its own bounding box, so the grid splits exactly at its
/// faces; a convex shape is not, so a face of it can cross a grid cell whose centre the grid
/// then certifies as free. Both are needed for the test to mean anything.
fn random_obstacles(rng: &mut Rng) -> Vec<Polygon> {
    let n_obs = 1 + rng.below(2) as usize;
    let mut obstacles = Vec::with_capacity(n_obs);
    for _ in 0..n_obs {
        let cx = rng.range(20.0, 80.0);
        let cy = rng.range(20.0, 80.0);
        let r = rng.range(3.0, 18.0);
        let poly = if rng.next() % 2 == 0 {
            box_poly([cx - r, cy - r, cx + r, cy + r])
        } else {
            random_convex(rng, cx, cy, r)
        };
        if WORKSPACE.contains(p(cx, cy)) {
            obstacles.push(poly);
        }
    }
    obstacles
}

/// A validated request over `WORKSPACE`, or `None` if the endpoints do not form a legal query.
///
/// Building and validating the request is the same step for both routers, and a query that does
/// not validate is not a finding, so it lives here rather than being repeated per test.
fn try_request(
    obstacles: Vec<Polygon>,
    start: Point2D,
    goal: Point2D,
    margin: f64,
) -> Option<RouteRequest> {
    let req = RouteRequest {
        workspace: WORKSPACE,
        obstacles,
        start: clearpath::PortConstraint::free(start),
        goal: clearpath::PortConstraint::free(goal),
        config: Config { margin, ..Config::default() },
    };
    req.validate().is_ok().then_some(req)
}

/// The first of `n+1` evenly spaced parameters at which `at` lands somewhere `clearance` rejects.
///
/// Both routes are checked the same way, so the check is written once: a leg is a straight
/// segment and a spline piece is a curve, and the only thing that differs is the parameterisation.
fn first_blocked_at(
    clearance: &Clearance<'_>,
    n: usize,
    at: impl Fn(f64) -> Point2D,
) -> Option<(f64, Point2D)> {
    (0..=n)
        .map(|k| k as f64 / n as f64)
        .map(|t| (t, at(t)))
        .find(|(_, q)| !clearance.is_free(*q))
}

/// The first sample of `leg` that `clearance` rejects, and its parameter along the leg.
fn first_blocked_sample(clearance: &Clearance<'_>, a: Point2D, b: Point2D) -> Option<(f64, Point2D)> {
    first_blocked_at(clearance, LEG_SAMPLES, |t| a.lerp(b, t))
}

/// The first sample of `segment` that `clearance` rejects, and its parameter along the segment.
fn first_blocked_sample_on_segment(
    clearance: &Clearance<'_>,
    segment: &CubicBezierSegment,
) -> Option<(f64, Point2D)> {
    first_blocked_at(clearance, SEGMENT_SAMPLES, |t| segment.evaluate(t))
}

/// The failing query in full, at full `f64` precision so that it replays exactly.
///
/// Rounded coordinates are useless for replay: whether a sample lands inside an obstacle is
/// decided by the last couple of digits, so every vertex is printed with `{:?}`.
fn describe(req: &RouteRequest) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "\nworkspace {:?}..{:?}\nstart {:?}\ngoal {:?}\nmargin {}\nobstacles:",
        req.workspace.min,
        req.workspace.max,
        req.start.point,
        req.goal.point,
        req.config.margin,
    ));
    for (i, o) in req.obstacles.iter().enumerate() {
        let verts: Vec<String> = o
            .vertices()
            .iter()
            .map(|v| format!("Point2D::new({:?}, {:?})", v.x, v.y))
            .collect();
        out.push_str(&format!("\n    {i}: Polygon::new(vec![{}]).unwrap(),", verts.join(", ")));
    }
    out
}

/// Where a check first failed: the parameter, the point there, and the interval it was on.
///
/// Grouped into one argument because the reporting helpers already sit at seven, and eight
/// separate parameters of which four are `Point2D` read as a data tuple rather than as inputs.
struct Hit {
    t: f64,
    at: Point2D,
    from: Point2D,
    to: Point2D,
}

/// Reports a blocked leg of `polyline` as a failure, with the whole query replayable.
///
/// The counterexample is built here rather than at the call site, so each test reads as the
/// property it asserts rather than as string formatting.
fn reject_blocked_leg(
    seed: Seed,
    req: &RouteRequest,
    polyline: &OrthogonalPolyline,
    leg: usize,
    hit: Hit,
) -> ! {
    let Hit { t, at: q, from: w0, to: w1 } = hit;
    panic!(
        "the orthogonal router returned a route that is not free\n\
         seed {:#x} ({})\n\
         leg {leg} of {}: {:?} -> {:?} (length {:?})\n\
         first blocked sample: t={t:?} point {q:?}\n\
         polyline: {polyline:?}\n{}",
        seed.value,
        seed.source,
        polyline.points.len().saturating_sub(1),
        w0,
        w1,
        w0.distance(w1),
        describe(req),
    )
}

/// Reports a blocked spline segment as a failure, with the whole query replayable.
fn reject_blocked_segment(
    seed: Seed,
    req: &RouteRequest,
    segments: &[CubicBezierSegment],
    i: usize,
    t: f64,
    q: Point2D,
) -> ! {
    let ends: Vec<String> = segments
        .iter()
        .map(|s| format!("{:?} -> {:?}", s.evaluate(0.0), s.evaluate(1.0)))
        .collect();
    panic!(
        "the smooth router returned a spline that is not free\n\
         seed {:#x} ({})\n\
         segment {i} of {}: {:?} -> {:?}\n\
         first blocked sample: t={t:?} point {q:?}\n\
         segment ends:\n    {}\n{}",
        seed.value,
        seed.source,
        segments.len(),
        segments[i].evaluate(0.0),
        segments[i].evaluate(1.0),
        ends.join("\n    "),
        describe(req),
    )
}

/// Every returned rectilinear route is entirely in the free space.
///
/// Ported from `probes/fuzz_router.rs`. Half the queries carry a margin, because a margin is
/// what turns a router's answer into a claim about clearance, and the obstacles must then be
/// rectilinear for the request to be answerable at all.
#[test]
fn every_orthogonal_route_is_free() {
    const ITERATIONS: usize = 2_000;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for _ in 0..ITERATIONS {
        let obstacles = random_obstacles(&mut rng);
        if obstacles.is_empty() {
            continue;
        }
        let start = p(rng.range(1.0, 99.0), rng.range(1.0, 99.0));
        let goal = p(rng.range(1.0, 99.0), rng.range(1.0, 99.0));
        let margin = if rng.next() % 2 == 0 { 0.0 } else { rng.range(0.5, 4.0) };
        let Some(req) = try_request(obstacles, start, goal, margin) else {
            continue;
        };

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
                reject_blocked_leg(seed, &req, &polyline, leg, Hit { t, at: q, from: w[0], to: w[1] });
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
/// Refusals are counted and printed, not asserted on. A funnel that refuses a route that exists
/// is a real defect, but a property test cannot assert the negative of it — there is no
/// in-library oracle for "a path exists here", and a crossing-number check written in the test
/// would be a second geometry engine whose disagreements would be indistinguishable from
/// router bugs.
#[test]
fn every_smooth_spline_is_free() {
    const ITERATIONS: usize = 2_000;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for _ in 0..ITERATIONS {
        // 1..=3 vertical walls.
        let n_obs = 1 + rng.below(3) as usize;
        let obstacles: Vec<Polygon> = (0..n_obs).map(|_| random_wall(&mut rng)).collect();
        let start = random_endpoint(&mut rng, &obstacles);
        let goal = random_endpoint(&mut rng, &obstacles);
        let Some(req) = try_request(obstacles, start, goal, 0.0) else {
            continue;
        };

        let Ok(segments) = planner.route_smooth(&req) else {
            refused += 1;
            continue;
        };
        routed += 1;

        let clearance = req.clearance();
        for (i, segment) in segments.iter().enumerate() {
            if let Some((t, q)) = first_blocked_sample_on_segment(&clearance, segment) {
                reject_blocked_segment(seed, &req, &segments, i, t, q);
            }
        }
    }

    println!(
        "smooth: {routed} routed, {refused} refused, {ITERATIONS} queries, seed {:#x} ({})",
        seed.value, seed.source
    );
}

/// A thin vertical slab: the shape that forces the cell graph to step sideways and
/// then step back, which is what makes a corridor non-monotone in x.
///
/// The two earlier searches use boxes and convex blobs, and both produce monotone
/// corridors, so the funnel's per-portal orientation branch is never taken. This is
/// the layout that reaches it: measured over 20000 queries, 51% of the corridors
/// built here carry a sign change, and none of them produced a collision.
///
/// A sign change is *not* itself a defect — `inside_corridor` catches a path that
/// leaves the corridor and reports `NoPathFound`. What is asserted here is the
/// weaker and load-bearing property: every spline that does come back is free.
#[test]
fn a_non_monotone_corridor_still_returns_a_free_spline() {
    const ITERATIONS: usize = 1_500;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for _ in 0..ITERATIONS {
        let n_obs = 1 + rng.below(3) as usize;
        let obstacles: Vec<Polygon> = (0..n_obs).map(|_| random_slab(&mut rng)).collect();
        let start = random_endpoint(&mut rng, &obstacles);
        let goal = random_endpoint(&mut rng, &obstacles);
        let Some(req) = try_request(obstacles, start, goal, 0.0) else {
            continue;
        };

        let Ok(segments) = planner.route_smooth(&req) else {
            refused += 1;
            continue;
        };
        routed += 1;

        let clearance = req.clearance();
        for (i, segment) in segments.iter().enumerate() {
            if let Some((t, q)) = first_blocked_sample_on_segment(&clearance, segment) {
                reject_blocked_segment(seed, &req, &segments, i, t, q);
            }
        }
    }

    assert!(routed > ITERATIONS / 2, "the search refused almost everything ({routed} of {ITERATIONS}), so it is not testing anything");
    println!(
        "slab: {routed} routed, {refused} refused, {ITERATIONS} queries, seed {:#x} ({})",
        seed.value, seed.source
    );
}

/// A directed port is honoured, or the direction is dropped for a stated reason.
///
/// `docs/SPEC.md` 607 makes a requested endpoint direction a *fixed* boundary
/// condition, and 337-338 say it is honoured only where the resulting tangent is
/// admissible at that knot, otherwise projected. So the contract has two halves and
/// a test that checks only the first will pass while the direction is ignored:
/// the tangent must either follow the request, or be absent because the request was
/// refused. A tangent pointing the *opposite* way satisfies neither.
///
/// The port sweep found 359 of 360 directions being discarded on a bent route.
#[test]
fn a_directed_port_is_honoured_or_the_route_is_refused() {
    const ITERATIONS: usize = 600;
    let seed = seed();
    let mut rng = Rng::new(seed.value);
    let mut planner = PathPlanner::new();
    let (mut routed, mut refused) = (0usize, 0usize);

    for _ in 0..ITERATIONS {
        let obstacles: Vec<Polygon> = (0..1 + rng.below(2) as usize)
            .map(|_| random_wall(&mut rng))
            .collect();
        let start = random_endpoint(&mut rng, &obstacles);
        let goal = random_endpoint(&mut rng, &obstacles);
        // A direction in any quadrant: a quarter of the requests are into the obstacle
        // and must be projected or refused, not silently reversed.
        let direction = p(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
        let Some(mut req) = try_request(obstacles, start, goal, 0.0) else {
            continue;
        };
        if direction.is_zero() {
            continue;
        }
        req.start = clearpath::PortConstraint::directed(start, direction);
        if req.validate().is_err() {
            continue;
        }

        let Ok(segments) = planner.route_smooth(&req) else {
            refused += 1;
            continue;
        };
        routed += 1;

        let wanted = direction.normalize().expect("a non-zero direction normalises");
        let actual = segments
            .first()
            .expect("a route with at least one segment has a first")
            .tangent(0.0);
        if actual.is_zero() {
            // A zero end derivative is permitted only at an undirected port
            // (`docs/SPEC.md` 262-264), so a directed one reaching here is a defect.
            panic!(
                "a directed port produced a zero start tangent\n\
                 requested direction {:?}\n{}\n",
                direction,
                describe(&req)
            );
        }
        let got = actual.normalize().expect("checked non-zero above");
        // A projection onto the convex admissible set cannot reverse the requested
        // direction, so a negative dot means the constraint was dropped, not adjusted.
        assert!(
            got.dot(wanted) >= -1e-9,
            "the start tangent {:?} points away from the requested direction {:?} (dot {:e})\n{}",
            got,
            wanted,
            got.dot(wanted),
            describe(&req)
        );

        let clearance = req.clearance();
        for (i, segment) in segments.iter().enumerate() {
            if let Some((t, q)) = first_blocked_sample_on_segment(&clearance, segment) {
                reject_blocked_segment(seed, &req, &segments, i, t, q);
            }
        }
    }

    println!(
        "directed: {routed} routed, {refused} refused, {ITERATIONS} queries, seed {:#x} ({})",
        seed.value, seed.source
    );
}

/// A thin vertical slab, the layout that makes a corridor non-monotone in x.
fn random_slab(rng: &mut Rng) -> Polygon {
    let ws = WORKSPACE.width().min(WORKSPACE.height());
    let x = rng.range(0.15 * ws, 0.85 * ws);
    let w = rng.range(1.5, 6.0);
    let lo = rng.range(0.05 * ws, 0.45 * ws);
    box_poly([x, lo, x + w, lo + rng.range(0.08 * ws, 0.4 * ws)])
}
