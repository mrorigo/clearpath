# pathplan

Deterministic 2D obstacle-avoidance path planning and smooth cubic Bézier fitting, in `no_std` +
`alloc`.

Given a workspace, polygonal obstacles, and two endpoints, `pathplan` returns either a `C¹` chain
of cubic Bézier segments or a rectilinear polyline. Every returned curve is collision-free with the
configured margin, and the whole pipeline is deterministic: the same input produces the same output
bits on the same target.

```rust
use pathplan::{BoundingBox, Point2D, Polygon, PathPlanner, RouteRequest};

let mut req = RouteRequest::new(
    BoundingBox::new(Point2D::new(0.0, 0.0), Point2D::new(10.0, 10.0)),
    Point2D::new(1.0, 5.0),
    Point2D::new(9.0, 5.0),
);
req.obstacles = vec![Polygon::new(vec![
    Point2D::new(4.0, 4.0),
    Point2D::new(6.0, 4.0),
    Point2D::new(6.0, 6.0),
    Point2D::new(4.0, 6.0),
])?];

let mut planner = PathPlanner::new();
let segments = planner.route_smooth(&req)?;
# Ok::<(), pathplan::PathPlanError>(())
```

## How it works

Four stages, each consuming plain data:

| Stage | Module | What it does |
| --- | --- | --- |
| 1 | `decomp` | A `+x` sweep decomposes the free space into convex trapezoids with vertical portals. |
| 2 | `funnel` | A\* over the cell adjacency finds a corridor; the Lee–Preparata funnel pulls the shortest taut string through it. |
| 3 | `corridor` | Each knot gets the set of tangents that keep the neighbouring control points in the free space. |
| 4 | `spline` | A minimum-curvature tangent solve, then a repair that clamps until every control hull is clear. |

Three properties are worth knowing about, because they are where the design earns its keep:

**The margin is applied once, exactly.** The forbidden geometry is the workspace eroded by
`margin` and every obstacle grown by it — a *square* inflation of a rectilinear partition, which is
exact and needs no arc approximation. A nonzero `margin` against a non-rectilinear obstacle is
reported as `MarginUnsupportedGeometry` rather than approximated, because approximating it
silently produces curves that pass within `margin` of the obstacle they were pushed off.

**Testing a cell's centre is exact, not a sample.** The sweep's event lines pass through every
obstacle face, so no face crosses a cell's interior and a cell is wholly in or wholly out of each
obstacle.

**Collisions are prevented, not filtered.** A cubic Bézier lies inside the convex hull of its four
control points, so a clear hull is a clear curve. Every hull is tested exactly against the
obstacles — vertices by distance, edges for both crossings *and* clearance, obstacle vertices for
containment — and a route that fails is refused rather than emitted.

## What it does not do

* No incremental decomposition: the free space is rebuilt per query. A warm query therefore
  allocates (see the performance note below).
* No global optimality across homotopy classes. The corridor search returns one corridor and the
  funnel is optimal *within* it.
* No self-intersecting obstacle polygons, and no polygons with holes. A hole is expressible as a
  reversed-wound simple polygon.
* No `simd` kernels yet. The feature exists and is the only place `unsafe` is permitted.

`docs/SPEC.md` is the normative specification, and it is kept honest: it records every ambiguity
that was resolved, including the several where the original design turned out to be unimplementable
and was replaced.

## Performance

Measured on an aarch64 core, release build, 1000×1000 workspace, axis-aligned boxes on a lattice.
`cargo bench` reproduces them.

| Obstacles | `route_smooth`, margin 0 | margin 1 |
| --- | --- | --- |
| 10 | 157 µs | 122 µs |
| 50 | 510 µs | 407 µs |
| 200 | 1.65 ms | 1.49 ms |

The specification's original target was 50 µs for 10–50 boxes; **it is not met**, and
`docs/SPEC.md` §8.2 carries the measurements rather than the aspiration.

The stage split is not where the design would suggest. At 10 boxes: the sweep is 9 µs (6%), the
corridor search 1 µs, the tangent solve 0.25 µs — and the *funnel* 77 µs and the *repair* 31 µs,
together three quarters of the query. Both are dominated by exact predicates in their verification
layers, not by the geometry they are verifying. That is where an optimisation pass should start.

A warm query makes 292 allocations at 10 boxes and 703 at 50, against 298 and 710 for a cold one:
the retained scratch is real but marginal, because the decomposition is rebuilt every time.

## Development

```
./scripts/ci.sh
```

builds, checks `no_std`, runs clippy with `-D warnings`, enforces that `unsafe` appears only in
`src/spline/simd.rs`, and runs the test suite **six times** — the property tests reseed every run,
and two of the bugs found during development reproduced on only some seeds, so one green run is
not evidence. `cargo bench --bench route` produces the tables above.

Licensed under MIT OR Apache-2.0.
