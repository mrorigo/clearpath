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
| 10 | 26 µs | 16 µs |
| 50 | 118 µs | 90 µs |
| 200 | 713 µs | 649 µs |

The specification's original target was 50 µs for 10–50 boxes. It was missed by 3x when the
optimisation pass started and is now **met at 10 boxes with a 2x margin**; 50 boxes is 2.4x over.
10 boxes went 157 µs → 26 µs and 50 boxes 510 µs → 118 µs.

The stage split is not where the design suggests, and that is the most useful thing the benchmarks
produced. At 10 boxes the sweep is 7.7 µs. The *funnel* was 77 µs — half the query — and none of
it was the funnel's geometry: it was the corridor-membership check, which sampled every segment
twice per unit of length and asked each sample whether any corridor cell contained it. A sample
standing in for a proof, two thousand exact predicates per long segment. Stating it exactly — a
cell is a trapezoid, so containment is a pair of half-interval intersections — took that stage to
0.4 µs.

The repair had the same shape of problem. Its obstacle loop was gated on the *hull's* bounding box,
but a control hull on a long route is long and thin, so its box is most of the workspace and every
obstacle passed it: 960 segment tests per query at 10 boxes, 6000 at 50. Replacing that proxy with
an exact local test — does this hull *edge* come within `margin` of this obstacle's box? — is a slab
test, and it took the repair from 25 µs to 7.6 µs.

**A spatial index was not built, because the measurements did not call for one.** The cost was never
a lack of locality in the data; it was a bad locality proxy in the query, and four comparisons
answer that exactly where a lookup would answer it approximately. A grid remains the right answer
for a query that genuinely spans the workspace, and is the next step if the obstacle count grows
enough for the linear scan to dominate again.

A warm query makes 89 allocations at 10 boxes (232 at 50), down from 292 and 657 at the start of the
optimisation work. The decomposition — the only stage that *scaled* with the input — now allocates
nothing when warm, because it is refilled in place and owns its transient sweep state. **The
zero-allocation goal of `docs/SPEC.md` §8.1 is still not met**: about 30 small `Vec`s remain across
the other five stages, worth under a microsecond, and removing them is only worth doing for a caller
who needs a genuinely allocation-free query — which argues for exposing a caller-owned scratch
rather than more plumbing inside the stages.

## Development

```
./scripts/ci.sh
```

builds, checks `no_std`, runs clippy with `-D warnings`, enforces that `unsafe` appears only in
`src/spline/simd.rs`, and runs the test suite **six times** — the property tests reseed every run,
and two of the bugs found during development reproduced on only some seeds, so one green run is
not evidence. `cargo bench --bench route` produces the tables above.

Licensed under MIT OR Apache-2.0.
