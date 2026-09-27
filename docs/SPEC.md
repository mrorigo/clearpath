# Specification: `pathplan`

**Status:** Normative. This document replaces `docs/SPEC_DRAFT.md`.
Where the draft was ambiguous or self-contradictory, the resolution is stated inline as
**Resolution:** and listed in Appendix A.

**Crate:** `pathplan`
**License:** MIT OR Apache-2.0
**Edition:** Rust 2024 (single edition; the draft's "2024 / 2021" is dropped — `Cargo.toml` pins 2024).
**Target:** `#![no_std]` + `alloc` core. Optional `std` and `simd` features.

**Terminology used throughout**

| Term | Meaning |
| --- | --- |
| Workspace | The bounding box supplied by the caller; the region the planner may route inside. |
| Obstacle | A closed simple polygon, wound counter-clockwise, with its interior forbidden. |
| Free space | `Workspace \ Union(obstacle interiors)`. |
| Margin | Non-negative clearance required between the routed curve and every obstacle boundary. |
| Portal | A vertical shared boundary between two adjacent free-space cells. Portal endpoints lie on obstacle edges or workspace walls. |
| Cell | A convex trapezoid of the vertical decomposition: bounded above and below by straight boundary edges and left/right by portals. Every cell is a subset of free space. |
| Corridor | A convex polygon that is a subset of free space and that contains one taut-string segment plus its slack space. |
| Knot | A vertex of the taut-string path; a Bezier segment is emitted per knot pair. |
| Eroded workspace | The workspace bounding box shrunk inward by `margin` on all four sides. All routing happens inside it. |
| Inflated domain | The workspace eroded by `margin` and every obstacle grown by `margin` (§3.3). All planning happens here. |

---

## 1. Executive Summary

`pathplan` is a deterministic 2D obstacle-avoidance path planner and smooth spline fitter for
2D polygonal free space. Given a bounding workspace, polygonal obstacles, start/goal points and
optional forced tangent directions, it produces either

* a `$C^1`-continuous (optionally `$C^2$`) chain of cubic Bezier segments, or
* a rectilinear (Manhattan) polyline.

Every returned curve is guaranteed collision-free with respect to the configured margin by
construction, not by post-hoc checking.

The library is intended as a modern, memory-safe replacement for Graphviz' `libpathplan`.

### 1.1 What "replacement" means, and what it does not

**v1 is a geometry library, not a Graphviz drop-in.** Its contract is: *given polygonal free
space and two points, return a collision-free curve of a stated continuity class, fast and
deterministically.* Route quality is defined by the stated invariants and the §8 latency target,
nothing else.

**Non-goals for v1** (each is a real difference from `libpathplan`, and each is deliberate):

| Not in v1 | Why |
| --- | --- |
| Node bounding boxes, ports, ranks, edge splines, `polyline`/`ortho` Graphviz option parity | These are Graphviz *graph layout* concerns layered on routing. v1 takes free points and tangents only. |
| Byte-identical routes to `libpathplan` | Depends on undocumented tension defaults and traversal order. Not a meaningful target for a rewrite. |
| Polygon obstacles with holes | Input is a list of simple polygons. A hole is expressible as a reversed-wound simple polygon, which the offset step makes unnecessary in practice. |
| Incremental / cached decomposition updates | §8.4. |
| Global optimum over *all* homotopy classes | The search returns one corridor; the funnel is optimal *within* it, not across corridors. |
| Curve queries, curvature bounds, arc-length parameterisation, obstacle-mutation callbacks | Out of scope; would be separate crates. |

**Milestone order for v1.** Each milestone is independently mergeable and green:

| M | Deliverable | Gate |
| --- | --- | --- |
| M1 | `geom`: `Point2D`, `Polygon`, exact `orient2d`/`incircle`, winding test | unit tests incl. degenerate/collinear/coincident |
| M2 | `geom::inflate` + `geom::clip`, §3.3 margin model | offset + clip property tests |
| M3 | `decomp`: vertical decomposition, cell graph, cell location | decomposition is a valid partition (M3 property tests) |
| M4 | `funnel`: A* cell search + Lee–Preparata string pull | funnel path is collision-free and in-corridor |
| M5 | `corridor`: per-segment convex corridor | corridors are convex, in free space, and overlap at knots |
| M6 | `spline` + `router::smooth` | §4.1–§4.3 invariants, §9.3 straight-line reduction |
| M7 | `router::orthogonal` | §9.8 |
| M8 | `proptest` suite, benches, `no_std` + `unsafe` CI gates | §10 |

M1–M5 are geometry with no public API beyond internal types; M6 is where `route_smooth` first
compiles, M7 adds `route_orthogonal`. If v1 is cut short, M1–M6 alone ship a usable library.


---

## 2. Pipeline

Four sequential, decoupled stages. Each stage consumes a plain data structure and produces
another; no stage reaches back into an earlier one.

```
[workspace bbox + obstacles + margin]
                |
                v
Stage 1a margin Erode the workspace by margin; partition each rectilinear
          |        obstacle into rectangles and grow each by margin
                v
Stage 1b decomp Vertical (trapezoidal) sweep decomposition of the eroded free space
          |  -> Vec<Cell> (convex) + cell adjacency + cell location
                v
Stage 2  funnel  A* cell search for a corridor, then Lee-Preparata funnel string pull
          |  -> TautPath (ordered knot list, on cell boundaries)
                v
Stage 3  corridor Per-segment convex corridor, built from the cells the taut segment crosses
          |  -> Vec<Corridor>, one per knot pair
                v
Stage 4  spline   Tridiagonal tangent solve + convex-hull containment enforcement
                |
                v
        [Vec<CubicBezierSegment> | OrthogonalPolyline]
```


**Stage independence rule.** Stage 2/3/4 operate on data, not on the mesh. This is what makes
`corridor` and `spline` testable against hand-built fixtures without constructing a
decomposition.

### 2.1 Resolution: the orthogonal path is a sibling, not a fifth stage

The draft's pipeline diagram ended in `[Vec or OrthogonalPolyline]` without defining the
orthogonal algorithm. This spec defines `route_orthogonal` as a **separate planner** that
reuses only Stage 1 (the workspace decomposition). See §7. It is *not* a post-processing of the
smooth spline and shares none of Stage 4.

---

## 3. Mathematical Foundations

### 3.1 Robust predicates

All topological decisions (orientation, in-circle, point-in-polygon winding, side-of-line
membership, hull construction ordering) use adaptive exact arithmetic in the style of
Shewchuk's `orient2d` / `incircle`: a floating-point filter evaluates the determinant, and only
when the result is within the filter's error bound is the exact expansion-arithmetic fallback
evaluated. Symbolic perturbation is **not** used.

Raw naive determinants are forbidden for any decision that can change topology.

**Distance and area arithmetic is not exact.** `sqrt`, division, and the containment/clearance
checks in §6 are ordinary `f64` operations and are subject to rounding. The guarantee in §4.1 is
a guarantee of *topological* correctness plus clearance to within floating-point rounding, not
exact arithmetic.

`sqrt` comes from `libm` when the `std` feature is disabled.

### 3.2 Precision constants

Exactly one internal tolerance exists, `geom::EPS: f64 = 1e-12`, used only to decide "is this
length/determinant effectively zero" for non-topological tests (duplicate-vertex detection,
degenerate-edge rejection). It is never used for orientation or containment decisions. It is not
part of the public API.

### 3.3 The margin model

The draft left margin to be re-applied downstream, and its two attempts disagreed (§4.1 said
corridors are free space "eroded by margin" while §6.4 inflated the taut segment along its normal
and re-clipped). **This spec applies `margin` exactly once, before anything else looks at the
geometry:**

> The forbidden region is the workspace eroded inward by `margin` and every obstacle **grown
> outward** by `margin`. The decomposition, the funnel, the corridors, and the spline all operate
> on that inflated domain, and no stage after this one refers to `margin` at all.

Consequences, all intentional:

* A curve is `>= margin` from every obstacle **and** from the workspace boundary. A caller that
  wants the curve to graze the workspace wall sets `margin = 0`.
* `margin` may exceed the workspace's slack. The eroded box then has no positive extent, the
  decomposition is empty, and the query reports `NoPathFound` — not a silent clamp, and not a
  silently ignored margin.
* Obstacles that no longer reach the eroded workspace are dropped, which is correct: they no longer
  constrain anything.
* The erosion does not usually *swallow* an obstacle inside the workspace. An obstacle at distance
  `$d$` from a wall survives iff `$d \ge 2 \cdot 	ext{margin}$`, so a large margin makes the
  workspace smaller, not the obstacles.

**How the growth is realised, and what it costs.** The growth is a **square inflation**, not a
circular one: an axis-aligned box grown by `margin` is a box, exactly. For a *rectilinear* obstacle
the model is exact, because

$$\bigl\{p : \mathrm{dist}(p,\textstyle\bigcup_i R_i) \ge m\bigr\}
= \bigcap_i \bigl\{p : \mathrm{dist}(p, R_i) \ge m\bigr\}
= \mathbb{R}^2 \setminus \bigcup_i \mathrm{outset}(R_i, m)$$

where `$R_i$` are disjoint rectangles partitioning the obstacle. The pieces need not compose — only
their complements do — so a rectilinear obstacle is partitioned into rectangles by a scanline
(`geom::polygon::rectilinear_rectangles`) and each is grown independently. Overlapping and nested
obstacles need no special handling at all, because the sweep's free-space test is a point-in-box
test over the grown boxes.

Square rather than round inflation is conservative in the corner directions only: a point at true
distance `margin` from a corner along the diagonal is excluded, so the realised clearance is
`>= margin` everywhere and the free space is at most `margin (\sqrt{2} - 1)` smaller than optimal
near each corner. The guarantee is one-sided by design, and the M3 gate asserts the one-sided form.

**A nonzero `margin` against a non-rectilinear obstacle is rejected**, with
[`PathPlanError::MarginUnsupportedGeometry`]. This is deliberate and is the one v1 limitation in the
margin model:

* Growing a general polygon's boundary requires the offset features to be **trimmed against each
  other** — where one feature's offset is cut off because another feature's offset covers it. An
  untrimmed offset is not a conservative approximation, it is simply wrong: a `margin` large
  relative to a notch produces a boundary that runs *through* the obstacle it came from, and a
  `margin` relative to two nearby features produces a ring that is not simple at all. Measured on
  a comb-shaped obstacle, the untrimmed construction put offset boundary points at zero distance
  from the obstacle on 854 sampled positions.
* The alternatives were considered and rejected. Eroding the *decomposition cells* instead does not
  work either, and for a less obvious reason worth recording: a cell's own faces do not account for
  an obstacle feature just beyond one of its portals, so a cell can sit arbitrarily close to an
  obstacle it does not touch, and the eroded region is then not even convex.
* Reporting the error is the only option that cannot produce a collision. Silently dropping the
  obstacle would route straight through it.

`margin == 0` is exact for every polygon and is the identity operation.

`margin` must be finite and `>= 0`; otherwise `InvalidConfig`.

---

## 4. Invariants

### 4.1 Safety

For every output curve and every `$t \in [0,1]$ on every segment, the point lies in the **offset
free space** (§3.3), i.e. its distance to every obstacle boundary *and* to the workspace boundary
is `>= margin`, up to floating-point rounding of the distance computation. The implementation
ensures this by construction, in two steps:

1. Every control point of every segment lies inside the segment's convex corridor (§6.5).
2. Every corridor is a subset of the offset free space.

By the convex-hull property of cubic Beziers, (1) implies the whole curve lies in the corridor;
by construction of the corridors from the decomposition cells, (2) implies the clearance. Since
`margin` is applied exactly once, before the decomposition, no stage can forget it. **Sampling is a test
strategy, not the guarantee** (see §9).

### 4.1a Required lemma (decomposition)

> **L1.** The cells produced by Stage 1b are convex, have pairwise disjoint interiors, have union
> equal to the offset free space, and each cell's closure is contained in it.

L1 is the load-bearing assumption of everything downstream: the funnel's portal sequence, the
corridors, and hence §4.1. M3's exit gate is a property test for L1, and §9.5 tests the
control-hull containment that relies on it. A cell that violates L1 is a bug, not a
degenerate input.


### 4.2 Continuity

For adjacent segments `$S_k`, `$S_{k+1}$`:

* `$C^1$` (default, always satisfied):
  `$S_k(1) = S_{k+1}(0)` exactly, and `$S_k'(1) = S_{k+1}'(0)$ exactly, as `f64` vectors.
* `$C^2$` (opt-in feature `c2`): additionally `$S_k''(1) = S_{k+1}''(0)` within `1e-9`.

**Resolution.** The draft wrote `$S_k'(1) = \alpha\,S_{k+1}'(0)$ for `alpha > 0` while calling
the result `$C^1$`. That is geometric ($G^1$) continuity, not parametric `$C^1$. This spec
requires true `$C^1$`: the derivative vectors are equal, and tangent *direction* continuity is
satisfied as a corollary. When a caller requests only a tangent *direction* at a port
(`PortConstraint::direction`), the routine scales the interior tangent to the exact length
required by `$C^1$` while preserving that direction.

Degenerate case: a segment whose end derivative is the zero vector is only permitted at a port
whose `direction` is `None` and whose local knot spacing is zero; in that case the zero is
propagated to the neighbouring segment so `$C^1$` still holds exactly.

### 4.3 Straight-line reduction

With `obstacles` empty (or all obstacles entirely outside the workspace interior), `start != goal`,
and no `direction` constraints, the output is exactly one segment whose four control points are
collinear and monotonically ordered along `start -> goal`. See §9.3 for the test.

**Resolution.** The draft titled this "Monotonic Convergence", which describes no iteration. It
is renamed *Straight-line reduction*.

### 4.4 Determinism

Two tiers, because the draft's "bitwise identical across all supported architectures" is not
achievable alongside the SIMD feature and is contradicted by it:

* **D1 — Semantic determinism (portable).** For identical inputs, the planner produces the same
  sequence of knots, the same cell assignment, the same classification of every query point, and
  the same discrete decisions (corridor count, repair iteration count, solver branch choices).
  This is guaranteed on all targets and is what the test suite asserts.
* **D2 — Bitwise determinism (same target).** On a fixed target triple (architecture, target
  feature set, compiler), the returned `f64` values are bit-identical across runs.

To make D2 hold:

* No reliance on iteration order of hash containers; all maps are sorted `Vec`s or `BTreeMap`s.
* No reliance on pointer or address values in any output or branch.
* No thread or task scheduling influence on query execution (one query runs on one thread).
* `simd`-gated kernels must produce results identical to the scalar path; they may only change
  throughput. A test asserts `cfg(feature = "simd")` outputs equal the scalar outputs bitwise.

**Cross-architecture bitwise equality of `f64` outputs is explicitly not a goal** and is not
tested. FMA contraction is disabled for the crate's numerics modules via
`#[allow]`-free opt-out where the toolchain permits; where it cannot be disabled, the affected
kernels are kept in the scalar path.

### 4.5 Unsafe code

`#![deny(unsafe_code)]` at the crate root, and `#![allow(unsafe_code)]` inside `spline::simd`
only. Note this is `deny`, not `forbid`: an inner `allow` **cannot** override an outer `forbid`,
so the draft's pairing was not compilable. `deny` plus the CI grep (section 9.7) gives the same
practical guarantee.

`geom::predicates` needs no `unsafe` at all. Its exact product is Dekker split-and-multiply over
`f64`, not `f64::mul_add`, because `mul_add` is only in `core` when the target has a hardware FMA
— relying on it would make the predicate silently target-dependent, which section 4.4 forbids.

### 4.6 Shared-tangent clamping (the multi-cell rule)

An interior knot `$k_i$ is shared by exactly two segments, `$S_{i-1}$` and `$S_i`, and their
corridors `$C_{i-1}$`, `$C_i$ are *different* convex polygons. A single tangent `$T_i` serves
both. This is a genuine constraint, not a detail, and the draft left it implicit.

> **Clamp rule.** For each interior knot, the admissible set is
> `$A_i = C_{i-1} \cap C_i$`. Because the corridors are built to share the portal at `$k_i`
> (§6.4), `$k_i \in A_i$ and `$A_i` is non-empty; it is a convex polygon. The tangent `$T_i`
> must satisfy `$k_i + T_i/3 \in A_i$ (from `$S_i$`) **and** `$k_i - T_i/3 \in A_i` (from
> `$S_{i-1}$`)`, which is exactly
> `$k_i \pm T_i/3 \in A_i \iff T_i/3 \in A_i - k_i$`, i.e. `$T_i \in 3\,(A_i - k_i)$` — a convex
> region symmetric about the origin.

Consequences:

* Clamping is done on the *tangent*, never on an individual control point. Repairs (§6.5) and
  the containment projection both operate on `$T_i$ in this region, so `$C^1$` survives every
  repair by construction. This is why repair cannot be per-segment.
* Projection onto `$A_i - k_i$ is a single well-defined operation: closest point of a convex
  polygon to a point, found by the standard "max over violated half-planes" walk using exact
  `orient2d`. No iteration, no tolerance ladder.
* End knots have a single incident corridor, so the admissible set is `$3\,(C_0 - k_0)$` and
  the same rule applies. This is what makes a `PortConstraint::direction` satisfiable: the
  direction is used if the resulting tangent is in the admissible set, otherwise the projection
  is taken and the requested direction is *not* honoured exactly. A caller that needs an exact
  direction must supply a `margin`/geometry where it is achievable. This is documented
  behaviour, not an error.
* If a knot's admissible set is `{0}` (both corridors pinch at the knot), `$T_i = 0` and the two
  segments are straight chords meeting at `$k_i` — a legitimate, if visually angular, result.
  The degenerate-derivative clause of §4.2 covers it.

### 4.7 Required lemma (funnel validity)

> **L2.** For any cell corridor returned by Stage 2a, the Lee–Preparata funnel over its portal
> sequence returns a polyline contained in the union of the corridor's cells, and every knot of
> that polyline lies on a cell boundary.

L2 holds because each cell is convex (L1) and consecutive cells meet in a portal, so the portal
sequence is a valid funnel input by construction. M4's exit gate tests L2 directly (every output
point is inside some cell of the corridor); §4.1 then rests on L1 and L2, not on hope.


---

## 5. Crate layout

```
pathplan/
  Cargo.toml
  src/
    lib.rs            # public facade, re-exports, crate-level lints
    error.rs          # PathPlanError
    geom/
      mod.rs
      point.rs        # Point2D, Vec2 ops, Segment
      polygon.rs      # Polygon, BoundingBox, robust winding test
      predicates.rs   # exact orient2d, incircle, expansions
      decomp/           # free-space decomposition (replaces the draft's `cdt/`)
      mod.rs
      cell.rs         # convex Cell, Portal, cell adjacency
      sweep.rs        # vertical decomposition sweep
      locate.rs       # point -> cell index
    funnel/
      mod.rs
      cell_search.rs  # A* corridor search over cell adjacency
      string_pull.rs  # Lee & Preparata funnel
    corridor/
      mod.rs
      convex_hull.rs  # corridor construction
    spline/
      mod.rs
      bezier.rs       # CubicBezierSegment, evaluation, flattening
      solver.rs       # tridiagonal solve (Thomas), C1/C2 tangent systems
      containment.rs  # corridor containment, repair, de Casteljau
      simd.rs         # feature = "simd"; the only unsafe in the crate
    router/
      mod.rs
      smooth.rs       # stages 1-4 orchestration
      orthogonal.rs   # rectilinear planner
      scratch.rs      # reusable buffers
  benches/
  tests/
```

`error.rs` and `scratch.rs` are additions to the draft; `pathplan::error` exists so that
`no_std` builds do not depend on the layout of `lib.rs`.

---

## 6. Module specifications

### 6.1 `geom`

* `Point2D { x: f64, y: f64 }` with `Add`, `Sub`, `Neg`, `Mul<f64>`, `dot`, `cross`,
  `norm_squared`, `length`, `normalize` (`normalize` returns `Option`, `None` for zero-length).
* `predicates::orient2d(a, b, c) -> Orientation` where
  `Orientation::{Clockwise, CounterClockwise, Collinear}`, exact per §3.1.
* `predicates::incircle(a, b, c, d) -> Incircle` (`In`, `Out`, `OnCircle`).
* `Polygon::contains(p)` is a robust crossing count; the `on-boundary` case returns `true` (a point
  exactly on the boundary is treated as inside, so endpoints may never sit on an obstacle).
* `Rect`, `is_rectilinear` and `rectilinear_rectangles` support the margin model (§3.3).

**Resolution.** `Polygon` is wound counter-clockwise and the header comment says so; the draft
also implied this in a type comment. Rings with repeated consecutive vertices, zero-length edges,
or fewer than 3 distinct vertices are rejected at construction (§7, `InvalidObstacle`). Collinear
sequences of vertices are *not* rejected — they are common and harmless.

### 6.2 `decomp`: vertical decomposition of the free space

**This replaces the draft's constrained Delaunay triangulation, and that is a deliberate
change of Stage 1.** Reasons, in order of weight:

1. A CDT requires *non-crossing constraint edges*. Overlapping obstacles — one box partly over
   another, which the test suite must accept (section 9.2) — have constraint edges that cross in
   their interiors. CDT then needs constraint *recovery*: split every crossed edge at its
   intersections, delete the crossed triangles, retriangulate each cavity. That is the single
   largest source of fragility in a naive implementation, and it is exactly the class of
   failure the draft's own "Graphviz crasher" test is about.
2. The funnel algorithm needs convex regions with *shared portals*. Triangles give convexity but
   a CDT's edge set is not organised for portal sequences, and the corridor then has to be
   recovered from triangles. A vertical decomposition produces the funnel's input directly.
3. Trapezoidal cells make the interior/exterior classification unnecessary: a cell is free space
   by construction, so there is no point-in-polygon labelling to get wrong.

The cost: cells are not Delaunay, and the "mesh" is a list of convex trapezoids rather than a
half-edge triangle mesh. Neither matters for the stated goals (section 1.1).

**Construction.** A sweep in `$+x`:

* Every edge of every offset obstacle is a *boundary edge*. The four edges of the eroded
  workspace are added as boundary edges with the region outside them marked forbidden, so the
  workspace bound is handled by the same machinery as obstacles.
* Events are the `$x`-coordinates of all boundary-edge endpoints, sorted by `orient2d`-based
  lexicographic `$(x, y)$` compare (not by a raw `f64` tuple compare, and not by
  `partial_cmp().unwrap()`).
* Between two consecutive event `$x$-values the active edge set is constant. Sorting the active
  edges by their `$y$` at the slab mid-`$x$` (exact `orient2d` of the two edges' y-intercepts)
  yields the vertical order; the slabs between consecutive active edges are candidate cells.
* A candidate cell is kept iff a point strictly inside it (slab mid-`$x$`, mid-`$y$ between its
  two boundary edges, nudged by `geom::EPS`) is outside every obstacle and inside the eroded
  workspace. The nudge exists so a cell whose interior touches a boundary is classified
  consistently.
* At each event `$x$, a **portal** is emitted in every gap whose lower or upper boundary edge
  changes, or which is adjacent to a new/removed edge. A portal is the vertical segment at that
  `$x$` between the two boundary edges forming the gap. Consecutive portals are clipped against
  the actual boundary edges at that `$x` (a portal never extends into an obstacle).
* Degenerate cells (zero width, zero height, or a cell whose two boundary edges are the same
  edge) are discarded, not stored.

**Invariants (L1, section 4.1a).** Cells are convex (two lines and two verticals), have disjoint
interiors, and their union is the offset free space. The M3 gate tests this exhaustively on
small random inputs.

**Cell location.** A point is located by binary-searching the sorted event `$x$-values to get the
slab, then binary-searching that slab's active edge list by `$y`. `$O(\log n)$`, no hash map, no
grid build cost. Points outside the eroded workspace, on a boundary, or in no cell return
`None` rather than guessing.

**Storage.** `Vec<Cell>` where `Cell { edges: [EdgeRef; 2], left: Option<PortalId>,
right: Option<PortalId> }` plus a `Vec<Portal>` and adjacency derived from portal sharing. No
`Rc`, no raw pointers, no reference counting; all indices are `u32` with a documented `$2^{32}`
vertex cap that is checked, not assumed.


### 6.3 `funnel`

Two distinct algorithms; the draft blurred them.

**(a) Corridor search — `funnel::cell_search`.** A* over the adjacency graph of free-space
cells, from the cell containing `start` to the cell containing `goal`. Two cells are adjacent
iff they share a portal; the edge cost is the portal length.

* Node cost: cell centroid. Portals are vertical, so centroids differ in `$x$` only within a
  slab; cost is well defined regardless.
* Heuristic: Euclidean distance from the node centroid to the **goal point** (not the goal
  cell centroid). Admissible: no path may be shorter than the straight-line distance, and every
  path between the two points lies in the offset free space, which the centroids also lie in.
* Tie-breaking: strictly by `(f, insertion_counter)`, never by pointer identity, to keep D2.
* The open list is a binary heap of `(f, counter, cell_id)` with a `BinaryHeap` and a
  `visited: Vec<bool>`; no hash containers anywhere.
* A* returns a **cell corridor**, not a path. It is a *guide*, not an optimality claim: the
  result is optimal only within the corridor it selects, and the spec does not claim global
  optimality over homotopy classes (section 1.1).
* If `start` or `goal` is in no cell, or the two cells are in different connected components,
  the result is `NoPathFound`. Connectivity is not precomputed; it falls out of the exhausted
  open list.

**(b) String pull — `funnel::string_pull`.** Lee–Preparata funnel over the portal sequence of
that corridor.

* `$O(N)$` in the number of portals, where `$N$` is the corridor length found in (a). The
  draft attributed `$O(N)$` to "the funnel algorithm" while listing A* in the same section; A*
  is `$O(E \log V)`. Both bounds are stated here separately and both are tested by counting
  operations.
* The first and last portals are replaced by the degenerate portals `{start}` and `{goal}`
  respectively, per Lee & Preparata; this is why the funnel cannot be handed a corridor whose
  ends are not the query points, and why the goal's own cell must terminate the corridor.
* **The orientation sign depends on the direction of travel.** Mirroring `x` flips the orientation
  of every triangle the wedge inequalities are written against, so a right-to-left corridor needs
  the opposite sign. The sign is carried as an explicit per-portal `s` rather than by swapping which
  endpoint is called `left`, because splitting the meaning across two conventions is how they drift
  apart; they were wrong together here, and the result was a path that looked plausible and left
  the corridor.
* **The scan restarts at the portal that produced the committed knot, plus one** — not at the length
  of the chain. The two differ by however many times each side has been retightened, and using the
  length restarts one portal too far: a knot at a slot's end is silently dropped and the next
  segment cuts through the obstacle it was meant to go around. `tests/funnel.rs` and the
  `the_funnel_touches_both_ends_of_a_slot` unit test both pin this.
* The committed corners are the *interior* knots; `start` and `goal` are prepended and appended
  here. Omitting the start is a silent bug that looks like a path beginning at a corner.
* The result is **verified against the corridor before it is returned**, and a path that leaves it
  is refused with `NoPathFound`. This is a safety net, not the mechanism: the funnel is proven for a
  monotone sleeve and the search is arranged to produce one, but the equal-abscissa fallback can
  still hand it a sleeve the algorithm was not written for. Refusing is the only safe answer — the
  alternative is a collision. The check is `O(segments x cells)` and is a candidate for removal
  once the fallback is either fixed or removed.
* Output is the taut polyline. Its interior knots lie on cell boundary edges (L2); knots that
  are redundant (collinear or duplicated within `geom::EPS`) are dropped before emission.
* Scratch buffers for the funnel state are passed in by the caller; the inner loop performs no
  allocation. The corridor from (a) is already owned by the caller, so (b) never allocates.

### 6.4 `corridor`

**The corridor is derived from the cells the taut segment passes through, not from a fresh
clip of the obstacles.** `margin` is not applied here again (section 3.3); the draft's
"inflate the segment by `margin` and re-clip" step is deleted.

For each taut segment `(k_i, k_{i+1})`:

1. Collect the cells of the corridor whose closure the segment meets, in order. The segment
   lies in the union of these cells by L2, so the set is well defined and contiguous.
2. For each such cell, take the sub-region of the cell on the taut segment's side of the cell's
   boundary edges, bounded by the taut segment itself. This "slab" is convex, being the
   intersection of a convex cell with two half-planes.
3. The corridor for the segment is the union of those slabs' **convex hull**, computed as the
   convex hull of the union of their vertex sets, not the union of the hulls. The union of two
   convex sets is not convex; the hull of the union is, and is still a subset of the offset free
   space only if the cells are, so the hull is taken over the cells' *common* region — namely
   the hull of the portal-to-portal envelope, which the slabs all contain.

   *Resolution:* this is the one place the draft's "hull of the free-space component" is not
   enough, because the per-cell slabs are not nested. The implementation computes the hull of
   the union of slab vertices; since the slabs form a chain joined at portals and each slab is
   contained in its cell, the hull of the chain is contained in the hull of the cells, and the
   cells' union is within the offset free space. The M5 gate asserts the containment directly
   rather than relying on this argument.
* Consecutive corridors share the portal at the common knot; this is what makes
  `$A_i = C_{i-1} \cap C_i \ne \emptyset` (section 4.6) hold.
* Degenerate case: a corridor of zero area. This is legal — it is what a taut string pinned
  tightly between two obstacles looks like. It forces `$T_i = 0` there and the segment
  degenerates to a straight line (section 4.6). It is not an error.
* Convex hull construction uses the robust `orient2d` comparator, so no `f64::atan2` sort key
  and no epsilon-based sort is permitted.

### 6.5 `spline`

**Control-point construction.** A Bezier segment for knots `$k_i, k_{i+1}$` with tangents
`$T_i, T_{i+1}$` is

```
p0 = k_i
p1 = k_i + T_i / 3
p2 = k_{i+1} - T_{i+1} / 3
p3 = k_{i+1}
```

which makes `$C^1$` exact and automatic at every interior knot by *sharing* `$T_i$.

**Tangent solve.** The unknowns are the 2-vectors `$T_0..T_{n-1}$` (tangent at each knot). The
solver minimises

```
sum_i  |T_{i+1} - 2 T_i + T_{i-1}|^2  +  lambda * sum_i |T_i - d_i|^2
```

where `$d_i$ is the unconstrained seed tangent (the chord direction), and `$lambda > 0` is
`Config::tangent_bias`. Minimising this separates into two independent tridiagonal systems, one
per coordinate, each solved with the Thomas algorithm. `$lambda = 0` yields the natural
(second-difference zero) solution; `$lambda \to \infty$ yields the chord. This is the standard
minimum-curvature variational formulation; the draft's phrase "enforcing C1 (or C2) continuity"
did not say what the system actually solved.

* Endpoint tangents: if `PortConstraint::direction` is `Some(d)`, `$T_0` is fixed to
  `$|k_0 - k_1| / 3 \cdot \hat{d}$` (and symmetrically for the goal). **Resolution:** the draft
  called `normal` a "forced unit tangent" in one place and "tangent direction and magnitude" in
  another; this spec has a direction only, and magnitude is derived from knot spacing as above,
  which is what makes `$C^1$` hold.
* `$C^2$` (feature `c2`): add the constraint `$T_i$ shared with curvature continuity by
  construction and solve for curvature-continuous joints; the resulting system is still
  tridiagonal, solved with the same routine. When `c2` is off, `$C^2$` is neither claimed nor
  tested.
* Tridiagonal solve is guarded: a zero or near-zero pivot falls back to the `$d_i$ seed rather
  than producing `NaN`/`inf`.
* With `n == 2` knots there are no interior unknowns and the result is one segment whose
  tangents are `$T_0, T_1` as above.

**Containment and repair.** For each segment, all four control points must lie in the segment's
corridor (section 4.1). Since `p0`, `p3` are knots, which lie in every incident corridor, and
`p1 = k_i + T_i/3`, `p2 = k_{i+1} - T_{i+1}/3`, containment is **equivalent** to `$T_i/3` and
`T_{i+1}/3` lying in the corridors of the two segments. That is precisely the admissible-set
condition of section 4.6, so the repair is a statement about tangents, not about control points.

**Resolution.** The draft offered "subdivision (de Casteljau) *or* dampening relaxation" with no
choice rule and no termination condition. Damping is the primary mechanism because de Casteljau
subdivision of a breaching segment does not by itself move a breaching control point. de
Casteljau is retained, but for a different, well-defined purpose: **flattening**, i.e.
`CubicBezierSegment::flatten(tolerance)` for the inspection and hit-testing paths, and for the
section 9 collision test sampling.

**Repair operates on shared tangents, never on individual control points.** This is forced by
section 4.6 and is the actual answer to the draft's unstated "what if two segments share a knot
and only one of them breaches?":

```
knots k[0..n]; tangents T[0..n]; admissible A[0..n]      # A[i] = 3*(corridor_i - k_i), section 4.6
for iter in 0..Config::max_repair_iters:
    if every T[i] lies in A[i]: done
    for i in 0..n:                      # simultaneous update; order-independent
        T[i] = damp * T[i]              # damp = Config::repair_dampening, default 0.5
    for i in 0..n:                      # then project back into the admissible set
        T[i] = project_onto_convex(A[i], T[i])
# after the loop, one final projection pass, so the result is always admissible:
for i in 0..n: T[i] = project_onto_convex(A[i], T[i])
```

* `project_onto_convex` is the standard closest-point-of-convex-polygon walk using exact
  `orient2d`; it is idempotent and needs no tolerance.
* Because damping shrinks `$|T|$` monotonically and `A[i]` contains `$0$`, the iteration
  converges; `max_repair_iters` bounds it and the final projection makes the bound irrelevant
  to correctness. **The routine cannot fail and cannot return a breaching segment.**
* `$C^1$ is preserved exactly at every step`, because `$T_i$ is a single shared value, not a
  per-segment quantity. The knots `$p_0, p_3` of each segment are the knots themselves and are
  never moved.
* The chord fallback (`$T_i = 0` everywhere) corresponds to a straight line between consecutive
  knots and is reachable by setting `Config::repair_dampening` to a small value; it is not a
  separate code path.
* Default `max_repair_iters = 16`, `repair_dampening = 0.5`. Both are `Config` fields;
  `repair_dampening` must lie in `(0, 1)` and `max_repair_iters > 0`, else `InvalidConfig`.


### 6.6 `spline::simd` (feature `simd`)

Optional. Contiguous-batch cubic evaluation and flattening over `&[CubicBezierSegment]`.
Must be bitwise-equal to the scalar path (§4.4). The module is the crate's only `unsafe`.

---

## 7. Public API

```rust
#![no_std]
extern crate alloc;

pub use alloc::vec::Vec;

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug)]
pub struct Polygon {
    /// Counter-clockwise, >= 3 distinct vertices, no repeated consecutive vertex.
    pub vertices: Vec<Point2D>,
}

impl Polygon {
    pub fn new(vertices: Vec<Point2D>) -> Result<Self, PathPlanError>;
    pub fn winding_number_contains(&self, p: Point2D) -> bool;
    pub fn bounds(&self) -> BoundingBox;
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingBox {
    pub min: Point2D,
    pub max: Point2D,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortConstraint {
    pub point: Point2D,
    /// Forced tangent *direction*; magnitude is derived from knot spacing (§6.5).
    /// Must be finite and non-zero if `Some`.
    pub direction: Option<Point2D>,
}

#[derive(Clone, Copy, Debug)]
pub struct Config {
    /// Non-negative, finite. Default 0.0.
    pub margin: f64,
    /// Tangent-solve bias lambda. Default 1.0. Must be >= 0 and finite.
    pub tangent_bias: f64,
    /// Outer repair iterations. Default 16. Must be > 0.
    pub max_repair_iters: u32,
    /// Per-iteration damping. Default 0.5. Must be in (0, 1).
    pub repair_dampening: f64,
    /// Arc-segment cap per offset corner. Default 64. Must be > 0.
    pub max_arc_segments: u32,
}

#[derive(Clone, Debug)]
pub struct RouteRequest {
    pub workspace: BoundingBox,
    pub start: PortConstraint,
    pub goal: PortConstraint,
    pub obstacles: Vec<Polygon>,
    pub config: Config,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicBezierSegment {
    pub p0: Point2D,
    pub p1: Point2D,
    pub p2: Point2D,
    pub p3: Point2D,
}

impl CubicBezierSegment {
    pub fn evaluate(&self, t: f64) -> Point2D;
    pub fn tangent(&self, t: f64) -> Point2D;
    pub fn bounding_box(&self) -> BoundingBox;
    /// De Casteljau subdivision; returns the two half segments.
    pub fn split(&self, t: f64) -> (Self, Self);
    /// Polyline approximation with max chord deviation <= tolerance.
    pub fn flatten(&self, tolerance: f64) -> Vec<Point2D>;
}

#[derive(Clone, Debug, PartialEq)]
pub enum Route {
    Smooth(Vec<CubicBezierSegment>),
    Orthogonal(OrthogonalPolyline),
}

#[derive(Clone, Debug, PartialEq)]
pub struct OrthogonalPolyline {
    /// Axis-aligned: consecutive points share x or y exactly.
    pub points: Vec<Point2D>,
}

pub struct PathPlanner { /* scratch arenas */ }

impl PathPlanner {
    pub fn new() -> Self;

    pub fn route_smooth(&mut self, req: &RouteRequest) -> Result<Vec<CubicBezierSegment>, PathPlanError>;
    pub fn route_orthogonal(&mut self, req: &RouteRequest) -> Result<OrthogonalPolyline, PathPlanError>;
    pub fn route(&mut self, req: &RouteRequest, kind: RouteKind) -> Result<Route, PathPlanError>;
}
```

`PathPlanner` is `Send`; it is deliberately **not** `Sync` (it holds mutable scratch).
`Default` is implemented.

### 7.1 `PathPlanError`

```rust
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum PathPlanError {
    #[error("start or goal point is not in free space (inside or on an obstacle)")]
    EndpointInObstacle { point: Point2D },

    #[error("start or goal point lies outside the workspace bounds")]
    EndpointOutsideWorkspace { point: Point2D },

    #[error("start and goal are the same point")]
    DegenerateEndpoints,

    #[error("no feasible path exists between start and goal")]
    NoPathFound,

    #[error("workspace bounds are empty or inverted")]
    DegenerateWorkspace,


    #[error("margin {margin} requires rectilinear obstacles (axis-aligned edges only)")]
    MarginUnsupportedGeometry { margin: f64 },

    #[error("invalid obstacle: {reason:?}")]
    InvalidObstacle { reason: InvalidObstacleReason },

    #[error("invalid configuration: {0}")]
    InvalidConfig(&'static str),
}
```

`InvalidObstacle` carries a `#[non_exhaustive]` fieldless enum:

```rust
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InvalidObstacleReason {
    TooFewVertices,     // < 3 distinct vertices
    RepeatedVertex,     // a vertex repeated anywhere in the ring, not only consecutively
    ZeroLengthEdge,     // consecutive vertices are bitwise equal
    NaNOrInfinite,      // any coordinate is not finite
    SelfIntersecting,   // a non-adjacent edge pair crosses
    DegenerateArea,     // zero total enclosed area
}
```

**Resolution.** The draft's `InvalidObstacle(String)` allocated in a `no_std` crate and was not
matchable in tests; the fieldless enum is not.

`WorkspaceEroded` and `MarginTooLarge` from my earlier draft are **removed**: the first
became unnecessary when Stage 1 stopped needing non-crossing constraints (section 6.2), and the
second was wrong — an obstacle may extend arbitrarily far outside the workspace, because
section 3.3 clips the *offset* obstacle to the eroded workspace. Only a non-finite, degenerate,
or self-intersecting ring is rejected.


### 7.2 `route_orthogonal`

`route_orthogonal` builds the Hanan grid induced by the x- and y-coordinates of the **offset**
obstacles' vertices (section 3.3), clipped to the eroded workspace; performs rectilinear A* over
the grid cells (4-neighbourhood, cost = Manhattan distance) with a free cell iff its centre is
outside every offset obstacle; and emits the resulting cell path as an axis-aligned polyline.

Guarantees: every emitted point is inside the eroded workspace and outside every offset obstacle;
consecutive points differ on exactly one axis; `points.first() == start.point` and
`points.last() == goal.point`; interior collinear points are merged. If no 4-connected grid path
exists the result is `NoPathFound` — **there is no `u`-shaped fallback**, because any such
fallback would be a guess about topology and would be wrong whenever the sealed region is not
rectangular. The grid is the approximation, and its limits are the error.

**This is a separate algorithm.** It does not use the decomposition, the funnel, the corridors,
or the spline. It shares only the margin model of section 3.3 and the ingest validation of
section 6.1. It is not a post-processing of the smooth route.


---

## 8. Performance targets

1. **Allocation.** After the first query on a given `PathPlanner` with a given
   `(workspace, obstacle-count, config)` shape, subsequent queries perform **zero heap
   allocations** on the success path, verified by a counting global allocator in the test
   suite. The first query of a shape, and any query whose result needs more capacity than the
   arena holds, may allocate. The draft's unqualified "zero persistent heap allocations" is
   otherwise untestable.
2. **Latency.** Warm query, 10–50 axis-aligned box obstacles, `margin = 1.0`, start/goal
   forcing 2–4 spline segments: **< 50 microseconds** median on a modern x86_64 or aarch64
   core, measured by Criterion with `--warm-up-time 1`.
   *Status: **not met, and the target is not currently achievable.** Measured decomposition-only
   cost (Stages 1a+1b, the whole of M3) on an aarch64 core, 100x100 workspace, boxes on a lattice:
   10 boxes 74 us, 50 boxes 227 us, at either `margin` 0 or 1. The margin is now nearly free
   because square inflation adds no vertices; the cost is the sweep, which is `O(events x edges)`.
   The end-to-end figure will be worse. The target is retained as the goal and the implementation
   is expected to need a sweep with an incrementally maintained active-edge list before it can be
   approached; `tests/scale.rs` prints the current figures so the gap stays visible.
3. **Benchmarks.**
   * End-to-end `route_smooth` at N = 10/50/200 boxes.
   * Stage 1 cost split: `inflate`+`clip` vs. `sweep` vs. cell location.
   * `string_pull` throughput in portals/`ns` with an operation counter assertion.
   * `solver` + `containment` (damping repair) throughput.
   * Scalar vs. `simd` batched evaluation, with a bitwise-equality assertion.
4. **Out of scope for v1.** Incremental or cached decomposition updates. The draft's "full CDT
   generation vs incremental updates" benchmark implies an incremental API that does not exist in
   the API surface. It is deferred; no incremental API is promised.

---

## 9. Verification

All tests run with `--features std` unless noted. The suite is a *test strategy layered under*
the section 4 invariants, not the source of the guarantee.

### 9.0 Milestone gates

| Gate | Property | Test |
| --- | --- | --- |
| M3 | L1: cells convex, disjoint interiors, union = the inflated free space; plus the one-sided margin property on cell boundaries and chord interiors | `tests/decomposition.rs`: 200x200 probe grid against an independent free-space predicate, for 15 hand-built fixtures and 24 random box fields |
| M4 | L2: every funnel output point is inside some corridor cell | 9.1's containment check, restricted to the returned corridor |
| M5 | each corridor is convex and contained in the offset free space | exact `orient2d` against every offset obstacle edge |
| M6 | section 4.1-4.3 | 9.1, 9.3, 9.4, 9.5 |

### 9.1 Randomized obstacle field (`proptest`)

* Generate `$N \in [1,100]` **overlapping or not** axis-aligned boxes inside the workspace, plus
  random `start`/`goal` (rejection-sampled against the offset domain, not the raw one), and
  random `margin \in [0, 2]`.
* Assert `route_smooth` returns `Ok` **or** one of `NoPathFound` / `EndpointInObstacle` /
  `EndpointOutsideWorkspace` / `InvalidObstacle` / `MarginUnsupportedGeometry`. Never a panic, never an
  undocumented variant.
* Collision check: flatten every segment with tolerance `1e-3` **and** independently sample at
  `dt = 0.005`; assert every sample point is outside every obstacle *by at least* `margin` minus
  `1e-9`, and inside the workspace eroded by `margin` minus `1e-9`. No percentage slack is needed
  now that the margin is realised by exact square inflation rather than an arc approximation; the
  `1e-9` absorbs `f64` distance rounding only. Inputs with a non-rectilinear obstacle and a
  nonzero `margin` are excluded, since that is `MarginUnsupportedGeometry` rather than a route.
* Additionally assert, for the same input, that **the exact containment condition of section
  4.1 holds**: `k_i \pm T_i/3` lies in the admissible set of section 4.6. This is the check that
  has no tolerance at all, and it is what makes the sampled check trustworthy.

### 9.2 Collinear / touching / overlapping obstacles

Boxes sharing edges, sharing corners, nested, duplicated, mutually overlapping, obstacle edges
collinear with the workspace bbox, and `margin` large enough to merge obstacles together. Assert:
no panic, and the result is `Ok` or a documented error. The draft's "without returning an error"
is not achievable — genuinely sealed boxes must produce `NoPathFound`, and asserting otherwise
would force the implementation to lie.

Also assert the *one-sided* margin property directly and with no tolerance ladder: for every point
of every cell boundary, and of every chord between consecutive boundary points, the distance to
every obstacle ring is `>= margin`.

### 9.3 Straight-line reduction

Empty obstacles, `start != goal`, no `direction` constraints, any `config`. Assert:

* exactly one segment returned,
* `orient2d(start, goal, p1) == Collinear` and `orient2d(start, goal, p2) == Collinear`,
* the chord-length parameterisation is monotone: `|p1-p0| == |p2-p1| == |p3-p2| / 2` within
  `1e-9` (the natural-spline solution on two knots),
* `max over 1000 samples of distance(sample, segment(start,goal)) <= 1e-9`.

### 9.4 Continuity

For every adjacent pair, assert `S_k(1) == S_{k+1}(0)` bitwise and
`normalize(S_k'(1)) ≈ normalize(S_{k+1}'(0))` within `1e-9`, and
`norm(S_k'(1) - S_{k+1}'(0)) <= 1e-12 * max(1, norm(S_k'(1)))`.

### 9.5 Control-hull containment

For every segment, assert all four control points lie in the segment's corridor, using exact
`orient2d` against the corridor's edges. This is the assertion that actually mirrors section 4.1;
the sampling in 9.1 is a cross-check against it.

### 9.6 Determinism

For a fixed corpus: 200 runs of the same request produce bit-identical output (D2), and the
`scalar` and `simd` builds produce bit-identical output (when `simd` is enabled).

### 9.7 `no_std` and unsafe

* `cargo build --no-default-features --target thumbv7em-none-eabihf` (or
  `aarch64-unknown-none`) must succeed.
* A CI job greps for `unsafe` and fails unless every hit is under `src/spline/simd.rs`.

### 9.8 Orthogonal router

Assert every consecutive point pair in `OrthogonalPolyline::points` differs on exactly one axis,
all points are inside the eroded workspace and outside every offset obstacle, and
`points[0] == start.point`, `last == goal.point`.


---

## 10. Definition of done

* `cargo build` succeeds with `--no-default-features`, with `std`, and with `simd`.
* `cargo clippy --all-features -- -D warnings` is clean.
* `cargo test --all-features` passes; `proptest` runs 256 cases per property in CI.
* `cargo bench` produces the five benchmark groups in section 8.3.
* Every claim in section 4 has a test in section 9 referencing it.

---

## Appendix A — Draft ambiguities and their resolutions

Part 1 is the review of `SPEC_DRAFT.md`. Part 2 is the second review, of this document's own
first draft, which found three remaining gaps (one internal contradiction, one unresolved
product question, one undecided algorithm) and closed them.

### A.1 From the original draft

| # | Draft text | Problem | Resolution |
| --- | --- | --- | --- |
| 1 | "Rust 2024 / 2021" | Edition cannot be both | Rust 2024 (§ header) |
| 2 | `margin`; safety `>= 0` | `>= 0` is vacuously true for any curve outside the polygon | §4.1: distance to obstacle boundary `>= margin`, guaranteed via corridors |
| 3 | `$C^1$` written as `$S_k'(1) = \alpha S_{k+1}'(0)$ | That is `$G^1$`; contradicts the stated `$C^1$` | §4.2: exact vector equality; direction-only ports are rescaled |
| 4 | "Monotonic Convergence" | No iteration exists to converge | Renamed *Straight-line reduction* (§4.3) |
| 5 | "bitwise identical across all architectures" | Contradicts SIMD and FMA | Split into D1/D2 (§4.4) |
| 6 | "zero unsafe outside SIMD" vs. exact predicates | Unclear whether predicates need `unsafe` | §4.5: they don't; `forbid` + one `allow` site |
| 7 | `PathPlanError::InvalidObstacle(String)` | Allocates in `no_std`; untestable | Fieldless reason enum (§7.1) |
| 8 | `Vec` used unqualified in a `no_std` API sketch | Missing `alloc` | `extern crate alloc` + re-export (§7) |
| 9 | A* "metric is centroids **or** edge midpoints" | Undefined choice | Centroids + admissibility proof (§6.3a) |
| 10 | A* and funnel both `$O(N)$ | A* is not `$O(N)$ | `$O(E log V)$` and `$O(N)$` stated separately (§6.3) |
| 11 | Bbox is "a bounding box" input, but free-space coverage unstated | Lookup can fail on uncovered regions | Workspace bound handled as boundary edges; full free-space partition (L1) |
| 12 | "tridiagonal system enforcing C1 (or C2)" | No objective stated | Minimum-curvature variational form, λ = `tangent_bias` (§6.5) |
| 13 | "hull property: *strictly* contained" | Hull property is non-strict | §4.1 states non-strict containment |
| 14 | Repair: "subdivision **or** dampening" | No rule, no termination, no failure mode | Damping on shared tangents + projection, bounded, provably cannot fail (§6.5, §4.6) |
| 15 | `PortConstraint::normal: Option` | Untyped; called both "unit" and "direction + magnitude"; "node bounding boxes" undefined | `direction: Option<Point2D>`, direction only, magnitude derived (§6.5, §7) |
| 16 | Pipeline output "Vec or OrthogonalPolyline" | Orthogonal algorithm never defined | §7.2, separate Hanan-grid planner |
| 17 | Test 2: "assert ... without returning an error" | Unreachable for sealed boxes | Assert no panic + documented-error set (§9.2) |
| 18 | Test 1: sampling `dt = 0.005` as the safety check | Sampling cannot prove safety | Containment is the guarantee; sampling is a cross-check (§9.1, §9.5) |
| 19 | Test 3: "equivalent to a straight line within ε" | Metric undefined | Collinearity + monotone params + 1e-9 max deviation (§9.3) |
| 20 | "under 50 µs typical" | Unpinned corpus; cold vs. warm | Corpus pinned, warm median, Criterion (§8.2) |
| 21 | Benchmark "full CDT vs incremental" | No incremental API exists | Deferred to v2 (§8.4) |
| 22 | Draft module tree | No home for `PathPlanError` or scratch buffers | `error.rs`, `router/scratch.rs` added (§5) |
| 23 | Missing error variants | Endpoint-outside-workspace, degenerate endpoints, invalid config unreported | Added (§7.1) |
| 24 | No termination/correctness criteria section | — | §10 added |

### A.2 From the second review

| # | Gap | Problem | Resolution |
| --- | --- | --- | --- |
| 25 | **Margin mechanism contradiction** | §4.1 said corridors are free space "eroded by `margin`" while §6.4 implemented margin by inflating the taut segment along its normal and re-clipping. Two different mechanisms; they disagree at corners, and a bug in either would silently produce a margin violation | **§3.3, new.** One mechanism, applied once, at Stage 1a: erode the workspace, offset every obstacle, clip to the eroded box, and never touch `margin` again. Arc approximation and its one-sided error bound are specified there |
| 26 | **Product framing unresolved** | "Replacement for `libpathplan`" implied Graphviz node/port/rank semantics, which §6.5 had quietly dropped. Acceptance criteria for "good" were unstated, and there was no scope cut-line | **§1.1, new.** v1 is a geometry library, not a drop-in router. An explicit non-goal table, an explicit "not global optimum over homotopy classes", and an 8-milestone order with per-milestone gates, including the statement that M1–M6 alone ship a usable library |
| 27a | **CDT vs. overlapping obstacles** | A CDT requires non-crossing constraint edges, so overlapping or edge-sharing obstacles need constraint *recovery* (split crossed edges, retriangulate cavities). §9.2 *requires* those inputs to work. The draft's chosen Stage 1 and its own test suite were in direct conflict, and the draft never acknowledged it | **§6.2, rewritten.** Stage 1 is a vertical (trapezoidal) sweep decomposition. Crossing edges are a non-event; cells come out convex, which is exactly the funnel's precondition; and interior/exterior labelling disappears. `cdt/` is renamed `decomp/`. This is the one place this spec departs from the draft's named algorithm, and it is argued rather than assumed |
| 27b | **Coincident-vertex snapping** | Would have needed a merge rule for a CDT. With a sweep it is still needed for the event ordering, and the rule is now stated | §6.2: exact `orient2d`-based lexicographic event sort, no `partial_cmp().unwrap()`, and `geom::EPS`-nudged cell classification |
| 27c | **Shared-tangent clamp** | A knot's tangent serves two segments with two different corridors; the draft had repair acting per segment, which cannot preserve `$C^1$` and was therefore impossible as written | **§4.6, new.** The admissible set is `$A_i = 3\,(C_{i-1} \cap C_i - k_i)$, a convex polygon containing the origin; repair is damping + closest-point projection on that set, which is idempotent, convergent, and `$C^1$`-preserving by construction |
| 27d | **Funnel validity** | The draft's funnel assumed the A* corridor was a valid funnel input without saying so | **§4.7 L2, new**, with an M4 gate that tests it directly rather than relying on the argument |
| 27e | **Orthogonal `u`-fallback** | My own first draft's fallback would have guessed topology and been wrong for non-rectangular sealed regions | Deleted (§7.2); `NoPathFound` is returned instead |
| 27f | **`OutsideWorkspace` / `ConstraintEdgesCross` errors** | Both became wrong or unnecessary once `margin` was applied at ingest and Stage 1 stopped needing crossing-free constraints | Removed from the error enum, with the reason stated (§7.1) |
| 28 | `#![forbid]` at the root plus `#![allow]` in a module | Not compilable: an inner `allow` cannot override an outer `forbid` | §4.5 uses `deny` plus the CI grep gate |
| 29 | Arc bound `2 arccos(1 - r)` | Guarantees the arc's *endpoints* clear the margin, not the interior of each chord, which passes at `R cos(step/2)` — up to `margin * r` too close at every corner. Found by the M2 chord-interior test | §3.3: `2 arccos(1/(1+r))`, and the ratio raised to `1e-2` now that the bound is exact rather than approximate |
| 30 | `InvalidObstacleReason` had no variant for "the offset ring self-intersects" | The only options were to drop the obstacle — which is a route through it — or to add a variant | `MarginTooLarge` added, then removed again once §3.3 stopped needing offsets (A.32, A.33) |
| 31 | §3.3 claimed the erosion "swallows" obstacles near a wall | It does not: an obstacle at distance `d` survives iff `d >= 2 * margin` | §3.3 states the actual condition; a regression test pins it |
| 32 | §3.3's original "offset each obstacle by `margin`, approximating convex corners with arcs" | **Unsound.** The offset boundary of a *reflex* corner is not the crossing of the two offset lines and not an arc about the vertex, and the failure is not a small inaccuracy: on a comb-shaped obstacle the construction put offset-boundary points at zero distance from the obstacle on 854 sampled positions, and whether it failed was non-monotone in `margin` (0.4, 0.5 and 0.6 failed while 0.8 and 1.0 succeeded). Silently violating the clearance is the worst failure mode in the whole crate | §3.3 rewritten: the growth is a **square inflation** applied to a rectilinear partition of each obstacle, which is exact. A nonzero `margin` against a non-rectilinear obstacle is now a reported error rather than a wrong answer. The alternatives (eroding cells; trimming offsets) are recorded in §3.3 with the reason each fails |
| 33 | §7.1's `MarginTooLarge` and `WorkspaceEroded` variants | Both existed only to paper over the offsetting failures; with exact square inflation neither condition can arise | Removed, with the reason stated (§7.1) |
| 34 | §8.2's "< 50 us" target | Not met, and the measurement is recorded rather than the target quietly dropped | §8.2 now carries the measured decomposition cost (74 us at 10 boxes, 227 us at 50) and names the optimisation the target depends on |
| 35 | A cell has one left portal and one right portal | A cell's side can overlap *several* cells on the other side, so the portal is overwritten and the adjacency is lost | §6.3a: both sides are lists, held in two compressed lists so there is no per-cell allocation |
| 36 | The funnel's portal `left`/`right` were assigned once, globally, from the corridor's overall direction | A corridor is not necessarily monotone in x — the cell graph fans out — so half the portals can be entered in the opposite direction to the other half, and the funnel then commits knots from the wrong chain. The result is a path that looks plausible and leaves the corridor | §6.3b: the orientation is carried per portal, and the corridor is made monotone by construction (§6.3a) so the case does not arise in the first place |
| 37 | The funnel's scan restarted at the length of the tightened chain | Off by however many times the sides were retightened, i.e. one portal too far. A knot at a slot's end is dropped and the next segment crosses the obstacle | §6.3b: the sides record the portal that *produced* them, and the scan restarts at that portal plus one. Pinned by `the_funnel_touches_both_ends_of_a_slot` |
| 38 | The M4 gate asserted the taut path is in the free space | Wrong assertion: a taut path hugs the obstacles it wraps around, so with `margin == 0` its knots lie exactly on obstacle boundaries. A containment test rejects every useful route | §9.0: a distance-based clearance test (`>= margin`), which is the actual guarantee. L2, the corridor-membership test, stays a containment test because a cell boundary is not an obstacle |

