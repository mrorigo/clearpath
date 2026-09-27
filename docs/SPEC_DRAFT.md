# Specification: `pathplan` (Obstacle-Avoidance Path Planning & Spline Engine)

**Crate:** `pathplan`

**License:** MIT OR Apache-2.0

**Language Edition:** Rust 2024 / 2021

**Target:** `#![no_std]` core with optional `std` and SIMD features. Zero unsafe code outside strictly audited SIMD intrinsics.

---

## 1. Executive Summary

`pathplan` is a high-performance, deterministic 2D obstacle-avoidance path planning and smooth spline fitting library. It serves as a modern, memory-safe, mathematically robust replacement for Graphviz's legacy `libpathplan`.

Given a set of polygonal obstacles, start/end endpoints, and optional boundary departure/arrival tangent constraints, `pathplan` calculates a collision-free, $C^1$-continuous cubic Bézier spline (or alternatively, an orthogonal Manhattan route) that navigates the free space without clipping obstacle boundaries.

---

## 2. Core Architecture Pipeline

The pipeline operates in four decoupled, sequential stages:

```
[Polygonal Obstacles + Waypoint Queries]
                   │
                   ▼
  Stage 1: Constrained Delaunay Triangulation (CDT) & Mesh Setup
                   │  (Computes free-space planar decomposition)
                   ▼
  Stage 2: Discrete Route & Funnel Algorithm
                   │  (Extracts minimum-turn taut-string polygonal path)
                   ▼
  Stage 3: Convex Corridor Generation
                   │  (Derives maximal open convex polygon channel)
                   ▼
  Stage 4: Constrained Cubic Bézier Spline Fitting
                   │  (Solves Thomas algorithm + iterative hull containment)
                   ▼
[Vec or OrthogonalPolyline]

```

---

## 3. Mathematical Foundations & Invariants

### 3.1. Robust Geometric Predicates

Floating-point rounding errors are catastrophic in computational geometry.

* All orientation tests ($\text{ccw}(a, b, c)$) and in-circle tests **must** use adaptive precision floating-point arithmetic (Shewchuk-style exact predicates) or symbolic perturbation to handle collinear points, coincident vertices, and zero-width corridors deterministically.
* Never use raw naive floating-point determinants for topological decisions.

### 3.2. Geometric Invariants

1. **Safety Guarantee:** For any parameter $t \in [0, 1]$ across all generated spline segments $S_i(t)$, the minimum Euclidean distance to any obstacle polygon boundary must be $\ge 0$ (or $\ge \text{margin}$ if margin $> 0$).
2. **Smoothness:** Spline joints must maintain $C^1$ continuity:

$$S_k(1) = S_{k+1}(0) \quad \text{and} \quad S_k'(1) = \alpha \cdot S_{k+1}'(0) \quad (\alpha > 0)$$

3. **Monotonic Convergence:** In the absence of obstacles, the routed path between two unconstrained points must reduce to a single straight-line segment.
4. **Determinism:** Given identical inputs, outputs must be bitwise identical across all supported architectures.

---

## 4. Crate Structure & Modules

```
pathplan/
├── Cargo.toml
├── src/
│   ├── lib.rs                  # Public facade, high-level routing API
│   ├── geom/                   # Low-level geometric primitives & predicates
│   │   ├── mod.rs
│   │   ├── point.rs            # Vec2 / Point2D definitions
│   │   ├── polygon.rs          # Polygon, BoundingBox, winding validation
│   │   └── predicates.rs       # Robust exact orientation & incircle tests
│   ├── cdt/                    # Constrained Delaunay Triangulation
│   │   ├── mod.rs
│   │   ├── mesh.rs             # Half-edge or quad-edge mesh representation
│   │   └── triangulate.rs      # CDT builder with constraint edge insertion
│   ├── funnel/                 # Shortest path & taut string
│   │   ├── mod.rs
│   │   ├── dual_graph.rs       # A* search through triangle adjacency
│   │   └── string_pull.rs      # Lee & Preparata Funnel algorithm
│   ├── corridor/               # Channel expansion
│   │   ├── mod.rs
│   │   └── convex_hull.rs      # Maximal convex corridor builder
│   ├── spline/                 # Spline generation & constraint solver
│   │   ├── mod.rs
│   │   ├── bezier.rs           # Cubic segment representation & evaluation
│   │   ├── solver.rs           # Tridiagonal matrix solver (Thomas algorithm)
│   │   └── containment.rs      # Convex hull & SAT obstacle collision check
│   └── router/                 # High-level orchestrators
│       ├── mod.rs
│       ├── smooth.rs           # Smooth spline pipeline
│       └── orthogonal.rs       # Manhattan / right-angle router

```

---

## 5. Module Technical Specifications

### 5.1. `geom`: Primitives & Predicates

* **`Point2D`**: Standard 2D vector with `Add`, `Sub`, `Mul`, dot product, cross product, and normalization.
* **`predicates`**:
* `orient2d(pa: Point2D, pb: Point2D, pc: Point2D) -> Orientation`
Returns `Clockwise`, `CounterClockwise`, or `Collinear`. Must perform fast-filter floating-point check first, escalating to exact extended-precision arithmetic on near-zero results.



### 5.2. `cdt`: Constrained Delaunay Triangulation

* Accepts a bounding box and a list of closed `Polygon` obstacles.
* Inserts obstacle edges as **hard constraints**.
* Marks triangles as `Interior` (inside obstacle) or `Traversable` (free space).
* Storage: Compact index-based half-edge structure using `Vec` instead of pointer-heavy nodes to ensure cache locality and zero reference counting.

### 5.3. `funnel`: Shortest Path & String Pulling

* **A* Search:** Traverses the dual graph of traversable triangles from the triangle containing `start` to the triangle containing `goal`. Distance metric is Euclidean distance between triangle centroids or edge midpoints.
* **Funnel Algorithm (Lee & Preparata, 1984):**
* Inputs: Sequence of portal edges $(L_i, R_i)$ connecting adjacent traversable triangles.
* Output: Minimal polygonal chain (taut string) connecting start to goal.
* Invariant: Runs in strictly $O(N)$ time where $N$ is the number of portal edges.
* Zero heap allocations inside the hot loop: pass in reusable scratch buffers.



### 5.4. `corridor`: Convex Channel Construction

* For each linear segment of the taut-string path, expand a sequence of overlapping convex polygons (corridors) enclosing the segment without intersecting any obstacle boundary.
* Each corridor shares a boundary portal with the next corridor, defining valid slack space for Bézier control points.

### 5.5. `spline`: Constrained Fitting & Smoothing

* **Initial Fitting:** Parameterize path points by chord length. Set up a tridiagonal system for internal Bézier control points enforcing $C^1$ (or $C^2$) continuity.
* **Boundary Conditions:**
* Support optional `PortConstraint` defining the tangent direction and magnitude at start and end points (e.g., perpendicular exit from node bounding boxes).


* **Convex Hull Property Validation:**
* A cubic Bézier curve is strictly contained within the convex hull of its four control points $P_0, P_1, P_2, P_3$.
* Verify that the control polygon of each segment lies entirely within its corresponding convex corridor.
* If a control point breaches corridor bounds: apply subdivision (de Casteljau's algorithm) or dampening relaxation until containment is met.



---

## 6. Public API Design

```rust
use std::fmt::Debug;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug)]
pub struct Polygon {
    pub vertices: Vec, // Counter-clockwise for outer boundaries
}

#[derive(Clone, Copy, Debug)]
pub struct PortConstraint {
    pub point: Point2D,
    pub normal: Option, // Optional forced unit tangent vector
}

#[derive(Clone, Debug)]
pub struct CubicBezierSegment {
    pub p0: Point2D,
    pub p1: Point2D,
    pub p2: Point2D,
    pub p3: Point2D,
}

impl CubicBezierSegment {
    pub fn evaluate(&self, t: f64) -> Point2D { ... }
    pub fn tangent(&self, t: f64) -> Point2D { ... }
    pub fn bounding_box(&self) -> (Point2D, Point2D) { ... }
}

#[derive(Clone, Debug)]
pub struct RouteRequest {
    pub start: PortConstraint,
    pub goal: PortConstraint,
    pub obstacles: Vec,
    pub margin: f64, // Clearance padding around obstacles
}

pub struct PathPlanner {
    // Reusable buffers to avoid allocations across sequential queries
    arena: ScratchBuffer,
}

impl PathPlanner {
    pub fn new() -> Self;

    /// Generates a smooth, collision-free C1 cubic spline route
    pub fn route_smooth(&mut self, req: &RouteRequest) -> Result, PathPlanError>;

    /// Generates an orthogonal (Manhattan / right-angle) route
    pub fn route_orthogonal(&mut self, req: &RouteRequest) -> Result, PathPlanError>;
}

#[derive(thiserror::Error, Debug)]
pub enum PathPlanError {
    #[error("Start or goal point resides inside an obstacle")]
    EndpointInsideObstacle(Point2D),
    #[error("No feasible path exists between start and goal")]
    NoPathFound,
    #[error("Degenerate polygon input: {0}")]
    InvalidObstacle(String),
}

```

---

## 7. Performance & Memory Targets

1. **Allocations:** Zero persistent heap allocations during query execution when using a pre-instantiated `PathPlanner` workspace. Scratch buffers must be recycled across invocations.
2. **Latency:** Under 50 microseconds for typical graph edge routes (10–50 obstacle boxes, 2–4 spline segments) on modern x86_64 / aarch64 cores.
3. **Benchmarking:** Criterion benchmarks comparing:
* Full CDT generation vs incremental updates.
* Funnel algorithm throughput.
* Thomas solver + de Casteljau subdivision throughput.



---

## 8. Verification & Property-Based Testing Directives

The implementing agent must establish the following tests using `proptest`:

1. **Randomized Obstacle Field:** Generate $N \in [1, 100]$ non-overlapping random axis-aligned bounding boxes and arbitrary start/end points in free space.
* Assert `route_smooth` returns `Ok(spline)`.
* Sample $t \in [0.0, 1.0]$ at step $\Delta t = 0.005$ across each segment; assert no point is within distance $< \text{margin}$ of any box.


2. **Collinear Edge Handling:** Generate adjacent boxes sharing common edges or corners (a notorious Graphviz crasher). Assert the planner executes without panicking or returning an error.
3. **Straight-Line Reduction:** When `obstacles` is empty and start/goal normals are collinear or `None`, assert the resulting spline is equivalent to a straight line within $\epsilon = 10^{-7}$.
4. **Continuity Verification:** Assert that for every adjacent segment pair $(S_k, S_{k+1})$, $S_k(1.0) == S_{k+1}(0.0)$ and their normalized tangent vectors match within $\epsilon = 10^{-6}$.
