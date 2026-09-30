# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `PathPlanner::route_smooth` no longer returns a curve that passes through an obstacle. A taut
  path runs along obstacle faces, and reconstructing a control point as `knot + t/3` lands a
  fraction of an ulp inside the face. `Clearance::is_free` accepted any point with a non-zero
  distance to the obstacle's ring, so such a point was reported free and the curve was shipped
  through the obstacle. The interior test now requires a point to be more than the boundary
  tolerance inside, which is the same tolerance the rest of the clearance tests use.
- `Clearance::is_free` and `Clearance::hull_is_free` derive their boundary tolerance from one
  place, so the test that certifies a route and the test that checks it cannot disagree.

## [0.1.0] — unreleased

First release. Named `clearpath`; the crate was developed as `pathplan`, after Graphviz's
`libpathplan`, and nothing was derived from it (see the README's "Origin and provenance").

### Added

- `PathPlanner::route_smooth`: a `C¹` chain of cubic Bézier segments around polygonal obstacles.
  Collision freedom is constructive at the segment ends — the tangent at each knot is clamped to
  the set that keeps the neighbouring control points inside a cell — and checked exactly over each
  segment's control hull.
- `PathPlanner::route_orthogonal`: a rectilinear polyline, over the Hanan grid of the forbidden
  faces. A cell's centre decides the cell, because no face crosses a cell's interior.
- A `margin` model applied exactly once, at ingest: the workspace eroded and every rectilinear
  obstacle grown, by a square inflation. A nonzero `margin` against a non-rectilinear obstacle is
  reported as `MarginUnsupportedGeometry` rather than approximated.
- `PortConstraint` for a forced departure or arrival direction, honoured when it is admissible and
  projected when it is not.
- `CubicBezierSegment::flatten` and `split` (de Casteljau), for sampling and hit-testing a curve.
- `no_std` + `alloc`, with no `unsafe` anywhere in the crate and no lint suppression to
  allow it.

### Known limitations

- Not allocation-free. A warm query allocates about 89 times at 10 obstacles, nearly all of it
  small `Vec`s spread across the stages; the decomposition itself allocates nothing when warm.
- Not incrementally updatable. The decomposition is rebuilt per query.
- Optimal only within the corridor the search selects, not across homotopy classes.
- The containment checks are `O(hull edges x nearby obstacle edges)`; a spatial index would be the
  next step for large obstacle counts.
