//! Scale and determinism of the decomposition. Not an assertion of the latency target: see
//! `docs/SPEC.md` section 8.2 and the note there about the measured figures.

use clearpath::decomp::sweep::decompose;
use clearpath::geom::{BoundingBox, FreeSpace, Point2D, Polygon};
use std::time::Instant;

#[test]
fn scale_and_determinism() {
    for n in [0usize, 10, 50] {
        for margin in [0.0f64, 1.0] {
            let obstacles: Vec<Polygon> = (0..n)
                .map(|i| {
                    let x = 5.0 + (i % 7) as f64 * 13.0;
                    let y = 5.0 + (i / 7) as f64 * 13.0;
                    Polygon::new(vec![
                        Point2D::new(x, y),
                        Point2D::new(x + 6.0, y),
                        Point2D::new(x + 6.0, y + 6.0),
                        Point2D::new(x, y + 6.0),
                    ])
                    .unwrap()
                })
                .collect();
            let ws = BoundingBox::new(Point2D::ZERO, Point2D::new(100.0, 100.0));
            let space = FreeSpace::new(ws, obstacles).unwrap();
            let t = Instant::now();
            let d = decompose(&space, margin).unwrap();
            let elapsed = t.elapsed();
            let again = decompose(&space, margin).unwrap();
            let same = d.cells().len() == again.cells().len()
                && d
                    .cells()
                    .iter()
                    .zip(again.cells())
                    .all(|(a, b)| a.corners() == b.corners());
            println!(
                "n={n:3} margin={margin}: {:>8.0?}  cells {:5} portals {:5}  deterministic={same}",
                elapsed,
                d.cell_count(),
                d.portal_count()
            );
            assert!(same, "non-deterministic decomposition");
        }
    }
}
