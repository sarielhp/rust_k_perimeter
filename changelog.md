# Changelog

All notable changes to this project will be documented in this file.

## [0.2.2] - 2026-07-19

### Added
- **Multi-page PDF Information Section**: PDFs generated for figures now contain all standard output log information (execution stages, timing metrics, Pick's theorem verification, and final statistics) starting on page 2.
- **Color & Polygon Legend**: Added an explanation header at the top of the PDF info pages detailing all rendered elements:
  - **Optimal Solution Polygon**: Blue fill with black boundary (line width 1.0).
  - **Baseline Disk Polygon (`ch_m`)**: Red/orange outlined baseline convex hull.
  - **Expanded Disk Polygon (`ch_m_exp`)**: Green outlined outer search boundary.
  - **Yellow Dot**: Origin $(0,0)$.
  - **Blue Dots**: "Good" grid points (candidate boundary vertices).
  - **Red Dots**: Interior ("bad set") and exterior grid points.
- **Stdout Logger (`logger.rs`)**: Introduced a global log collector module (`logger.rs` and `log_println!` macro) to capture standard output during execution for PDF generation.
- **Multi-page Pagination**: Dynamic page creation and text wrapping in Cairo (`render_text_pages`) to ensure log output flowing across multiple pages is not truncated.

### Changed
- **Optimal Polygon Styling**: Optimal solution polygon boundary is now rendered in black with a line width of `1.0` (previously blue with width `2.0`).
- **Rendering Order**: Optimal polygon is drawn first, followed by estimated bound polygons (`ch_m` and `ch_m_exp`) on top so they remain visible, followed by grid points on top of all polygons.
- **Red Grid Point Filtering**: Red grid points are now rendered in the PDF only if they are within Euclidean distance $\le 5$ of a differently colored grid point (yellow origin or blue good set point) or within distance $\le 1$ of the grid bounding box boundary.
- **Cleaned-up Polygon Storage & Linear $O(V)$ Algorithm**: Implemented a single-pass $O(V)$ stack-based cleanup algorithm (`polygon_rm_redundant_vertices` in `geom.rs`) that pops redundant collinear vertices as points are inserted and resolves wraparound boundary collinearity. Retained the legacy multi-pass algorithm (`polygon_rm_redundant_vertices_old`) and added unit tests asserting exact equivalence across all test cases.
