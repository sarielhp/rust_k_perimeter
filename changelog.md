# Changelog

All notable changes to this project will be documented in this file.

## [0.2.2] - 2026-07-19

### Added
- **Multi-page PDF Information Section**: PDFs generated for figures now contain the execution summary metrics starting on page 2, filtered to include only lines starting with `#`.
- **Color & Polygon Legend**: Added an explanation header at the top of the PDF info pages detailing all rendered elements:
  - **Optimal Solution Polygon**: Blue fill with black boundary (line width 1.0).
  - **Baseline Disk Polygon (`ch_m`)**: Red/orange outlined baseline convex hull.
  - **Expanded Disk Polygon (`ch_m_exp`)**: Green outlined outer search boundary.
  - **Yellow Dot**: Origin $(0,0)$.
  - **Blue Dots**: "Good" grid points (candidate boundary vertices).
  - **Red Dots**: Interior ("bad set") and exterior grid points.
- **Stdout Logger (`logger.rs`)**: Introduced a global log collector module (`logger.rs` and `log_println!` macro) to capture standard output during execution for PDF generation.
- **Comma-Separated Integer Formatting**: All integer metrics and progress counts printed to standard output and rendered in PDF info pages (such as `k`, `vertices`, `boundary grid points`, `Configs computed`, and `UB primitive edge len`) are now formatted with commas (e.g. `1,633,661` and `10,000`) for enhanced readability.
- **Multi-page Pagination**: Dynamic page creation and text wrapping in Cairo (`render_text_pages`) to ensure log output flowing across multiple pages is not truncated.

### Changed
- **Optimal Polygon Styling**: Optimal solution polygon boundary is now rendered in black with a line width of `1.0` (previously blue with width `2.0`).
- **Rendering Order**: Optimal polygon is drawn first, followed by estimated bound polygons (`ch_m` and `ch_m_exp`) on top so they remain visible, followed by grid points on top of all polygons.
- **Red Grid Point Filtering**: Red grid points are now rendered in the PDF only if they are within Euclidean distance $\le 5$ of a differently colored grid point (yellow origin or blue good set point) or within distance $\le 1$ of the grid bounding box boundary.
- **Cleaned-up Polygon Storage & Linear $O(V)$ Algorithm**: Implemented a single-pass $O(V)$ stack-based cleanup algorithm (`polygon_rm_redundant_vertices` in `geom.rs`) that pops redundant collinear vertices as points are inserted and resolves wraparound boundary collinearity. Retained the legacy multi-pass algorithm (`polygon_rm_redundant_vertices_old`), added unit tests, and added runtime assertions in `main.rs` asserting exact equivalence across all test cases.
- **Output Deduplication**: Removed redundant unformatted configuration count output line from `dp.rs`, keeping the single comma-formatted `# Configs computed` entry in the summary metrics block.
