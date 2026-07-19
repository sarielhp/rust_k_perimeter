//! Visualization utilities using Cairo.
//!
//! Provides functions to render the resulting polygon and the grid of points
//! to a PDF file for inspection.

use crate::geom::{bound, GridSet};
use crate::point::*;
use cairo::{Context, FontSlant, FontWeight, PdfSurface};
use num_format::{Locale, ToFormattedString};

/// Helper to wrap a single line of text if it exceeds `max_width` points.
fn wrap_line(cr: &Context, line: &str, max_width: f64) -> Vec<String> {
    if line.is_empty() {
        return vec!["".to_string()];
    }
    let extents = match cr.text_extents(line) {
        Ok(ext) => ext,
        Err(_) => return vec![line.to_string()],
    };
    if extents.x_advance() <= max_width {
        return vec![line.to_string()];
    }

    let mut result = Vec::new();
    let mut current_subline = String::new();

    for word in line.split(' ') {
        let test_line = if current_subline.is_empty() {
            word.to_string()
        } else {
            format!("{} {}", current_subline, word)
        };

        let advance = cr.text_extents(&test_line).map(|e| e.x_advance()).unwrap_or(0.0);
        if advance <= max_width {
            current_subline = test_line;
        } else {
            if !current_subline.is_empty() {
                result.push(current_subline);
            }
            current_subline = word.to_string();
        }
    }

    if !current_subline.is_empty() {
        result.push(current_subline);
    }

    if result.is_empty() {
        vec![line.to_string()]
    } else {
        result
    }
}

/// Renders multi-page text information onto the PDF surface.
fn render_text_pages(cr: &Context, width: f64, height: f64, text_content: &str) {
    let margin_x = 35.0;
    let margin_y = 35.0;
    let font_size = 9.0;
    let line_height = 13.0;

    let page_width = width.max(600.0);
    let page_height = height.max(750.0);
    let max_text_width = page_width - 2.0 * margin_x;
    let max_y = page_height - margin_y;

    let prepare_page = |cr: &Context| {
        cr.set_source_rgb(1.0, 1.0, 1.0);
        cr.rectangle(0.0, 0.0, page_width, page_height);
        cr.fill().unwrap();
        cr.set_source_rgb(0.0, 0.0, 0.0);
    };

    prepare_page(cr);
    let mut current_y = margin_y + font_size;

    for line in text_content.lines() {
        let apply_style = |cr: &Context, _subline: &str| {
            if line.starts_with("===") || line.starts_with("---") {
                cr.select_font_face("Monospace", FontSlant::Normal, FontWeight::Bold);
                cr.set_font_size(font_size);
                cr.set_source_rgb(0.3, 0.3, 0.3);
            } else if line.starts_with("FIRST PAGE")
                || line.starts_with("STANDARD OUTPUT LOG")
                || line.ends_with(':')
            {
                cr.select_font_face("Sans", FontSlant::Normal, FontWeight::Bold);
                cr.set_font_size(font_size + 1.0);
                cr.set_source_rgb(0.0, 0.2, 0.6);
            } else if line.starts_with("  •") {
                cr.select_font_face("Sans", FontSlant::Normal, FontWeight::Bold);
                cr.set_font_size(font_size);
                cr.set_source_rgb(0.1, 0.1, 0.1);
            } else {
                cr.select_font_face("Monospace", FontSlant::Normal, FontWeight::Normal);
                cr.set_font_size(font_size);
                cr.set_source_rgb(0.0, 0.0, 0.0);
            }
        };

        let sublines = wrap_line(cr, line, max_text_width);
        for subline in sublines {
            if current_y > max_y {
                cr.show_page().unwrap();
                prepare_page(cr);
                current_y = margin_y + font_size;
            }

            apply_style(cr, &subline);

            cr.move_to(margin_x, current_y);
            cr.show_text(&subline).unwrap();
            current_y += line_height;
        }
    }
}

/// Helper to draw a filled and outlined polygon.
fn draw_polygon(cr: &Context, poly: &[(f64, f64)], rgba: (f64, f64, f64, f64)) {
    if poly.is_empty() {
        return;
    }
    let pt = poly[0];
    cr.move_to(pt.0, pt.1);
    for i in 1..poly.len() {
        cr.line_to(poly[i].0, poly[i].1);
    }
    cr.close_path();
    cr.stroke_preserve().unwrap();
    cr.set_source_rgba(rgba.0, rgba.1, rgba.2, rgba.3);
    cr.fill().unwrap();
}

/// Simple wrapper for perimeter calculation.
pub fn compute_perimeter(poly: &[Point2D]) -> f64 {
    crate::geom::euclidean_length(poly)
}

/// Renders the optimal polygon, estimated circle-like boundaries, grid points,
/// and multi-page execution log summary on subsequent pages.
pub fn draw_polygon_with_grid(
    dir_pdfs: &str,
    poly: &[Point2D],
    poly_circ: &[Point2D],
    poly_circ_exp: &[Point2D],
    k: usize,
    _ub_circle: f64,
    good: &GridSet,
    stdout_info: &str,
) {
    let filename: String = format!("{}/{:06}.pdf", dir_pdfs, k);
    let (min_x, max_x, _, max_y) = bound(&[poly, poly_circ, poly_circ_exp], 10);
    let min_y = 0;

    let margin = 50.0;
    // Dynamic scale based on k to keep the output size reasonable.
    let scale = (100.0 / (k as f64).sqrt()).max(2.0).round() as i64;

    let width = (scale * (max_x as i64 - min_x as i64)) as f64 + 2.0 * margin;
    let height = (scale * (max_y as i64 - min_y as i64)) as f64 + 2.0 * margin;

    let translate = |p: Point2D| -> (f64, f64) {
        (
            (scale * (p.x as i64 - min_x as i64)) as f64 + margin,
            (scale * (p.y as i64 - min_y as i64)) as f64 + margin,
        )
    };

    let surface = PdfSurface::new((width).round(), (height).round(), filename.clone()).unwrap();
    let cr = Context::new(&surface).unwrap();

    // Fill Page 1 background with solid white
    cr.set_source_rgb(1.0, 1.0, 1.0);
    cr.rectangle(0.0, 0.0, width, height);
    cr.fill().unwrap();

    let trans_poly: Vec<(f64, f64)> = poly.iter().map(|&p| translate(p)).collect();
    let trans_poly_circ: Vec<(f64, f64)> = poly_circ.iter().map(|&p| translate(p)).collect();
    let trans_poly_circ_exp: Vec<(f64, f64)> =
        poly_circ_exp.iter().map(|&p| translate(p)).collect();

    // Draw the optimal polygon FIRST (with black boundary width 1.0 and blue fill).
    cr.set_line_width(1.0);
    cr.set_source_rgb(0.0, 0.0, 0.0);
    draw_polygon(&cr, &trans_poly, (0.2, 0.5, 0.8, 1.0));

    // Draw the circle-like estimated bounds on top so they remain visible.
    cr.set_line_width(1.0);
    cr.set_source_rgb(0.9, 0.5, 0.8);
    draw_polygon(&cr, &trans_poly_circ, (0.9, 0.5, 0.2, 0.01));
    cr.set_source_rgb(0.1, 1.0, 0.1);
    draw_polygon(&cr, &trans_poly_circ_exp, (0.1, 0.9, 0.1, 0.1));

    let rad: f64 = width / (20.0 * (max_x as f64));
    // Draw grid points: yellow origin, blue good points, and filtered red points.
    for x in min_x..=max_x {
        for y in min_y..=max_y {
            let p = Point2D::new(x as CoordType, y as CoordType);
            let (tx, ty) = translate(p);

            if x == 0 && y == 0 {
                cr.set_source_rgb(1.0, 1.0, 0.0);
            }
            // Origin
            else if good.contains(&p) {
                cr.set_source_rgb(0.0, 0.0, 1.0);
            }
            // Search set
            else {
                // Red point: draw only if distance <= 5 to a non-red point or distance <= 1 to grid boundary
                let near_grid_boundary = (x - min_x) <= 1
                    || (max_x - x) <= 1
                    || (y - min_y) <= 1
                    || (max_y - y) <= 1;

                let mut near_non_red = false;
                if !near_grid_boundary {
                    'search: for dx in -5..=5 {
                        for dy in -5..=5 {
                            if dx * dx + dy * dy <= 25 {
                                let q = Point2D::new(p.x + dx as CoordType, p.y + dy as CoordType);
                                if (q.x == 0 && q.y == 0) || good.contains(&q) {
                                    near_non_red = true;
                                    break 'search;
                                }
                            }
                        }
                    }
                }

                if !(near_grid_boundary || near_non_red) {
                    continue;
                }

                cr.set_source_rgb(0.8, 0.0, 0.0);
            } // Filtered Interior/Exterior

            cr.arc(tx, ty, rad, 0.0, 2.0 * std::f64::consts::PI);
            cr.fill().unwrap();
        }
    }

    // End Page 1 (Figure)
    cr.show_page().unwrap();

    let text_page_w = width.max(600.0);
    let text_page_h = height.max(750.0);
    surface.set_size(text_page_w, text_page_h).unwrap();

    let hash_info = stdout_info
        .lines()
        .filter(|line| line.trim_start().starts_with('#'))
        .collect::<Vec<&str>>()
        .join("\n");

    // Construct full text info for Page 2 and subsequent pages
    let full_info = format!(
r#"================================================================================
FIRST PAGE COLOR & POLYGON EXPLANATION
================================================================================

Polygons:
  • Blue Polygon with Black Boundary (line width 1.0):
    Optimal Solution Polygon — Minimum Euclidean perimeter polygon enclosing
    exactly k={} grid points found by the DP solver.

  • Red/Orange Outlined Polygon (line width 1.0):
    Baseline Disk Polygon (ch_m) — Convex hull of ~k points around origin
    shifted by delta, providing the baseline geometric estimate.

  • Green Outlined Polygon (line width 1.0):
    Expanded Disk Polygon (ch_m_exp) — Convex hull of expanded disk estimate,
    defining outer search boundary for grid partitioning.

Grid Points (Dots):
  • Yellow Dot:
    Origin (0,0) — Reference origin point for Pick's Theorem & DP search.

  • Blue Dots:
    Good Grid Points — Points within search width l of estimated boundary;
    these form candidate vertices in the visibility graph.

  • Red Dots:
    Interior / Exterior Points — Points excluded from candidate vertices
    (interior bad-set points deep inside or exterior points outside search set).

================================================================================
SUMMARY METRICS
================================================================================
{}"#,
        k.to_formatted_string(&Locale::en), hash_info
    );

    render_text_pages(&cr, text_page_w, text_page_h, &full_info);

    surface.finish();
    println!("Polygon saved to {}", filename);
}
