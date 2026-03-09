use crate::geom::{bound, CoordType, GridSet, Point2D};
use cairo::{Context, FontSlant, FontWeight, PdfSurface};

fn cr_string_pdf(cr: &Context, w: f64, h: f64, text_content: &str) {
    cr.set_source_rgb(1.0, 1.0, 1.0);
    cr.rectangle(0.0, 0.0, w, h);
    cr.fill().unwrap();

    cr.select_font_face("Sans", FontSlant::Normal, FontWeight::Bold);
    cr.set_font_size(8.0);
    let fh = 16.0;

    let lines: Vec<&str> = text_content.lines().collect();
    let mut current_y = fh * 1.5;

    cr.set_source_rgb(0.0, 0.0, 0.0);
    for line in lines {
        cr.move_to(10.0, current_y);
        cr.show_text(line).unwrap();
        current_y += fh;
    }
}

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
pub fn compute_perimeter(poly: &[Point2D]) -> f64 {
    crate::geom::euclidean_length(poly)
}

pub fn draw_polygon_with_grid(
    poly: &[Point2D],
    poly_circ: &[Point2D],
    poly_circ_exp: &[Point2D],
    k: usize,
    ub_circle: f64,
    bad: &GridSet,
    good: &GridSet,
) {
    std::fs::create_dir_all("output").unwrap_or_default();
    let filename: String = format!("output/{}.pdf", k);

    let (min_x, max_x, min_y, max_y) = bound(&[poly, poly_circ, poly_circ_exp], 1);

    let margin = 50.0;
    let scale = (100.0 / (k as f64)).max(10.0).round() as i64;

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

    let trans_poly: Vec<(f64, f64)> = poly.iter().map(|&p| translate(p)).collect();
    let trans_poly_circ: Vec<(f64, f64)> = poly_circ.iter().map(|&p| translate(p)).collect();
    let trans_poly_circ_exp: Vec<(f64, f64)> =
        poly_circ_exp.iter().map(|&p| translate(p)).collect();

    cr.set_line_width(2.0);
    cr.set_source_rgb(0.2, 0.5, 0.8);
    draw_polygon(&cr, &trans_poly, (0.2, 0.5, 0.8, 1.0));

    cr.set_line_width(1.0);
    cr.set_source_rgb(0.9, 0.5, 0.8);
    draw_polygon(&cr, &trans_poly_circ, (0.9, 0.5, 0.2, 0.01));
    cr.set_source_rgb(0.1, 1.0, 0.1);
    draw_polygon(&cr, &trans_poly_circ_exp, (0.1, 0.9, 0.1, 0.1));

    for x in min_x..=max_x {
        for y in min_y..=max_y {
            let (tx, ty) = translate(Point2D::new(x as CoordType, y as CoordType));
            let p = Point2D::new(x as CoordType, y as CoordType);

            if x == 0 && y == 0 {
                cr.set_source_rgb(1.0, 1.0, 0.0);
            } else if bad.contains(&p) {
                cr.set_source_rgb(0.8, 0.0, 0.0);
            } else if good.contains(&p) {
                cr.set_source_rgb(0.0, 0.0, 1.0);
            } else {
                cr.set_source_rgb(0.2, 0.2, 0.2);
            }

            cr.arc(tx, ty, 1.5, 0.0, 2.0 * std::f64::consts::PI);
            cr.fill().unwrap();
        }
    }

    cr.show_page().unwrap();

    let str_main = format!("k = {}\nCircle ub on perimeter: {}\n", k, ub_circle);
    let str_b = format!(
        "solution parimeter: {}\n# vertices: {}\n",
        compute_perimeter(poly),
        poly.len()
    );

    let text = format!("{}{}", str_main, str_b);
    cr_string_pdf(&cr, width, height, &text);

    surface.finish();
    println!("Polygon saved to {}", filename);
}
