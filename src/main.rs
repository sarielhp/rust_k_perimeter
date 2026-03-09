mod dp;
mod draw;
mod geom;

use dp::minimize_perimeter_dp;
use draw::{compute_perimeter, draw_polygon_with_grid};
use geom::{
    ch_disk_origin, compute_good_bad_sets, compute_max_turn_angle, len_longest_edge, vtrans,
    CoordType, Point2D,
};
use std::env;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k]");
        std::process::exit(1);
    }

    let k = match args[1].parse::<usize>() {
        Ok(val) => val,
        Err(_) => {
            println!("Error: '{}' is not a valid integer.", args[1]);
            std::process::exit(1);
        }
    };

    let ch_z = ch_disk_origin(k, false);
    let ch_z_exp = ch_disk_origin(k, true);

    let mut delta: i32 = ((k as f64).powf(0.25) / 4.0).ceil() as i32;
    if k > 400 {
        delta = 0;
    }
    let ch_m = vtrans(&ch_z, Point2D::new(-delta as CoordType, 0));
    let ch_m_exp = vtrans(&ch_z_exp, Point2D::new(-delta as CoordType, 0));

    let l_f = ((k as f64).powf(0.25) / 3.0).round() as i64;
    let l = if l_f > 3 { l_f } else { 3 };

    let (good, bad, bad_ch) = compute_good_bad_sets(&ch_m_exp, l as f64);

    println!("# Good points: {}", good.length());
    println!("# bad  points: {}", bad.length());

    let (sol, ub_circle) = minimize_perimeter_dp(k as u32, &good, &bad, &bad_ch);

    println!("k: {}", k);
    let perimeter = compute_perimeter(&sol);
    draw_polygon_with_grid(&sol, &ch_m, &ch_m_exp, k, ub_circle, &bad, &good);
    let ch_m_perimeter = compute_perimeter(&ch_m);
    println!("Perimeter        : {}", perimeter);
    println!("circle perimeter : {}", ch_m_perimeter);
    println!("Naive perimeter  : {}", ub_circle);
    let c_max_angle = compute_max_turn_angle(&sol);

    assert!(
        perimeter <= ch_m_perimeter * 1.00001,
        "perimeter not computed correctly A"
    );
    assert!(
        perimeter <= ub_circle * 1.00001,
        "perimeter not computed correctly B"
    );
    println!("Max angle of sol : {}", c_max_angle);
    println!("Longest edge len : {} ", len_longest_edge(&sol));
    let duration = start.elapsed();
    println!("Running time in seconds: {}", duration.as_secs())
}
