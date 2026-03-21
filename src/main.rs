mod dp;
mod draw;
mod geom;
mod point;
mod polygon;

use dp::{
    max_edge_length, minimize_perimeter_dp, NgDtoStrategy, NgPerimDtoStrategy, NgThenIdxStrategy,
    NgThenPerimStrategy, PerimNgDtogStrategy, PerimThenIdxStrategy, PerimThenNgStrategy,
};
use draw::{compute_perimeter, draw_polygon_with_grid};
use geom::{
    build_visibility_graph, ch_disk_origin, compute_good_set, compute_max_turn_angle,
    len_longest_edge, vtrans,
};
use point::*;
use polygon::*;

use std::env;
use std::fmt::Write as StrWrite;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use crate::geom::{boundary_grid_points, polygon_area};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k] [--queue=strategy]");
        println!();
        println!("Arguments:");
        println!("  [k]               The number of grid points the polygon should enclose.");
        println!(
            "  --no-cache        Disable the caching for the `is_all_left_turns` function call."
        );
        println!("  --queue=strategy  Select the priority queue ordering strategy.");
        println!();
        println!("Available queue strategies:");
        println!("  ng_perim_dto  (default) Sort by n_g (ascending), then (perimeter + DTO).");
        println!("  ng_perim                Sort by n_g (ascending), then perimeter (ascending).");
        println!("  perim_ng                Sort by perimeter (ascending), then n_g (ascending).");
        println!("  perim_ng_dtog           Sort by perimeter (ascending), then (n_g + DTO_G).");
        println!("  ng_dto                  Sort by n_g (ascending), then DTO (ascending).");
        println!("  ng_idx                  Sort by n_g (ascending), then queue insertion index.");
        println!(
            "  perim_idx               Sort by perimeter (ascending), then queue insertion index."
        );
        std::process::exit(1);
    }

    let dir_pdfs = "output/pdfs";
    let dir_summary = "output/summary";
    let dir_polys = "output/polys";

    std::fs::create_dir_all("output").unwrap_or_default();
    std::fs::create_dir_all(dir_pdfs).unwrap_or_default();
    std::fs::create_dir_all(dir_polys).unwrap_or_default();
    std::fs::create_dir_all(dir_summary).unwrap_or_default();

    let mut queue_strategy = "ng_perim_dto".to_string();
    for arg in args.iter().skip(1) {
        if arg.starts_with("--queue=") {
            queue_strategy = arg.split('=').nth(1).unwrap_or("ng_perim_dto").to_string();
        }
    }

    let k = match args[1].replace('_', "").parse::<usize>() {
        Ok(val) => val,
        Err(_) => {
            eprintln!("Error: Could not parse '{}' as a number.", args[1]);
            std::process::exit(1);
        }
    };

    println!("Computing convex-hull of disk around origin with k points... Kind of...");
    let ch_z = ch_disk_origin(k, false);
    let ch_z_exp = ch_disk_origin(k, true);

    let mut delta: i32 = ((k as f64).powf(0.25) / 4.0).ceil() as i32;
    if k > 400 {
        delta = 0;
    }
    let ch_m = vtrans(&ch_z, Point2D::new(-delta as CoordType, 0));
    let ch_m_exp = vtrans(&ch_z_exp, Point2D::new(-delta as CoordType, 0));

    //let l_f = ((k as f64).powf(0.25) / 3.0).round() as i64;
    //let l = if l_f > 3 { l_f } else { 3 };
    let l_f = 2 + ((k as f64).powf(0.25) / 4.0).round() as i64;
    let l = if l_f > 3 { l_f } else { 3 };

    println!("Computing good and bad sets with width = {}...", l);
    let (mut good, bad_ch) = compute_good_set(&ch_m_exp, l as f64);

    println!("Computing distance of good points to the origin...");
    good.fill_dist_to_origin(&bad_ch);
    println!("   ...done");

    let max_edge_l = max_edge_length(k);
    let dirs = crate::geom::generate_primitive_vectors(max_edge_l);

    println!("Building visibility graph...");
    let vg = build_visibility_graph(&good, &bad_ch, &dirs, k);

    let mut log = String::new();

    write!(log, "# k: {}\n", k)?;
    writeln!(log, "# Good points: {}", good.length())?;
    println!("# Good points: {}", good.length());

    let Ok((sol, ub_circle, conf_count)) = (if queue_strategy == "perim_ng" {
        minimize_perimeter_dp::<PerimThenNgStrategy>(k, &good, &vg)
    } else if queue_strategy == "ng_idx" {
        minimize_perimeter_dp::<NgThenIdxStrategy>(k, &good, &vg)
    } else if queue_strategy == "perim_idx" {
        minimize_perimeter_dp::<PerimThenIdxStrategy>(k, &good, &vg)
    } else if queue_strategy == "perim_ng_dtog" {
        minimize_perimeter_dp::<PerimNgDtogStrategy>(k, &good, &vg)
    } else if queue_strategy == "ng_dto" {
        minimize_perimeter_dp::<NgDtoStrategy>(k, &good, &vg)
    } else if queue_strategy == "ng_perim" {
        minimize_perimeter_dp::<NgThenPerimStrategy>(k, &good, &vg)
    } else {
        minimize_perimeter_dp::<NgPerimDtoStrategy>(k, &good, &vg)
    }) else {
        eprintln!("Error: The DP solver failed to initialize or write to disk.");
        std::process::exit(1);
    };

    println!("k: {}", k);

    let perimeter = compute_perimeter(&sol);
    let area = polygon_area(&sol);
    let b_n = boundary_grid_points(&sol);

    let v_n = sol.len();

    // Pick's theorem:  A= I + B/2 -1
    // As k = I + B, we have I = k - B, so A = k - B + B/2 - 1 = k - B/2 - 1
    // A = k - B/2 - 1
    if (area - (k as f64) + (b_n as f64) / 2.0 + 1.0).abs() > 1e-6 {
        panic!("Area not computed correctly");
    }

    draw_polygon_with_grid(dir_pdfs, &sol, &ch_m, &ch_m_exp, k, ub_circle, &good);
    let ch_m_perimeter = compute_perimeter(&ch_m);
    writeln!(log, "# Perimeter            : {}", perimeter)?;
    writeln!(log, "# circle perimeter     : {}", ch_m_perimeter)?;
    writeln!(log, "# Naive perimeter      : {}", ub_circle)?;
    writeln!(log, "# Area                 : {}", area)?;
    writeln!(log, "# vertices             : {}", v_n)?;
    writeln!(log, "# boundary grid points : {}", b_n)?;
    writeln!(log, "# Configs computed     : {}", conf_count)?;

    println!("# Perimeter            : {}", perimeter);
    println!("# circle perimeter     : {}", ch_m_perimeter);
    println!("# Naive perimeter      : {}", ub_circle);
    println!("# Area                 : {}", area);
    println!("# vertices             : {}", v_n);
    println!("# boundary grid points : {}", b_n);
    println!("# Configs computed     : {}", conf_count);

    let c_max_angle = compute_max_turn_angle(&sol);

    assert!(
        perimeter <= ch_m_perimeter * 1.00001,
        "perimeter not computed correctly A"
    );
    assert!(
        perimeter <= ub_circle * 1.00001,
        "perimeter not computed correctly B"
    );
    writeln!(log, "# Max angle of sol : {}", c_max_angle)?;
    writeln!(log, "# Longest edge len : {}", len_longest_edge(&sol))?;

    println!("Max angle of sol : {}", c_max_angle);
    println!("Longest edge len : {}", len_longest_edge(&sol));
    let duration = start.elapsed();

    writeln!(log, "# Running time in seconds: {}", duration.as_secs_f64())?;
    println!("Running time in seconds: {}", duration.as_secs_f64());

    let filename_poly: String = format!("{}/{:06}_poly.txt", dir_polys, k);
    save_polygon(&filename_poly, &sol, Some(&log))?;

    let filename_log: String = format!("{}/{:06}_s.txt", dir_summary, k);
    let mut log_fl = File::create(filename_log)?;
    writeln!(log_fl, "{}", &log)?;

    Ok(())
}
