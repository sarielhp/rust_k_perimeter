mod dp;
mod draw;
mod point;
mod polygon;
mod geom;


use dp::{minimize_perimeter_dp, NgThenPerimStrategy, PerimThenNgStrategy, NgThenIdxStrategy, PerimThenIdxStrategy, NgPerimDtoStrategy, PerimNgDtogStrategy, NgDtoStrategy};
use draw::{compute_perimeter, draw_polygon_with_grid};
use point::*;
use geom::{
    ch_disk_origin, compute_good_bad_sets, compute_max_turn_angle, len_longest_edge, vtrans,
};
use polygon::*;

use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use std::fmt::Write as StrWrite;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k] [--no-cache] [--queue=strategy]");
        println!();
        println!("Arguments:");
        println!("  [k]               The number of grid points the polygon should enclose.");
        println!("  --no-cache        Disable the caching for the `is_all_left_turns` function call.");
        println!("  --queue=strategy  Select the priority queue ordering strategy.");
        println!();
        println!("Available queue strategies:");
        println!("  ng_perim_dto  (default) Sort by n_g (ascending), then (perimeter + DTO).");
        println!("  ng_perim                Sort by n_g (ascending), then perimeter (ascending).");
        println!("  perim_ng                Sort by perimeter (ascending), then n_g (ascending).");
        println!("  perim_ng_dtog           Sort by perimeter (ascending), then (n_g + DTO_G).");
        println!("  ng_dto                  Sort by n_g (ascending), then DTO (ascending).");
        println!("  ng_idx                  Sort by n_g (ascending), then queue insertion index.");
        println!("  perim_idx               Sort by perimeter (ascending), then queue insertion index.");
        std::process::exit(1);
    }

    let  dir_pdfs = "output/pdfs";
    let  dir_summary = "output/summary";
    let  dir_polys = "output/polys";
    
    std::fs::create_dir_all("output").unwrap_or_default();
    std::fs::create_dir_all( dir_pdfs ).unwrap_or_default();
    std::fs::create_dir_all( dir_polys ).unwrap_or_default();
    std::fs::create_dir_all( dir_summary ).unwrap_or_default();
    

    let use_cache = !args.contains(&String::from("--no-cache"));

    let mut queue_strategy = "ng_perim_dto".to_string();
    for arg in args.iter().skip(1) {
        if arg.starts_with("--queue=") {
            queue_strategy = arg.split('=').nth(1).unwrap_or("ng_perim_dto").to_string();
        }
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

    //let l_f = ((k as f64).powf(0.25) / 3.0).round() as i64;
    //let l = if l_f > 3 { l_f } else { 3 };
    let l_f = 2 + ((k as f64).powf(0.25) / 4.0).round() as i64;
    let l = if l_f > 3 { l_f } else { 3 };

    let (mut good, bad, bad_ch) = compute_good_bad_sets(&ch_m_exp, l as f64);

    println!("Computing distance of good points to the origin...");
    good.fill_dist_to_origin(&bad_ch);
    println!("   ...done");

    let mut log = String::new();

    write!(log, "# k: {}\n", k)?;
    writeln!(log, "# Good points: {}", good.length())?;
    writeln!(log, "# bad  points: {}", bad.length())?;
    println!("# Good points: {}", good.length());
    println!("# bad  points: {}", bad.length());
    
    let (sol, ub_circle) = if queue_strategy == "perim_ng" {
        minimize_perimeter_dp::<PerimThenNgStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else if queue_strategy == "ng_idx" {
        minimize_perimeter_dp::<NgThenIdxStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else if queue_strategy == "perim_idx" {
        minimize_perimeter_dp::<PerimThenIdxStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else if queue_strategy == "perim_ng_dtog" {
        minimize_perimeter_dp::<PerimNgDtogStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else if queue_strategy == "ng_dto" {
        minimize_perimeter_dp::<NgDtoStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else if queue_strategy == "ng_perim" {
        minimize_perimeter_dp::<NgThenPerimStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    } else {
        minimize_perimeter_dp::<NgPerimDtoStrategy>(k as u32, &good, &bad, &bad_ch, use_cache)
    };

    println!("k: {}", k);

    let perimeter = compute_perimeter(&sol);
    draw_polygon_with_grid( dir_pdfs, &sol, &ch_m, &ch_m_exp, k, ub_circle, &bad, &good);
    let ch_m_perimeter = compute_perimeter(&ch_m);
    writeln!(log, "# Perimeter        : {}", perimeter)?;
    writeln!(log, "# circle perimeter : {}", ch_m_perimeter)?;
    writeln!(log, "# Naive perimeter  : {}", ub_circle)?;
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
    writeln!(log, "# Max angle of sol : {}", c_max_angle)?;
    writeln!(log, "# Longest edge len : {}", len_longest_edge(&sol))?;

    println!("Max angle of sol : {}", c_max_angle);
    println!("Longest edge len : {}", len_longest_edge(&sol));
    let duration = start.elapsed();

    writeln!(log, "# Running time in seconds: {}", duration.as_secs())?;
    println!("Running time in seconds: {}", duration.as_secs());

    let filename_poly: String = format!("{}/{:05}_poly.txt", dir_polys, k);
    save_polygon( &filename_poly, &sol, Some( &log ) )?;

    let filename_log: String = format!("{}/{:05}_summmary.txt", dir_summary, k);
    let mut log_fl = File::create(filename_log)?;
    writeln!( log_fl, "{}", &log )?;

    Ok(())
}
