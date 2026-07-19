//! Entry point for the k-perimeter solver.
//!
//! This program finds a polygon with minimum Euclidean perimeter that encloses exactly k grid points.
//! The algorithm follows several stages:
//! 1.  Estimate a search area based on a disk-like boundary.
//! 2.  Partition grid points into a "good set" (search area) and a "bad set" (forbidden interior).
//! 3.  Construct a Visibility Graph of valid segments between good points.
//! 4.  Perform a Topological Sort to ensure the graph is a DAG and establish a processing order.
//! 5.  Run a Dynamic Programming (DP) solver using a priority queue and memory-mapped storage.
//! 6.  Reconstruct and visualize the optimal polygon.

mod dp;
mod draw;
mod geom;
mod kd_tree;
mod logger;
mod point;
mod polygon;
mod v_graph;

use dp::{max_edge_length, minimize_perimeter_dp};
use draw::{compute_perimeter, draw_polygon_with_grid};
use geom::{ch_disk_origin, compute_good_set, vtrans};
use point::*;
use num_format::{Locale, ToFormattedString};
use polygon::*;
use v_graph::build_visibility_graph;

use std::fs::File;
use std::io::Write;
use std::time::Instant;
use std::{env, fmt::Write as StrWrite};

use crate::geom::{
    boundary_grid_points, compute_max_turn_angle, len_longest_edge, len_longest_primitive_edge,
    polygon_area, polygon_boundary_distance, polygon_rm_redundant_vertices,
    polygon_rm_redundant_vertices_old,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    assert!(
        std::mem::size_of::<usize>() >= 8,
        "This program requires a 64-bit architecture (usize must be at least 64 bits)."
    );
    logger::clear_log();

    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k] [--vg-only] [--retain factor]");
        println!();
        println!("Arguments:");
        println!("  [k]               The number of grid points the polygon should enclose.");
        println!("  --vg-only         Stop after visibility graph construction.");
        println!("  --retain factor   The factor by which the dynamic threshold is increased after cleanup (default: 1.5).");
        println!();
        std::process::exit(1);
    }

    let dir_pdfs = "output/pdfs";
    let dir_summary = "output/summary";
    let dir_polys = "output/polys";

    // Ensure all output directories exist.
    std::fs::create_dir_all("output").unwrap_or_default();
    std::fs::create_dir_all(dir_pdfs).unwrap_or_default();
    std::fs::create_dir_all(dir_polys).unwrap_or_default();
    std::fs::create_dir_all(dir_summary).unwrap_or_default();

    let mut vg_only = false;
    let mut retain_factor = 1.1;

    let mut i = 1;
    while i < args.len() {
        if args[i] == "--vg-only" {
            vg_only = true;
        } else if args[i] == "--retain" && i + 1 < args.len() {
            retain_factor = args[i + 1].parse::<f64>().unwrap_or(1.5);
            i += 1;
        } else if args[i].starts_with("--retain=") {
            retain_factor = args[i]["--retain=".len()..].parse::<f64>().unwrap_or(1.5);
        }
        i += 1;
    }

    let k = match args[1].replace('_', "").parse::<usize>() {
        Ok(val) => val,
        Err(_) => {
            eprintln!("Error: Could not parse '{}' as a number.", args[1]);
            std::process::exit(1);
        }
    };

    // Stage 1: Initial Bound Estimation.
    // Generates a convex hull of ~k points around origin to establish a search baseline.
    log_println!("Computing convex-hull of disk around origin with k points... Kind of...");
    let ch_z = ch_disk_origin(k, false);
    let ch_z_exp = ch_disk_origin(k, true);

    let mut delta: i32 = ((k as f64).powf(0.25) / 4.0).ceil() as i32;
    if k > 400 {
        delta = 0;
    }
    let ch_m = vtrans(&ch_z, Point2D::new(-delta as CoordType, 0));
    let ch_m_exp = vtrans(&ch_z_exp, Point2D::new(-delta as CoordType, 0));

    // Stage 2: Grid Partitioning.
    // Identifies "good" points (near the boundary) and "bad" points (deep inside).
    let l_f = 1 + (1.1 * (k as f64).powf(0.25) / 4.0).round() as i64;
    let l = if l_f > 3 { l_f } else { 3 };

    log_println!("Computing good and bad sets with width = {}...", l);
    let (mut good, bad_ch, so_so, bad_out) = compute_good_set(&ch_m_exp, l as f64);

    log_println!("Computing distance of good points to the origin...");
    good.fill_dist_to_origin(&bad_ch);
    log_println!("   ...done");

    // Stage 3: Visibility Graph Construction.
    // Build edges between good points while respecting turn angles and bad set obstacles.
    let max_edge_l = max_edge_length(k);
    let dirs = geom::generate_primitive_vectors(max_edge_l);
    let max_turn_angle = 3.0 * std::f64::consts::PI / (k as f64).powf(1.0 / 3.0);
    log_println!("max_turn_angle: {}", max_turn_angle);

    log_println!("Building visibility graph...");
    let mut vg = build_visibility_graph(&good, &bad_ch, &so_so, &dirs, k, max_turn_angle);

    if vg.num_edges() == 0 {
        panic!("Visibility graph is empty. No valid edges found between grid points for k={}. Consider increasing the search area or checking convexity constraints.", k);
    }

    if vg_only {
        log_println!("Stopping after visibility graph construction as requested.");
        return Ok(());
    }

    // Stage 4: Topological Sorting.
    // Ensures the graph is acyclic and provides an optimal processing order for the DP.
    log_println!("Performing topological sort...");
    let origin = Point2D::new(0, 0);
    if good.contains(&origin) {
        let origin_id = good.get_point_id(origin);
        for edges in vg.adjacency_list.iter_mut() {
            edges.retain(|e| e.target_id != origin_id);
        }
    }

    // Now the sorting...
    if let Some(topo_order) = v_graph::topological_sort(&vg) {
        for (idx, &id) in topo_order.iter().enumerate() {
            good.set_topo_idx(id, idx as u32);
        }
    } else {
        eprintln!("Error: Visibility graph is not a DAG.");
        std::process::exit(1);
    }

    // Stage 5: DP Solver.
    // Priority-queue based exploration of the configuration space.
    let start_dp = Instant::now();
    let dp_res = minimize_perimeter_dp(k, &good, &vg, retain_factor);
    let dp_duration = start_dp.elapsed();

    let Ok((sol, ub_circle, conf_count)) = dp_res else {
        eprintln!("Error: The DP solver failed.");
        std::process::exit(1);
    };

    // Stage 6: Validation and Output.
    // Verify results with Pick's Theorem and save outputs.
    let rt = start.elapsed().as_secs_f64();

    let t_start_linear = Instant::now();
    let sol_c = polygon_rm_redundant_vertices(&sol);
    let cleanup_linear_dur = t_start_linear.elapsed();

    let t_start_old = Instant::now();
    let sol_c_old = polygon_rm_redundant_vertices_old(&sol);
    let cleanup_old_dur = t_start_old.elapsed();

    assert_eq!(
        sol_c, sol_c_old,
        "Mismatch between linear and old polygon_rm_redundant_vertices!"
    );
    let perimeter = compute_perimeter(&sol_c);
    let area = polygon_area(&sol_c);
    let b_n = boundary_grid_points(&sol_c);
    let v_n = sol_c.len();

    let mut log = String::new();

    // Save before checking if the solution is valid, to have a record of the output even if Pick's theorem fails.
    save_polygon(
        &format!("{}/{:07}_poly.txt", dir_polys, k),
        &sol_c,
        Some(&log),
    )?;

    // Pick's theorem verification: A = I + B/2 - 1 => A = k - B/2 - 1
    let k_by_pick = area + (b_n as f64) / 2.0 + 1.0;
    if (k as f64 - k_by_pick).abs() > 1e-3 {
        log_println!("Error: Pick's Theorem verification failed.");
        log_println!("solution # of edges  : {}", sol_c.len());
        log_println!("Computed area        : {}", area);
        log_println!("Boundary points (B)  : {}", b_n);
        log_println!("k by Pick's          : {}", k_by_pick);
        log_println!("k                    : {}", k);

        panic!("Area/Point-count mismatch via Pick's Theorem");
    } else {
        log_println!("Pick's Theorem verification passed.");
    }

    let ch_m_perimeter = compute_perimeter(&ch_m);

    // compute the minimum distance from bad_out points to the solution polygon, to verify that the solution
    // is not too close to bad points.
    let mut min_bad_d: f64 = k as f64;
    for p in &bad_out {
        let d = polygon_boundary_distance(&sol_c, *p);
        if d < min_bad_d {
            min_bad_d = d;
        }
    }

    let mut log_and_print =
        |label: &str, value: &dyn std::fmt::Display| -> Result<(), Box<dyn std::error::Error>> {
            writeln!(log, "# {:26} : {}", label, value)?;
            log_println!("# {:26} : {}", label, value);
            Ok(())
        };

    log_and_print("k", &k.to_formatted_string(&Locale::en))?;
    log_and_print("Perimeter", &perimeter)?;
    log_and_print("circle perimeter", &ch_m_perimeter)?;
    log_and_print("Naive perimeter", &ub_circle)?;
    log_and_print("Area", &area)?;
    log_and_print("vertices", &v_n.to_formatted_string(&Locale::en))?;
    log_and_print("boundary grid points", &b_n.to_formatted_string(&Locale::en))?;
    log_and_print("Configs computed", &conf_count.to_formatted_string(&Locale::en))?;

    let c_max_angle = compute_max_turn_angle(&sol_c);
    log_and_print("Max angle of sol", &c_max_angle)?;
    log_and_print("UB max turn_angle", &max_turn_angle)?;

    let len_p_longest = len_longest_primitive_edge(&sol_c);
    let len_longest = len_longest_edge(&sol_c);
    log_and_print("Longest edge len", &len_longest)?;

    log_and_print("Longest primitive edge len", &len_p_longest)?;
    log_and_print("UB primitive edge len", &max_edge_l.to_formatted_string(&Locale::en))?;

    log_and_print("Min dist bad pnt to sol", &min_bad_d)?;

    log_and_print("DP duration", &dp_duration.as_secs_f64())?;
    log_and_print("Cleanup linear duration", &cleanup_linear_dur.as_secs_f64())?;
    log_and_print("Cleanup old duration", &cleanup_old_dur.as_secs_f64())?;
    log_and_print("Running time in seconds", &rt)?;

    // Draw PDF with page 1 figure and page 2+ color explanations and standard output log.
    draw_polygon_with_grid(
        dir_pdfs,
        &sol_c,
        &ch_m,
        &ch_m_exp,
        k,
        ub_circle,
        &good,
        &logger::get_log(),
    );

    // Save polygon with complete log comments.
    save_polygon(
        &format!("{}/{:07}_poly.txt", dir_polys, k),
        &sol_c,
        Some(&log),
    )?;

    // Save summary log for this run.
    let mut log_fl = File::create(format!("{}/{:07}_s.txt", dir_summary, k))?;
    writeln!(log_fl, "{}", &log)?;

    Ok(())
}
