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
mod point;
mod polygon;
mod v_graph;

use dp::{max_edge_length, minimize_perimeter_dp};
use draw::{compute_perimeter, draw_polygon_with_grid};
use geom::{ch_disk_origin, compute_good_set, vtrans};
use v_graph::build_visibility_graph;
use point::*;
use polygon::*;

use std::env;
use std::fmt::Write as StrWrite;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use crate::geom::{
    boundary_grid_points, polygon_area,
    polygon_rm_redundant_vertices,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k] [--topo] [--vg-only] [--retain factor]");
        println!();
        println!("Arguments:");
        println!("  [k]               The number of grid points the polygon should enclose.");
        println!("  --topo            Build visibility graph, perform topological sort, output vertices to output/topo.txt and exit.");
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

    let mut topo_mode = false;
    let mut vg_only = false;
    let mut retain_factor = 1.5;

    let mut i = 1;
    while i < args.len() {
        if args[i] == "--topo" {
            topo_mode = true;
        } else if args[i] == "--vg-only" {
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
    println!("Computing convex-hull of disk around origin with k points... Kind of...");
    let ch_z = ch_disk_origin(k, false);
    let ch_z_exp = ch_disk_origin(k, true);

    let mut delta: i32 = ((k as f64).powf(0.25) / 4.0).ceil() as i32;
    if k > 400 { delta = 0; }
    let ch_m = vtrans(&ch_z, Point2D::new(-delta as CoordType, 0));
    let ch_m_exp = vtrans(&ch_z_exp, Point2D::new(-delta as CoordType, 0));

    // Stage 2: Grid Partitioning.
    // Identifies "good" points (near the boundary) and "bad" points (deep inside).
    let l_f = 2 + ((k as f64).powf(0.25) / 4.0).round() as i64;
    let l = if l_f > 3 { l_f } else { 3 };

    println!("Computing good and bad sets with width = {}...", l);
    let (mut good, bad_ch, so_so) = compute_good_set(&ch_m_exp, l as f64);

    println!("Computing distance of good points to the origin...");
    good.fill_dist_to_origin(&bad_ch);
    println!("   ...done");

    // Stage 3: Visibility Graph Construction.
    // Build edges between good points while respecting turn angles and bad set obstacles.
    let max_edge_l = max_edge_length(k);
    let dirs = geom::generate_primitive_vectors(max_edge_l);
    let max_turn_angle = 3.0 * std::f64::consts::PI / (k as f64).powf(1.0 / 3.0);
    println!("max_turn_angle: {}", max_turn_angle);

    println!("Building visibility graph...");
    let mut vg = build_visibility_graph(&good, &bad_ch, &so_so, &dirs, k, max_turn_angle);

    if vg_only {
        println!("Stopping after visibility graph construction as requested.");
        return Ok(());
    }

    // Stage 4: Topological Sorting.
    // Ensures the graph is acyclic and provides an optimal processing order for the DP.
    println!("Performing topological sort...");
    let origin = Point2D::new(0, 0);
    if good.contains(&origin) {
        let origin_id = good.get_point_id(origin);
        for edges in vg.adjacency_list.iter_mut() {
            edges.retain(|e| e.target_id != origin_id);
        }
    }

    if let Some(topo_order) = v_graph::topological_sort(&vg) {
        for (idx, &id) in topo_order.iter().enumerate() {
            good.set_topo_idx(id, idx as u32);
        }

        if topo_mode {
            let mut fl = File::create("output/topo.txt")?;
            for id in &topo_order {
                let p = good.get_point_by_id(*id);
                writeln!(fl, "{} {}", p.x, p.y)?;
            }
            println!("Topological order written to output/topo.txt.");
            return Ok(());
        }
    } else {
        eprintln!("Error: Visibility graph is not a DAG.");
        std::process::exit(1);
    }

    // Stage 5: DP Solver.
    // Priority-queue based exploration of the configuration space.
    println!("Starting DP solver with retain factor {}...", retain_factor);
    let start_dp = Instant::now();
    let dp_res = minimize_perimeter_dp(k, &good, &vg, retain_factor);
    let dp_duration = start_dp.elapsed();

    let Ok((sol, ub_circle, conf_count)) = dp_res else {
        eprintln!("Error: The DP solver failed.");
        std::process::exit(1);
    };

    // Stage 6: Validation and Output.
    // Verify results with Pick's Theorem and save outputs.
    let sol_c = polygon_rm_redundant_vertices(&sol);
    let perimeter = compute_perimeter(&sol);
    let area = polygon_area(&sol);
    let b_n = boundary_grid_points(&sol_c);

    // Pick's theorem verification: A = I + B/2 - 1
    if (area - (k as f64) + (b_n as f64) / 2.0 + 1.0).abs() > 1e-6 {
        panic!("Area/Point-count mismatch via Pick's Theorem");
    }

    draw_polygon_with_grid(dir_pdfs, &sol, &ch_m, &ch_m_exp, k, ub_circle, &good);
    
    let mut log = String::new();
    let mut log_and_print = |label: &str, value: &dyn std::fmt::Display| -> Result<(), Box<dyn std::error::Error>> {
        writeln!(log, "# {:26} : {}", label, value)?;
        println!("# {:26} : {}", label, value);
        Ok(())
    };

    log_and_print("DP duration", &dp_duration.as_secs_f64())?;
    log_and_print("Perimeter", &perimeter)?;
    log_and_print("Area", &area)?;
    log_and_print("Configs computed", &(conf_count as f64))?;
    log_and_print("Running time in seconds", &start.elapsed().as_secs_f64())?;

    save_polygon(&format!("{}/{:06}_poly.txt", dir_polys, k), &sol, Some(&log))?;
    let mut log_fl = File::create(format!("{}/{:06}_s.txt", dir_summary, k))?;
    writeln!(log_fl, "{}", &log)?;

    Ok(())
}
