mod dp;
mod draw;
mod geom;
mod kd_tree;
mod point;
mod polygon;
mod v_graph;

use dp::{
    max_edge_length, minimize_perimeter_dp, NgDtoStrategy, NgPerimDtoStrategy, NgThenIdxStrategy,
    NgThenPerimStrategy, PerimNgDtogStrategy, PerimThenIdxStrategy, PerimThenNgStrategy,
    TopoThenNgStrategy,
};
use draw::{compute_perimeter, draw_polygon_with_grid};
use geom::{ch_disk_origin, compute_good_set, compute_max_turn_angle, vtrans};
use v_graph::build_visibility_graph;
use point::*;
use polygon::*;

use std::env;
use std::fmt::Write as StrWrite;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use crate::geom::{
    boundary_grid_points, len_longest_edge, len_longest_primitive_edge, polygon_area,
    polygon_rm_redundant_vertices,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Usage: k_perimeter [k] [--queue=strategy] [--topo] [--vg-only]");
        println!();
        println!("Arguments:");
        println!("  [k]               The number of grid points the polygon should enclose.");
        println!(
            "  --no-cache        Disable the caching for the `is_all_left_turns` function call."
        );
        println!("  --queue=strategy  Select the priority queue ordering strategy.");
        println!("  --topo            Build visibility graph, perform topological sort, output vertices to output/topo.txt and exit.");
        println!("  --vg-only         Stop after visibility graph construction.");
        println!();
        println!("Available queue strategies:");
        println!("  topo_ng       (default) Sort by topological index (ascending), then n_g (ascending).");
        println!("  ng_idx        Sort by n_g (ascending), then queue insertion index.");
        println!("  ng_perim_dto            Sort by n_g (ascending), then (perimeter + DTO).");
        println!("  ng_perim                Sort by n_g (ascending), then perimeter (ascending).");
        println!("  perim_ng                Sort by perimeter (ascending), then n_g (ascending).");
        println!("  perim_ng_dtog           Sort by perimeter (ascending), then (n_g + DTO_G).");
        println!("  ng_dto                  Sort by n_g (ascending), then DTO (ascending).");
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

    let mut queue_strategy = "topo_ng".to_string();
    let mut topo_mode = false;
    let mut vg_only = false;
    for arg in args.iter().skip(1) {
        if arg.starts_with("--queue=") {
            queue_strategy = arg.split('=').nth(1).unwrap_or("topo_ng").to_string();
        }
        if arg == "--topo" {
            topo_mode = true;
        }
        if arg == "--vg-only" {
            vg_only = true;
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
    let (mut good, bad_ch, so_so) = compute_good_set(&ch_m_exp, l as f64);

    println!("Computing distance of good points to the origin...");
    good.fill_dist_to_origin(&bad_ch);
    println!("   ...done");

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

    if topo_mode || queue_strategy == "topo_ng" {
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
                println!(
                    "Topological order written to output/topo.txt. Points in topo order: {} / {}",
                    topo_order.len(),
                    good.num_points()
                );
                return Ok(());
            }
        } else {
            if queue_strategy == "topo_ng" {
                eprintln!("Error: Visibility graph is not a DAG, cannot use topo_ng strategy.");
                std::process::exit(1);
            }
        }
    }

    let mut log = String::new();

    write!(log, "# k: {}\n", k)?;
    writeln!(log, "# Good points: {}", good.length())?;
    println!("# Good points: {}", good.length());

    println!("Starting DP solver...");
    let start_dp = Instant::now();
    let dp_res = if queue_strategy == "perim_ng" {
        minimize_perimeter_dp::<PerimThenNgStrategy>(k, &good, &vg)
    } else if queue_strategy == "topo_ng" {
        minimize_perimeter_dp::<TopoThenNgStrategy>(k, &good, &vg)
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
    };
    let dp_duration = start_dp.elapsed();

    let Ok((sol, ub_circle, conf_count)) = dp_res else {
        eprintln!("Error: The DP solver failed to initialize or write to disk.");
        std::process::exit(1);
    };

    println!("k: {}", k);

    let sol_c = polygon_rm_redundant_vertices(&sol);

    let perimeter = compute_perimeter(&sol);
    let area = polygon_area(&sol);
    let perimeter_c = compute_perimeter(&sol_c);
    let area_c = polygon_area(&sol_c);
    let b_n = boundary_grid_points(&sol_c);

    assert!((area - area_c).abs() < 1e-6, "Something is wrong with area");
    assert!(
        (perimeter - perimeter_c).abs() < 1e-6,
        "Something is wrong with perimeter"
    );

    let v_n = sol_c.len();
    //let v_n_o = sol.len();

    /*if v_n != v_n_o {
        println!("Removed Redundant vertices");
    }*/
    // Pick's theorem:  A= I + B/2 -1
    // As k = I + B, we have I = k - B, so A = k - B + B/2 - 1 = k - B/2 - 1
    // A = k - B/2 - 1
    if (area - (k as f64) + (b_n as f64) / 2.0 + 1.0).abs() > 1e-6 {
        panic!("Area not computed correctly");
    }

    draw_polygon_with_grid(dir_pdfs, &sol, &ch_m, &ch_m_exp, k, ub_circle, &good);
    let ch_m_perimeter = compute_perimeter(&ch_m);

    // Explicitly tell the closure it returns a compatible Result
    let mut log_and_print =
        |label: &str, value: &dyn std::fmt::Display| -> Result<(), Box<dyn std::error::Error>> {
            writeln!(log, "# {:26} : {}", label, value)?;
            println!("# {:26} : {}", label, value);

            Ok(())
        };

    log_and_print("DP duration", &dp_duration.as_secs_f64())?;
    log_and_print("Perimeter", &perimeter)?;
    log_and_print("circle perimeter", &ch_m_perimeter)?;
    log_and_print("Naive perimeter", &ub_circle)?;
    log_and_print("Area", &area)?;
    let v_n_f64 = v_n as f64;
    log_and_print("vertices", &v_n_f64)?;
    let b_n_f64 = b_n as f64;
    log_and_print("boundary grid points", &b_n_f64)?;
    let conf_count_f64 = conf_count as f64;
    log_and_print("Configs computed", &conf_count_f64)?;

    let c_max_angle = compute_max_turn_angle(&sol);

    assert!(
        perimeter <= ch_m_perimeter * 1.00001,
        "perimeter not computed correctly A"
    );
    assert!(
        perimeter <= ub_circle * 1.00001,
        "perimeter not computed correctly B"
    );
    log_and_print("Max angle of sol", &c_max_angle)?;
    log_and_print("UB max turn_angle", &max_turn_angle)?;

    let len_p_longest = len_longest_primitive_edge(&sol);
    let len_longest = len_longest_edge(&sol);
    log_and_print("Longest edge len", &len_longest)?;
    log_and_print("Longest primitive edge len", &len_p_longest)?;
    log_and_print("UB primitive edge len", &max_edge_l)?;
    let duration = start.elapsed();

    let secs = duration.as_secs_f64();
    log_and_print("Running time in seconds", &secs)?;

    let filename_poly: String = format!("{}/{:06}_poly.txt", dir_polys, k);
    save_polygon(&filename_poly, &sol, Some(&log))?;

    let filename_log: String = format!("{}/{:06}_s.txt", dir_summary, k);
    let mut log_fl = File::create(filename_log)?;
    writeln!(log_fl, "{}", &log)?;

    Ok(())
}
