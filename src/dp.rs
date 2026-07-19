// dp.rs
//! This module implements the Dynamic Programming (DP) solver for the k-perimeter problem.
//! It uses a priority queue to explore configurations in a topological order,
//! ensuring that we find the minimum perimeter for exactly k points.
//!
//! Large state spaces are handled using memory-mapped files via `mmap_vec`.

use crate::point::*;
use bytemuck::AnyBitPattern;
use std::cmp::max;

use mmap_vec::MmapVec;

use crate::geom::GridSet;
use crate::v_graph::VisibilityGraph;
use num_format::{Locale, ToFormattedString};
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Key for the DP state.
/// Represents a configuration by the current point (`loc_id`) and the number of grid points enclosed (`n_g`).
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, AnyBitPattern)]
pub struct DPStateKey {
    /// Index of the current point in the GridSet.
    pub loc_id: u32,
    /// Total number of grid points enclosed by the polygon so far.
    pub n_g: u32,
}

impl PartialOrd for DPStateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DPStateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        // Sort primarily by n_g to process configurations in a natural order for the DP.
        self.n_g
            .cmp(&other.n_g)
            .then(self.loc_id.cmp(&other.loc_id))
    }
}

/// Value stored for each DP state.
/// Tracks the optimal perimeter and the path for reconstruction.
///
/// This struct is stored in a memory-mapped file, so it must be `#[repr(C)]`
/// and implement `AnyBitPattern` for zero-copy I/O.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, AnyBitPattern)]
pub struct DPStateValue {
    /// The key corresponding to this state.
    pub cfg: DPStateKey,
    /// Minimum perimeter found to reach this configuration.
    pub perimeter_so_far: f64,
    /// Index of the previous state in the `dp_vals` vector for solution reconstruction.
    /// This allows us to trace back the vertices of the optimal polygon.
    pub prev_idx: u64,
}

/// Reconstructs the polygon vertices by following the `prev_idx` pointers.
///
/// Starting from the best complete solution, it jumps back through the `dp_vals`
/// vector until it reaches the origin (where `prev_idx` points to itself).
pub fn extract_solution(
    dp_vals: &MmapVec<DPStateValue>,
    best_sol_idx: usize,
    good: &GridSet,
) -> Vec<Point2D> {
    let mut out = Vec::new();
    let mut curr_idx = best_sol_idx;

    loop {
        let val = &dp_vals[curr_idx];
        out.push(good.get_point_by_id(val.cfg.loc_id as usize).clone());
        if val.prev_idx as usize == curr_idx {
            break;
        }
        curr_idx = val.prev_idx as usize;
    }

    out
}

/// Item stored in the priority queue.
///
/// The `key` is the `Reverse(topo_idx)`, which ensures that the `BinaryHeap`
/// (which is a max-heap) behaves as a min-heap for the topological index.
/// This processes points in a strict topological order.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct QueueItem {
    /// Sorting key: Reverse(topo_idx). Smaller topo indices are popped first.
    pub key: std::cmp::Reverse<u32>,
    /// Index into the `dp_vals` vector where the full state details are stored.
    pub idx: usize,
    /// Location ID of the current point.
    pub loc_id: u32,
}

/// Context passed around during the DP execution.
///
/// Encapsulates all shared state to avoid passing multiple mutable references.
pub struct DPContext<'a> {
    /// Total number of configurations explored.
    pub conf_count: &'a mut i64,
    /// Hash table for deduplication: maps (loc_id, n_g) to an index in `dp_vals`.
    /// Uses a packed u64 key (loc_id << 32 | n_g) for performance.
    pub d_all: &'a mut FxHashMap<u64, u64>,
    /// Persistent storage for DP values, backed by a memory-mapped file for large k.
    pub dp_vals: &'a mut MmapVec<DPStateValue>,
    /// Priority queue of configurations to explore.
    pub pq: &'a mut BinaryHeap<QueueItem>,
    /// Global upper bound on the optimal perimeter found so far.
    pub opt_perim: &'a mut f64,
    /// Index of the best complete solution found.
    pub best_sol_idx: &'a mut Option<usize>,
    /// Target number of grid points to enclose.
    pub k: usize,
    /// The set of valid grid points.
    pub good: &'a GridSet,
    /// Bitmask for periodic status updates (e.g., print info every 2^19 configs).
    pub mask: u32,
    /// The visibility graph precalculated for the grid points.
    pub vg: &'a VisibilityGraph,
    /// The factor by which the dynamic threshold is increased after cleanup.
    pub retain_factor: f64,
    pub start_dp: std::time::Instant,
}

/// Prints current progress of the DP solver, including memory usage of mmap storage.
fn print_info(ctx: &DPContext, l_i: u32) {
    let used_bytes = ctx.dp_vals.len() * std::mem::size_of::<DPStateValue>();
    let used_mb = used_bytes / 1_048_576;
    crate::log_println!(
        "k: {} c: {:>12}  {:3}%  dpA: {:>6}MB  hash: {}  {:.2}s",
        ctx.k,
        (*ctx.conf_count).to_formatted_string(&Locale::en),
        (ctx.good.get_topo_idx(l_i as usize) as f64 / ctx.good.num_points() as f64 * 100.0).round(),
        used_mb.to_formatted_string(&Locale::en),
        ctx.d_all.len().to_formatted_string(&Locale::en),
        ctx.start_dp.elapsed().as_secs_f64()
    );
}

/// Packs loc_id and n_g into a single u64 for efficient hashing.
#[inline(always)]
fn make_key(loc_id: u32, n_g: u32) -> u64 {
    ((loc_id as u64) << 32) | (n_g as u64)
}

/// Processes a set of configurations with the same loc_id by attempting to extend them with all visible neighbors.
///
/// Optimization: All ids in the slice share the same `loc_id`. We fetch the neighbor list once.
/// Optimization: The ids are sorted by their index in `dp_vals` to improve cache locality during sequential reads.
fn process_configurations(ctx: &mut DPContext, mut ids: Vec<usize>) {
    if ids.is_empty() {
        return;
    }

    // Sort by cfg_idx to improve memory locality when accessing ctx.dp_vals.
    // Sequential access is critical for memory-mapped storage performance.
    ids.sort_unstable();

    // Since all ids share the same loc_id, we fetch neighbors once from the adjacency list.
    let loc_id = ctx.dp_vals[ids[0]].cfg.loc_id as usize;
    let nbrs = &ctx.vg.adjacency_list[loc_id];

    for cfg_idx in ids {
        *ctx.conf_count += 1;

        let (cfg, perimeter_so_far, prev_idx) = {
            let val = &ctx.dp_vals[cfg_idx];
            (val.cfg, val.perimeter_so_far, val.prev_idx as usize)
        };

        let (start_idx, end_idx) = if prev_idx == cfg_idx {
            // Start state
            (0, nbrs.len())
        } else {
            // Find the edge that led to loc_id from prev_loc_id
            let prev_loc_id = ctx.dp_vals[prev_idx].cfg.loc_id as usize;
            let prev_nbrs = &ctx.vg.adjacency_list[prev_loc_id];

            // Search for the edge (prev_loc_id, loc_id)
            let mut found_range = None;
            for edge in prev_nbrs {
                if edge.target_id == loc_id {
                    found_range = Some((
                        edge.next_edge_start_idx as usize,
                        edge.next_edge_end_idx as usize,
                    ));
                    break;
                }
            }
            found_range.expect("Incoming edge not found in visibility graph")
        };

        // Explore outgoing edges (pre-filtered by turn angle and convexity).
        for edge in nbrs[start_idx..end_idx].iter() {
            let next_n_g = cfg.n_g + edge.n_g_delta;

            // Early pruning: don't exceed target k.
            if next_n_g as usize > ctx.k {
                continue;
            }
            // Early pruning: ensure we can still reach k from here.
            if (next_n_g + edge.max_addion_g) < (ctx.k as u32) {
                continue;
            }

            let next_perim = perimeter_so_far + edge.edge_len;
            // Pruning: skip if the current path already exceeds the global best perimeter.
            if next_perim > *ctx.opt_perim {
                continue;
            }

            // Admissibility heuristic: total_perim = perimeter so far + min distance back to origin.
            let total_perim = next_perim + edge.target_dto;
            if total_perim > *ctx.opt_perim {
                continue;
            }

            let key = make_key(edge.target_id as u32, next_n_g);
            let mut f_queued = false;
            let mut existing_val_idx = None;

            // Deduplication: check if this (location, n_g) state was already reached with a better perimeter.
            if let Some(&idx) = ctx.d_all.get(&key) {
                let idx_u = idx as usize;
                existing_val_idx = Some(idx_u);
                f_queued = true;
                if ctx.dp_vals[idx_u].perimeter_so_far <= next_perim {
                    continue;
                }
            }

            let next_cfg = DPStateKey {
                loc_id: edge.target_id as u32,
                n_g: next_n_g,
            };

            let next_val = DPStateValue {
                cfg: next_cfg,
                perimeter_so_far: next_perim,
                prev_idx: cfg_idx as u64,
            };

            let push_idx: usize;
            if let Some(idx_u) = existing_val_idx {
                // Update the existing state with the new, better perimeter.
                ctx.dp_vals[idx_u] = next_val;
                push_idx = idx_u;
            } else {
                // Register a newly discovered state.
                push_idx = ctx.dp_vals.len();
                ctx.dp_vals
                    .push(next_val)
                    .expect("Failed to push to DP state MmapVec (possible disk/memory I/O error)");
                ctx.d_all.insert(key, push_idx as u64);
            }

            // Add to priority queue only if it's a new state or not currently queued.
            if !f_queued {
                ctx.pq.push(QueueItem {
                    key: std::cmp::Reverse(ctx.good.get_topo_idx(edge.target_id)),
                    idx: push_idx,
                    loc_id: edge.target_id as u32,
                });
            }

            // If we've enclosed exactly k points, check if we've found a new global optimum.
            if next_n_g as usize == ctx.k {
                if total_perim < *ctx.opt_perim {
                    *ctx.opt_perim = total_perim;
                    *ctx.best_sol_idx = Some(push_idx);
                }
            }
        }
    }
}

/// Heuristic for the maximum allowed length of a single polygon segment based on k.
pub fn max_edge_length(k: usize) -> u32 {
    //    ((k as f64).powf(1.0 / 3.0) / 2.0).round() as u32 + 2
    // Safe:
    // (5.0 + 2.0 * (k as f64).powf(1.0 / 6.0)).round() as u32 + 2
    (2.0 + 1.5 * (k as f64).powf(1.0 / 6.0)).round() as u32 + 2
}

/// Entry point for the DP solver.
///
/// Finds the minimum perimeter polygon enclosing exactly k grid points.
/// Returns the vertices, the naive perimeter upper bound, and the total configuration count.
pub fn minimize_perimeter_dp(
    k: usize,
    good: &GridSet,
    vg: &VisibilityGraph,
    retain_factor: f64,
) -> anyhow::Result<(Vec<Point2D>, f64, i64)> {
    let sqrt_k = (k as f64).sqrt().ceil() as u32 + 1;
    let ub_circle = 2.0 * (std::f64::consts::PI * (k as f64)).sqrt();
    let sq_perim = (4 * sqrt_k + 6) as f64;
    let power = 20;
    let mask = (1 << power) - 1;

    // Start with a naive upper bound on the perimeter.
    let mut opt_perim = 0.2 + ub_circle.min(sq_perim);

    // Initial state: origin (0,0) which encloses 1 point by definition.
    let start_loc_id = good.get_point_id(Point2D::new(0, 0));
    let start_key = DPStateKey {
        loc_id: start_loc_id as u32,
        n_g: 1,
    };

    let mut best_sol_idx = None;

    let mut pq = BinaryHeap::new();
    pq.push(QueueItem {
        key: std::cmp::Reverse(good.get_topo_idx(start_loc_id)),
        idx: 0,
        loc_id: start_loc_id as u32,
    });

    // Estimate total capacity for the mmap vector.
    let hint_size = if k > 100000 {
        10_000_000
    } else {
        std::cmp::min(2 * k * k, 10_000_000)
    };
    let mut threshold = std::cmp::min(k * 100, 10_000_000);

    let mut d_all = FxHashMap::default();
    d_all.reserve(std::cmp::min(threshold, 1_000_000));
    d_all.insert(make_key(start_loc_id as u32, 1), 0);

    // MmapVec automatically uses a temporary file on disk if the capacity is large.
    let mut dp_vals = MmapVec::<DPStateValue>::with_capacity(hint_size)?;

    let start_val = DPStateValue {
        cfg: start_key,
        perimeter_so_far: 0.0,
        prev_idx: 0,
    };
    dp_vals.push(start_val)?;

    let mut conf_count: i64 = 0;
    let mut ctx = DPContext {
        conf_count: &mut conf_count,
        d_all: &mut d_all,
        dp_vals: &mut dp_vals,
        pq: &mut pq,
        opt_perim: &mut opt_perim,
        best_sol_idx: &mut best_sol_idx,
        k,
        good,
        mask,
        vg,
        retain_factor,
        start_dp: std::time::Instant::now(),
    };

    let mut last_loc_popped = 0;
    // Main DP loop.
    loop {
        let item = ctx.pq.pop();
        if item.is_none() {
            print_info(&ctx, last_loc_popped);
            break;
        }
        let QueueItem {
            idx: popped_idx,
            loc_id: popped_loc_id,
            ..
        } = item.unwrap();

        last_loc_popped = popped_loc_id;
        // Group together all items in the queue that have the same loc_id.
        // Since we sort by topo_idx, all items with the same loc_id are adjacent in the queue.
        let mut ids = vec![popped_idx];
        while let Some(next_item) = ctx.pq.peek() {
            if next_item.loc_id == popped_loc_id {
                if let Some(matched_item) = ctx.pq.pop() {
                    ids.push(matched_item.idx);
                }
            } else {
                break;
            }
        }

        // Periodic pruning of the deduplication map.
        // We can safely remove any state with a topo_idx smaller than the current one,
        // because in a DAG, those states can never be improved or reached again.
        if ctx.d_all.len() > threshold {
            let min_topo_idx = ctx.good.get_topo_idx(popped_loc_id as usize);
            ctx.d_all.retain(|&key, _| {
                let loc_id = (key >> 32) as u32;
                ctx.good.get_topo_idx(loc_id as usize) >= min_topo_idx
            });
            threshold = max(
                threshold,
                (ctx.retain_factor * ctx.d_all.len() as f64) as usize,
            );
        }

        let old_mask = *ctx.conf_count & (ctx.mask as i64);
        // Process the grouped configurations as a single batch.
        process_configurations(&mut ctx, ids);

        if *ctx.conf_count & (ctx.mask as i64) < old_mask {
            print_info(&ctx, popped_loc_id);
        }
    }

    assert!(
        ctx.best_sol_idx.is_some(),
        "No solution found by DP solver. This should never happen since the origin itself is a valid solution."
    );
    let sol = extract_solution(&ctx.dp_vals, ctx.best_sol_idx.unwrap(), ctx.good);
    Ok((sol, ub_circle, *ctx.conf_count))
}
