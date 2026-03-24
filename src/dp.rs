// dp.rs
use crate::point::*;
#[warn(unused_imports)]
use bytemuck::AnyBitPattern;

use mmap_vec::MmapVec;

use crate::geom::{GridSet, VisibilityGraph};
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
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, AnyBitPattern)]
pub struct DPStateValue {
    /// The key corresponding to this state.
    pub cfg: DPStateKey,
    /// Minimum perimeter found to reach this configuration.
    pub perimeter_so_far: f64,
    /// Index of the previous state in the `dp_vals` vector for solution reconstruction.
    pub prev_idx: u32,
}

/// Reconstructs the polygon vertices by following the `prev_idx` pointers.
pub fn extract_solution(
    dp_vals: &MmapVec<DPStateValue>,
    best_sol_idx: u32,
    good: &GridSet,
) -> Vec<Point2D> {
    let mut out = Vec::new();
    let mut curr_idx = best_sol_idx as usize;

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

// Wrapper for f64 to implement Ord
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrderedFloat(pub f64);

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

pub trait QueueStrategy {
    type Key: Ord;
    fn compute_key(
        n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        loc_id: u32,
        good: &GridSet,
        total_perim_contribution: f64,
        target_dtog: u32,
    ) -> Self::Key;
}

pub struct NgThenPerimStrategy;
impl QueueStrategy for NgThenPerimStrategy {
    type Key = (
        std::cmp::Reverse<u32>,
        std::cmp::Reverse<OrderedFloat>,
        std::cmp::Reverse<usize>,
    );

    fn compute_key(
        n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        _total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        (
            std::cmp::Reverse(n_g),
            std::cmp::Reverse(OrderedFloat(perimeter_so_far)),
            std::cmp::Reverse(idx),
        )
    }
}

pub struct PerimThenNgStrategy;
impl QueueStrategy for PerimThenNgStrategy {
    type Key = (
        std::cmp::Reverse<OrderedFloat>,
        std::cmp::Reverse<u32>,
        std::cmp::Reverse<usize>,
    );

    fn compute_key(
        n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        _total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        (
            std::cmp::Reverse(OrderedFloat(perimeter_so_far)),
            std::cmp::Reverse(n_g),
            std::cmp::Reverse(idx),
        )
    }
}

pub struct NgThenIdxStrategy;
impl QueueStrategy for NgThenIdxStrategy {
    type Key = (std::cmp::Reverse<u32>, std::cmp::Reverse<usize>);

    fn compute_key(
        n_g: u32,
        _perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        _total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        (std::cmp::Reverse(n_g), std::cmp::Reverse(idx))
    }
}

pub struct PerimThenIdxStrategy;
impl QueueStrategy for PerimThenIdxStrategy {
    type Key = (std::cmp::Reverse<OrderedFloat>, std::cmp::Reverse<usize>);

    fn compute_key(
        _n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        _total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        (
            std::cmp::Reverse(OrderedFloat(perimeter_so_far)),
            std::cmp::Reverse(idx),
        )
    }
}

pub struct NgPerimDtoStrategy;
impl QueueStrategy for NgPerimDtoStrategy {
    type Key = (
        std::cmp::Reverse<u32>,
        std::cmp::Reverse<OrderedFloat>,
        std::cmp::Reverse<usize>,
    );

    fn compute_key(
        n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        let total_perim = perimeter_so_far + total_perim_contribution;
        (
            std::cmp::Reverse(n_g),
            std::cmp::Reverse(OrderedFloat(total_perim)),
            std::cmp::Reverse(idx),
        )
    }
}

pub struct PerimNgDtogStrategy;
impl QueueStrategy for PerimNgDtogStrategy {
    type Key = (
        std::cmp::Reverse<OrderedFloat>,
        std::cmp::Reverse<i64>,
        std::cmp::Reverse<usize>,
    );

    fn compute_key(
        n_g: u32,
        perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        _total_perim_contribution: f64,
        target_dtog: u32,
    ) -> Self::Key {
        let total_g = n_g as i64 + target_dtog as i64;
        (
            std::cmp::Reverse(OrderedFloat(perimeter_so_far)),
            std::cmp::Reverse(total_g),
            std::cmp::Reverse(idx),
        )
    }
}

pub struct NgDtoStrategy;
impl QueueStrategy for NgDtoStrategy {
    type Key = (
        std::cmp::Reverse<u32>,
        std::cmp::Reverse<OrderedFloat>,
        std::cmp::Reverse<usize>,
    );

    fn compute_key(
        n_g: u32,
        _perimeter_so_far: f64,
        idx: usize,
        _loc_id: u32,
        _good: &GridSet,
        total_perim_contribution: f64,
        _target_dtog: u32,
    ) -> Self::Key {
        (
            std::cmp::Reverse(n_g),
            std::cmp::Reverse(OrderedFloat(total_perim_contribution)),
            std::cmp::Reverse(idx),
        )
    }
}

/// Item stored in the priority queue.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct QueueItem<K: Ord> {
    /// Sorting key determined by the selected `QueueStrategy`.
    key: K,
    /// Current number of enclosed points.
    n_g: u32,
    /// Index into the `dp_vals` vector.
    idx: usize,
}

/// Window size for the circular buffer of maps.
/// Must be larger than the maximum possible `n_g_delta`.
const D_ALL_WINDOW_SIZE: usize = 65536;

/// Context passed around during the DP execution.
pub struct DPContext<'a, K: Ord> {
    /// Total number of configurations explored.
    pub conf_count: &'a mut i64,
    /// Circular buffer of hash tables for deduplication: maps loc_id to an index in `dp_vals` for a specific n_g.
    pub d_all_window: &'a mut Vec<FxHashMap<u32, u32>>,
    /// Persistent storage for DP values, potentially backed by a memory-mapped file.
    pub dp_vals: &'a mut MmapVec<DPStateValue>,
    /// Priority queue of configurations to explore.
    pub pq: &'a mut BinaryHeap<QueueItem<K>>,
    /// Global upper bound on the optimal perimeter found so far.
    pub opt_perim: &'a mut f64,
    /// Index of the best complete solution found.
    pub best_sol_idx: &'a mut u32,
    /// Target number of grid points to enclose.
    pub k: usize,
    /// The set of valid grid points.
    pub good: &'a GridSet,
    /// Bitmask for periodic status updates.
    pub mask: u32,
    /// The visibility graph precalculated for the grid points.
    pub vg: &'a VisibilityGraph,
}

/// Prints current progress of the DP solver.
fn print_info<S: QueueStrategy>(ctx: &DPContext<S::Key>, n_g: u32) {
    let used_bytes = ctx.dp_vals.len() * std::mem::size_of::<DPStateValue>();
    let cap_bytes = ctx.dp_vals.capacity() * std::mem::size_of::<DPStateValue>();
    let used_mb = used_bytes / 1_048_576;
    let cap_mb = cap_bytes / 1_048_576;
    println!(
        "c: {}  n_g: {}  dp_vals (u/t): {} MB / {} MB",
        (*ctx.conf_count).to_formatted_string(&Locale::en),
        n_g,
        used_mb.to_formatted_string(&Locale::en),
        cap_mb.to_formatted_string(&Locale::en),
    );
}

/// Processes a single configuration by attempting to extend it with all visible neighbors.
fn process_configuration<S: QueueStrategy>(ctx: &mut DPContext<S::Key>, cfg_idx: usize) {
    *ctx.conf_count += 1;

    let (cfg_n_g, cfg_loc_id, perimeter_so_far) = {
        let val = &ctx.dp_vals[cfg_idx];
        (val.cfg.n_g, val.cfg.loc_id, val.perimeter_so_far)
    };

    if *ctx.conf_count & (ctx.mask as i64) == 0 {
        print_info::<S>(ctx, cfg_n_g);
    }

    // Capture current opt_perim to a local variable for faster access in the loop.
    let mut local_opt_perim = *ctx.opt_perim;

    // Neighbors are pre-calculated and sorted by total_perim_contribution in the visibility graph.
    let nbrs = &ctx.vg.adjacency_list[cfg_loc_id as usize];

    for edge in nbrs.iter() {
        // Pruning: skip if the current path already exceeds the global best perimeter.
        // Since nbrs are sorted by total_perim_contribution, we can break early.
        let total_perim = perimeter_so_far + edge.total_perim_contribution;
        if total_perim > local_opt_perim {
            break;
        }

        let next_n_g = cfg_n_g + edge.n_g_delta;
        // Optimization: early exit if too many grid points are enclosed.
        if next_n_g as usize > ctx.k {
            continue;
        }

        // Tighter grid point pruning: check if the return path already exceeds k.
        if next_n_g as usize + edge.target_dtog as usize > ctx.k {
            continue;
        }

        let next_perim = perimeter_so_far + edge.edge_len;

        let map_idx = (next_n_g as usize) % D_ALL_WINDOW_SIZE;
        let layer_map = &mut ctx.d_all_window[map_idx];
        
        let mut f_queued = false;
        let mut existing_val_idx = u32::MAX;

        // Deduplication: check if this location was already reached at this n_g with a better perimeter.
        if let Some(&idx) = layer_map.get(&edge.target_id) {
            existing_val_idx = idx;
            f_queued = true;
            if ctx.dp_vals[idx as usize].perimeter_so_far <= next_perim {
                continue;
            }
        }

        let next_cfg = DPStateKey {
            loc_id: edge.target_id,
            n_g: next_n_g,
        };

        let next_val = DPStateValue {
            cfg: next_cfg,
            perimeter_so_far: next_perim,
            prev_idx: cfg_idx as u32,
        };

        let push_idx;
        if existing_val_idx != u32::MAX {
            // Update the existing state with the new, better perimeter.
            ctx.dp_vals[existing_val_idx as usize] = next_val;
            push_idx = existing_val_idx as usize;
        } else {
            // Register a newly discovered state.
            push_idx = ctx.dp_vals.len();
            ctx.dp_vals.push(next_val).expect("push failed");
            layer_map.insert(edge.target_id, push_idx as u32);
        }

        // Add to priority queue only if it's a new state or not currently queued.
        if !f_queued {
            ctx.pq.push(QueueItem {
                key: S::compute_key(
                    next_n_g,
                    next_perim,
                    push_idx,
                    edge.target_id,
                    ctx.good,
                    edge.total_perim_contribution,
                    edge.target_dtog,
                ),
                n_g: next_n_g,
                idx: push_idx,
            });
        }

        // If we've enclosed exactly k points, check if we've found a new global optimum.
        if next_n_g as usize == ctx.k {
            if total_perim < local_opt_perim {
                local_opt_perim = total_perim;
                *ctx.best_sol_idx = push_idx as u32;
            }
        }
    }
    
    // Write back updated opt_perim.
    *ctx.opt_perim = local_opt_perim;
}

pub fn max_edge_length(k: usize) -> u32 {
    (k as f64).powf(1.0 / 3.0).round() as u32 + 1
}

pub fn minimize_perimeter_dp<S: QueueStrategy>(
    k: usize,
    good: &GridSet,
    vg: &VisibilityGraph,
) -> anyhow::Result<(Vec<Point2D>, f64, i64)> {
    let sqrt_k = (k as f64).sqrt().ceil() as u32 + 1;
    let ub_circle = 2.0 * (std::f64::consts::PI * (k as f64)).sqrt();
    let sq_perim = (4 * sqrt_k + 6) as f64;
    let power = 19;
    let mask = (1 << power) - 1;

    let mut opt_perim = ub_circle.min(sq_perim);

    let start_loc_id = good.get_point_id(Point2D::new(0, 0));
    let (start_dto, start_dtog) = good.get_dto(Point2D::new(0, 0));
    let start_key = DPStateKey {
        loc_id: start_loc_id as u32,
        n_g: 1,
    };

    let mut best_sol_idx = 0;

    let mut pq = BinaryHeap::new();
    pq.push(QueueItem {
        key: S::compute_key(
            1,
            0.0,
            0,
            start_loc_id as u32,
            good,
            start_dto,
            start_dtog as u32,
        ),
        n_g: 1,
        idx: 0,
    });

    let hint_size = if k > 100000 {
        100_000_000
    } else {
        std::cmp::min(2 * k * k, 10_000_000)
    };

    // Initialize circular buffer of maps.
    let mut d_all_window = Vec::with_capacity(D_ALL_WINDOW_SIZE);
    for _ in 0..D_ALL_WINDOW_SIZE {
        d_all_window.push(FxHashMap::default());
    }
    d_all_window[1 % D_ALL_WINDOW_SIZE].insert(start_loc_id as u32, 0);

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
        d_all_window: &mut d_all_window,
        dp_vals: &mut dp_vals,
        pq: &mut pq,
        opt_perim: &mut opt_perim,
        best_sol_idx: &mut best_sol_idx,
        k,
        good,
        mask,
        vg,
    };

    let mut ids = Vec::with_capacity(1024);
    let mut last_cleared_n_g = 0;

    loop {
        let item = ctx.pq.pop();
        if item.is_none() {
            break;
        }
        let QueueItem {
            n_g: curr_n_g,
            idx: popped_idx,
            ..
        } = item.unwrap();

        ids.clear();
        ids.push(popped_idx);
        while let Some(next_item) = ctx.pq.peek() {
            if next_item.n_g == curr_n_g {
                if let Some(matched_item) = ctx.pq.pop() {
                    ids.push(matched_item.idx);
                }
            } else {
                break;
            }
        }

        // Clear maps in the circular buffer that are no longer reachable.
        // We clear from last_cleared_n_g up to curr_n_g - 1.
        for n_g_to_clear in last_cleared_n_g..curr_n_g {
            ctx.d_all_window[n_g_to_clear as usize % D_ALL_WINDOW_SIZE].clear();
        }
        last_cleared_n_g = curr_n_g;

        for &id in &ids {
            process_configuration::<S>(&mut ctx, id);
        }
    }

    println!("# of configurations generated: {}", *ctx.conf_count);
    let sol = extract_solution(&ctx.dp_vals, *ctx.best_sol_idx, ctx.good);
    Ok((sol, ub_circle, *ctx.conf_count))
}
