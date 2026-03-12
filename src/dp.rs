// dp.rs
//mod point;
use  crate::point::*;

use crate::geom::{
    is_all_left_turns, is_colinear, is_lefteq_turn, is_right_turn, triangle_count_new_points,
    GridSet,
};
use num_format::{Locale, ToFormattedString};
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DPStateKey {
    pub dir_index: u16,
    pub loc: Point2D,
    pub n_g: u32,
}

impl PartialOrd for DPStateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DPStateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.n_g
            .cmp(&other.n_g)
            .then(self.loc.cmp(&other.loc))
            .then(self.dir_index.cmp(&other.dir_index))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DPStateValue {
    pub cfg: DPStateKey,
    pub perimeter_so_far: f64,
    pub prev_idx: usize,
    pub handled: bool,
}

pub fn comp_next_conf(
    cfg: &DPStateKey,
    perimeter_so_far: f64,
    i_v: u16,
    k: u32,
    sqrt_k: u32,
    opt_perim: f64,
    dirs: &Vec<Point2D>,
) -> (bool, DPStateKey, f64) {
    let a = cfg.loc - dirs[cfg.dir_index as usize];
    let b = cfg.loc;
    let c = cfg.loc + dirs[i_v as usize];

    let origin = Point2D::new(0, 0);

    if !is_lefteq_turn(a, b, c)
        || (c.y < 0)
        || (a == c)
        || c.is_zero()
        || (c.x as i64) > (sqrt_k as i64)
        || (c.y as i64) > ((2 * sqrt_k) as i64)
        || (-c.x as i64) > (sqrt_k as i64)
        || (-c.y as i64) > (sqrt_k as i64)
        || perimeter_so_far > opt_perim
        || is_right_turn(b, c, origin)
    {
        return (false, *cfg, -1.0);
    }

    if is_colinear(b, c, origin) {
        if c.norm_sq() < b.norm_sq() {
            return (false, *cfg, -1.0);
        }
    }

    let &v = &dirs[i_v as usize];
    assert!(b != c);
    assert!(v.x != 0 || v.y != 0);

    let (tri_i_new, tri_b_new) = triangle_count_new_points(origin, b, c);

    //assert!(tri_i_new >= 0 && tri_b_new >= 0);
    let n_g: u32 = cfg.n_g + tri_i_new + tri_b_new;

    if n_g > 0 && c.is_zero() {
        return (false, *cfg, -1.0);
    }

    if n_g > k {
        return (false, *cfg, -1.0);
    }

    let new_cfg = DPStateKey {
        dir_index: i_v,
        loc: c,
        n_g,
    };

    let new_perim = perimeter_so_far + (b - c).norm();
    (true, new_cfg, new_perim)
}

pub fn extract_solution(
    d_all: &FxHashMap<DPStateKey, usize>,
    dp_vals: &[DPStateValue],
    best_sol: DPStateKey,
) -> Vec<Point2D> {
    let mut out = Vec::new();
    let mut curr_idx = *d_all.get(&best_sol).unwrap_or_else(|| {
        eprintln!(
            "Error in extract_solution: Key '{:#?}' not found in map.",
            best_sol
        );
        std::process::exit(1);
    });

    loop {
        let val = &dp_vals[curr_idx];
        out.push(val.cfg.loc);
        if val.prev_idx == curr_idx {
            break;
        }
        curr_idx = val.prev_idx;
    }

    // The Julia code pushes in reverse chronological order (end -> start)
    // We reverse it to match counter-clockwise tracing from origin if needed.
    // Wait, julia implementation extracts backwards but does not reverse:
    // `extract_solution` returns `out` directly.
    out
}

pub fn is_store(
    d_all: &FxHashMap<DPStateKey, usize>,
    dp_vals: &[DPStateValue],
    next_cfg: &DPStateKey,
    new_perim: f64,
    opt_perim: f64,
    k: u32,
) -> bool {
    if new_perim > opt_perim {
        return false;
    }
    if next_cfg.n_g > k {
        println!("FLOGI");
        return false;
    }
    if let Some(&idx) = d_all.get(next_cfg) {
        if dp_vals[idx].perimeter_so_far < new_perim {
            return false;
        }
    }
    true
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
        loc: Point2D,
        good: &GridSet,
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
        _loc: Point2D,
        _good: &GridSet,
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
        _loc: Point2D,
        _good: &GridSet,
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
        _loc: Point2D,
        _good: &GridSet,
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
        _loc: Point2D,
        _good: &GridSet,
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
        loc: Point2D,
        good: &GridSet,
    ) -> Self::Key {
        let total_perim = perimeter_so_far + good.get_dto(loc).0;
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
        loc: Point2D,
        good: &GridSet,
    ) -> Self::Key {
        let total_g = n_g as i64 + good.get_dto(loc).1;
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
        loc: Point2D,
        good: &GridSet,
    ) -> Self::Key {
        let dto = good.get_dto(loc).0;
        (
            std::cmp::Reverse(n_g),
            std::cmp::Reverse(OrderedFloat(dto)),
            std::cmp::Reverse(idx),
        )
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub struct QueueItem<K: Ord> {
    key: K,
    n_g: u32,
    idx: usize,
}

pub fn minimize_perimeter_dp<S: QueueStrategy>(
    k: u32,
    good: &GridSet,
    bad: &GridSet,
    bad_in_ch: &Vec<Point2D>,
    use_cache: bool,
) -> (Vec<Point2D>, f64) {
    let sqrt_k = (k as f64).sqrt().ceil() as u32 + 1;
    // max_angle = 3 * pi / k^(1/3)
    let max_angle = 5.0 / (k as f64).powf(1.0 / 3.0);
    let ub_circle = 2.0 * (std::f64::consts::PI * (k as f64)).sqrt();
    let sq_perim = (4 * sqrt_k + 6) as f64;
    let max_edge_l: u32 = (k as f64).powf(1.0 / 3.0).round() as u32 + 1;
    let power = 19;
    let mask = (1 << power) - 1;

    println!("k           : {}", k);
    println!("sqrt_k      : {}", sqrt_k);
    println!("max_edge_l  : {}", max_edge_l);
    println!("max angle   : {}", max_angle);

    let mut opt_perim = ub_circle.min(sq_perim);

    let dirs = crate::geom::generate_primitive_vectors(max_edge_l);
    let stops = crate::geom::comp_stop_indexes(&dirs, max_angle);

    let start_key = DPStateKey {
        dir_index: 0,
        loc: Point2D::new(0, 0),
        n_g: 1,
    };

    let mut best_sol = start_key;

    let mut pq = BinaryHeap::new();
    pq.push(QueueItem {
        key: S::compute_key(1, 0.0, 0, start_key.loc, good),
        n_g: 1,
        idx: 0,
    });

    let hint_size = 2 * (k as usize) * (k as usize);

    println!(
        "hint_size   : {}",
        (hint_size as i64).to_formatted_string(&Locale::en)
    );

    let mut d_all = FxHashMap::default();
    //let mut d_all = FxHashMap::with_capacity( hint_size );
    let mut all_left_turns_cache: FxHashMap<(Point2D, Point2D), bool> = FxHashMap::default();
    d_all.reserve(hint_size);
    let mut dp_vals = Vec::with_capacity(hint_size);
    let start_val = DPStateValue {
        cfg: start_key,
        perimeter_so_far: 0.0,
        prev_idx: 0,
        handled: false,
    };
    if let Err(_) = d_all.try_reserve(1) {
        let msg = "Error: Out of memory when attempting to insert into hash table";
        println!("{}", msg);
        eprintln!("{}", msg);
        std::process::exit(1);
    }
    d_all.insert(start_key, 0);
    dp_vals.push(start_val);

    let mut conf_count: i64 = 0;
    let mut conf_useless_count: i64 = 0;
    while let Some(QueueItem {
        key: _,
        n_g: _,
        idx: popped_idx,
    }) = pq.pop()
    {
        let mut store_count = 0;
        conf_count += 1;

        let cfg_idx = popped_idx;
        let (cfg, perimeter_so_far) = {
            let val = &mut dp_vals[cfg_idx];
            if val.handled {
                println!("HANDLED???");
                continue;
            }
            val.handled = true;
            (val.cfg, val.perimeter_so_far)
        };

        if conf_count & mask == 0 {
            let used_bytes = dp_vals.len() * std::mem::size_of::<DPStateValue>();
            let cap_bytes = dp_vals.capacity() * std::mem::size_of::<DPStateValue>();
            let used_mb = used_bytes / 1_048_576;
            let cap_mb = cap_bytes / 1_048_576;
            println!(
                "c: {}  n_g: {}  dp_vals (u/t): {} MB / {} MB  d_all {}",
                conf_count.to_formatted_string(&Locale::en),
                cfg.n_g,
                used_mb.to_formatted_string(&Locale::en),
                cap_mb.to_formatted_string(&Locale::en),
                d_all.len().to_formatted_string(&Locale::en)
            );
        }
        let dir_start = cfg.dir_index;

        let dir_end = stops[dir_start as usize];
        let _prev_dir = dirs[dir_start as usize];

        for dir_i in (dir_start as usize)..=dir_end {
            let v = dirs[dir_i];

            if cfg.loc.y == 0 && v.y <= 0 {
                continue;
            }

            let p_next = cfg.loc + v;
            if bad.contains(&p_next) {
                continue;
            }

            let cache_key = (cfg.loc, p_next);
            let is_valid_turn = if use_cache {
                if let Some(&res) = all_left_turns_cache.get(&cache_key) {
                    res
                } else {
                    let res = is_all_left_turns(cfg.loc, p_next, bad_in_ch);
                    if let Err(_) = all_left_turns_cache.try_reserve(1) {
                        let msg = "Error: Out of memory when attempting to insert into hash table";
                        println!("{}", msg);
                        eprintln!("{}", msg);
                        std::process::exit(1);
                    }
                    all_left_turns_cache.insert(cache_key, res);
                    res
                }
            } else {
                is_all_left_turns(cfg.loc, p_next, bad_in_ch)
            };

            if !is_valid_turn {
                //                println!("BANGO!");
                continue;
            }

            let (f_valid, next_cfg, new_perim) = comp_next_conf(
                &cfg,
                perimeter_so_far,
                dir_i as u16,
                k,
                sqrt_k,
                opt_perim,
                &dirs,
            );
            if !f_valid {
                continue;
            }

            // let total_perim = new_perim + good.dto( next_cfg.loc.norm();
            let total_perim = new_perim + good.get_dto(next_cfg.loc).0;
            if total_perim > opt_perim {
                continue;
            }

            let mut f_queued = false;
            let mut existing_idx = None;
            if let Some(&idx) = d_all.get(&next_cfg) {
                f_queued = true;
                existing_idx = Some(idx);
            }

            if !is_store(&d_all, &dp_vals, &next_cfg, new_perim, opt_perim, k) {
                continue;
            }

            store_count += 1;
            let next_val = DPStateValue {
                cfg: next_cfg,
                perimeter_so_far: new_perim,
                prev_idx: cfg_idx,
                handled: false,
            };
            let push_idx;
            if let Some(idx) = existing_idx {
                dp_vals[idx] = next_val;
                push_idx = idx;
            } else {
                let idx = dp_vals.len();
                dp_vals.push(next_val);
                if let Err(_) = d_all.try_reserve(1) {
                    let msg = "Error: Out of memory when attempting to insert into hash table";
                    println!("{}", msg);
                    eprintln!("{}", msg);
                    std::process::exit(1);
                }
                d_all.insert(next_cfg, idx);
                push_idx = idx;
            }

            if !f_queued {
                pq.push(QueueItem {
                    key: S::compute_key(next_cfg.n_g, new_perim, push_idx, next_cfg.loc, good),
                    n_g: next_cfg.n_g,
                    idx: push_idx,
                });
            } else {
            }

            if next_cfg.n_g == k {
                if total_perim < opt_perim {
                    opt_perim = total_perim;
                    best_sol = next_cfg;
                }
            }
        }
        if store_count == 0 {
            //d_all.remove(&cfg);
            conf_useless_count += 1;
            //            println!("Woga?");
        }
    }
    println!("# of configurations generated: {}", conf_count);
    println!("# of dead-end  configurations: {}", conf_useless_count);
    let sol = extract_solution(&d_all, &dp_vals, best_sol);
    (sol, ub_circle)
}

// End of file
