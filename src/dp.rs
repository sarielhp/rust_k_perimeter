// dp.rs
//mod point;
use crate::point::*;
#[warn(unused_imports)]
use bytemuck::{AnyBitPattern};
use std::cmp::max;
use std::mem;

use mmap_vec::{MmapVec};

use crate::geom::{
    is_all_left_turns, is_colinear, is_lefteq_turn, is_right_turn, triangle_count_new_points,
    GridSet,
};
use num_format::{Locale, ToFormattedString};
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, AnyBitPattern)]
pub struct DPStateKey {
    pub loc: Point2D,
    pub n_g: u32,
    pub dir_index: u16,
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

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, AnyBitPattern)]
//#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DPStateValue {
    pub cfg: DPStateKey,
    pub perimeter_so_far: f64,
    pub prev_idx: usize
    //pub handled: bool,
}

pub fn comp_next_conf<K: Ord>(
    ctx: &DPContext<K>,
    cfg: &DPStateKey,
    perimeter_so_far: f64,
    i_v: u16,
) -> (bool, DPStateKey, f64) {
    let a = cfg.loc - ctx.dirs[cfg.dir_index as usize];
    let b = cfg.loc;
    let c = cfg.loc + ctx.dirs[i_v as usize];

    let origin = Point2D::new(0, 0);

    if !is_lefteq_turn(a, b, c)
        || (c.y < 0)
        || (a == c)
        || c.is_zero()
        || (c.x as i64) > (ctx.sqrt_k as i64)
        || (c.y as i64) > ((2 * ctx.sqrt_k) as i64)
        || (-c.x as i64) > (ctx.sqrt_k as i64)
        || (-c.y as i64) > (ctx.sqrt_k as i64)
        || perimeter_so_far > *ctx.opt_perim
        || is_right_turn(b, c, origin)
    {
        return (false, *cfg, -1.0);
    }

    if is_colinear(b, c, origin) {
        if c.norm_sq() < b.norm_sq() {
            return (false, *cfg, -1.0);
        }
    }

    let &v = &ctx.dirs[i_v as usize];
    assert!(b != c);
    assert!(v.x != 0 || v.y != 0);

    let (tri_i_new, tri_b_new) = triangle_count_new_points(origin, b, c);

    //assert!(tri_i_new >= 0 && tri_b_new >= 0);
    let n_g: u32 = cfg.n_g + tri_i_new + tri_b_new;

    if n_g > 0 && c.is_zero() {
        return (false, *cfg, -1.0);
    }

    if n_g > ctx.k {
        return (false, *cfg, -1.0);
    }

    let new_cfg = DPStateKey {
        dir_index: i_v,
        loc: c,
        n_g
        };

    let new_perim = perimeter_so_far + (b - c).norm();
    (true, new_cfg, new_perim)
}

pub fn extract_solution<K: Ord>(ctx: &DPContext<K>) -> Vec<Point2D> {
    let mut out = Vec::new();
    let mut curr_idx = *ctx.d_all.get(&*ctx.best_sol).unwrap_or_else(|| {
        eprintln!(
            "Error in extract_solution: Key '{:#?}' not found in map.",
            *ctx.best_sol
        );
        std::process::exit(1);
    });

    loop {
        let val = &ctx.dp_vals[curr_idx];
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

pub fn is_store<K: Ord>(
    ctx: &DPContext<K>,
    next_cfg: &DPStateKey,
    new_perim: f64,
) -> bool {
    if new_perim > *ctx.opt_perim {
        return false;
    }
    if next_cfg.n_g > ctx.k {
        println!("FLOGI");
        return false;
    }
    if let Some(&idx) = ctx.d_all.get(next_cfg) {
        if ctx.dp_vals[idx].perimeter_so_far < new_perim {
            return false;
        }
    }
    true
}

#[warn(unused)]
pub fn filter_d_all_by_n_g<K: Ord>(ctx: &mut DPContext<K>, min_n_g: u32) {
    let original_count = ctx.d_all.len();
    ctx.d_all.retain(|key, _| key.n_g >= min_n_g);
    let filtered_count = ctx.d_all.len();
    println!("Filtered {} entries out of {} in the hash table.", filtered_count, original_count);
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

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct QueueItem<K: Ord> {
    key: K,
    n_g: u32,
    idx: usize,
}

pub struct DPContext<'a, K: Ord> {
    pub conf_count: &'a mut i64,
    pub conf_useless_count: &'a mut i64,
    pub d_all: &'a mut FxHashMap<DPStateKey, usize>,
    pub dp_vals: &'a mut MmapVec<DPStateValue>,
    pub pq: &'a mut BinaryHeap<QueueItem<K>>,
    pub opt_perim: &'a mut f64,
    pub best_sol: &'a mut DPStateKey,
    pub all_left_turns_cache: &'a mut FxHashMap<(Point2D, Point2D), bool>,
    pub dirs: &'a Vec<Point2D>,
    pub stops: &'a Vec<usize>,
    pub bad: &'a GridSet,
    pub bad_in_ch: &'a Vec<Point2D>,
    pub use_cache: bool,
    pub k: u32,
    pub sqrt_k: u32,
    pub good: &'a GridSet,
    pub mask: u32,
}

fn process_configuration<S: QueueStrategy>(ctx: &mut DPContext<S::Key>, cfg_idx: usize ) {
    let mut store_count = 0;
    *ctx.conf_count += 1;

    let (cfg, perimeter_so_far) = {
        let val = &mut ctx.dp_vals[cfg_idx];
        (val.cfg, val.perimeter_so_far)
    };

    if *ctx.conf_count & (ctx.mask as i64) == 0 {
        let used_bytes = ctx.dp_vals.len() * std::mem::size_of::<DPStateValue>();
        let cap_bytes = ctx.dp_vals.capacity() * std::mem::size_of::<DPStateValue>();
        let used_mb = used_bytes / 1_048_576;
        let cap_mb = cap_bytes / 1_048_576;
        println!(
            "c: {}  n_g: {}  dp_vals (u/t): {} MB / {} MB  d_all {}",
            (*ctx.conf_count).to_formatted_string(&Locale::en),
            cfg.n_g,
            used_mb.to_formatted_string(&Locale::en),
            cap_mb.to_formatted_string(&Locale::en),
            ctx.d_all.len().to_formatted_string(&Locale::en)
        );
    }
    let dir_start = cfg.dir_index;

    let dir_end = ctx.stops[dir_start as usize];
    let _prev_dir = ctx.dirs[dir_start as usize];

    for dir_i in (dir_start as usize)..=dir_end {
        let v = ctx.dirs[dir_i];

        if cfg.loc.y == 0 && v.y <= 0 {
            continue;
        }

        let p_next = cfg.loc + v;
        if ctx.bad.contains(&p_next) {
            continue;
        }

        let cache_key = (cfg.loc, p_next);
        let is_valid_turn = if ctx.use_cache {
            if let Some(&res) = ctx.all_left_turns_cache.get(&cache_key) {
                res
            } else {
                let res = is_all_left_turns(cfg.loc, p_next, ctx.bad_in_ch);
                if let Err(_) = ctx.all_left_turns_cache.try_reserve(1) {
                    let msg = "Error: Out of memory when attempting to insert into hash table";
                    println!("{}", msg);
                    eprintln!("{}", msg);
                    std::process::exit(1);
                }
                ctx.all_left_turns_cache.insert(cache_key, res);
                res
            }
        } else {
            is_all_left_turns(cfg.loc, p_next, ctx.bad_in_ch)
        };

        if !is_valid_turn {
            continue;
        }

        let (f_valid, next_cfg, new_perim) = comp_next_conf(
            ctx,
            &cfg,
            perimeter_so_far,
            dir_i as u16,
        );
        if !f_valid {
            continue;
        }

        let total_perim = new_perim + ctx.good.get_dto(next_cfg.loc).0;
        if total_perim > *ctx.opt_perim {
            continue;
        }

        let mut f_queued = false;
        let mut existing_idx = None;
        if let Some(&idx) = ctx.d_all.get(&next_cfg) {
            f_queued = true;
            existing_idx = Some(idx);
        }

        if !is_store(ctx, &next_cfg, new_perim) {
            continue;
        }

        store_count += 1;
        let next_val = DPStateValue {
            cfg: next_cfg,
            perimeter_so_far: new_perim,
            prev_idx: cfg_idx,
        };
        let push_idx;
        if let Some(idx) = existing_idx {
            ctx.dp_vals[idx] = next_val;
            push_idx = idx;
        } else {
            let idx = ctx.dp_vals.len();
            ctx.dp_vals.push(next_val).expect("push failed");
            if let Err(_) = ctx.d_all.try_reserve(1) {
                let msg = "Error: Out of memory when attempting to insert into hash table";
                println!("{}", msg);
                eprintln!("{}", msg);
                std::process::exit(1);
            }
            ctx.d_all.insert(next_cfg, idx);
            push_idx = idx;
        }

        if !f_queued {
            ctx.pq.push(QueueItem {
                key: S::compute_key(next_cfg.n_g, new_perim, push_idx, next_cfg.loc, ctx.good),
                n_g: next_cfg.n_g,
                idx: push_idx,
            });
        }

        if next_cfg.n_g == ctx.k {
            if total_perim < *ctx.opt_perim {
                *ctx.opt_perim = total_perim;
                *ctx.best_sol = next_cfg;
            }
        }
    }
    if store_count == 0 {
        *ctx.conf_useless_count += 1;
    }
}


pub fn minimize_perimeter_dp<S: QueueStrategy>(
    k: u32,
    good: &GridSet,
    bad: &GridSet,
    bad_in_ch: &Vec<Point2D>,
    use_cache: bool,
) -> anyhow::Result<(Vec<Point2D>, f64)> { //
    let sqrt_k = (k as f64).sqrt().ceil() as u32 + 1;
    // max_angle = 3 * pi / k^(1/3)
    let max_angle = 5.0 / (k as f64).powf(1.0 / 3.0);
    let ub_circle = 2.0 * (std::f64::consts::PI * (k as f64)).sqrt();
    let sq_perim = (4 * sqrt_k + 6) as f64;
    let max_edge_l: u32 = (k as f64).powf(1.0 / 3.0).round() as u32 + 1;
    let power = 19;
    let mask = (1 << power) - 1;

    println!("Size of DPStateValue: {} bytes",
        mem::size_of::<DPStateValue>() );
    println!("Size of DPStateKey: {} bytes",
        mem::size_of::<DPStateKey>() );
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
        n_g: 1
    };

    let mut best_sol = start_key;

    let mut pq = BinaryHeap::new();
    pq.push(QueueItem {
        key: S::compute_key(1, 0.0, 0, start_key.loc, good),
        n_g: 1,
        idx: 0,
    });

    let hint_size =  2*(k as usize) * (k as usize) ;
    //let hint_size_d = 2 * (k as usize) * (k as usize) + 20 * (k as usize) + 1000;
    //let hint_size_s = (k as usize) + (k as usize) * (k as usize) / 8;
    let mut threshold = (k as usize) * 100;

    println!(
        "hint_size   : {}",
        (hint_size as i64).to_formatted_string(&Locale::en)
    );

    let mut d_all = FxHashMap::default();
    //let mut d_all = FxHashMap::with_capacity( hint_size );
    let mut all_left_turns_cache: FxHashMap<(Point2D, Point2D), bool> = FxHashMap::default();
    println!( "Tring to resize d_all.. {}", 2*threshold );
    d_all.reserve( 2 * threshold );


    //let byte_size = hint_size * std::mem::size_of::<DPStateValue>();
    // 1. Open/Create the file and set its size
    
    println!( "Trying to allocate dp_vals... " );
    //let mut dp_vals = MmapVec::<DPStateValue>::new();
    let mut dp_vals = MmapVec::<DPStateValue>::with_capacity( hint_size )?;
    println!( "Path of file: {}", dp_vals.path().display() );
    //let mut dp_vals: MmapVec<> = config.create()?;

    //let mut dp_vals = Vec::with_capacity(hint_size);
    let start_val = DPStateValue {
        cfg: start_key,
        perimeter_so_far: 0.0,
        prev_idx: 0
        //handled: false,
    };
    if let Err(_) = d_all.try_reserve(1) {
        let msg = "Error: Out of memory when attempting to insert into hash table";
        println!("{}", msg);
        eprintln!("{}", msg);
        std::process::exit(1);
    }
    d_all.insert(start_key, 0);
    dp_vals.push(start_val)?;

    let mut conf_count: i64 = 0;
    let mut conf_useless_count: i64 = 0;
        let mut ctx = DPContext {
            conf_count: &mut conf_count,
            conf_useless_count: &mut conf_useless_count,
            d_all: &mut d_all,
            dp_vals: &mut dp_vals,
            pq: &mut pq,
            opt_perim: &mut opt_perim,
            best_sol: &mut best_sol,
            all_left_turns_cache: &mut all_left_turns_cache,
            dirs: &dirs,
            stops: &stops,
            bad,
            bad_in_ch,
            use_cache,
            k,
            sqrt_k,
            good,
            mask,
        };
        loop {
            let item = ctx.pq.pop();
            if item.is_none() { break; }
            let QueueItem {
                key: _,
                n_g: curr_n_g,
                idx: popped_idx,
            } = item.unwrap();

            let mut ids: Vec<usize> = Vec::new();
            ids.push( popped_idx );

            // Continuously peek at the next item

            while let Some(next_item) = ctx.pq.peek() {
                // If the next item's n_g matches our current one, pop and store it
                if next_item.n_g == curr_n_g {
                    if let Some(matched_item) = ctx.pq.pop() {
                        ids.push(matched_item.idx);
                    }
                } else {
                    // Since BinaryHeap is sorted, once we find a different n_g, 
                    // we know no other matches exist.
                    break;
                }   
            }

            //println!( "n_g: {}, batch size: {}", curr_n_g, ids.len() );
            if ctx.d_all.len() > threshold {
                let min_n_g = curr_n_g;
                filter_d_all_by_n_g(&mut ctx, min_n_g);
                threshold = max( threshold, 3 * ctx.d_all.len() / 2 );
            }
            for  id  in ids {
                process_configuration::<S>(&mut ctx, id );
            }
        }
        println!("# of configurations generated: {}", *ctx.conf_count);
        println!("# of dead-end  configurations: {}", *ctx.conf_useless_count);
    let sol = extract_solution(&ctx);
    Ok( (sol, ub_circle) )
}

// End of file
