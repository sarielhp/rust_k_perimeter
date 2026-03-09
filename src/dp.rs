use crate::geom::{
    is_all_left_turns, is_colinear, is_lefteq_turn, is_right_turn, triangle_count_new_points,
    GridSet, Point2D,
};
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DPStateKey {
    pub loc_prev: Point2D,
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
            .then(self.loc_prev.cmp(&other.loc_prev))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DPStateValue {
    pub perimeter_so_far: f64,
    pub prev_key: DPStateKey,
    pub dir_index: u32,
    pub handled: bool,
}

pub fn comp_next_conf(
    cfg: &DPStateKey,
    perimeter_so_far: f64,
    v: Point2D,
    k: u32,
    sqrt_k: u32,
    opt_perim: f64,
) -> (bool, DPStateKey, f64) {
    let a = cfg.loc_prev;
    let b = cfg.loc;
    let c = cfg.loc + v;

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

    assert!(b != c);
    assert!(v.x != 0 || v.y != 0);

    let (tri_i_new, tri_b_new) = triangle_count_new_points(origin, b, c);

    //assert!(tri_i_new >= 0 && tri_b_new >= 0);
    let n_g = cfg.n_g + tri_i_new + tri_b_new;

    if n_g > 0 && c.is_zero() {
        return (false, *cfg, -1.0);
    }

    if n_g > k {
        return (false, *cfg, -1.0);
    }

    let new_cfg = DPStateKey {
        loc_prev: b,
        loc: c,
        n_g,
    };

    let new_perim = perimeter_so_far + (b - c).norm();
    (true, new_cfg, new_perim)
}

pub fn extract_solution(
    d_all: &FxHashMap<DPStateKey, DPStateValue>,
    best_sol: DPStateKey,
) -> Vec<Point2D> {
    let mut out = Vec::new();
    let mut curr = best_sol;

    loop {
        let val = d_all.get(&curr).unwrap();
        out.push(curr.loc);
        if val.prev_key == curr {
            break;
        }
        curr = val.prev_key;
    }

    // The Julia code pushes in reverse chronological order (end -> start)
    // We reverse it to match counter-clockwise tracing from origin if needed.
    // Wait, julia implementation extracts backwards but does not reverse:
    // `extract_solution` returns `out` directly.
    out
}

pub fn is_store(
    d_all: &FxHashMap<DPStateKey, DPStateValue>,
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
    if let Some(val) = d_all.get(next_cfg) {
        if val.perimeter_so_far < new_perim {
            return false;
        }
    }
    true
}

#[derive(PartialEq)]
struct QueueItem {
    n_g: u32,
    cfg: DPStateKey,
}

impl Eq for QueueItem {}

impl PartialOrd for QueueItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for QueueItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap is a max-heap, so we reverse the ordering
        // based on n_g and then cfg to emulate a min-heap
        other
            .n_g
            .cmp(&self.n_g)
            .then_with(|| other.cfg.cmp(&self.cfg))
    }
}

pub fn minimize_perimeter_dp(
    k: u32,
    _good: &GridSet,
    bad: &GridSet,
    bad_in_ch: &Vec<Point2D>,
) -> (Vec<Point2D>, f64) {
    let sqrt_k = (k as f64).sqrt().ceil() as u32 + 1;
    // max_angle = 3 * pi / k^(1/3)
    let max_angle = 5.0 / (k as f64).powf(1.0 / 3.0);
    let ub_circle = 2.0 * (std::f64::consts::PI * (k as f64)).sqrt();
    let sq_perim = (4 * sqrt_k + 6) as f64;
    let max_edge_l: u32 = (k as f64).powf(1.0 / 3.0).round() as u32 + 1;

    println!("sqrt_k      : {}", sqrt_k);
    println!("max_edge_l  : {}", max_edge_l);
    println!("max angle   : {}", max_angle);

    let mut opt_perim = ub_circle.min(sq_perim);

    let start_key = DPStateKey {
        loc_prev: Point2D::new(0, 0),
        loc: Point2D::new(0, 0),
        n_g: 1,
    };

    let mut best_sol = start_key;

    let mut pq = BinaryHeap::new();
    pq.push(QueueItem {
        n_g: 1,
        cfg: start_key,
    });

    let mut d_all = FxHashMap::default();
    let start_val = DPStateValue {
        perimeter_so_far: 0.0,
        prev_key: start_key,
        dir_index: 0,
        handled: false,
    };
    d_all.insert(start_key, start_val);

    let v_vec = crate::geom::generate_primitive_vectors(max_edge_l);
    let stops = crate::geom::comp_stop_indexes(&v_vec, max_angle);

    while let Some(QueueItem { n_g: _, cfg }) = pq.pop() {
        let (dir_start, perimeter_so_far) = {
            let val = d_all.get_mut(&cfg).unwrap();
            if val.handled {
                continue;
            }
            val.handled = true;
            (val.dir_index, val.perimeter_so_far)
        };

        let dir_end = stops[dir_start as usize];
        let _prev_dir = v_vec[dir_start as usize];

        for dir_i in (dir_start as usize)..=dir_end {
            let v = v_vec[dir_i];

            if cfg.loc.y == 0 && v.y <= 0 {
                continue;
            }

            let p_next = cfg.loc + v;
            if bad.contains(&p_next) {
                continue;
            }

            if !is_all_left_turns(cfg.loc, p_next, bad_in_ch) {
                //                println!("BANGO!");
                continue;
            }

            let (f_valid, next_cfg, new_perim) =
                comp_next_conf(&cfg, perimeter_so_far, v, k, sqrt_k, opt_perim);
            if !f_valid {
                continue;
            }

            if !is_store(&d_all, &next_cfg, new_perim, opt_perim, k) {
                continue;
            }

            d_all.insert(
                next_cfg,
                DPStateValue {
                    perimeter_so_far: new_perim,
                    prev_key: cfg,
                    dir_index: (dir_i as u32),
                    handled: false,
                },
            );
            pq.push(QueueItem {
                n_g: next_cfg.n_g,
                cfg: next_cfg,
            });

            if next_cfg.n_g == k {
                let total_perim = new_perim + next_cfg.loc.norm();
                if total_perim < opt_perim {
                    opt_perim = total_perim;
                    best_sol = next_cfg;
                }
            }
        }
    }

    let sol = extract_solution(&d_all, best_sol);
    (sol, ub_circle)
}
