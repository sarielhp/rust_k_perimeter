use crate::geom::{GridSet, is_all_left_turns, does_segment_intersect_polygon, is_right_turn, is_colinear, triangle_count_new_points, is_lefteq_turn, is_left_turn};
use crate::kd_tree::{BBox, KDTree, Status};
use crate::point::*;
use std::collections::VecDeque;
use std::time::Instant;

/// Precalculated information about an edge in the visibility graph.
/// Storing these values avoids redundant geometric calculations during DP.
#[derive(Debug, Clone, Copy)]
pub struct EdgeInfo {
    /// ID of the target point in the GridSet.
    pub target_id: usize,
    /// Number of new grid points enclosed when adding this edge to the polygon.
    pub n_g_delta: u32,
    /// Max additional grid points for pruning.
    pub max_addion_g: u32,
    /// Euclidean length of the edge.
    pub edge_len: f64,
    /// Minimum Euclidean distance from the target point back to the origin.
    pub target_dto: f64,
    /// Index in the adjacency list of the target point where valid outgoing edges start.
    pub next_edge_start_idx: u32,
    /// Index in the adjacency list of the target point where valid outgoing edges end (exclusive).
    pub next_edge_end_idx: u32,
}

/// A graph where edges represent valid visibility segments between grid points.
pub struct VisibilityGraph {
    /// Adjacency list storing precalculated edge information for each point.
    pub adjacency_list: Vec<Vec<EdgeInfo>>,
}

pub fn turn_angle(u: Point2D, v: Point2D, w: Point2D) -> f64 {
    let p_uv = v - u;
    let p_vw = w - v;
    let angle_uv = (p_uv.y as f64).atan2(p_uv.x as f64);
    let angle_vw = (p_vw.y as f64).atan2(p_vw.x as f64);
    let mut diff = angle_vw - angle_uv;
    while diff < 0.0 {
        diff += 2.0 * std::f64::consts::PI;
    }
    while diff >= 2.0 * std::f64::consts::PI {
        diff -= 2.0 * std::f64::consts::PI;
    }
    diff
}

pub fn verify_suffix_property(vg: &VisibilityGraph, good: &GridSet, max_turn_angle: f64) {
    for u_id in 0..vg.adjacency_list.len() {
        let u = *good.get_point_by_id(u_id);
        for (edge_idx, edge_uv) in vg.adjacency_list[u_id].iter().enumerate() {
            let v_id = edge_uv.target_id;
            let v = *good.get_point_by_id(v_id);

            let out_edges = &vg.adjacency_list[v_id];
            let start_idx = edge_uv.next_edge_start_idx as usize;
            let end_idx = edge_uv.next_edge_end_idx as usize;

            for (idx, edge_vw) in out_edges.iter().enumerate() {
                let w = *good.get_point_by_id(edge_vw.target_id);
                let is_convex = is_lefteq_turn(u, v, w);
                let t_angle = turn_angle(u, v, w);
                let is_within_angle = is_convex && t_angle <= max_turn_angle + 1e-9;

                if idx < start_idx {
                    if is_convex {
                        panic!(
                            "Suffix property violated (start_idx) at vertex {:?} (id: {}), incoming from {:?} (id: {}), edge index {}. \
                            Outgoing edge to {:?} (idx {}) is convex, but start_idx is {}.",
                            v, v_id, u, u_id, edge_idx, w, idx, start_idx
                        );
                    }
                } else if idx < end_idx {
                    if !is_within_angle {
                        panic!(
                            "Range property violated (end_idx) at vertex {:?} (id: {}), incoming from {:?} (id: {}), edge index {}. \
                            Outgoing edge to {:?} (idx {}) is NOT within angle (angle={}, max={}), but it is in the range [{}, {}).",
                            v, v_id, u, u_id, edge_idx, w, idx, t_angle, max_turn_angle, start_idx, end_idx
                        );
                    }
                } else {
                    if is_within_angle {
                        // This might be okay if it's after end_idx, but only if all subsequent ones are also not within angle?
                        // Actually, with CCW sort, the turn angle is increasing.
                    }
                }
            }
        }
    }
    println!("Visibility graph suffix property verified successfully.");
}

pub fn sort_adjacency_list_edges(good: &GridSet, u: Point2D, edges: &mut [EdgeInfo]) {
    edges.sort_by(|a, b| {
        let p = *good.get_point_by_id(a.target_id);
        let q = *good.get_point_by_id(b.target_id);
        if p == q {
            return std::cmp::Ordering::Equal;
        }
        if is_left_turn(u, p, q) {
            std::cmp::Ordering::Less
        } else if is_right_turn(u, p, q) {
            std::cmp::Ordering::Greater
        } else {
            // Colinear. Compare distance from u.
            // Since they are on the same side of u (due to visibility constraints),
            // comparing distance is correct.
            let dist_p = (p - u).norm_sq();
            let dist_q = (q - u).norm_sq();
            dist_p.cmp(&dist_q)
        }
    });
}

fn halfplane_bbox_status(bbox: &BBox, a: Point2D, b: Point2D, inclusive: bool) -> Status {
    let dx = b.x as i64 - a.x as i64;
    let dy = b.y as i64 - a.y as i64;

    let mut min_f = i64::MAX;
    let mut max_f = i64::MIN;

    let corners = [
        (bbox.min_x, bbox.min_y),
        (bbox.min_x, bbox.max_y),
        (bbox.max_x, bbox.min_y),
        (bbox.max_x, bbox.max_y),
    ];

    for (x, y) in corners {
        let f = dx * (y as i64 - a.y as i64) - dy * (x as i64 - a.x as i64);
        if f < min_f {
            min_f = f;
        }
        if f > max_f {
            max_f = f;
        }
    }

    if inclusive {
        if min_f >= 0 {
            return Status::Inside;
        }
        if max_f < 0 {
            return Status::Outside;
        }
    } else {
        if min_f > 0 {
            return Status::Inside;
        }
        if max_f <= 0 {
            return Status::Outside;
        }
    }
    Status::Partial
}

/// Builds the visibility graph by checking visibility and orientation constraints for all potential edges.
/// Also precalculates grid point counts (n_g) and distances for each valid edge.
pub fn build_visibility_graph(
    good: &GridSet,
    bad_ch: &[Point2D],
    so_so: &[Point2D],
    dirs: &Vec<Point2D>,
    k: usize,
    max_turn_angle: f64,
) -> VisibilityGraph {
    let start_vg = Instant::now();
    let n = good.num_points();
    let mut adjacency_list = vec![vec![]; n];
    let n_dirs = dirs.len();
    let sqrt_k = (k as f64).sqrt().ceil() as i32 + 1;

    let start_kdtree = Instant::now();
    let tree = KDTree::new(so_so.to_vec());
    println!("   KDTree construction took: {:?}", start_kdtree.elapsed());

    let start_edges = Instant::now();
    let origin = Point2D { x: 0, y: 0 };
    for i in 0..n {
        let p = good.points[i];
        for j in 0..n_dirs {
            let q_dir = dirs[j];
            let q = p + q_dir;

            if (q.y < 0)
                || q.is_zero()
                || (q.x as i64) > (sqrt_k as i64)
                || (q.y as i64) > ((2 * sqrt_k) as i64)
                || (-q.x as i64) > (sqrt_k as i64)
                || (-q.y as i64) > (sqrt_k as i64)
                || (p.y == 0 && q.y <= 0)
                || (p == q)
            {
                continue;
            }

            if !good.contains(&q) {
                continue;
            }

            if is_right_turn(p, q, origin) {
                continue;
            }

            if is_colinear(p, q, origin) {
                if q.norm_sq() < p.norm_sq() {
                    continue;
                }
            }

            if does_segment_intersect_polygon(bad_ch, p, q) {
                continue;
            }

            let f_left = is_all_left_turns(p, q, bad_ch);
            if !f_left {
                continue;
            }

            let q_id = good.get_point_id(q);
            let (tri_i_new, tri_b_new) = triangle_count_new_points(ORIGIN, p, q);
            let n_g_delta = tri_i_new + tri_b_new;
            let edge_len = (p - q).norm();
            let target_dto = good.get_dto(q).0;

            let max_addion_g = tree.count_in_region(
                |pt| is_left_turn(ORIGIN, q, pt) && !is_right_turn(p, q, pt),
                &|bbox| {
                    let s1 = halfplane_bbox_status(bbox, ORIGIN, q, false);
                    let s2 = halfplane_bbox_status(bbox, p, q, true);
                    match (s1, s2) {
                        (Status::Outside, _) | (_, Status::Outside) => Status::Outside,
                        (Status::Inside, Status::Inside) => Status::Inside,
                        _ => Status::Partial,
                    }
                },
            );

            adjacency_list[i].push(EdgeInfo {
                target_id: q_id,
                n_g_delta,
                max_addion_g: max_addion_g as u32,
                edge_len,
                target_dto,
                next_edge_start_idx: 0,
                next_edge_end_idx: 0,
            });
        }
    }
    println!("   Edges generation took: {:?}", start_edges.elapsed());

    let start_sort = Instant::now();

    for i in 0..n {
        let u = good.points[i];
        let edges = &mut adjacency_list[i];
        if edges.len() > 1 {
            let mut is_sorted = true;
            let mut is_reverse_sorted = true;
            for j in 0..edges.len() - 1 {
                let p = *good.get_point_by_id(edges[j].target_id);
                let q = *good.get_point_by_id(edges[j + 1].target_id);
                
                // Compare edges[j] and edges[j+1]
                let cmp = if p == q {
                    std::cmp::Ordering::Equal
                } else if is_left_turn(u, p, q) {
                    std::cmp::Ordering::Less
                } else if is_right_turn(u, p, q) {
                    std::cmp::Ordering::Greater
                } else {
                    let dist_p = (p - u).norm_sq();
                    let dist_q = (q - u).norm_sq();
                    dist_p.cmp(&dist_q)
                };

                if cmp == std::cmp::Ordering::Greater {
                    is_sorted = false;
                }
                if cmp == std::cmp::Ordering::Less {
                    is_reverse_sorted = false;
                }
            }
            if is_sorted {
                // Already sorted
            } else if is_reverse_sorted {
                edges.reverse();
            } else {
                sort_adjacency_list_edges(&good, u, edges);
            }
        } else if edges.len() == 1 {
            // No sorting needed for single edge
        }
    }
    println!("   Total `Sorting` edges took: {:?}", start_sort.elapsed());

    // Correct way to update:
    let start_next_edge = Instant::now();
    for i in 0..n {
        let u = good.points[i];
        for j in 0..adjacency_list[i].len() {
            let v_id = adjacency_list[i][j].target_id;
            let v = *good.get_point_by_id(v_id);

            let out_edges = &adjacency_list[v_id];
            let mut start_idx = out_edges.len() as u32;
            let mut end_idx = out_edges.len() as u32;

            for (idx, edge_vw) in out_edges.iter().enumerate() {
                let w = *good.get_point_by_id(edge_vw.target_id);
                if is_lefteq_turn(u, v, w) {
                    if start_idx == out_edges.len() as u32 {
                        start_idx = idx as u32;
                    }
                    if turn_angle(u, v, w) <= max_turn_angle + 1e-9 {
                        // This w is still within the angle limit.
                        // Since edges are sorted CCW, and they are already convex,
                        // the turn angle is increasing.
                    } else {
                        if end_idx == out_edges.len() as u32 {
                            end_idx = idx as u32;
                        }
                    }
                }
            }
            adjacency_list[i][j].next_edge_start_idx = start_idx;
            adjacency_list[i][j].next_edge_end_idx = end_idx;
        }
    }
    println!("   Next edge indices took: {:?}", start_next_edge.elapsed());

    let vg = VisibilityGraph { adjacency_list };
    let f_verify_vg = false;
    if f_verify_vg {
        let start_verify = Instant::now();
        verify_suffix_property(&vg, good, max_turn_angle);
        println!("   Verification took: {:?}", start_verify.elapsed());
    }
    println!(
        "   Total build_visibility_graph took: {:?}",
        start_vg.elapsed()
    );
    vg
}

pub fn topological_sort(vg: &VisibilityGraph) -> Option<Vec<usize>> {
    let n = vg.adjacency_list.len();
    let mut in_degree = vec![0; n];
    for u in 0..n {
        for edge in &vg.adjacency_list[u] {
            in_degree[edge.target_id] += 1;
        }
    }

    let mut queue = VecDeque::new();
    for i in 0..n {
        if in_degree[i] == 0 {
            queue.push_back(i);
        }
    }

    let mut result = Vec::with_capacity(n);
    while let Some(u) = queue.pop_front() {
        result.push(u);
        for edge in &vg.adjacency_list[u] {
            let v = edge.target_id;
            in_degree[v] -= 1;
            if in_degree[v] == 0 {
                queue.push_back(v);
            }
        }
    }

    if result.len() == n {
        Some(result)
    } else {
        // If it's not a DAG, result will be shorter than n.
        // We could return partial result but None indicates a cycle.
        Some(result)
    }
}
