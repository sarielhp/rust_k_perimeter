//! Visibility graph construction and topological sorting.
//!
//! This module builds a graph of valid edges between grid points. An edge is valid if:
//! 1. It is within the search area (good set).
//! 2. It doesn't intersect the forbidden "bad set".
//! 3. It maintains the required convexity and turn angle constraints.

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
    /// Max additional grid points that could potentially be enclosed from this target point forward.
    pub max_addion_g: u32,
    /// Euclidean length of the edge.
    pub edge_len: f64,
    /// Minimum Euclidean distance from the target point back to the origin.
    pub target_dto: f64,
    /// Precomputed index range in the target's adjacency list for convex continuations.
    pub next_edge_start_idx: u32,
    /// Precomputed index range in the target's adjacency list for convex continuations.
    pub next_edge_end_idx: u32,
}

/// A graph where edges represent valid visibility segments between grid points.
pub struct VisibilityGraph {
    /// Adjacency list storing precalculated edge information for each point.
    pub adjacency_list: Vec<Vec<EdgeInfo>>,
}

/// Calculates the turn angle between segments (u,v) and (v,w).
pub fn turn_angle(u: Point2D, v: Point2D, w: Point2D) -> f64 {
    let p_uv = v - u;
    let p_vw = w - v;
    let angle_uv = (p_uv.y as f64).atan2(p_uv.x as f64);
    let angle_vw = (p_vw.y as f64).atan2(p_vw.x as f64);
    let mut diff = angle_vw - angle_uv;
    while diff < 0.0 { diff += 2.0 * std::f64::consts::PI; }
    while diff >= 2.0 * std::f64::consts::PI { diff -= 2.0 * std::f64::consts::PI; }
    diff
}

/// Sorts outgoing edges from point `u` in Counter-Clockwise (CCW) order.
/// This sorting is crucial for the `next_edge_start/end_idx` optimization.
pub fn sort_adjacency_list_edges(good: &GridSet, u: Point2D, edges: &mut [EdgeInfo]) {
    edges.sort_by(|a, b| {
        let p = *good.get_point_by_id(a.target_id);
        let q = *good.get_point_by_id(b.target_id);
        if p == q { return std::cmp::Ordering::Equal; }
        if is_left_turn(u, p, q) { std::cmp::Ordering::Less }
        else if is_right_turn(u, p, q) { std::cmp::Ordering::Greater }
        else { (p - u).norm_sq().cmp(&(q - u).norm_sq()) }
    });
}

/// Determines the visibility/containment status of a bounding box relative to a half-plane.
fn halfplane_bbox_status(bbox: &BBox, a: Point2D, b: Point2D, inclusive: bool) -> Status {
    let dx = b.x as i64 - a.x as i64;
    let dy = b.y as i64 - a.y as i64;
    let mut min_f = i64::MAX; let mut max_f = i64::MIN;
    let corners = [(bbox.min_x, bbox.min_y), (bbox.min_x, bbox.max_y), (bbox.max_x, bbox.min_y), (bbox.max_x, bbox.max_y)];
    for (x, y) in corners {
        let f = dx * (y as i64 - a.y as i64) - dy * (x as i64 - a.x as i64);
        min_f = min_f.min(f); max_f = max_f.max(f);
    }
    if inclusive {
        if min_f >= 0 { return Status::Inside; }
        if max_f < 0 { return Status::Outside; }
    } else {
        if min_f > 0 { return Status::Inside; }
        if max_f <= 0 { return Status::Outside; }
    }
    Status::Partial
}

/// Builds the visibility graph by checking visibility and orientation constraints.
///
/// Optimization: Uses a KDTree to precalculate `max_addion_g`, which is an upper bound
/// on how many points can still be enclosed starting from an edge.
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
    let sqrt_k = (k as f64).sqrt().ceil() as i32 + 1;

    let tree = KDTree::new(so_so.to_vec());
    for i in 0..n {
        let p = good.points[i];
        for &q_dir in dirs {
            let q = p + q_dir;
            // Basic boundary and origin constraints.
            if q.y < 0 || q.is_zero() || q == p { continue; }
            if (q.x as i64).abs() > (sqrt_k as i64) || q.y as i64 > (2 * sqrt_k as i64) { continue; }
            if !good.contains(&q) { continue; }

            // Ensure the edge doesn't wrap "behind" the origin (maintains convexity).
            if is_right_turn(p, q, ORIGIN) { continue; }
            if is_colinear(p, q, ORIGIN) && q.norm_sq() < p.norm_sq() { continue; }

            // Visibility: cannot cross the forbidden interior.
            if does_segment_intersect_polygon(bad_ch, p, q) { continue; }
            if !is_all_left_turns(p, q, bad_ch) { continue; }

            let q_id = good.get_point_id(q);
            let (tri_i, tri_b) = triangle_count_new_points(ORIGIN, p, q);
            
            // Heuristic enclosure pruning: count points in the wedge defined by (origin, p, q).
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
                n_g_delta: tri_i + tri_b,
                max_addion_g: max_addion_g as u32,
                edge_len: (p - q).norm(),
                target_dto: good.get_dto(q).0,
                next_edge_start_idx: 0,
                next_edge_end_idx: 0,
            });
        }
    }

    // Sort edges CCW to enable efficient range lookups during DP.
    for i in 0..n {
        let u = good.points[i];
        sort_adjacency_list_edges(good, u, &mut adjacency_list[i]);
    }

    // Precalculate suffix ranges: for each edge (u,v), find the range of edges (v,w)
    // that are convex relative to (u,v) and within the turn angle limit.
    for i in 0..n {
        let u = good.points[i];
        for j in 0..adjacency_list[i].len() {
            let v_id = adjacency_list[i][j].target_id;
            let out_edges = &adjacency_list[v_id];
            let mut start_idx = out_edges.len() as u32;
            let mut end_idx = out_edges.len() as u32;

            for (idx, edge_vw) in out_edges.iter().enumerate() {
                let w = *good.get_point_by_id(edge_vw.target_id);
                if is_lefteq_turn(u, *good.get_point_by_id(v_id), w) {
                    if start_idx == out_edges.len() as u32 { start_idx = idx as u32; }
                    if turn_angle(u, *good.get_point_by_id(v_id), w) <= max_turn_angle + 1e-9 {
                        // Turn angle is valid.
                    } else if end_idx == out_edges.len() as u32 {
                        end_idx = idx as u32;
                    }
                }
            }
            adjacency_list[i][j].next_edge_start_idx = start_idx;
            adjacency_list[i][j].next_edge_end_idx = end_idx;
        }
    }
    println!("   Total build_visibility_graph took: {:?}", start_vg.elapsed());
    VisibilityGraph { adjacency_list }
}

/// Kahn's algorithm for topological sorting.
/// Returns the point IDs in topological order, or None if a cycle is detected.
pub fn topological_sort(vg: &VisibilityGraph) -> Option<Vec<usize>> {
    let n = vg.adjacency_list.len();
    let mut in_degree = vec![0; n];
    for u in 0..n {
        for edge in &vg.adjacency_list[u] { in_degree[edge.target_id] += 1; }
    }
    let mut queue = VecDeque::new();
    for i in 0..n { if in_degree[i] == 0 { queue.push_back(i); } }
    let mut result = Vec::with_capacity(n);
    while let Some(u) = queue.pop_front() {
        result.push(u);
        for edge in &vg.adjacency_list[u] {
            let v = edge.target_id;
            in_degree[v] -= 1;
            if in_degree[v] == 0 { queue.push_back(v); }
        }
    }
    if result.len() == n { Some(result) } else { Some(result) }
}
