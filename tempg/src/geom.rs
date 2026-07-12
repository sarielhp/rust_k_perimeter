//! Geometric primitives and algorithms for the k-perimeter problem.
//!
//! This module provides the foundational tools for working with 2D grid points,
//! calculating areas using Pick's Theorem and the Shoelace formula, computing
//! convex hulls, and managing the "good set" of points allowed for DP exploration.

use crate::point::*;
use num::integer::gcd;
use std::cmp::{max, min};

/// Translates a set of points by a given vector `v`.
pub fn vtrans(p: &[Point2D], v: Point2D) -> Vec<Point2D> {
    p.iter().map(|&pt| pt + v).collect()
}

/// Checks if three points (a, b, c) form a left turn or are collinear.
/// Calculated using the 2D cross product.
#[allow(dead_code)]
pub fn is_lefteq_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x as i64 - a.x as i64) * (c.y as i64 - a.y as i64)
        >= (b.y as i64 - a.y as i64) * (c.x as i64 - a.x as i64)
}

/// Euclidean distance between two points.
pub fn d_y(a: Point2D, b: Point2D) -> f64 {
    (a - b).norm()
}

/// Checks if three points (a, b, c) form a strict left turn (counter-clockwise).
pub fn is_left_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x as i64 - a.x as i64) * (c.y as i64 - a.y as i64)
        > (b.y as i64 - a.y as i64) * (c.x as i64 - a.x as i64)
}

/// Checks if three points (a, b, c) form a strict right turn (clockwise).
pub fn is_right_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x as i64 - a.x as i64) * (c.y as i64 - a.y as i64)
        < (b.y as i64 - a.y as i64) * (c.x as i64 - a.x as i64)
}

/// Checks if three points (a, b, c) form a right turn or are collinear.
#[allow(dead_code)]
pub fn is_righteq_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x as i64 - a.x as i64) * (c.y as i64 - a.y as i64)
        <= (b.y as i64 - a.y as i64) * (c.x as i64 - a.x as i64)
}

/// Calculates the shortest distance from point `p` to the line segment `ab`.
/// Projects `p` onto the line and clamps the projection to the segment's interval [0, 1].
pub fn distance_to_segment(p: Point2D, a: Point2D, b: Point2D) -> f64 {
    let ab = b - a;
    let ap = p - a;
    let ab_len_sq = dot(&ab, &ab);

    if ab_len_sq == 0 {
        return (p - a).norm();
    }

    let t = dot(&ap, &ab) as f64 / ab_len_sq as f64;
    let t_clamped = t.clamp(0.0, 1.0);

    let a_f = Point2DF::new(a.x as f64, a.y as f64);
    let ab_f = Point2DF::new(ab.x as f64, ab.y as f64);
    let p_f = Point2DF::new(p.x as f64, p.y as f64);

    let closest_point = a_f + ab_f * t_clamped;
    (p_f - closest_point).norm()
}

/// Calculates the total Euclidean perimeter of a closed polygon defined by a sequence of points.
pub fn euclidean_length(sol: &[Point2D]) -> f64 {
    if sol.is_empty() {
        return 0.0;
    }
    let mut l = (sol.first().unwrap().clone() - sol.last().unwrap().clone()).norm();
    for i in 1..sol.len() {
        l += (sol[i - 1] - sol[i]).norm();
    }
    l
}

/// Finds the minimum distance from `query_point` to any point on the boundary of the `poly`.
pub fn polygon_boundary_distance(poly: &[Point2D], query_point: Point2D) -> f64 {
    let n = poly.len();
    if n < 2 {
        eprintln!("Error: Polygon must have at least 2 vertices");
        std::process::exit(1);
    }

    let mut min_dist = f64::INFINITY;
    for i in 0..n {
        let p1 = poly[i];
        let p2 = poly[(i + 1) % n];
        let d = distance_to_segment(query_point, p1, p2);
        if d < min_dist {
            min_dist = d;
        }
    }
    min_dist
}

/// Ray-casting algorithm to determine if a point `p` is inside the polygon `poly`.
pub fn is_point_in_polygon(poly: &[Point2D], p: Point2D) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;

    for i in 0..n {
        let v1 = poly[i];
        let v2 = poly[j];

        if ((v1.y > p.y) != (v2.y > p.y))
            && ((p.x as f64)
                < ((v2.x - v1.x) as f64) * ((p.y - v1.y) as f64) / ((v2.y - v1.y) as f64)
                    + (v1.x as f64))
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Recursively generates Farey-sequence-like primitive vectors between `u` and `v` up to `max_d`.
pub fn fary(vec: &mut Vec<Point2D>, u: Point2D, v: Point2D, max_d: i64) {
    if v.x as i64 > max_d {
        return;
    }
    let mid = Point2D::new(u.x + v.x, u.y + v.y);
    if mid.x as i64 > max_d {
        return;
    }
    fary(vec, u, mid, max_d);
    vec.push(mid);
    fary(vec, mid, v, max_d);
}

/// Twice the signed area of a triangle (a, b, c).
pub fn double_triangle_area(a: Point2D, b: Point2D, c: Point2D) -> i64 {
    (a.x as i64 * (b.y as i64 - c.y as i64)
        + (b.x as i64) * (c.y as i64 - a.y as i64)
        + (c.x as i64) * (a.y as i64 - b.y as i64))
        .abs()
}

/// Signed area of a triangle (a, b, c).
pub fn triangle_area(a: Point2D, b: Point2D, c: Point2D) -> f64 {
    double_triangle_area(a, b, c) as f64 / 2.0
}

/// Returns the number of grid points strictly on the interior of segment (u, v).
/// Calculated as gcd(|dx|, |dy|) - 1.
pub fn grid_points_inside_edge(u: Point2D, v: Point2D) -> u32 {
    if u == v {
        return 0;
    }
    let dx = (u.x as i64 - v.x as i64).abs();
    let dy = (u.y as i64 - v.y as i64).abs();
    let t = gcd(dx, dy);
    (t - 1) as u32
}

/// Total number of grid points on the boundary of a polygon.
pub fn boundary_grid_points(poly: &[Point2D]) -> u32 {
    let n = poly.len();
    if n < 3 {
        return n as u32;
    }
    let mut total = n as u32;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        total += grid_points_inside_edge(a, b);
    }
    total
}

/// Andrew's Monotone Chain algorithm for computing the convex hull of a set of points.
pub fn convex_hull(points: &[Point2D]) -> Vec<Point2D> {
    let n = points.len();
    if n <= 2 {
        return points.to_vec();
    }

    let mut sorted_points = points.to_vec();
    sorted_points.sort_by(|a, b| a.x.cmp(&b.x).then(a.y.cmp(&b.y)));

    let mut lower = Vec::new();
    for &p in &sorted_points {
        while lower.len() >= 2 && !is_left_turn(lower[lower.len() - 2], lower[lower.len() - 1], p) {
            lower.pop();
        }
        lower.push(p);
    }

    let mut upper = Vec::new();
    for &p in sorted_points.iter().rev() {
        while upper.len() >= 2 && !is_left_turn(upper[upper.len() - 2], upper[upper.len() - 1], p) {
            upper.pop();
        }
        upper.push(p);
    }

    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

fn get_y_bounds(points: &[Point2D]) -> Option<(CoordType, CoordType)> {
    if points.is_empty() {
        return None;
    }
    let initial_y = points[0].y;
    let bounds = points
        .iter()
        .skip(1)
        .fold((initial_y, initial_y), |(min, max), p| {
            (min.min(p.y), max.max(p.y))
        });
    Some(bounds)
}

fn add_points(points: &mut Vec<Point2D>, y_old: CoordType, y_new: CoordType) -> &mut Vec<Point2D> {
    let mut new_points: Vec<Point2D> = points
        .iter()
        .filter(|p| p.y == y_old)
        .map(|p| Point2D { x: p.x, y: y_new })
        .collect();
    points.append(&mut new_points);
    points.sort();
    points.dedup();
    points
}

fn average_of_min_y(points: &[Point2D]) -> Option<Point2D> {
    if points.is_empty() {
        return None;
    }
    let min_y = points.iter().map(|p| p.y).min()?;
    let mut sum_x: i64 = 0;
    let mut sum_y: i64 = 0;
    let mut count: i64 = 0;
    for p in points.iter().filter(|p| p.y == min_y) {
        sum_x += p.x as i64;
        sum_y += p.y as i64;
        count += 1;
    }
    let avg_x = (sum_x / count) as CoordType;
    let avg_y = ((sum_y + count) / count) as CoordType;
    Some(Point2D { x: avg_x, y: avg_y })
}

fn translate_v_points(points: Vec<Point2D>, avg: Point2D) -> Vec<Point2D> {
    points
        .iter()
        .map(|p| Point2D {
            x: p.x - avg.x,
            y: p.y - avg.y,
        })
        .collect()
}

/// Checks if point `b` lies on the line segment `ac`.
pub fn is_point_on_segment(a: Point2D, b: Point2D, c: Point2D) -> bool {
    let cross = (b.x as i64 - a.x as i64) * (c.y as i64 - a.y as i64)
        - (b.y as i64 - a.y as i64) * (c.x as i64 - a.x as i64);
    if cross != 0 {
        return false;
    }
    if a == c {
        return b == a;
    }
    let min_x = a.x.min(c.x);
    let max_x = a.x.max(c.x);
    let min_y = a.y.min(c.y);
    let max_y = a.y.max(c.y);
    b.x >= min_x && b.x <= max_x && b.y >= min_y && b.y <= max_y
}

/// Simplifies a polygon by removing vertices that are collinear with their adjacent neighbors.
pub fn polygon_rm_redundant_vertices(poly: &[Point2D]) -> Vec<Point2D> {
    let n = poly.len();
    if n <= 2 {
        return poly.to_vec();
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let prev = poly[(i + n - 1) % n];
        let cur = poly[i];
        let next = poly[(i + 1) % n];
        if !is_point_on_segment(prev, cur, next) {
            out.push(cur);
        }
    }
    if out.is_empty() && !poly.is_empty() {
        out.push(poly[0]);
    }
    out
}

/// Approximates a disk of ~k grid points around origin.
pub fn ch_disk_origin(k: usize, f_expand: bool) -> Vec<Point2D> {
    let r = ((k as f64) / std::f64::consts::PI).sqrt() + 2.0;
    let mut v = Vec::new();
    let l_bound = r.ceil() as i32 + 1;
    for x in -l_bound..=l_bound {
        for y in -l_bound..=l_bound {
            let p = Point2D::new(x as CoordType, y as CoordType);
            if p.norm() <= r {
                v.push(p);
            }
        }
    }
    assert!(v.len() > k);
    v.sort_by(|a, b| a.norm_sq().partial_cmp(&b.norm_sq()).unwrap());
    v.truncate(k);

    let Some((min_y, max_y)) = get_y_bounds(&v) else {
        std::process::exit(1);
    };
    if f_expand {
        add_points(&mut v, max_y - 1, max_y);
        add_points(&mut v, min_y + 1, min_y);
    }
    v.sort_by(|a, b| a.x.cmp(&b.x).then(a.y.cmp(&b.y)));

    let mut out = Vec::new();
    let len = v.len();
    for i in 0..len {
        if i == 0 || i == len - 1 {
            out.push(v[i]);
            continue;
        }
        if v[i].x == v[i - 1].x && v[i].x == v[i + 1].x {
            continue;
        }
        out.push(v[i]);
    }

    let ch = convex_hull(&out);
    let min_y = ch.iter().map(|p| p.y).min().unwrap();
    let max_x = ch
        .iter()
        .filter(|p| p.y == min_y)
        .map(|p| p.x)
        .max()
        .unwrap();
    let Some(avg) = average_of_min_y(&ch) else {
        std::process::exit(1);
    };

    let l = ((max_x - avg.x) / 2).abs();
    let mut mv = Point2D {
        x: max_x - 1,
        y: min_y,
    };
    if f_expand {
        mv = Point2D {
            x: max_x - l,
            y: min_y,
        };
    }
    translate_v_points(ch, mv)
}

/// Generates all primitive integer vectors (dx, dy) with length <= `max_d`.
pub fn generate_primitive_vectors(max_d: u32) -> Vec<Point2D> {
    let mut vec = Vec::new();
    let mut tvec = Vec::new();
    fary(
        &mut tvec,
        Point2D::new(1, 0),
        Point2D::new(1, 1),
        max_d as i64,
    );
    vec.push(Point2D::new(1, 0));
    vec.extend(&tvec);
    vec.push(Point2D::new(1, 1));
    for v in tvec.iter().rev() {
        vec.push(Point2D::new(v.y, v.x));
    }
    let len = vec.len();
    for i in 0..len {
        let p = vec[i];
        vec.push(Point2D::new(-p.y, p.x));
    }
    let len2 = vec.len();
    for i in 0..len2 {
        let p = vec[i];
        vec.push(Point2D::new(-p.x, -p.y));
    }
    vec
}

/// Checks if three points are collinear.
pub fn is_colinear(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x - a.x) * (c.y - a.y) == (b.y - a.y) * (c.x - a.x)
}

fn count_distinct(a: Point2D, b: Point2D, c: Point2D) -> u64 {
    if a == b && b == c {
        1
    } else if a == b || a == c || b == c {
        2
    } else {
        3
    }
}

/// Counts new internal and boundary grid points added by forming a triangle (origin, a, b).
/// Uses Pick's Theorem ($A = I + B/2 - 1$) to solve for $I$.
pub fn triangle_count_new_points(a: Point2D, b: Point2D, c: Point2D) -> (u32, u32) {
    let area2 = double_triangle_area(a, b, c);
    let ab_g_n = grid_points_inside_edge(a, b);
    let bc_g_n = grid_points_inside_edge(b, c);
    let ac_g_n = grid_points_inside_edge(a, c);
    let boundary_n: u64 =
        ab_g_n as u64 + bc_g_n as u64 + ac_g_n as u64 + count_distinct(a, b, c) as u64;

    let tri_i_new = if area2 > 0 {
        //let area_i64 = area2 as i64;
        //let boundary_i64 = boundary_n as i64;
        let result = (area2 - boundary_n as i64 + 2) / 2;
        result as u32
    } else {
        0
    };

    let tri_b_new: i64 = if a == c {
        std::process::exit(1);
    } else if a == b {
        ac_g_n as i64 + 1
    } else if area2 == 0 {
        bc_g_n as i64 + 1
    } else {
        ac_g_n as i64 + bc_g_n as i64 + 1
    };

    (tri_i_new, tri_b_new as u32)
}

/// Calculates the bounding box of a collection of polygons.
pub fn bound(polys: &[&[Point2D]], expand: i32) -> (CoordType, CoordType, CoordType, CoordType) {
    let mut min_x = 0;
    let mut max_x = 0;
    let mut min_y = 0;
    let mut max_y = 0;
    let mut f_init = false;
    for &poly in polys {
        if poly.is_empty() {
            continue;
        }
        let px_min = poly.iter().map(|p| p.x).min().unwrap() - 1;
        let px_max = poly.iter().map(|p| p.x).max().unwrap() + 1;
        let py_min = poly.iter().map(|p| p.y).min().unwrap() - 1;
        let py_max = poly.iter().map(|p| p.y).max().unwrap() + 1;
        if !f_init {
            min_x = px_min;
            max_x = px_max;
            min_y = py_min;
            max_y = py_max;
            f_init = true;
        } else {
            min_x = min_x.min(px_min);
            max_x = max_x.max(px_max);
            min_y = min_y.min(py_min);
            max_y = max_y.max(py_max);
        }
    }
    (
        min_x - expand as CoordType,
        max_x + expand as CoordType,
        min_y - expand as CoordType,
        max_y + expand as CoordType,
    )
}

/// A structure for managing a grid of points and associated properties.
#[derive(Debug)]
pub struct GridSet {
    pub min_x: CoordType,
    pub max_x: CoordType,
    pub min_y: CoordType,
    pub max_y: CoordType,
    _size: usize,
    width: usize,
    data: Vec<bool>,
    dto: Vec<f64>,
    point_id: Vec<usize>,
    pub points: Vec<Point2D>,
    dto_g: Vec<i64>,
    topo_idx: Vec<u32>,
}

impl GridSet {
    pub fn num_points(&self) -> usize {
        self.points.len()
    }
    pub fn new(min_x: CoordType, max_x: CoordType, min_y: CoordType, max_y: CoordType) -> Self {
        let width = (max_x - min_x + 1).max(0) as usize;
        let height = (max_y - min_y + 1).max(0) as usize;
        let _size = width * height;
        Self {
            min_x,
            max_x,
            min_y,
            max_y,
            _size,
            width,
            points: Vec::new(),
            data: vec![false; _size],
            dto: vec![-1.0; _size],
            point_id: vec![usize::MAX; _size],
            dto_g: vec![0; _size],
            topo_idx: Vec::new(),
        }
    }
    pub fn get_topo_idx(&self, id: usize) -> u32 {
        self.topo_idx[id]
    }
    pub fn set_topo_idx(&mut self, id: usize, idx: u32) {
        if self.topo_idx.is_empty() {
            self.topo_idx = vec![0; self.points.len()];
        }
        self.topo_idx[id] = idx;
    }
    #[allow(dead_code)]
    pub fn length(&self) -> usize {
        self.data.iter().filter(|&&x| x).count()
    }
    pub fn insert(&mut self, p: Point2D) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            if self.data[idx] {
                return;
            }
            self.data[idx] = true;
        }
    }
    pub fn delete(&mut self, p: Point2D) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.data[idx] = false;
        }
    }
    fn get_index(&self, p: Point2D) -> usize {
        (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize
    }
    pub fn insert_val(&mut self, p: Point2D, val: f64, g: i64) {
        let idx = self.get_index(p);
        self.dto[idx] = val;
        self.dto_g[idx] = g;
    }
    pub fn get_dto(&self, p: Point2D) -> (f64, i64) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = self.get_index(p);
            (self.dto[idx], self.dto_g[idx])
        } else {
            (self.min_x as f64 * self.min_x as f64, 0)
        }
    }
    pub fn contains(&self, p: &Point2D) -> bool {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            self.data[self.get_index(*p)]
        } else {
            false
        }
    }
    pub fn compute_points(&mut self) {
        self.points.clear();
        self.topo_idx.clear();
        for y in self.min_y..=self.max_y {
            for x in self.min_x..=self.max_x {
                let p = Point2D::new(x, y);
                if self.contains(&p) {
                    let id = self.points.len();
                    self.points.push(p);
                    self.topo_idx.push(0);
                    let index = self.get_index(p);
                    self.point_id[index] = id;
                }
            }
        }
    }
    pub fn get_point_by_id(&self, id: usize) -> &Point2D {
        &self.points[id]
    }
    pub fn get_point_id(&self, p: Point2D) -> usize {
        self.point_id[self.get_index(p)]
    }
    pub fn fill_dist_to_origin(&mut self, bad_ch: &[Point2D]) {
        for y in self.min_y..=self.max_y {
            if y < 0 {
                continue;
            }
            for x in self.min_x..=self.max_x {
                let p = Point2D::new(x, y);
                if self.contains(&p) {
                    let (l, g) = distance_to_origin(bad_ch, p);
                    self.insert_val(p, l, g);
                }
            }
        }
    }
}

fn segments_intersect(p1: Point2D, q1: Point2D, p2: Point2D, q2: Point2D) -> bool {
    fn orientation(p: Point2D, q: Point2D, r: Point2D) -> i32 {
        let val = (q.y as i64 - p.y as i64) * (r.x as i64 - q.x as i64)
            - (q.x as i64 - p.x as i64) * (r.y as i64 - q.y as i64);
        if val == 0 {
            0
        } else if val > 0 {
            1
        } else {
            2
        }
    }
    fn on_segment(p: Point2D, q: Point2D, r: Point2D) -> bool {
        q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
    }
    let o1 = orientation(p1, q1, p2);
    let o2 = orientation(p1, q1, q2);
    let o3 = orientation(p2, q2, p1);
    let o4 = orientation(p2, q2, q1);
    if o1 != o2 && o3 != o4 {
        return true;
    }
    if o1 == 0 && on_segment(p1, p2, q1) {
        return true;
    }
    if o2 == 0 && on_segment(p1, q2, q1) {
        return true;
    }
    if o3 == 0 && on_segment(p2, p1, q2) {
        return true;
    }
    if o4 == 0 && on_segment(p2, q1, q2) {
        return true;
    }
    false
}

pub fn does_segment_intersect_polygon(polygon: &[Point2D], s1: Point2D, s2: Point2D) -> bool {
    let n = polygon.len();
    for i in 0..n {
        if segments_intersect(s1, s2, polygon[i], polygon[(i + 1) % n]) {
            return true;
        }
    }
    false
}

pub fn is_all_left_turns(p: Point2D, p_next: Point2D, v: &[Point2D]) -> bool {
    let ax = p_next.x - p.x;
    let ay = p_next.y - p.y;
    for q in v {
        if (ax * (q.y - p_next.y) - ay * (q.x - p_next.x)) <= 0 {
            return false;
        }
    }
    true
}

/// Calculates the shortest distance and enclosed points from origin to point `p`, respecting obstacles.
pub fn distance_to_origin(bad_ch: &[Point2D], p: Point2D) -> (f64, i64) {
    let origin = Point2D { x: 0, y: 0 };
    if !does_segment_intersect_polygon(bad_ch, origin, p) {
        return (d_y(origin, p), 0);
    }
    let mut pnts = bad_ch.to_vec();
    pnts.push(origin);
    pnts.push(p);
    let poly = convex_hull(&pnts);
    let i_o = poly.iter().position(|&x| x == origin).unwrap();
    let i_p = poly.iter().position(|&x| x == p).unwrap();
    let (s, t) = (min(i_o, i_p), max(i_o, i_p));
    let mut t_l = 0.0;
    let mut t_b = 1;
    let mut area = 0.0;
    for k in s..t {
        t_l += d_y(poly[k], poly[k + 1]);
        t_b += grid_points_inside_edge(poly[k], poly[k + 1]) + 1;
        area += triangle_area(origin, poly[k], poly[k + 1]);
    }
    let b_edge = grid_points_inside_edge(poly[s], poly[t]);
    t_b += b_edge;
    let i_verts = (area - t_b as f64 / 2.0 + 1.0) as i64;
    let g_overall: i64 = i_verts + t_b as i64 - b_edge as i64 - 2;
    if i_p < i_o {
        return (t_l, g_overall);
    }
    let total_area = polygon_area(&poly);
    let mut total_b = 0;
    for k in 0..poly.len() {
        total_b += grid_points_inside_edge(poly[k], poly[(k + 1) % poly.len()]) + 1;
    }
    let comp_area = total_area - area;
    let comp_b = total_b + b_edge + 2 - t_b;
    let g_comp: i64 = (comp_area - comp_b as f64 / 2.0 + 1.0) as i64;
    (euclidean_length(&poly) - t_l, g_comp)
}

pub fn polygon_area(poly: &[Point2D]) -> f64 {
    let mut area: f64 = 0.0;
    for i in 0..poly.len() {
        let (p1, p2) = (poly[i], poly[(i + 1) % poly.len()]);
        area += (p1.x as f64 * p2.y as f64) - (p2.x as f64 * p1.y as f64);
    }
    (area / 2.0).abs()
}

pub fn compute_good_set(
    ch_m: &[Point2D],
    d_bad: f64,
) -> (GridSet, Vec<Point2D>, Vec<Point2D>, Vec<Point2D>) {
    let expand = (3.5 * d_bad + 5.0).ceil() as i32;
    let (min_x, max_x, _, max_y) = bound(&[ch_m], expand);
    let mut good = GridSet::new(min_x, max_x, 0, max_y);
    let (mut bad_in, mut so_so) = (Vec::new(), Vec::new());
    let mut bad_out: Vec<Point2D> = Vec::new();
    for y in 0..=max_y {
        for x in min_x..=max_x {
            let p = Point2D::new(x, y);
            let (f_in, d) = (
                is_point_in_polygon(ch_m, p),
                polygon_boundary_distance(ch_m, p),
            );
            if d > d_bad {
                if f_in {
                    bad_in.push(p);
                    so_so.push(p);
                } else {
                    bad_out.push(p);
                }
                continue;
            }
            if f_in || d <= d_bad {
                good.insert(p);
                so_so.push(p);
            }
        }
    }
    bad_in.push(Point2D { x: 0, y: 0 });
    let bad_ch_ext = convex_hull(&bad_in);
    bad_in.pop();
    for y in 0..=max_y {
        for x in min_x..=max_x {
            let p = Point2D::new(x, y);
            if (x != 0 || y != 0) && good.contains(&p) && is_point_in_polygon(&bad_ch_ext, p) {
                good.delete(p);
                bad_in.push(p);
            }
        }
    }
    let bad_ch = convex_hull(&bad_in);
    good.compute_points();
    (good, bad_ch, so_so, bad_out)
}

#[allow(dead_code)]
pub fn len_longest_edge(poly: &[Point2D]) -> f64 {
    let mut max_sq = 0;
    for i in 0..poly.len() {
        let (p1, p2) = (poly[i], poly[(i + 1) % poly.len()]);
        max_sq = max_sq.max(dot(&(p2 - p1), &(p2 - p1)));
    }
    (max_sq as f64).sqrt()
}

#[allow(dead_code)]
pub fn len_longest_primitive_edge(poly: &[Point2D]) -> f64 {
    let mut max_sq = 0;
    for i in 0..poly.len() {
        let (p1, p2) = (poly[i], poly[(i + 1) % poly.len()]);
        let mut d = p2 - p1;
        let g = gcd(d.x as i64, d.y as i64).abs();
        if g > 0 {
            d.x /= g as CoordType;
            d.y /= g as CoordType;
        }
        max_sq = max_sq.max(dot(&d, &d));
    }
    (max_sq as f64).sqrt()
}

#[allow(dead_code)]
pub fn compute_max_turn_angle(poly: &[Point2D]) -> f64 {
    let mut max_angle: f64 = 0.0;
    for i in 0..poly.len() {
        let (p1, p2, p3) = (
            poly[i],
            poly[(i + 1) % poly.len()],
            poly[(i + 2) % poly.len()],
        );
        let (v1, v2) = (p2 - p1, p3 - p2);
        let (m1, m2) = (v1.norm(), v2.norm());
        if m1 > 0.0 && m2 > 0.0 {
            max_angle = max_angle.max((dot(&v1, &v2) as f64 / (m1 * m2)).clamp(-1.0, 1.0).acos());
        }
    }
    max_angle
}
