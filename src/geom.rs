use std::{
    cmp::{max, min},
    ops::{Add, Div, Mul, Sub},
};

pub type CoordType = i16;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Point2D {
    pub x: CoordType,
    pub y: CoordType,
}

impl Point2D {
    pub fn new(x: CoordType, y: CoordType) -> Self {
        Self { x, y }
    }

    pub fn is_zero(&self) -> bool {
        self.x == 0 && self.y == 0
    }

    pub fn norm_sq(&self) -> f64 {
        ((self.x * self.x) + (self.y * self.y)) as f64
    }

    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }
}

impl Add for Point2D {
    type Output = Point2D;
    fn add(self, rhs: Self) -> Self::Output {
        Point2D::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Point2D {
    type Output = Point2D;
    fn sub(self, rhs: Self) -> Self::Output {
        Point2D::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<i32> for Point2D {
    type Output = Point2D;
    fn mul(self, rhs: i32) -> Self::Output {
        Point2D::new(self.x * rhs as CoordType, self.y * rhs as CoordType)
    }
}

impl Div<i32> for Point2D {
    type Output = Point2D;
    fn div(self, rhs: i32) -> Self::Output {
        Point2D::new(self.x / rhs as CoordType, self.y / rhs as CoordType)
    }
}

pub fn dot(a: &Point2D, b: &Point2D) -> i64 {
    a.x as i64 * b.x as i64 + a.y as i64 * b.y as i64
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Point2DF {
    pub x: f64,
    pub y: f64,
}

impl Point2DF {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

impl Add for Point2DF {
    type Output = Point2DF;
    fn add(self, rhs: Self) -> Self::Output {
        Point2DF::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Point2DF {
    type Output = Point2DF;
    fn sub(self, rhs: Self) -> Self::Output {
        Point2DF::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Point2DF {
    type Output = Point2DF;
    fn mul(self, rhs: f64) -> Self::Output {
        Point2DF::new(self.x * rhs, self.y * rhs)
    }
}

pub fn vtrans(p: &[Point2D], v: Point2D) -> Vec<Point2D> {
    p.iter().map(|&pt| pt + v).collect()
}

pub fn is_lefteq_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x - a.x) * (c.y - a.y) >= (b.y - a.y) * (c.x - a.x)
}

pub fn d_y(a: Point2D, b: Point2D) -> f64 {
    (a - b).norm()
}

pub fn is_left_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x - a.x) * (c.y - a.y) > (b.y - a.y) * (c.x - a.x)
}

pub fn is_right_turn(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x - a.x) * (c.y - a.y) < (b.y - a.y) * (c.x - a.x)
}

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

pub fn euclidean_length(sol: &[Point2D]) -> f64 {
    if sol.is_empty() {
        return 0.0;
    }
    let mut l = (sol.first().unwrap_or_else(|| {
        eprintln!("Error: sol is empty in euclidean_length");
        std::process::exit(1);
    }).clone() - sol.last().unwrap_or_else(|| {
        eprintln!("Error: sol is empty in euclidean_length");
        std::process::exit(1);
    }).clone()).norm();
    for i in 1..sol.len() {
        l += (sol[i - 1] - sol[i]).norm();
    }
    l
}

pub fn polygon_boundary_distance(poly: &[Point2D], query_point: Point2D) -> f64 {
    let n = poly.len();
    if n < 2 {
        eprintln!("Error: Polygon must have at least 2 vertices in polygon_boundary_distance");
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

pub fn double_triangle_area(a: Point2D, b: Point2D, c: Point2D) -> i32 {
    (a.x as i64 * (b.y as i64 - c.y as i64)
        + (b.x as i64) * (c.y as i64 - a.y as i64)
        + (c.x as i64) * (a.y as i64 - b.y as i64))
        .abs() as i32
}

pub fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

pub fn grid_points_inside_edge(u: Point2D, v: Point2D) -> u32 {
    if u == v {
        return 0;
    }
    let e = u - v;
    let t = gcd(e.x as i64, e.y as i64);
    assert!(t >= 1);
    (t - 1) as u32
}

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

    // Initialize min and max with the first element's y
    let initial_y = points[0].y;

    let bounds = points
        .iter()
        .skip(1)
        .fold((initial_y, initial_y), |(min, max), p| {
            (
                min.min(p.y), // Returns the minimum of the two
                max.max(p.y), // Returns the maximum of the two
            )
        });

    Some(bounds)
}

fn add_points(points: &mut Vec<Point2D>, y_old: CoordType, y_new: CoordType) -> &mut Vec<Point2D> {
    // 1. Find matches and create new copies
    let mut new_points: Vec<Point2D> = points
        .iter()
        .filter(|p| p.y == y_old)
        .map(|p| Point2D { x: p.x, y: y_new })
        .collect();

    // 2. Add them to the existing vector
    points.append(&mut new_points);

    // 3. Sort the vector
    // Since Point2D derives Ord, it sorts by x then y automatically
    points.sort();

    // 4. Remove consecutive duplicates
    points.dedup();

    points
}

#[allow(dead_code)]
fn average_of_min_y(points: &[Point2D]) -> Option<Point2D> {
    if points.is_empty() {
        return None;
    }

    // 1. Find the minimum y-coordinate
    let min_y = points.iter().map(|p| p.y).min()?;

    // 2. Filter points that have this min_y and sum their coordinates
    let mut sum_x: i64 = 0;
    let mut sum_y: i64 = 0;
    let mut count: i64 = 0;

    for p in points.iter().filter(|p| p.y == min_y) {
        //println!("p.x: {}", p.x);
        sum_x += p.x as i64;
        sum_y += p.y as i64;
        count += 1;
    }

    // 3. Calculate averages with rounding up
    // Formula for ceiling division: (a + b - 1) / b
    //println!("sum: {}", sum_x);
    //println!("avg: {}", sum_x / count);
    //println!("count: {}", count);
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

pub fn ch_disk_origin(k: usize, f_expand: bool) -> Vec<Point2D> {
    let r = ((k as f64) / std::f64::consts::PI).sqrt() + 2.0;

    //println!("expand: {}", f_expand);
    let mut v = Vec::new();
    //    let mut v_exp = Vec::new();
    let l_bound = r.ceil() as i32 + 1;
    for x in -l_bound..=l_bound {
        for y in -l_bound..=l_bound {
            let p = Point2D::new(x as CoordType, y as CoordType);
            /*
            let mut p_alt = p;
            if f_expand {
                if y >= (-l_bound + 1) {
                    //println!("yes");
                    p_alt = Point2D::new(x, y + 1);
                }
                if y >= (l_bound - 1) {
                    //println!("no");
                    p_alt = Point2D::new(x, y - 1);
                }
            }*/
            if p.norm() <= r {
                v.push(p);
            }
        }
    }
    assert!(v.len() > k);
    v.sort_by(|a, b| a.norm_sq().partial_cmp(&b.norm_sq()).unwrap_or_else(|| {
        eprintln!("Error: NaN encountered in partial_cmp in ch_disk_origin");
        std::process::exit(1);
    }));
    v.truncate(k);
    //println!("v_len: {}", v.len());

    let Some((min_y, max_y)) = get_y_bounds(&v) else {
        eprintln!("Error: The vector was empty in ch_disk_origin");
        std::process::exit(1);
    };

    if f_expand {
        add_points(&mut v, max_y - 1, max_y);
        add_points(&mut v, min_y + 1, min_y);
    }

    //println!("{}...{}", min_y, max_y);
    // Equivalent of sort!( V, by = p -> (p.x, p.y) )
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
    let min_y = ch.iter().map(|p| p.y).min().unwrap_or_else(|| {
        eprintln!("Error: Convex hull was empty when getting min_y in ch_disk_origin");
        std::process::exit(1);
    });
    let max_x = ch
        .iter()
        .filter(|p| p.y == min_y)
        .map(|p| p.x)
        .max()
        .unwrap_or_else(|| {
            eprintln!("Error: Convex hull was empty when getting max_x in ch_disk_origin");
            std::process::exit(1);
        });

    let Some(avg) = average_of_min_y(&ch) else {
        eprintln!("Error: Vector empty in average_of_min_y called by ch_disk_origin");
        std::process::exit(1);
    };

    let l = ((max_x - avg.x) / 2).abs();

    //println!("AVERGAE: ({},{})", avg.x, avg.y);
    //println!("max_x: {}", max_x);
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
    //println!("MV: ( {}, {} )", mv.x, mv.y);
    let ch_n = translate_v_points(ch, mv);

    //println!("{:#?}", ch_n);

    ch_n
    //ch.iter().map(|&p| p - Point2D::new(avg.x, avg.y)).collect();
}

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

pub fn is_colinear(a: Point2D, b: Point2D, c: Point2D) -> bool {
    (b.x - a.x) * (c.y - a.y) == (b.y - a.y) * (c.x - a.x)
}

pub fn count_distinct(a: Point2D, b: Point2D, c: Point2D) -> u32 {
    if a == b && b == c {
        1
    } else if a == b || a == c || b == c {
        2
    } else {
        3
    }
}

pub fn triangle_count_new_points(a: Point2D, b: Point2D, c: Point2D) -> (u32, u32) {
    let area2 = double_triangle_area(a, b, c);
    let ab_g_n = grid_points_inside_edge(a, b);
    let bc_g_n = grid_points_inside_edge(b, c);
    let ac_g_n = grid_points_inside_edge(a, c);

    //    assert!(ab_g_n >= 0);
    //    assert!(bc_g_n >= 0);
    //    assert!(ac_g_n >= 0);

    let boundary_n = ab_g_n + bc_g_n + ac_g_n + count_distinct(a, b, c);

    let tri_i_new = if area2 > 0 {
        (area2 - (boundary_n as i32) + 2) / 2
    } else {
        0
    };

    let tri_b_new: i32 = if a == c {
        eprintln!("Error: a == c inside triangle_count_new_points");
        std::process::exit(1);
    } else if a == b {
        ac_g_n as i32 + 1
    } else if area2 == 0 {
        bc_g_n as i32 + 1
    } else {
        ac_g_n as i32 + bc_g_n as i32 + 1
    };

    assert!(tri_i_new >= 0);
    assert!(tri_b_new >= 0);

    (tri_i_new as u32, tri_b_new as u32)
}

pub fn angle_between(p1: Point2D, p2: Point2D) -> f64 {
    let dot_prod = (p1.x * p2.x + p1.y * p2.y) as f64;
    let det_prod = (p1.x * p2.y - p1.y * p2.x) as f64;
    det_prod.atan2(dot_prod)
}

pub fn comp_stop_indexes(v: &[Point2D], max_angle: f64) -> Vec<usize> {
    let mut stops = Vec::new();
    let dir_end = v.len() - 1;
    for i in 0..v.len() {
        let prev_dir = v[i];
        let mut f_break = false;

        for dir_i in i..v.len() {
            let vec_dir = v[dir_i];
            let alpha = angle_between(prev_dir, vec_dir);
            if alpha > max_angle {
                if dir_i > 0 {
                    stops.push(dir_i - 1);
                } else {
                    stops.push(0);
                }
                f_break = true;
                break;
            }
        }
        if !f_break {
            stops.push(dir_end);
        }
    }
    assert_eq!(stops.len(), v.len());
    stops
}

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
        let px_min = poly.iter().map(|p| p.x).min().unwrap_or_else(|| {
            eprintln!("Error: poly is empty in bound");
            std::process::exit(1);
        }) - 1;
        let px_max = poly.iter().map(|p| p.x).max().unwrap_or_else(|| {
            eprintln!("Error: poly is empty in bound");
            std::process::exit(1);
        }) + 1;
        let py_min = poly.iter().map(|p| p.y).min().unwrap_or_else(|| {
            eprintln!("Error: poly is empty in bound");
            std::process::exit(1);
        }) - 1;
        let py_max = poly.iter().map(|p| p.y).max().unwrap_or_else(|| {
            eprintln!("Error: poly is empty in bound");
            std::process::exit(1);
        }) + 1;

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
    min_x = min_x - expand as CoordType;
    max_x = max_x + expand as CoordType;
    min_y = min_y - expand as CoordType;
    max_y = max_y + expand as CoordType;
    (
        min_x - expand as CoordType,
        max_x + expand as CoordType,
        min_y - expand as CoordType,
        max_y + expand as CoordType,
    )
}

pub struct GridSet {
    pub min_x: CoordType,
    pub max_x: CoordType,
    pub min_y: CoordType,
    pub max_y: CoordType,
    _size: usize,
    width: usize,
    data: Vec<bool>,
    dto: Vec<f64>, // Distance to origin
}

impl GridSet {
    pub fn new(min_x: CoordType, max_x: CoordType, min_y: CoordType, max_y: CoordType) -> Self {
        let width = (max_x - min_x + 1).max(0) as usize;
        let height = (max_y - min_y + 1).max(0) as usize;
        
        let _size = width * height;
        println!( "Size: {}", _size );
        Self {
            min_x,
            max_x,
            min_y,
            max_y,
            _size, 
            width,
            data: vec![false; _size],
            dto: vec![-1.0; _size],
        }
    }

    pub fn length(&self) -> usize {
        let count = self.data.iter().filter(|&&x| x).count();

        count
    }

    pub fn insert(&mut self, p: Point2D) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.data[idx] = true;
        }
    }

    pub fn delete(&mut self, p: Point2D) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.data[idx] = false;
        }
    }

    pub fn insert_val(&mut self, p: Point2D, val: f64) {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.dto[idx] = val;
        }
    }

    #[allow(dead_code)]
    pub fn get_dto(&self, p: Point2D) -> f64 {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.dto[idx]
        } else {
            self.min_x as f64 * self.min_x as f64
        }
    }

    pub fn contains(&self, p: &Point2D) -> bool {
        if p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y {
            let idx = (p.y - self.min_y) as usize * self.width + (p.x - self.min_x) as usize;
            self.data[idx]
        } else {
            false
        }
    }

    pub fn fill_dist_to_origin(&mut self, bad_ch: &Vec<Point2D>) {
        for y in self.min_y..=self.max_y {
            if y < 0 {
                continue;
            }
            for x in self.min_x..=self.max_x {
                let p = Point2D::new(x, y);
                if !self.contains(&p) {
                    continue;
                }
                let l = distance_to_origin(bad_ch, p);
                self.insert_val(p, l);
            }
        }
    }
}

/// Determines the orientation of ordered triplet (p, q, r).
/// Returns:
///  0 -> p, q and r are collinear
///  1 -> Clockwise
///  2 -> Counterclockwise
fn orientation(p: Point2D, q: Point2D, r: Point2D) -> i32 {
    let val = (q.y as i64 - p.y as i64) * (r.x as i64 - q.x as i64)
        - (q.x as i64 - p.x as i64) * (r.y as i64 - q.y as i64);
    if val == 0 {
        return 0;
    }
    if val > 0 {
        1
    } else {
        2
    }
}

/// Checks if point 'q' lies on line segment 'pr'
fn on_segment(p: Point2D, q: Point2D, r: Point2D) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
}

/// Returns true if segment p1q1 and p2q2 intersect.
fn segments_intersect(p1: Point2D, q1: Point2D, p2: Point2D, q2: Point2D) -> bool {
    let o1 = orientation(p1, q1, p2);
    let o2 = orientation(p1, q1, q2);
    let o3 = orientation(p2, q2, p1);
    let o4 = orientation(p2, q2, q1);

    // General case: segments straddle each other
    if o1 != o2 && o3 != o4 {
        return true;
    }

    // Special Cases: Collinear points
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

/// Main function: Checks if segment (s1, s2) intersects the polygon boundary
pub fn does_segment_intersect_polygon(polygon: &[Point2D], s1: Point2D, s2: Point2D) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }

    for i in 0..n {
        let v1 = polygon[i];
        let v2 = polygon[(i + 1) % n]; // Wraps around to close the polygon

        if segments_intersect(s1, s2, v1, v2) {
            return true;
        }
    }

    false
}

pub fn is_all_left_turns(p: Point2D, p_next: Point2D, v: &[Point2D]) -> bool {
    // Vector A = p_next - p
    let ax = p_next.x - p.x;
    let ay = p_next.y - p.y;

    for q in v {
        // Vector B = q - p_next
        let bx = q.x - p_next.x;
        let by = q.y - p_next.y;

        // The 2D cross product (Z-component of 3D cross product)
        // formula: ax * by - ay * bx
        let cross_product = ax * by - ay * bx;

        // If cross_product > 0, it's a left turn.
        // If cross_product < 0, it's a right turn.
        // If cross_product == 0, the points are collinear.
        if cross_product <= 0 {
            return false;
        }
    }

    true
}

pub fn distance_to_origin(bad_ch: &[Point2D], p: Point2D) -> f64 {
    let origin = Point2D { x: 0, y: 0 };
    if !does_segment_intersect_polygon(bad_ch, origin, p) {
        return d_y(origin, p);
    }
    let mut pnts = bad_ch.to_vec();
    pnts.push(origin);
    pnts.push(p);
    let poly = convex_hull(&pnts);

    let n = poly.len();
    assert!(n > 2);
    assert!(is_left_turn(poly[0], poly[1], poly[2]));

    // This will return the index or panic with the message provided
    let i_o = poly
        .iter()
        .position(|&x| x == origin)
        .unwrap_or_else(|| {
            eprintln!("Error: Origin not found in CH! in distance_to_origin");
            std::process::exit(1);
        });
    let i_p = poly
        .iter()
        .position(|&x| x == p)
        .unwrap_or_else(|| {
            eprintln!("Error: p not found in CH! in distance_to_origin");
            std::process::exit(1);
        });

    let s = min(i_o, i_p);
    let t = max(i_o, i_p);
    assert!(s < t);
    let mut t_l = 0.0;
    for k in s..t {
        t_l += d_y(poly[k], poly[k + 1]);
    }
    if i_p < i_o {
        return t_l;
    }

    let total_l = euclidean_length(&poly);
    total_l - t_l
}

pub fn compute_good_bad_sets(ch_m: &[Point2D], l: f64) -> (GridSet, GridSet, Vec<Point2D>) {
    let expand = (3.0 * l + 3.0).ceil() as i32;
    let (min_x, max_x, _, max_y) = bound(&[ch_m], expand);
    let min_y = 0;

    let mut bad = GridSet::new(min_x, max_x, min_y, max_y);
    let mut good = GridSet::new(min_x, max_x, min_y, max_y);
    let mut bad_in: Vec<Point2D> = Vec::new();
    for y in min_y..=max_y {
        if y < 0 {
            continue;
        }
        for x in min_x..=max_x {
            let p = Point2D::new(x, y);
            let f_in = is_point_in_polygon(ch_m, p);
            let d = polygon_boundary_distance(ch_m, p);

            if d > l {
                bad.insert(p);
                if f_in {
                    bad_in.push(p);
                }

                continue;
            }

            if f_in || d <= l {
                good.insert(p);
            }
        }
    }

    bad_in.push(Point2D { x: 0, y: 0 });
    let bad_ch_ext = convex_hull(&bad_in);
    bad_in.pop();

    for y in min_y..=max_y {
        for x in min_x..=max_x {
            if x == 0 && y == 0 {
                continue;
            }
            let p = Point2D::new(x, y);
            if !good.contains(&p) {
                continue;
            }
            if !is_point_in_polygon(&bad_ch_ext, p) {
                continue;
            }
            good.delete(p);
            bad.insert(p);
            bad_in.push(p);
        }
    }
    let bad_ch = convex_hull(&bad_in);

    (good, bad, bad_ch)
}

pub fn len_longest_edge(poly: &Vec<Point2D>) -> f64 {
    let n = poly.len();
    if n < 2 {
        return 0.0;
    }

    let mut max_dist_sq: i64 = 0;

    for i in 0..n {
        let p1 = poly[i];
        let p2 = poly[(i + 1) % n]; // Wraps around to the first point

        let dx = p2.x - p1.x;
        let dy = p2.y - p1.y;

        // Calculate squared distance: dx^2 + dy^2
        let dist_sq = dx as i64 * dx as i64 + dy as i64 * dy as i64;

        if dist_sq > max_dist_sq {
            max_dist_sq = dist_sq;
        }
    }

    // Convert to f64 and take the square root once
    (max_dist_sq as f64).sqrt()
}

/// Computes the maximum turn angle in degrees for a convex polygon.
pub fn compute_max_turn_angle(poly: &Vec<Point2D>) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }

    let mut max_angle_rad: f64 = 0.0;

    for i in 0..n {
        // Points: current, next, and the one after that
        let p1 = poly[i];
        let p2 = poly[(i + 1) % n];
        let p3 = poly[(i + 2) % n];

        // Vector A (p1 -> p2)
        let ax = (p2.x - p1.x) as f64;
        let ay = (p2.y - p1.y) as f64;

        // Vector B (p2 -> p3)
        let bx = (p3.x - p2.x) as f64;
        let by = (p3.y - p2.y) as f64;

        let mag_a = (ax * ax + ay * ay).sqrt();
        let mag_b = (bx * bx + by * by).sqrt();

        if mag_a > 0.0 && mag_b > 0.0 {
            // Dot product: A · B = |A||B| cos(theta)
            let dot = ax * bx + ay * by;
            let cos_theta = (dot / (mag_a * mag_b)).clamp(-1.0, 1.0);

            let angle = cos_theta.acos();
            if angle > max_angle_rad {
                max_angle_rad = angle;
            }
        }
    }

    // Convert radians to degrees
    max_angle_rad
}
