use crate::point::{CoordType, Point2D};
use std::cmp::{max, min};

#[derive(Debug, Clone, Copy)]
pub struct BBox {
    pub min_x: CoordType,
    pub max_x: CoordType,
    pub min_y: CoordType,
    pub max_y: CoordType,
}

impl BBox {
    pub fn new(p: Point2D) -> Self {
        Self {
            min_x: p.x,
            max_x: p.x,
            min_y: p.y,
            max_y: p.y,
        }
    }

    pub fn extend(&mut self, p: Point2D) {
        self.min_x = min(self.min_x, p.x);
        self.max_x = max(self.max_x, p.x);
        self.min_y = min(self.min_y, p.y);
        self.max_y = max(self.max_y, p.y);
    }

    pub fn merge(&mut self, other: &BBox) {
        self.min_x = min(self.min_x, other.min_x);
        self.max_x = max(self.max_x, other.max_x);
        self.min_y = min(self.min_y, other.min_y);
        self.max_y = max(self.max_y, other.max_y);
    }

    #[allow(dead_code)]
    pub fn contains(&self, p: Point2D) -> bool {
        p.x >= self.min_x && p.x <= self.max_x && p.y >= self.min_y && p.y <= self.max_y
    }
}

#[allow(dead_code)]
pub enum KDNode {
    Leaf {
        points: Vec<Point2D>,
        bbox: BBox,
    },
    Internal {
        split_val: CoordType,
        split_dim: usize, // 0 for x, 1 for y
        left: Box<KDNode>,
        right: Box<KDNode>,
        bbox: BBox,
        count: usize,
    },
}

impl KDNode {
    pub fn bbox(&self) -> &BBox {
        match self {
            KDNode::Leaf { bbox, .. } => bbox,
            KDNode::Internal { bbox, .. } => bbox,
        }
    }

    pub fn count(&self) -> usize {
        match self {
            KDNode::Leaf { points, .. } => points.len(),
            KDNode::Internal { count, .. } => *count,
        }
    }

    pub fn build(mut points: Vec<Point2D>, depth: usize) -> Self {
        if points.len() <= 8 {
            let mut bbox = BBox::new(points[0]);
            for &p in &points[1..] {
                bbox.extend(p);
            }
            return KDNode::Leaf { points, bbox };
        }

        let split_dim = depth % 2;
        if split_dim == 0 {
            points.sort_by_key(|p| p.x);
        } else {
            points.sort_by_key(|p| p.y);
        }

        let mid = points.len() / 2;
        let split_val = if split_dim == 0 {
            points[mid].x
        } else {
            points[mid].y
        };

        let right_points = points.split_off(mid);
        let left = KDNode::build(points, depth + 1);
        let right = KDNode::build(right_points, depth + 1);

        let mut bbox = *left.bbox();
        bbox.merge(right.bbox());
        let count = left.count() + right.count();

        KDNode::Internal {
            split_val,
            split_dim,
            left: Box::new(left),
            right: Box::new(right),
            bbox,
            count,
        }
    }

    pub fn count_in_region<F>(&self, is_inside: F, bbox_status: &dyn Fn(&BBox) -> Status) -> usize
    where
        F: Fn(Point2D) -> bool + Copy,
    {
        match bbox_status(self.bbox()) {
            Status::Inside => self.count(),
            Status::Outside => 0,
            Status::Partial => match self {
                KDNode::Leaf { points, .. } => points.iter().filter(|&&p| is_inside(p)).count(),
                KDNode::Internal { left, right, .. } => {
                    left.count_in_region(is_inside, bbox_status)
                        + right.count_in_region(is_inside, bbox_status)
                }
            },
        }
    }
}

pub enum Status {
    Inside,
    Outside,
    Partial,
}

pub struct KDTree {
    root: Option<KDNode>,
}

impl KDTree {
    pub fn new(points: Vec<Point2D>) -> Self {
        if points.is_empty() {
            return Self { root: None };
        }
        Self {
            root: Some(KDNode::build(points, 0)),
        }
    }

    pub fn count_in_region<F>(&self, is_inside: F, bbox_status: &dyn Fn(&BBox) -> Status) -> usize
    where
        F: Fn(Point2D) -> bool + Copy,
    {
        if let Some(root) = &self.root {
            root.count_in_region(is_inside, bbox_status)
        } else {
            0
        }
    }
}
