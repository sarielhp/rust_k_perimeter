use bytemuck::{AnyBitPattern};

use std::{
    //cmp::{max, min},
    ops::{Add, Div, Mul, Sub},
};

pub type CoordType = i16;

#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, AnyBitPattern)]
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
