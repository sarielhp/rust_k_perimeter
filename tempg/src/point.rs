//! 2D point structures and basic vector operations.
//!
//! Provides `Point2D` for integer coordinates (i16) and `Point2DF` for 
//! floating-point coordinates (f64), along with operator overloads.

use bytemuck::AnyBitPattern;
use std::ops::{Add, Div, Mul, Sub};

/// Integer coordinate type. i16 is used to keep the `Point2D` struct small (4 bytes).
pub type CoordType = i16;

/// A 2D point with integer coordinates.
/// `#[repr(C)]` and `AnyBitPattern` allow for efficient storage in memory-mapped files.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, AnyBitPattern)]
pub struct Point2D {
    pub x: CoordType,
    pub y: CoordType,
}

impl Point2D {
    pub const fn new(x: CoordType, y: CoordType) -> Self {
        Self { x, y }
    }

    pub fn is_zero(&self) -> bool {
        self.x == 0 && self.y == 0
    }

    /// Squared Euclidean norm (distance from origin squared).
    pub fn norm_sq(&self) -> i64 {
        ((self.x as i64 * self.x as i64) + (self.y as i64 * self.y as i64)) as i64
    }

    /// Euclidean norm (distance from origin).
    pub fn norm(&self) -> f64 {
        (self.norm_sq() as f64).sqrt()
    }
}

pub const ORIGIN: Point2D = Point2D::new(0, 0);

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

/// 2D Dot product of two integer vectors.
pub fn dot(a: &Point2D, b: &Point2D) -> i64 {
    a.x as i64 * b.x as i64 + a.y as i64 * b.y as i64
}

/// A 2D point with floating-point coordinates.
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
