//! Polygon I/O utilities.
//!
//! Provides functions to save and load polygon vertices to/from text files.

use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};
use crate::point::*;

/// Saves a sequence of points to a text file.
/// Each line contains "x y" coordinates.
/// An optional comment can be added at the top of the file.
#[allow(dead_code)]
pub fn save_polygon(path: &str, points: &[Point2D], comment: Option<&str>) -> io::Result<()> {
    let mut file = File::create(path)?;
    if let Some(text) = comment { writeln!(file, "# {}", text)?; }
    for point in points { writeln!(file, "{} {}", point.x, point.y)?; }
    Ok(())
}

/// Loads a polygon from a text file where each line is "x y".
/// Skips comments (lines starting with '#') and empty lines.
#[allow(dead_code)]
pub fn load_polygon(path: &str) -> io::Result<Vec<Point2D>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut points = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.starts_with('#') || trimmed.is_empty() { continue; }
        let coords: Vec<&str> = trimmed.split_whitespace().collect();
        if coords.len() == 2 {
            let x = coords[0].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let y = coords[1].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            points.push(Point2D { x: (x as CoordType), y: (y as CoordType) });
        }
    }
    Ok(points)
}
