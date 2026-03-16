use std::fs::File;
use std::io::{self, Write, BufRead, BufReader};

use crate::point::*;
//mod geom;

#[allow(dead_code)]
pub fn save_polygon(path: &str, points: &[Point2D], comment: Option<&str>) -> io::Result<()> {
    let mut file = File::create(path)?;

    // Write the optional header comment
    if let Some(text) = comment {
        writeln!(file, "# {}", text)?;
    }

    // Write each point as "x y"
    for point in points {
        writeln!(file, "{} {}", point.x, point.y)?;
    }

    Ok(())
}

#[allow(dead_code)]
pub fn load_polygon(path: &str) -> io::Result<Vec<Point2D>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut points = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        // Parse coordinates
        let coords: Vec<&str> = trimmed.split_whitespace().collect();
        if coords.len() == 2 {
            let x = coords[0].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let y = coords[1].parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            points.push(Point2D { x: (x as CoordType), y: (y as CoordType) });
        }
    }

    Ok(points)
}