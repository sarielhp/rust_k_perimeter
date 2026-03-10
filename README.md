# Smallest k-Perimeter Polygon

Joint work with Elfarouk Harb.

## rust_k_perimeter

Rust code to compute the minimum perimeter grid polygon (i.e., polygon
whose vertices are integral points) containing k points. It is not
hard to come up with O(k^3) running time, but this solution does
slightly better, by aggressively restricting the configurations being
inspected

This project provides a Julia based algorithmic solver to find the
smallest perimeter polygon with vertices on the integer grid that
encapsulates exactly $k$ grid points, subject to the condition that it
strictly encloses a set of "good" points and strictly excludes a set
of "bad" points.

## Project Structure

- `src/*`: Contains the core algorithmic logic and geometric functions.
- `output/`: (Created upon running) Holds the PDF visualizations of the generated polygons.

## Core Algorithm

The algorithm leverages an **Even-Odd Dynamic Programming (DP)** approach:

1. First, an estimated bounding shape (`ch_m`) is created based on a disk of radius determined by $k$, to separate necessary "good" internal points and "bad" external points.

2. We then traverse the state space searching for a sequence of coordinate choices `(loc_prev, loc, n_g)` where `n_g` is the accumulated number of enclosed grid points, optimizing for minimal boundary length (`perimeter_so_far`).

3. State transitions (`comp_next_conf`) check for valid polygon properties: counter-clockwise (left) turns, strictly enclosing valid grid points, keeping the angle within a specific bound, and honoring point sets `good` and `bad`.

## Usage Instructions

To run the solver and generate a PDF output for a given $k$:

```bash
# Example for k = 10
cargo run --release -- 44
```

This will run the full computation, output the theoretical circle
bounds, and the exact coordinates of the generated polygon. It will
then generate a PDF visualization of the resulting optimal polygon
inside `output/44.pdf`.

## Requirements
- Rust

## Disclaimer

The program was originally written in Julia and then translated into
rust using (google) antigravity, and then a lot of fine tuning. 

## Hacks

- Using an annulus to identify good/bad points.

- Computing the convex-hull of the bad points together with the origin, to identify
  good points that are undercover bad points, and register them as bad points.

  + bad_ch: Compute the convex-hull of all the bad points (old and newly minted).

- For each good point p, we precompute distance from p to the origin, going around the
  bad_ch. This is dto (distance to origin) in the code. Now, whenever considering a 
  configuration, we compute a lower bound of completing this solution into a full 
  solution by adding the dto to the current perimeter. If this is bigger than our 
  current best solution, then we eliminate the configuration from the dp (by not 
  inserting it to the heap).

- Making the CoordType to be i16 reduced memory usage and enabled
  running this program on larger inputs (currently k=19_999, but bigger
  inputs should be possible). Note, that approximate solutions for this
  problem are easy to find, so this is useful only for exact algorithm
  fanatics, or people interested in the dynamic programming solution
  (which is admittedly quite interesting here).

