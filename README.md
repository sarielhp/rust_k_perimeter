# Smallest k-Perimeter Polygon

This project provides a Julia based algorithmic solver to find the smallest perimeter polygon with vertices on the integer grid that encapsulates exactly $k$ grid points, subject to the condition that it strictly encloses a set of "good" points and strictly excludes a set of "bad" points.

## Project Structure

- `src/k_perimeter.jl`: Contains the core algorithmic logic and geometric functions.
- `bench2.jl`: A benchmarking script to test the computation time for various values of $k$.
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
julia --project src/k_perimeter.jl 10
```

This will run the full computation, output the theoretical circle bounds, and the exact coordinates of the generated polygon. It will then generate a PDF visualization of the resulting optimal polygon inside `output/10.pdf`.

To run benchmarks:
```bash
julia --project bench2.jl
```

## Requirements
- Julia 1.x
- `DataStructures.jl`
- `Cairo.jl` (For PDF visualizations)
- `Printf` and `LinearAlgebra` (Part of standard library)
