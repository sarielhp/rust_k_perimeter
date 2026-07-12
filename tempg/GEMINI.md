I will now create the `GEMINI.md` file for the `rust_k_perimeter` project, detailing the architecture, data structures, and performance characteristics as requested.

```markdown
# rust_k_perimeter Project Overview

This project implements a solver for the **k-perimeter problem**: finding a polygon with minimum Euclidean perimeter that encloses exactly $k$ grid points. The implementation utilizes a Dynamic Programming (DP) approach accelerated by a precomputed Visibility Graph and memory-mapped storage for massive state spaces.

## High-Level Logic

1.  **Initial Bound Estimation**: The algorithm starts by generating a convex hull of a disk-like set of $k$ points around the origin to establish a baseline and a search area.
2.  **Grid Pruning**: A "good set" of grid points is identified within a specific width ($l$) of the estimated boundary. Points too far inside (the "bad set") are excluded to ensure the DP only explores potential boundary vertices.
3.  **Visibility Graph Construction**: A graph is built where edges represent valid segments between grid points. For each edge, the algorithm calculates the number of grid points enclosed in the triangle formed by the origin and the edge using **Pick's Theorem** ($A = I + B/2 - 1$).
4.  **DP Solver**: A priority-queue-based DP explores configurations. A state is defined by `(current_location, points_enclosed)`. The goal is to reach exactly $k$ points with minimum perimeter.
5.  **Solution Reconstruction**: Once the target $k$ is reached, the optimal path is reconstructed using back-pointers stored in the state values.

## Core Data Structures

### VisibilityGraph
Defined in `geom.rs`, the `VisibilityGraph` is the backbone of the search space. It stores an adjacency list of `EdgeInfo` structs:
- `target_id`: Index of the destination point.
- `n_g_delta`: The number of grid points added to the total count by choosing this edge (calculated via Pick's Theorem).
- `edge_len`: Precomputed Euclidean distance.
- `target_dto`: Distance from the target point back to the origin (used as a heuristic for pruning).

### Memory-Mapped DP States
To handle large $k$ values where the state space exceeds available RAM, the project uses `mmap_vec::MmapVec` in `dp.rs`.
- **`DPStateKey`**: A 64-bit packed representation of `loc_id` and `n_g`. It implements `AnyBitPattern` from the `bytemuck` crate.
- **`DPStateValue`**: Contains the `perimeter_so_far`, `prev_idx` (for reconstruction), and the key. 
- **`DPContext`**: Encapsulates the mutable state of the solver, including the `FxHashMap` (for state deduplication) and the `MmapVec`.

## Performance Notes

### Bytemuck Implementation
The use of `bytemuck::AnyBitPattern` combined with `#[repr(C)]` on DP structs is a critical performance optimization. It allows the `MmapVec` to treat the underlying file as a raw array of structs. This enables:
- **Zero-copy I/O**: The OS handles paging data in and out of memory.
- **Efficient Hashing**: `DPStateKey` can be cast and compared at the bit level safely.

### Deduplication and Pruning
- **`d_all` (FxHashMap)**: Uses the fast `rustc_hash` to track reached states.
- **Topological Cleanup (Default)**: The solver uses a topological ordering of the visibility graph to prune `d_all`. Since the DP explores points in a topological order, states with a topological index smaller than the current frontier can be safely removed.
- **`filter_d_all`**: Periodically prunes the hash map based on the selected `QueueStrategy` (e.g., by `n_g` or topological index).

## Known Issues and Constraints

- **Refutable Patterns**: The code utilizes `if let Ok(...)` patterns for DP initialization, which may hide disk I/O or mapping failures if not handled with explicit logging.
- **Integer Overflow**: While `CoordType` is `i16`, geometric calculations (like cross products) are promoted to `i64` to prevent overflow. However, extremely large $k$ values might approach the limits of `u32` for grid point counts in specific edge calculations.
- **Floating Point Comparison**: The DP uses `total_cmp` via an `OrderedFloat` wrapper to allow `f64` values in a `BinaryHeap`.
- **Strict Convexity**: The `is_all_left_turns` check ensures the generated polygon remains convex relative to the origin/bad-set, but may require careful tuning of the `l` parameter for non-standard $k$ distributions.
```
