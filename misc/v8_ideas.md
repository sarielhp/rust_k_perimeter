# Core Ideas and Optimizations in k_perimeter (v8_branch)

This document details the architectural and algorithmic optimizations implemented in the `v8_branch` of the `rust_k_perimeter` project to solve the minimum perimeter $k$-enclosure problem.

## 1. Geometric Search Space Pruning
The problem of finding a minimum perimeter polygon enclosing $k$ grid points is computationally expensive. The program uses several geometric strategies to prune the search space before and during the DP:

*   **Good Set / Bad Set Partitioning:** Instead of searching the entire grid, the algorithm estimates a disk-like boundary for $k$ points. It identifies a "good set" of grid points within a narrow width ($l$) of this boundary. Points significantly inside the boundary (the "bad set") are treated as an opaque interior, and the search is restricted to the "good set" to find the optimal boundary.
*   **Visibility Graph with DAG constraints:** A visibility graph is constructed where edges represent potential polygon segments. To ensure the DP progresses and the resulting polygon is convex, edges are filtered based on a `max_turn_angle` and must satisfy "all left turns" relative to the origin/bad-set.
*   **Pick's Theorem:** To efficiently count the number of grid points enclosed by a segment $(u, v)$, the algorithm uses Pick's Theorem: $Area = I + \frac{B}{2} - 1$. This allows calculating the number of *new* interior points added by an edge without explicitly scanning the grid.

## 2. Advanced DP Data Structures
To handle $k$ values where the state space (number of configurations) exceeds available RAM, several systems-level optimizations are employed:

*   **Memory-Mapped Storage (`MmapVec`):** The primary storage for DP configurations (`DPStateValue`) is backed by a memory-mapped file. This allows the OS to manage paging data in and out of RAM, enabling the solver to explore hundreds of millions of states (gigabytes of data) while keeping a low resident memory footprint.
*   **Zero-Copy I/O with `bytemuck`:** DP structs are marked with `#[repr(C)]` and implement `bytemuck::AnyBitPattern`. This allows the `MmapVec` to treat the file as a raw array of structs, eliminating serialization/deserialization overhead.
*   **Fast Deduplication with `rustc-hash`:** An `FxHashMap` is used to track reached states `(loc_id, n_g)`. It uses the fast `rustc_hash` algorithm (preferred by the Rust compiler) and stores packed `u64` keys for maximum efficiency.

## 3. Algorithmic Optimizations
*   **Mandatory Topological Sort:** The visibility graph is treated as a Directed Acyclic Graph (DAG). A mandatory topological sort is performed, and configurations are processed in topological order. This ensures that when a point is processed, all potential "previous" points have already been explored.
*   **Topological Pruning:** Because of the topological ordering, once the DP frontier moves past a certain topological index, any state in the deduplication map (`d_all`) with a lower index can never be reached again. The algorithm periodically "cleans" the hash map by retaining only states with a topological index $\ge$ the current frontier, preventing the map from growing indefinitely.
*   **Grouped and Sorted Processing:**
    *   **Grouping:** Multiple configurations reaching the same `loc_id` are popped from the priority queue together and processed as a batch. This allows fetching the neighbors list for that location from the visibility graph exactly once for the entire batch.
    *   **Sorting for Locality:** Inside the batch processing, configuration IDs are sorted by their index in the `MmapVec` (`ids.sort_unstable()`). This ensures that memory accesses to the large, potentially paged-out DP state vector are sequential rather than random, significantly improving cache hit rates and reducing disk I/O.

## 4. Heuristics and Early Exit
*   **Distance to Origin (DTO):** Each edge carries a precomputed `target_dto` (Euclidean distance from the destination point back to the origin/bad-set). This acts as a lower-bound heuristic: if `perimeter_so_far + target_dto` already exceeds the global `opt_perim`, the branch is pruned immediately.
*   **Enclosure Pruning:** If a segment adds more points than remaining to reach $k$, or if it's impossible to reach $k$ even by adding the maximum possible points from that location, the configuration is discarded.
*   **Global Upper Bound:** The solver initializes with a baseline perimeter (e.g., from a circle-like polygon) to enable aggressive pruning from the very first iterations.

## 5. Performance Monitoring
*   **Dynamic Thresholding:** The hash map pruning frequency is controlled by a dynamic threshold that scales with the size of the data. The `--retain` command-line flag allows tuning the factor by which this threshold increases after each cleanup pass, balancing the cost of the cleanup against the memory savings.
*   **Formatted Status Updates:** Periodic updates display configuration counts, `n_g` progress, and memory usage (used vs. total mapped), allowing for observation of the solver's health during long runs.
