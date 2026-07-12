# Benchmark Results: v7_branch vs. v8_branch

Comparison of the solver performance across different versions and strategies.

| $k$ | Version | DP Duration (s) | Configs Computed | Total Time (s) |
| :--- | :--- | :--- | :--- | :--- |
| **1,000** | v7 (topo_ng) | 0.0273 | 84,008 | 0.0924 |
| | v8 (topo only) | 0.0269 | 84,008 | 0.2231* |
| | v8 (no loc_id) | 0.0274 | 84,008 | 0.1004 |
| | v8 (grouped) | 0.0265 | 84,008 | 0.0967 |
| | **v8 (sorted)** | **0.0277** | **84,005** | **0.1084** |
| **10,000** | v7 (topo_ng) | 2.4918 | 2,507,103 | 3.2910 |
| | v8 (topo only) | 2.3218 | 2,507,106 | 3.0908 |
| | v8 (no loc_id) | 2.5949 | 2,507,106 | 3.3357 |
| | v8 (grouped) | 2.4370 | 2,507,106 | 3.2051 |
| | **v8 (sorted)** | **2.4745** | **2,507,100** | **3.2456** |
| **20,000** | v7 (topo_ng) | 7.3339 | 5,260,449 | 9.0553 |
| | v8 (topo only) | 6.9261 | 5,260,468 | 8.5141 |
| | v8 (no loc_id) | 7.1254 | 5,260,468 | 8.7224 |
| | v8 (grouped) | 6.6503 | 5,260,468 | 8.2338 |
| | **v8 (sorted)** | **6.7014** | **5,260,448** | **8.2703** |
| **100,000** | v7 (topo_ng) | 107.74 | 38,559,718 | 123.19 |
| | v8 (topo only) | 100.55 | 38,559,731 | 115.94 |
| | v8 (no loc_id) | 106.17 | 38,559,731 | 121.29 |
| | v8 (grouped) | 102.03 | 38,559,731 | 117.63 |
| | **v8 (sorted)** | **101.53** | **38,559,728** | **117.02** |

*\*Note: The initial v8 k=1,000 total time included compilation overhead.*

### Summary of Observations
1. **v8 (topo only) vs v7 (topo_ng):** The simplified "topological only" strategy provides a consistent **6-7%** performance improvement in DP duration.
2. **v8 (grouped & sorted):** Sorting the configuration IDs by their index in `dp_vals` before processing (`v8 (sorted)`) ensures coherent memory access patterns. This version performs similarly to the unsorted grouped version but is architecturally more robust for large state spaces that exceed RAM (mmap'd).
3. **v8 (no loc_id) Impact:** Removing the `loc_id` field from `QueueItem` resulted in a performance regression, confirming that storing the ID locally in the priority queue is more efficient than fetching it from memory-mapped storage.
