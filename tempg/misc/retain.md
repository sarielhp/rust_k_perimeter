# Impact of Retain Factor on DP Performance (k=300,000)

The `retain_factor` controls the growth of the dynamic threshold for hash map pruning. A smaller factor (e.g., 1.1) prunes more frequently, keeping memory usage low but potentially increasing the number of cleanup passes. A larger factor (e.g., 4.0) allows the hash map to grow larger between cleanups.

| Retain Factor | DP Duration (s) | Configs Computed | Total Time (s) |
| :--- | :--- | :--- | :--- |
| 1.1 | 603.48 | 209,121,614 | 673.67 |
| 1.2 | | | |
| 1.5 | | | |
| 1.8 | | | |
| 2.0 | | | |
| 2.5 | | | |
| 3.0 | | | |
| 4.0 | | | |

## Observations and Recommendation
(Pending completion of benchmarks)
