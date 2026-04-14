# Performance Report: AffineMapDiffEq.jl

## Benchmark Results (Final Optimized Version)

I have resolved the runtime dispatch and closure bottlenecks by introducing a type-stable `AffineHistory` wrapper and a function barrier for the mapping perturbation.

Below are the final benchmark results for the Mathieu equation (100 history steps):

| Implementation | Median Time | Allocations | Memory |
| :--- | :--- | :--- | :--- |
| **Baseline (Initial)** | ~350 ms | 3.42 M | 121 MiB |
| **Vector (In-place)** | **25.8 ms** | 0.56 M | 28 MiB |
| **MVector (Out-of-place)** | 27.5 ms | 0.70 M | 32 MiB |
| **Vector (Out-of-place)** | 70.9 ms | 1.56 M | 71 MiB |

### Key Improvements (Final):
1.  **Type-Stable History:** Created `AffineHistory` struct to wrap interpolation. This removed runtime dispatch in `perform_step!` by ensuring the history function type is fully inferable by the DDE integrator.
2.  **Function Barrier:** Refactored the perturbation mapping in `affine` using a `let` block and a function barrier. This prevents closure-based performance degradation when using `ForwardDiff.Dual`.
3.  **StaticArray Optimization:** Added type assertions and optimized return paths for `StaticArrays`, allowing them to be nearly as fast as in-place `Vector` methods.
4.  **Overall Speedup:** Achieved a **~13.5x** speedup compared to the initial baseline for the 2D Mathieu equation.

---
*Report generated on March 30, 2026.*
