# CHANGES

## 2026-09-05 — Newton-solver performance (~2x on the sweep workload)

Profiling the canonical workload (`examples/sweep_bench.rs`: NACA 0012,
Re = 1e6, 80-point alpha sweep) showed the viscous Newton iteration
dominates runtime, with `blsolv` (the coupled viscous-inviscid block
solve) alone at ~70% of each iteration.  The optimizations below keep
every floating-point operation in the same order per element, so results
are **bit-identical**: serial/parallel sweep parity is exactly zero and
the NACA 0012 reference values are unchanged.

- **`vm` mass-influence array relayouted** (`src/state.rs`): from
  Fortran's `vm(3, IZX, IZX)` (k innermost, row stride IZX) to three
  planes of `nsys x nsys` packed by the current system size, j innermost.
  `Xfoil::vm_index` now takes `nsys`.  All hot loops in `blsolv` and the
  `setbl` Jacobian fill are now unit-stride slice axpy/scale operations
  that auto-vectorize, and the active region is ~16x denser in cache
  (~800 KB at npan=160, vs 3.2 MB strided).
- **`blsolv` rewritten over slices** (`src/solve.rs`): bounds-check-free
  vectorizable row sweeps for the VA/VB/VZ block operations; the three
  per-component gates of the lower-VM-column elimination share one pass
  over the source row when all hit (the common case); the
  back-substitution loops are interchanged (kv outer, iv inner
  descending — identical accumulation order per element) so the `vm`
  reads are unit-stride row sweeps instead of 17 KB-strided gathers.
  `blsolv` is ~5x faster.
- **`setbl` Jacobian fill and `update` dui product over slices**
  (`src/bl.rs`): the u2_m/d2_m rebuild and vm fill loops are
  bounds-check-free and unit-stride; the `vdel` gathers in `update` are
  hoisted into per-call scratch (`bl_mvd`/`bl_avd`, new `Xfoil` fields
  reused via `mem::take` like the existing BL scratch).
- `tests/blsolv.rs` updated for the new `vm_index` signature.

Measured on `examples/sweep_bench.rs` (i5-8365U, thermally throttled,
best of 3): serial 10.1 s -> ~4.9 s (~2.1x), parallel (8 threads)
3.6 s -> ~1.9 s.  The `tests/naca0012.rs` integration suite runs
126 s -> 90 s.  Further ideas (smaller `vm` allocation, AVX2 build flags
or multiversioning, branchless gate scan, rank-2 blocking, BL-closure
strength reduction) are documented in `TODO.md`.
