# rust-foil performance notes and further optimization ideas

## State of play (2026-09)

Profiling the canonical workload (`examples/sweep_bench.rs`: NACA 0012,
Re = 1e6, 80-point alpha sweep) showed the viscous Newton iteration
(`oper::viscal`) dominates. Per-iteration phase costs before optimization:

- `blsolv` (block Newton solve): ~70% of iteration time
- `setbl` Jacobian fill (vm_fill + u2_m/d2_m loops): ~15%
- `mrchdu` (mixed-mode BL march): ~8%
- `update` (Newton update, dui matrix-vector product): ~4%

Optimizations landed (all bit-identical: same FP operations in the same order
per element; serial/parallel sweep parity remains exact):

- `vm` mass-influence array relayouted from Fortran's `(i, j, k)` (k
  innermost, stride IZX) to three packed `nsys x nsys` planes (`(k, i, j)`,
  j innermost).  Every hot loop in `blsolv` and the `setbl` fill is now a
  unit-stride slice axpy/scale that auto-vectorizes, and the active region
  (~3*nsys^2, ~800 KB at npan=160) is ~16x denser in cache.  `blsolv` alone
  got ~5x faster.
- `blsolv` back-substitution loops interchanged (kv outer, iv inner
  descending — same accumulation order per element) for unit-stride reads.
- `blsolv` kv-elimination: the three per-component gate branches share one
  pass over the source row when all gates hit (the common case).
- `setbl` Jacobian fill and `update`'s dui product rewritten over slices
  (bounds-check-free, unit-stride); the `vdel` gathers in `update` hoisted
  into per-call `bl_mvd`/`bl_avd` scratch buffers.
- Net effect on the sweep workload: serial ~2x, parallel ~1.9x.

## Ideas not yet pursued

Measurement caution: these were profiled on a thermally-throttled laptop
(i5-8365U); identical binaries varied +/-30% run to run.  Re-measure on a
quiet machine before/after any of these.

1. **Shrink `vm` to the live system.** `Xfoil::new()` still allocates
   `3*IZX*IZX` f64 = 12.5 MB per engine instance, but the packed layout only
   touches `3*nsys^2` (~800 KB at npan=160).  Allocating `3*nsys^2` when nsys
   is set (in `iblsys`) would cut per-instance memory ~16x and improve cache
   behaviour of `aseq_par`, where 8 engine clones stream 100 MB of buffers.

2. **Build flags.** `RUSTFLAGS="-C target-cpu=native"` (AVX2 + FMA, 4-wide)
   gave another ~15% on the sweep.  Alternatively ship runtime
   multiversioning for `blsolv`'s inner loops via
   `is_x86_feature_detected` + `#[target_feature(enable = "avx2,fma")]` —
   bit-identical (per-element ops unchanged), at the cost of `unsafe`.

3. **`blsolv` kv-column scan.** The gate scan does 3 strided loads per
   (iv, kv) pair with data-dependent branches (~45% irregular hit rate), and
   measured ~1.8 ms/call of scan overhead in isolation.  A blocked,
   branchless gate precompute (scan 8-16 kv rows into a bitmask, then
   eliminate) was tried and landed within measurement noise on the throttled
   laptop; worth re-testing on a quiet machine.  Software prefetch of the
   next column entries is another option.

4. **Rank-2 blocking of the trailing update.** Eliminating two adjacent
   columns per pass over a destination row would halve `vm` traffic in
   `blsolv`, but changes FP association (`d -= a*x + b*y` instead of two
   separate subtractions), so results are no longer bit-identical — only do
   this if the tolerance-based tests are deemed sufficient.

5. **BL closure arithmetic** (`blvar`/`blmid`/`blkin`/`trchek*`, the bulk of
   `mrchdu`/`mrchue`): chains of `powf`/`exp`/`log`/`sqrt`.  Some could be
   strength-reduced (e.g. `powf(x, 1.5)`), but not bit-identically; and
   `blkin`'s compressibility terms could be skipped when the local Mach is
   unchanged.  Moderate win, highest regression risk of the list.

6. **`dij_t` gathers.** `setbl`'s u2_m loop and `update`'s dui product read
   a gathered row of `dij_t` (4.2 MB array) per station.  Rows are hot only
   by accident of the marching order; a tiling/reordering pass could cut L3
   traffic, at significant complexity.

## Not worth it (measured/considered)

- Parallelizing `blsolv`'s kv loop with rayon: per-iv phases are ~30-100 us,
  spawn overhead dominates.  `aseq_par` scaling (~3x) is hardware-limited on
  the 4-core test laptop, not code-limited.
- The final back-substitution and the VA/VB block row ops: already
  vectorized and small vs the kv elimination.
