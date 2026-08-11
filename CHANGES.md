# CHANGES

## 2026-08-12 — Parallel alpha sweeps in rust-foil

### rust-foil

- `aseq_par`: parallel alpha sweep.  The sweep is split into one contiguous
  chunk per rayon thread; each chunk warm-starts point-to-point and only
  the chunk's first point cold-starts (the same initialization the serial
  sweep's first point does).  Measured on NACA 0012, Re = 1e6, 80 points
  (Intel Core i5-8365U, 8 threads): **2.98s -> 1.43s (2.1x)**, matching the
  serial sweep bit-for-bit (max |dCL| = 0, zero convergence-flag
  differences).
- **Re guard**: below Re = 1e6 the boundary-layer convergence is
  history-sensitive — the cold-started chunk boundaries can land on
  different (non-)converged states (at Re = 1.4e5 the values diverge by up
  to |dCL| ~ 0.35; at Re = 5e4 the convergence flags flip and parallel is
  slower than serial).  `aseq_par` therefore falls back to the serial path
  below 1e6, so it never changes results.  Verified by two tests:
  `aseq_par_matches_aseq` (Re = 1e6 parity) and
  `aseq_par_falls_back_below_1e6` (bit-identical below the threshold).
- First dependency: rayon.

### proply-rs

No functional changes: the design loop's polar simulations run at
Re < 1e6, where the alpha sweep cannot be split internally (see the Re
guard above), and the station shapes are only discovered during the
sequential design, so there is nothing independent to parallelize ahead
of time.  `aseq_par` is available to callers analyzing airfoils at
Re >= 1e6 (e.g. wind-turbine cases).

### Reverted during development (all changed the design result)

Three "optimizations" were implemented, measured, and reverted because
each changed the design (verified against a 2.58 N cold-run baseline):

1. **Polar preheat** (simulate the polars each station's optimizer could
   reach, in parallel, before the design loop).  Reverted: it over-
   simulated (169 polars vs the ~34 the design actually needs), making
   the cold run *slower* (2m48s vs 2m08s), and the hub stations' chord
   bound `2πr/(B+2)/cos(twist)` is twist-dependent, so the twist-0
   preheat shapes collided with the design's shapes via the 2-decimal
   cache hash and the design used the wrong airfoils' polars (2.48 N).
2. **Shared engine template** (clone a per-shape rust-foil engine instead
   of rebuilding per polar).  Reverted: the same 2-decimal hash rounding
   groups near-identical-but-not-identical shapes, so a cached engine
   simulated the wrong airfoil for later stations in the bucket (2.54 N).
3. **Mach-free polar cache key**.  Reverted: the Mach number is not
   applied to the XFOIL run, but it was *accidentally disambiguating* the
   hash collisions — different stations at the same (hash, Re) with
   different shapes have different velocities, hence different Mach
   values, hence separate cache keys and each station used its own polar.
   Dropping Mach merged the colliding shapes (2.55 N).  This is
   pre-existing behavior shared with the Python (its cache keyed the same
   way); fixing the 2-decimal hash precision is a separate decision.

The three attempts were caught by comparing the design result against the
baseline (2.58 N) — the design number is a sensitive oracle for polar-path
changes.

### Measured

- Single 80-point sweep, NACA 0012, Re = 1e6: 2.98s -> 1.43s (2.1x).
- Cold-cache design run (`test_prop.json`, resolution 30): 2m08s before
  and after — unchanged, by design (see above).  Reruns against a warm
  cache are fast (seconds) either way.

### Note on design numbers

A cold-cache design run designs to Total Thrust 2.58 N, not the 2.39 N
seen when running against the pre-existing `foil_cache.json`: that warm
cache holds polars simulated by an earlier state of the pipeline (the
closed-TE polar preparation change does not alter the cache key), so a
fresh cache gives a different design.  Delete `foil_cache.json` to
re-verify acceptance numbers against the current code.
