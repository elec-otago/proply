# CHANGES

## 2026-08-20 — CST (Kulfan) parametrization in rust-foil

### rust-foil

- **`KulfanParams`** — new canonical geometry representation (CST "with
  LEM", the exact 18-parameter flavor AeroSandbox/NeuralFoil use):
  `Default` (AeroSandbox defaults, 8 weights/side), `coordinates()`
  (closed TE→LE→TE loop at cosine stations), `upper_y()`/`lower_y()`,
  `fit_from_coordinates()` (linear least squares via the existing `gauss`,
  with the negative-TE re-solve), `fit_from_coordinates_n()`, and
  `from_naca()` (NACA 4/5-digit → parameters).  No new dependencies.
- **`XFoil::cst(&KulfanParams)`** — new public API; `KulfanParams` is
  re-exported at the crate root.
- **`XFoil::naca()` now generates through CST**: the designation is fitted
  to 8 weights/side and the buffer is produced by the CST path (cosine
  stations instead of TE-bunched).  The NACA shape is not exactly
  CST-representable; the fit error is ~1.1–1.5e-4 in y/c (worst case),
  which the reference polars in `tests/naca0012.rs` stay within — **the
  whole existing suite passes unchanged** (no re-baselining needed).
- The naca/cst buffer-setting tail is factored into a private
  `set_buffer` helper (scalc + segspl×2 + geopar + pangen).
- `tests/cst.rs`: forward sanity, fit round-trip (~1e-15 on CST points),
  NACA symmetry (lower ≈ −upper, LEM ≈ 0, TE ≈ 2.5e-3), geometric
  fidelity (thickness/camber within 1e-3), naca↔cst buffer identity, and
  end-to-end viscous convergence.  (Note: the `Default` section has a
  perfectly sharp TE, on which XFOIL's viscous solve does not converge;
  inviscid and blunt-TE solves are fine.)

### proply-rs

- **CST foil family**: `--cst` (or `"cst": true` in the design JSON)
  switches every station foil from the NACA 4-series to a CST (Kulfan)
  section: rust-foil's canonical 18-parameter `KulfanParams` shape (default
  AeroSandbox section), re-thicknessed (weights scaled linearly — thickness
  is linear in the weights) and cambered (LEM weight set from the same
  camber candidates) to the design's radial laws.  The trailing-edge gap
  maps onto the CST TE term.
- New `foil::Cst` (FoilLike: shape points, hash, bounding box, …) and the
  `foil::FoilFamily` enum the design loop now dispatches on (`Prop` holds
  `BladeElement<FoilFamily>` instead of `BladeElement<Naca4>`); the NACA4
  path is bit-identical (golden tests unchanged).  `Cst::from_naca(code)`
  exposes NACA sections as CST parameters.
- The simulator is unchanged: CST station shapes feed the same
  coordinate path (TE closed for the polar solve) as NACA4.
- Verified end-to-end: `--cst --bem` and `--cst --lifting-line` designs
  converge and write valid STEP output.

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
