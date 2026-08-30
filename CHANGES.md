# CHANGES

## 2026-08-30 — Expose the mechanical-law modulus on the web demo

### proply-rs web demo

The mechanical thickness law's `modulus` (material stiffness) is now
editable in the browser: the **Propeller Specifications** tab gains an
**Elastic modulus (GPa)** field (with a hover hint — nylon/ABS ≈ 3,
aluminium ≈ 70, carbon fibre ≈ 100 GPa), and the default design carries
`"modulus": "3 GPa"` so the key is visible in the JSON.  A stiffer
material sizes much thinner foils — the thickness goes as `E^(−1/3)`:
on the web-default design the root `t/c` drops from 0.28 at 3 GPa to
0.09 at 100 GPa (the outer sections then ride the `thickness_floor`).
JS-only change; the committed wasm already parses `modulus`.

## 2026-08-30 — Document the mechanical thickness law

### docs

- New [MECHANICAL.md](MECHANICAL.md) explains the mechanical strength
  calculation behind `mech_thickness`: the cantilever-beam model (thrust
  per unit span, bending moment, `I = c t³/12` stiffness, constant-
  curvature sizing `t = (6 M L² / (E c δ))^(1/3)`, the chord/twist and
  floor effects), the JSON/CLI parameters, a worked example with the
  sized sections, and the model's assumptions and limitations.
- The root [README.md](README.md) links the document from the repository
  layout and the usage section.

## 2026-08-30 — Mechanical airfoil-thickness law (beam deflection)

### proply-rs design loop

The TODO's second item: a mechanical mechanism for airfoil thickness —
"treating the blade as a beam and approximating deflection … the
thickness should be chosen to keep the blade shape from deforming too
much … the hub thickness should not be involved".  Opt-in via
`mech_thickness: true` in the design JSON (or `--mech-thickness`); the
geometric power law stays the default.

- **Beam model** (new `thickness.rs`): the blade is a cantilever beam
  anchored at the hub; the converged design's station loads (annular
  element thrusts, divided by panel width and blade count) form the load
  per unit span `q(r)` in the z direction, and the bending moment
  `M(r) = ∫ q(s)(s − r) ds` is integrated from the tip inwards.  The
  thickness is laid out for a constant-curvature bend with the tip
  deflection closed at the allowed value — closed form
  `t(r) = (6 M(r) L² / (E c(r) δ))^(1/3)` (`L = R − r_hub`,
  `δ = deflection_fraction · R`, rectangular section `I = c t³/12`) —
  floored at a minimum fraction of the local chord (default 0.06, the
  thinnest ARA-D table: the tip, where `M → 0`, would otherwise taper to
  a knife edge).  The chord enters the stiffness (`I ∝ c t³`, so a wider
  blade is stiffer), and the twist enters through the design's chord
  distribution.  The **hub thickness (`hub_depth`) is deliberately not
  involved** — the airfoil thickness follows the load, not the hub
  geometry.
- **Wiring** (`prop.rs`, `design_parameters.rs`, `pipeline.rs`,
  `main.rs`): `DesignParameters` gains `mech_thickness`, `modulus`
  (elastic modulus, Pa, unit-suffixed strings `"3 GPa"`; default
  3 GPa = moulded nylon/ABS), `deflection_fraction` (allowed tip
  deflection as a fraction of R, default 0.05) and `thickness_floor`
  (fraction of local chord, default 0.06).  The pipeline runs the design
  once for loads, sizes the thickness law (radius → t/c, so the sized
  section scales exactly with the final chord), and re-runs the design
  on that law so the reported operating point matches the mechanically
  sized blade.  The law is passed read-only into the worker threads' foil
  construction; lifting-line designs take their station loads from one
  extra circulation solve of the final geometry.  The YAML summary names
  the law (`thickness_law: geometric|mechanical`) and reports the
  predicted tip deflection (`tip_deflection_mm`).
- **Web demo** (`forms.js`, `index.html`): the Propeller Specifications
  tab gains a **Mechanical thickness (beam sizing)** checkbox; boolean
  fields are a new form type, composing `mech_thickness: true` into the
  JSON (unchecked drops the key, so the Rust default applies).

On the web-default design (68 mm, 3 blades, plate polars) the mechanical
law sizes root t/c ≈ 0.28 (vs 0.71 by the power law — by the beam model
the old law is over-thick at the root), holding the predicted tip
deflection at 3.38 mm of the 3.40 mm allowed (5% of R), versus 4.9 mm a
blade built to the geometric law is predicted to deflect.

## 2026-08-30 — ARA-D foils selectable in the web demo

### proply-rs web demo

The TODO's first item, "Add the ARA-D foils as an option in the web
interface. Make sure that foil family choice can be chosen though the
JSON input file as well".  The ARA-D family itself (and the `cst`/
`arad` JSON keys and `--cst`/`--arad` CLI flags) were already ported;
what was missing was any way to pick the family from the page.

- **Foil family select** (`forms.js`, `index.html`): the
  **Propeller Specifications** tab gains a **Foil family** dropdown —
  NACA 4-series (the default), CST (Kulfan) or ARA-D.  `forms.js` gains
  a new "select" pseudo-field type: `composeDesign` writes the choice as
  the JSON's `cst`/`arad` booleans (dropping both keys for the NACA
  default, so the Rust default applies), `syncForm` derives the dropdown
  position from those keys when a design is loaded (boot, localStorage
  or a hand edit in the textarea), and `buildForm` renders a proper
  `<select>`.
- **JSON input file**: the family remains entirely JSON-driven — the
  dropdown is just an editor for the `cst`/`arad` keys, which the CLI,
  the wasm and the JSON input files already parse (`arad.json` golden
  test).  No Rust or wasm changes were needed; the committed wasm
  package already contains the family.

## 2026-08-30 — Debug the browser-demo plate/full-polar thrust discrepancy

### proply-rs (design loop)

The TODO's `browser_demo` design (160 mm radius, 30 N target, pinned
motor_torque 1.6 Nm @ 5500 RPM) converged onto the operating point with
plate polars (Q = 1.587 Nm, T = 30.2 N) but collapsed to ~10% of the
thrust with simulated polars (Q = 0.107 Nm, T = 2.48 N).  Debugging showed
this is specific to infeasible operating points, not a general plate-vs-
polar bug: the web demo's default design (68 mm, 3 N) converges either way
(T = 1.19 N plate vs 1.41 N full polars, both at Q ≈ 0.016 Nm).  The
pinned point demands ~920 W from a 320 mm prop whose NACA sections are
t/c 0.26–0.59 at Re 30–150k — beyond what the simulated polars can
deliver (cached sweeps: cl_max ≈ 0.2–1.2).  The no-stall plate model
"absorbs" the torque in deep stall (alpha 18–38°) and reports ~30 N; with
real polars the station BEM solves fail, the torque-match loop stalls, and
the tool silently reported a broken design as if it had converged.

- **Polar continuity across stall** (`simulator.rs`): the flat-plate
  fallback switched at |alpha| > 30°, so a fitted cl of ~1 just below 30°
  jumped to 2π·alpha ≈ 3.3 just above — a discontinuity the optimizers
  exploited as a free lift source.  `get_cl`/`get_cd` now blend smoothly
  (cubic smoothstep) from the fitted polar to the analytic flat-plate
  model over |alpha| ∈ [18°, 33°]: pure polar below, pure flat plate
  above, C1-continuous throughout.  No polar-cache key change (evaluation-
  time only).
- **Honest failure reporting** (`prop.rs`, `pipeline.rs`, `yaml_out.rs`,
  `main.rs`): when the geometry cannot absorb the design torque,
  `design_for_torque` now returns an explicit `warning` ("design torque
  1.6000 Nm at 5500 rpm not achievable: the closest design absorbs 0.0977
  Nm …; 5/11 BEM stations converged") instead of silently proceeding.
  Non-converged stations keep their last induction state and are flagged
  (per-element `converged`, reported in the YAML per section and as a
  station-coverage count), rather than being silently reset to
  (dv=1, a_prime=0).  The warning travels through the CLI printout, the
  YAML summary (`warning:` key) and the browser demo status line.

The plate-polar optimum is **not** a good starting point for non-plate
designs (as the TODO asked): the plate trajectory operates at alpha
18–38°, deep stall for the real sections, so it is not near any real-
polar equilibrium; the failure is an infeasibility failure, not an
optimizer stop.  On the fixed code the browser_demo design reports its
closest achievable design explicitly (T = 1.66 N, Q = 0.098 Nm, warning
emitted; 5/11 stations converged), and the achievable web-default design
is unchanged in kind (T = 1.54 N, Q = 0.016 Nm at the 0.016 Nm target,
all stations converged, no warning).

## 2026-08-29 — Link to the live web demo from the README

### README

- The root README now links the deployed browser demo
  (<https://proply-zeta.vercel.app/>) right after the intro, so the web
  version is reachable from the project page without cloning the repo
  or running a local server.

## 2026-08-29 — Tabbed editors for the design JSON

### proply-rs web demo

- The design parameters can now be edited from three tabs above the
  JSON textarea — **Propeller Specifications** (name, blades, radius,
  thrust, tip/ hub geometry, trailing edge, scimitar, airspeed,
  altitude), **Electric Motor** (Kv, volts, no-load current, winding
  resistance) and **Other Motor** (torque and RPM, the direct operating
  point).  Each tab holds a field per JSON descriptor, in the same
  unit style as the defaults (`68 mm`, `3 N`); suffixed strings pass
  through to the wasm unchanged, and bare numbers stored in the JSON
  (SI) are converted into the display unit when the forms sync.
- The tabs compose into the design JSON shown in the textarea and
  persist it in localStorage (`proply-design-params`), so a design
  survives page reloads.  Blank fields drop their key from the JSON so
  the Rust defaults apply; the **Other Motor** pair only takes effect
  when both torque and RPM are filled — the wasm ignores one without
  the other and uses the electric motor model instead.  Hand-written
  run options (`bem`, `resolution`, `plate`, ...) are preserved by the
  merge, and hand edits in the textarea re-sync the tabs.
- New `web/forms.js` holds the per-tab field definitions and the
  compose/sync logic (pure functions, testable without a browser);
  `main.js` wires the tab and textarea changes to it; `index.html`
  gains the tabs container and its styles.  The textarea remains the
  source of truth for running a design — no changes to the design
  worker or the wasm.

## 2026-08-28 — 3D STEP preview of the finished design

### proply-rs web demo

- A "Show 3D preview" button next to the STEP download link renders the
  finished design in a 3D window: the STEP text from the design worker
  is tessellated in the tab by occt-import-js (an OpenCASCADE WASM
  build, 0.0.23) and displayed with three.js — the same stack as
  devonkcopeland/step-viewer.  The preview is opt-in and lazy: nothing
  3D loads until the button is clicked, so a design run no longer pulls
  three.js and the ~2 MB occt WASM from the CDN on its own.  The
  window shows a progress note ("tessellating the STEP model…") until
  the model is actually in the scene, supports orbit/zoom (OrbitControls)
  with lighting and a ground grid, and keeps its drawing buffer so the
  rendered model can always be captured or exported.  If tessellation
  or parsing fails, the window reports the reason and the download link
  still works.
- `viewer.js` holds the occt + three pipeline, imported by `main.js` on
  the first preview click (three.js resolves through the import map in
  `index.html`; occt-import-js and its .wasm come from the CDN);
  `index.html` gains the `#viewer` window and the button.  The button
  is enabled once a design completes and disabled while one runs;
  starting a new design clears the on-screen model (`clearStep` in
  viewer.js, disposing geometry and material) so a stale preview never
  lingers next to a fresh download link.

## 2026-08-28 — Design spinner with progress and completion info

### proply-rs web demo

- While a design runs, the status line shows an animated spinner and a
  live elapsed-time ticker (`designing… — 4.2 s elapsed`); the page
  stays interactive because the designer runs in its worker.  When the
  design completes, the spinner stops and the line reports the result:
  wall time, thrust, torque and power at the operating point, and how
  many polars were simulated and cached.

## 2026-08-28 — Attribution on the demo page

### proply-rs web demo

- The main screen now carries the author line: *proply by Tim Molteno
  (tim@elec.ac.nz)*, under the page heading.

## 2026-08-28 — The wasm designer runs in a web worker

### proply-rs web demo

- The browser demo's WebAssembly designer now runs in a dedicated worker
  (`proply-rs/web/designer.js`) instead of on the page's main thread.
  This builds on the WebAssembly port (merged today): a cold XFOIL
  design takes minutes single-threaded, and previously froze the tab for
  the whole run — now the page stays interactive (verified by reading
  page state mid-design, which used to time out).
- The worker owns the wasm module, the `PropSession` polar cache and its
  IndexedDB persistence; the page (`proply-rs/web/main.js`) only posts
  design requests and renders the results that come back.

## 2026-08-21 — Specified motor operating point

### proply-rs

- A `motor_torque` + `motor_RPM` pair in the design JSON now sets the
  design's operating point directly (engine-style: `rotax_912_uls`'
  334.8 N m at 1950 RPM — 68.4 kW shaft), overriding the electric motor
  model's maximum-efficiency point derived from `motor_Kv` /
  `motor_volts` & co.  One field without the other is ignored.  The
  resolution lives in `DesignParameters::motor_operating_point` (unit
  tested); the YAML summary's motor section, the `--auto` goal torque
  and the startup banner follow automatically.
- `props/rotax_912_uls.json` is the only shipped design using the pair;
  its committed summary and gallery numbers predate this change and
  refresh on the next `make summaries`.

## 2026-08-20 — Progress bars for the serial lifting-line phases

### proply-rs

The two long silent stretches of the lifting-line design now show
progress (indicatif bars, like the warm-up bars; hidden when stderr is
not a TTY):

- After the seeding warm-up, the per-station best-L/D seeding scan has
  a `station seeds` bar (one item per candidate x station, labelled with
  the camber candidate being scanned).  On a cold cache these scans
  simulate their own polars and the phase runs minutes with no output.
- After the composed-camber warm-up, the parallel candidate design
  passes share an evaluation counter (`design evaluations`), with the
  latest pass result as the message — the Nelder-Mead passes run full
  circulation solves (and, on a cold cache, fresh polars for the moved
  chords' Reynolds buckets), so minutes pass between printed
  improvements.

All six TODO.md items are now done; the file is empty.

## 2026-08-20 — Polar warm-up pool actually parallel

### proply-rs

The lifting-line design's polar warm-up pool (`warm_polar_pool`, used by
the seeding and composed-camber warm-ups) ran fully serial despite
spawning a worker per core — the whole pool contended on one mutex:

- **Queue-lock guard lifetime (the real serialization)**: the worker loop
  was `while let Some(task) = queue.lock().unwrap().pop() { simulate }` —
  the `MutexGuard` temporary lives to the end of the loop *body*, so the
  first worker held the queue lock for its entire rust-foil sweep and
  every other worker blocked on `pop()` (diagnosed from per-thread
  backtraces: one worker in `viscal`, seven in `lock_contended` with zero
  CPU time ever).  The loop now pops through a `match` so the guard drops
  before the simulation starts.
- **Bucket-level tasks**: the work list is flattened to one task per
  distinct `(foil, Reynolds-grid bucket, Mach)` warm target
  (`FoilSimulator::warm_plan` / `warm_bucket`, split out of
  `warm_polars`).  Adjacent stations share bracketing buckets, so
  per-station tasks would otherwise pile every worker onto the same first
  bucket behind the per-key claim gate.

Measured on a cold-cache `dys_1806_2300kv` lifting-line design
(resolution 30, 8 cores): warm-up now runs at ~7.25 cores, and the whole
cold run drops **944 s -> 622 s** with a bit-identical design result
(same winning camber, attack offset, thrust and torque).  The remaining
serial phases (the best-L/D seeding scan and the post-warm-up design
passes) are unaffected — that is what the remaining TODO progress-bar
items track.

## 2026-08-20 — YAML design summaries and gallery integration

### proply-rs

- Every design run now writes a YAML summary beside the STEP output
  (`<name>.yml` in the output directory, or next to an explicit
  `--step-file`): the design parameters, the motor operating point
  (optimum RPM / torque / power), the performance totals (thrust, torque,
  shaft power, propulsive efficiency `T u_0 / P`, hover figure of merit
  `T^{3/2}/(sqrt(2 rho A) P)`, tip speed) and the hub->tip per-station
  section list (radius, r/R, chord, twist, camber, thickness — plus the
  converged induced velocity and element loads for BEM designs; the
  lifting-line loop reports loads only in the totals, so those fields are
  omitted there instead of written as zeros).
- New `yaml_out` module (`serde_yaml`, new dependency); `Cst::camber()`
  added as the exact inverse of `set_camber` for per-station camber of
  CST designs.

### Gallery

- The Makefile design rule is now a grouped target producing the STEP
  model and the YAML summary together (one design run writes both;
  regenerates if either is missing), `make summaries` regenerates just
  the summaries, and `make clean` removes them.
- `build/out/*.yml` are committed — everything else under `build/` stays
  ignored — and every prop entry in [GALLERY.md](GALLERY.md) links its
  summary next to the parameter JSON.

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
