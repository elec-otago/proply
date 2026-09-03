# CHANGES
## 2026-09-03 — Reject unattached flow states (extreme induced inflow)

### proply-rs

A vision check of the gallery renders found blades with extreme wiggles
whose exported twists (phi + alpha) ballooned to 66-76 deg at mid span.
All three cases trace to the discrete wake model finding non-physical
roots outside its small-perturbation validity:

- flywoo_robo_rb1202.5 (20 mm hover rotor, 41 krpm, u0 = 0): induced
  axial velocity several times the local rotational speed mid-blade;
- dys_2814_910kv: a single-station circulation spike (gamma 7.6 vs ~0.8
  on its neighbours);
- turnigy_CA_120's mechanical re-design: a fully diverged solve (ui
  ratio ~8600 at one station).

A converged circulation is now usable only if, outside the inner hub
zone, the induced axial velocity stays within ~1.0 of u0 + w.r (checked
at every adoption point in the pass and in the cold export re-match)
and no station carries a circulation spike (largest within 4x of the
second largest).  Rejected states cannot win the candidate competition
or become the exported blade.  When the mechanical-thickness re-design
finds no usable state at all, the pipeline keeps the first design and
says so in the warning instead of reporting a blade that was never
built.  flywoo now exports a smooth design (T = 0.590 N at Q = 0.0023
Nm, cold-verified, err 0.00%).
df6c7a3
## 2026-09-03 — Reject spurious single-station circulation spikes

### proply-rs

A `make gallery` review found several blades rendered non-smooth.  The
cold re-solve that verifies the exported blade can settle on a
discrete-solve root where one station carries a circulation spike —
measured at gamma 7.6 vs ~0.8 on its neighbours (a 9x spike, on
low-aspect-ratio blades at 30+ stations) — with ui = 45 m/s there.
The spike is not a physical loading: it distorts the induced inflow
over several stations, so the exported twist (phi + alpha) develops a
large local excursion (dys_2814's max twist curvature went from 6.2 deg
in the old warm-exported design to 29.7 deg).

Spike states are now rejected wherever the brake branches are: a sample
is physical only if every station lifts (gamma >= 0) *and* the largest
circulation is within 4x of the second largest (smooth loadings peak at
~1.5x).  When the cold re-match at the target torque cannot find a
smooth physical state, the export falls back to the warm state the
optimizer matched (its twists and its operating point) instead of
exporting the spiky cold branch.  dys_2814 now exports the warm state:
T = 12.95 N at Q = 0.1304 Nm (err 0.01%), no spike.
de7dc35
## 2026-09-03 — Journal the polar cache (per-polar durability without whole-file rewrites)

### proply-rs

The per-polar checkpoint — every freshly simulated polar persisted the
moment it exists — rewrote the entire pretty-printed `foil_cache.json`
on each insert.  Once the cache reaches tens of megabytes that rewrite
dominated the design wall clock: a `make gallery` run issued ~950 MB of
disk writeback per 30 s (one ~86 MB rewrite per polar) against a device
sustaining ~32 MB/s — ~2.5 s of serialized writeback per polar, on the
global store lock, and worsening as the cache grows prop by prop.

Path-backed stores now append each new polar (a good sweep or a
degenerate failure marker) as one NDJSON line to `foil_cache.json.journal`
— still written at calculation time, but O(one polar) — and rewrite the
whole file only every 200 inserts or at save/exit, via a temp file +
atomic rename, then dropping the journal.  Loading replays the file then
the journal, skipping a torn tail and tolerating duplicate records.  The
in-memory wasm store is unchanged (its per-polar IndexedDB hook already
writes O(one polar)).

A full dji_phantom3 design (element-count 30, base + mechanical phases)
that was on track for ~25+ min completed in 191 s with 45 sweeps and no
writeback stall.
6549325
## 2026-09-03 — Design directly onto the motor operating point (torque, RPM)

### proply-rs

The operating point was already (torque, RPM) — the motor's torque at the
design RPM — with the JSON `thrust` used only as an iteration seed for an
outer fixed point that ran a *full redesign* per sample.  That response
was discontinuous: consecutive re-optimisations could land on different
flow branches (honda@30 sampled Q = 1.63 and Q = 0.067 N m at the same
thrust on successive matches — the observed torque cliff), so the damped
iteration could not converge and needed ever more machinery
(cross-match seeding, adaptive damping, restart pruning).

The lifting-line path now reaches the operating point in a single design
pass: each geometry evaluation matches a common attack offset `da` so
the absorbed torque equals the target (a monotone bisection, the same
shape as the old thrust match) and the chord/camber search maximises the
resulting thrust — maximum efficiency, since the shaft power Q·ω is
fixed by (torque, RPM).

- Best-L/D seed angles are floored at zero lift: the polar fits produce
  spurious negative-alpha optima on thin low-Re outboard sections (the
  m=0.04 tip came out at −5.4°), prescribing a negative-lift tip brake
  that the torque-matched objective would exploit.
- Circulation samples with any negative station gamma are rejected
  everywhere — brake-tip flow branches never seed the search.
- The exported blade is cold-verified: `da` is re-matched with cold
  solves on the winning geometry, so the reported operating point is
  exactly what an independent solve of the exported blade measures.
- The BEM path keeps its damped iteration (per-station momentum design
  has no global offset); it now seeds itself from a momentum-theory
  estimate of the thrust the demanded torque can produce.

Results at element-count 30 (which used to cliff): honda_gx35 converges
in one pass to Q = 1.6000 N m cold-verified, T = 39.9 N (the old loop
never converged; its fallback sat 3.9% off target at 35.0 N);
multistar_2209 and dys_2814_910kv likewise converge to their torque
targets within 0.03-0.00% in one pass.
30b6e1b

## 2026-09-03 — Remove the thrust target from the JSON input format and the web UI

### proply-rs / web

The achieved thrust is an output of the design, so the JSON `thrust` key
is no longer a design input: the field is removed from the parameters
(the N / kgf / gf unit parsing goes with it), the web demo's Thrust
field and default are gone, and the prop JSON files and docs follow the
new schema.  A legacy `thrust` key in an old design file is still
tolerated — it is ignored — and covered by a test.
881d222

## 2026-09-03 — Reject garbage circulation solves in the design competition

### proply-rs

Found auditing a 7-hour `make gallery` run: some props' torque matches
blew up with errors of 784% and 6144%, stalling at **negative absorbed
torque** (Q = −10.95 and −96.7 N m).  A warm-started Newton solve can
land on a garbage branch with negative torque whose thrust still happens
to match the target, and two adoption points trusted it:

- the warm da-acceptance (`err <= 3%`) returned any match, and
- the tight da-refinement accepted any bisection sample with a smaller
  thrust error.

Because the pass objective is `f = q + 50·err`, the garbage *won* the
candidate competition (negative q minimises f) and derailed the outer
torque match.  All adoption points — the warm acceptance, the failed-warm
best seeding, and the refinement — now require a physical result
(`q > 0`, finite).  The failing honda config no longer produces garbage
or stalls, and its torque match proceeds on sane designs.

## 2026-09-02 — Cache failed polar sweeps (the seeding warm-up)

### proply-rs

- The seeding warm-up (and any other first-touch sweep) on a fresh prop
  used to *fail* on roughly half its calculations — rust-foil returns
  zero converged points for thick mechanical-law root sections at moderate
  Reynolds numbers — and a failed sweep stored nothing, so every fresh
  run re-attempted each doomed bucket once (and an interrupted run lost
  even its successes before per-polar persistence).  A failed sweep now
  stores a **degenerate marker** (an empty polar) through the normal
  insert path, so it is persisted the moment it happens; markers never
  pass the polar checks, and any later fetch — the same run or a future
  one — treats the key as degenerate (flat-plate fallback) without
  re-sweeping.  The seeding warm-up now caches *every* outcome: run A of
  a cold prop re-swept only the 31 previously-failed keys (storing 31
  markers), and run B's seeding warm-up completed with **zero** sweeps.

## 2026-09-02 — Persist every polar the moment it is calculated

### proply-rs

- The polar cache is now written to disk with **each** completed
  rust-foil sweep (`PolarStore::insert` saves immediately), not in
  batches: an interrupted run loses at most the sweep in flight.  The
  whole-file JSON write takes milliseconds next to the seconds each sweep
  takes.  Verified by killing a discovery run mid-design: each finished
  calculation was already on disk (every polar appears the moment it is
  calculated).

## 2026-09-02 — Checkpoint the polar cache during a run

### proply-rs

- The polar cache used to be written to disk only once, at process exit —
  a killed or interrupted run lost every polar it had simulated (which is
  why interrupted `make` runs re-paid their 10-40 minute discovery cost
  each time).  `PolarStore::insert` now auto-saves after every 25 freshly
  simulated polars (`CHECKPOINT_EVERY`); a whole-file JSON write takes
  milliseconds next to the seconds each rust-foil sweep takes, so the
  steady-state cost is invisible.  Verified by killing a discovery run
  mid-design: the on-disk cache kept exactly the checkpointed batches
  (2500 → 2550 keys, two batches of 25), and a re-run reused them.

## 2026-09-02 — Rename `resolution` to `element_count`

### proply-rs CLI / JSON schema

- The `--resolution` option was described as a "radial resolution in
  millimetres", but it is actually the *number* of radial blade elements
  (spanwise stations) the design is solved on — each element covers
  `(radius - hub_radius)/N` of the span, and the station count is
  `N * radius / (radius - hub_radius)` (≈ N for small hubs).  It is
  renamed **`--element-count <N>`** (JSON key **`element_count`**), with
  the old `--resolution` flag and `resolution` JSON key kept as accepted
  aliases so existing design files and scripts keep working.
- The startup banner now reports what it means ("Blade elements: 10
  (each covers 25.20 mm of span)") and the design summary's
  `spanwise_resolution` field is `element_count`.  The CLI help for `--n`
  now states it is the STEP loft's chordwise sampling (geometry only).

## 2026-09-02 — Stop Nelder-Mead at the objective's noise floor

### proply-rs

- `optimize::NelderMead` now stops when the simplex shrinks into a value
  plateau: each time the simplex's extent halves, the search checks
  whether the function spread fell with it.  On a smooth objective the
  spread falls with the size (quadratically near a minimum), so nothing
  changes; on an objective with a deterministic roughness floor — the
  design objective's da-matching branches and its 50x error amplification
  leave the spread pinned at ~1e-3 once the simplex is small — the spread
  stops falling while the size keeps halving towards the fine `xatol`,
  and the search stops after four such halving events instead of grinding
  the noise floor.  Configurable via `plateau_halvings` (default 4).

Measured on the honda design (real polars, warm cache): the same design
bit-for-bit, with the Nelder-Mead evaluations dropping from ~19,000 to
~9,250 per design (and the repeated-evaluation tail from ~6,000 to
~1,800).  The Rosenbrock and BEM-solver tests still converge to their
tight tolerances — the detector only fires on true plateaus.

## 2026-09-02 — Converge the torque match (lifting-line design)

### proply-rs

Analysis of the design logs showed three convergence inefficiencies —
measured on the 268 mm honda design (real polars), the torque-match
iteration oscillated through 16 matches with ~46,700 circulation solves
(~half of them repeating the previous evaluation), and the plate-mode
run through 18:

- **Candidates now compete at the exact target thrust.**  The pass's
  warm `da` acceptance (thrust within 3%) let a high-torque design that
  happened to match thrust exactly outrank a low-torque one sitting a
  percent or two off target — the mis-ranking that made the outer torque
  iteration oscillate (both code versions show 15-20% error jumps when
  such a candidate won).  The final measured candidate refines `da` with
  a bounded bisection around the warm value (branch-safe: the matched
  thrust is monotone in `da`, every solve warm-starts from the previous
  circulation), so candidates are ranked on the torque at the target.
- **The previous match's design re-competes warm.**  Its chord controls,
  camber distribution and attack angles are kept on the `Prop` and the
  next match runs them as the "prev" incumbent, warm-started from its
  own controls (with the full-chord run kept as the reachability anchor,
  and warm-chained scaled restarts that stop as soon as they cannot beat
  the running best).  Consecutive matches stay on one geometry branch
  instead of re-optimizing from scratch and landing on arbitrary local
  optima.
- **The torque-match update adapts to the measured response.**  The
  fixed 0.8 damping exponent assumed an elasticity of ~1.25; the design
  logs show a thrust→torque elasticity of 2-5 with the mechanical
  thickness law (linearly unstable under fixed damping), so the exponent
  is now reduced towards 1/(measured elasticity), and every step is
  capped at ±25%.

Result (honda, real polars, warm cache): 16 → 6 torque matches,
~46,700 → ~19,000 circulation solves, 4.5 s → 1.8 s, with a monotone
error decrease (no 15-20% jumps) — and a genuinely better design, since
the old winner ranking was distorted: figure of merit 0.366 → 0.370,
propulsive efficiency 0.864 → 0.870, torque target met to 0.10% (was
0.45%).  Same-state runs remain bit-identical.

## 2026-09-02 — Speed up lifting-line designs (degenerate polars, hot path)

### proply-rs

Measured on the 268 mm honda design (real polars, lifting line): the pure
design math is ~1 s, everything else was the polar machinery.  Three
changes:

- **Don't re-simulate degenerate polars.**  rust-foil's viscous solve
  fails *deterministically* on some (foil, Reynolds) targets — too few
  converged points, or an unphysical near-zero-drag sweep — yet every
  fresh design pass (new simulators, empty fit caches) re-ran the doomed
  sweep for each stubborn bucket: 29 of a 348-key warm honda cache were
  degenerate and a re-run spent ~90% of its time re-simulating them.  A
  key proven degenerate is now remembered for the session and falls
  straight back to the flat-plate model (the value the doomed sweep would
  have produced), and a *stored* degenerate polar is trusted without
  re-sweeping.
- **Persist the degenerate-key markers.**  The CLI saves them to a sidecar
  (`foil_cache.json.bad.json`) at exit and loads it at startup, so keys
  proven degenerate in one run — or one prop of a `make` sweep — are not
  swept again by the next (rust-foil buckets that store *nothing* would
  otherwise be re-swept on every run).  Delete the sidecar to retry them.
- **Allocation-free polar lookup path.**  The hot `get_cl`/`get_cd` path
  now blends fixed-size `[f64; 10]` fit arrays and evaluates them with no
  per-call heap allocation: the per-simulator fit cache is keyed by the
  Reynolds bucket alone (the foil is fixed per simulator — no string
  hashing), the Reynolds grid is computed once instead of per call, and
  the bucket bracket is a single logarithm instead of a `powf` scan.
  Same coefficient-blend arithmetic throughout, so polar values are
  unchanged.  A cached lookup drops from ~1.3 µs to ~130 ns (10x).

The fully-warm honda design goes from ~28 min (cold, original code; a
warm re-run never finished) to ~4.5 s, with bit-identical outputs
verified across the code versions.

## 2026-09-02 — Design trace log (`--log <file>`)

### proply-rs CLI

- New `--log <FILE>` output option (equally settable as the `log` key in
  the JSON design file): every design-loop line the console prints is
  mirrored to FILE as it prints — the operating-point torque match,
  per-station BEM state, the lifting-line camber scan, warnings and the
  final totals — giving a detailed, durable record for debugging
  convergence issues (a failed run's trace ends with the failure reason).
- New `design_log` module holds the thread-safe trace sink and the
  `dprintln!` / `deprintln!` tee macros; the per-station BEM trace now
  also reports the solver residual (`bem_err`).  Nothing is written when
  no log is requested, so the WebAssembly build (no filesystem) and the
  library tests are unaffected.

## 2026-09-02 — Fix the 3D preview for large props (far-plane clip)

### proply-rs web demo

- `web/viewer.js` now calls `camera.updateProjectionMatrix()` after
  framing the STEP model.  The near/far planes were set per-model
  (`radius / 100` .. `radius * 100`) but the projection was never
  rebuilt, so the camera kept its initial `far = 1000`: any prop whose
  framing distance exceeded 1000 mm-units (roughly a 200 mm+ radius) was
  clipped at the far plane and the view rendered empty/dark.  The 3D
  preview now frames and lights a 280 mm prop correctly.  JS-only
  change; no Rust or wasm change.

## 2026-08-30 — Build label on the web demo

### proply-rs web demo

The demo page now shows a build label in the form `yyyy-mm-dd.xx` — the
date of the last commit and its **per-day build number** (the count of
commits on that date, so each build on a day increments `xx` and a new
day starts at `.01`):

- new `make build-date` (Makefile) writes `proply-rs/web/build.js`
  (`export const BUILD = "2026-08-30.08"`), derived from git
  (`git log -1 --format=%cs` for the date, the matching commit count for
  `xx`);
- `web/index.html` gains a `#build-line` and `web/main.js` imports the
  stamp and renders `build 2026-08-30.08` under the intro;
- AGENTS.md notes the stamp must be regenerated (`make build-date`)
  before committing web changes, so the deployed label matches the
  deployed sources.  JS/docs only; no Rust or wasm change.

## 2026-08-30 — Mechanical tab in the web demo; accuracy pass on MECHANICAL.md

### proply-rs web demo

- The web demo gains a **Mechanical** tab holding *all four* mechanical
  parameters — the **Mechanical thickness (beam sizing)** checkbox
  (`mech_thickness`), **Elastic modulus (GPa)** (`modulus`), **Deflection
  limit (fraction of R)** (`deflection_fraction`) and **Minimum thickness
  (fraction of chord)** (`thickness_floor`) — replacing the two fields
  that used to live in Propeller Specifications.  The default design now
  carries all four keys so the tab is populated.  JS-only change (the
  committed wasm already parses the keys).

### docs

- [MECHANICAL.md](MECHANICAL.md) accuracy pass against the current code:
  the geometric-law comparison is corrected (a blade built to the power
  law deflects ≈ 2.5 mm under the same twist/shape-aware beam model —
  its root `t/c ≈ 0.71` overshoots the deflection budget at roughly four
  times the mechanical design's root thickness — not 4.9 mm from the
  early rectangle model), the 100 GPa modulus example is corrected to
  land on the `thickness_floor`, and the web-demo instructions describe
  the new Mechanical tab.

## 2026-08-30 — Camber in the mechanical thickness law (real section inertia)

### proply-rs design loop

The mechanical thickness law now sizes with each section's **real
bending inertia** instead of the enclosing rectangle: `i_flat` and
`i_edge` are the foil's actual second moments about the chord line and
the thickness axis (relative to `c t³/12` and `t c³/12`), computed from
the shape points of every family (NACA 4-series, CST, ARA-D) via
polygon mass moments (new `FoilLike::section_shape_factors`, `foil.rs`).
The twisted-section inertia becomes
`I = i_flat·c t³/12·cos²θ + i_edge·c³ t/12·sin²θ` and the sizing cubic
`i_flat·cos²θ·t³ + i_edge·c²·sin²θ·t = 6ML²/(Ecδ)`.

- **A cambered (curved) section is stiffer than a flat one** (the user's
  point): the chord-line moment contains the camber — `i_flat ≈ 0.47`
  symmetric, ≈ 0.62 at 2% camber, ≈ 1.35 at 6% camber (12–15%
  thickness), so cambered sections size ≈ 10–30% thinner.  A real
  symmetric section, however, is only ≈ 0.47× as stiff as its
  rectangle, so the absolute sections also thicken ≈ 29% versus the
  rectangle model.
- On the web-default design (BEM, plate, 3 GPa): root `t/c ≈ 0.17` at
  zero camber (0.71 geometric, 0.08 rectangle-mechanical); with 6%
  camber the inboard sections drop to the `thickness_floor` (root
  `t/c = 0.06`, predicted deflection 2.77 mm of the 3.4 mm allowed).
- Docs (`MECHANICAL.md`, `README.md`) describe the real-section model,
  numbers and caveats; the wasm package is rebuilt.

## 2026-08-30 — Twist-aware mechanical thickness law

### proply-rs design loop

The mechanical thickness law now takes the local twist into account.
A foil section twisted by `θ` about the spanwise axis resists the
z-directed thrust with the second moment of the *rotated* rectangle,
`I = c³t/12·sin²θ + ct³/12·cos²θ` — the summed squares of the
z-projections of the section's chord (`c·sin θ`) and thickness
(`t·cos θ`).  Since the chord is much longer than the thickness, the
chord's projection dominates at high twist: the sized thickness "tends
to the chord" (the section's z-extent `c·sin θ + t·cos θ`), exactly as
the TODO's "take the twist into account, as the chord can make the beam
stiffer" intended.  The sizing becomes the cubic
`t³·cos²θ + t·c²·sin²θ = 6ML²/(Ecδ)` (closed form `t = (K)^(1/3)` only
for untwisted sections), solved per station with the element's own
twist (`thickness.rs`, `prop.rs`).  On the web-default design this
halves the root requirement again: root `t/c` ≈ 0.08 at 3 GPa (0.28
untwisted, 0.71 geometric), most inboard sections ride the 0.06 floor,
and the predicted tip deflection closes at 3.00 mm of the 3.4 mm
allowed.  Docs (`MECHANICAL.md`, `README.md`) describe the
twist-aware model; the wasm package is rebuilt.

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
