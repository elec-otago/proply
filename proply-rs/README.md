# proply-rs

A Rust port of the Python [proply](https://github.com/elec-otago/proply)
propeller design package.  Blade element momentum theory sizes each radial
blade element, airfoil polars come from [rust-foil](../rust-foil) (a Rust
port of XFOIL), and the finished propeller is written out as a single
STEP (AP242) file with true NURBS surfaces — no OpenSCAD, no STL.

## Status

The NACA path (`--naca --bem`) with STEP output is the shipped design
solver.  A CST (Kulfan) foil family (`--cst` / `"cst": true` in the JSON)
is available alongside it: every station foil is rust-foil's canonical
18-parameter `KulfanParams` section (default AeroSandbox shape),
re-thicknessed and cambered to the same radial laws.  So is the ARA-D
family (`--arad` / `"arad": true`), ported from the legacy proply: the
table-driven ARA-D sections (6, 10, 13 and 20% thick), smoothed and
blended over the radial thickness law (the family carries its own
camber).  Deferred: GMSH meshing (`--mesh`).

## Usage

```sh
cargo run --release -p proply-rs -- --naca --bem --n 40 --resolution 30 \
    --dir=build/out --param='props/test_prop.json'
```

This designs the propeller at the motor's maximum-efficiency operating
point, then writes `build/out/<name>.step` — an assembly containing the
hub and `blades` copies of the blade solid — and `build/out/<name>.yml`,
a YAML summary of the design (excerpt):

```yaml
performance:
  rpm: 13062.659355
  thrust_n: 3.774991
  torque_nm: 0.074784
  shaft_power_w: 102.298489
  figure_of_merit: 0.413492
sections:                        # hub -> tip
- r_m: 0.005
  r_over_R: 0.08
  chord_m: 0.005444
  twist_deg: 62.958388
  camber: 0.0
  thickness_m: 0.002613
```

The summary records the design parameters, the motor operating point
(RPM, torque, power), the performance totals (thrust, torque, shaft
power, propulsive efficiency, hover figure of merit, tip speed) and the
per-station section list (radius, chord, twist, camber, thickness, and
— for BEM designs — the converged induced velocity and element loads).

The design converges onto the motor's operating point: the thrust target
(the JSON `thrust` key) is iterated until the blade absorbs the design
torque at the design RPM.  The inner loops maximise efficiency at a
matched thrust, and the shaft power is fixed once (torque, RPM) is, so
the torque-matched design is the maximum-efficiency design at that
operating point — the achieved thrust is an output.  This replaces the
old `--auto` torque ceiling (`--auto` is still accepted, and implied).

Quantity keys in the design JSON may carry unit suffixes as quoted
strings — lengths in `m`, `cm` or `mm`, thrust in `N`, `kg` or `g`
(kilogram-force):

```json
    "radius": "6.8cm",
    "tip_chord": "5mm",
    "center_hole": "1.5mm",
    "hub_radius": "6 mm",
    "hub_depth": "6mm",
    "trailing_edge": "0.25mm",
    "thrust": "500g"
```

A bare number keeps its historical unit — metres for lengths, newtons
for thrust, millimetres for `trailing_edge` — so every existing file
parses exactly as before.

Design the same propeller with CST (Kulfan) station foils instead of the
NACA 4-series:

```sh
cargo run --release -p proply-rs -- --cst --bem --n 40 --resolution 30 \
    --dir=build/out --param='props/test_prop.json'
```

Design the same propeller with ARA-D station foils (the legacy proply's
table-driven family, blended over the thickness law):

```sh
cargo run --release -p proply-rs -- --arad --bem --n 40 --resolution 30 \
    --dir=build/out --param='props/test_prop.json'
```

The hub includes the `center_hole` mounting bore as a single closed shell:
outer and bore cylinders plus two annular caps (ring faces with an inner
bound), written so every circular edge is split into two semicircular arcs
shared by exactly two faces.  When the JSON omits `center_hole` it defaults
to half the `hub_radius`.  FreeCAD's OCCT importer reads the solid back
as watertight (`freecadcmd` + `Part.read` reports valid, closed solids).

The first run simulates a polar (an 80-point alpha sweep) for each blade
station; polars are cached in `foil_cache.json` in the working directory,
so reruns are fast.

## WebAssembly (browser build)

The whole design pipeline — BEM / lifting-line, rust-foil polars, the STEP
and YAML writers — also compiles to `wasm32-unknown-unknown` and runs
entirely in a browser tab; nothing is uploaded anywhere.  The wasm-only
bindings (`src/wasm.rs`, compiled out of native builds) wrap the same
`pipeline::run_design` the CLI uses.

```sh
# from the workspace root
make wasm            # wasm-pack build proply-rs --target web --release
                     # (the package lands in proply-rs/web/pkg/)
make check-wasm      # plain cargo check for the wasm32 target

# serve the demo page (the ES-module import needs http://, not file://)
python3 -m http.server -d proply-rs/web 8000
# then open http://localhost:8000/
```

The demo page edits the same JSON design parameters as the CLI's
`--param` file.  "Plate polars" (checked by default) uses the analytic
flat-plate model for a quick run; unchecking it computes real XFOIL
polars in-wasm, which is much slower on the single-threaded main thread.

### Deploying the demo (GitHub + Vercel)

`web/` is a self-contained static site once the wasm package is built
into it, so Vercel can serve it straight from the GitHub repository with
no build step:

1. Build the package and commit it (the one generated artifact that is
   meant to be committed):

   ```sh
   make wasm
   git add proply-rs/web/pkg
   git commit -m "Rebuild the browser wasm package"
   ```

2. Push the repository to GitHub.
3. In Vercel: *Add New → Project → Import* the repository, then set

   | Setting | Value |
   |---|---|
   | Root Directory | `proply-rs/web` |
   | Framework Preset | Other |
   | Build Command | *(leave empty)* |
   | Output Directory | *(default)* |

   and deploy.

Each push deploys whatever `web/pkg` is committed, so publishing a new
build is just `make wasm` plus a commit.  To keep build artifacts out of
the branch, move the build into CI instead — a GitHub Actions workflow
that runs `make wasm` and force-pushes `proply-rs/web/` to a `deploy`
branch, with the Vercel project pointed at that branch, works the same
way.

### Polar cache in the browser

There is no filesystem in the browser, so the polar cache lives in a
long-lived WASM session (`PropSession`) hydrated from IndexedDB at page
load; after each design, only the freshly simulated polars are handed
back to the page and written to IndexedDB.  The cache is a pure
performance artifact — a blocked or cleared database just means a cold
run.  (`PropSession::cache_to_json` / `hydrate_json` provide a whole-cache
export/import escape hatch.)

## Lifting-line analysis (`--lifting-line`)

The default `--bem` path is an **annular blade-element momentum theory**: each
radial station is sized independently, and blade finiteness enters only
through an empirical Prandtl tip/hub-loss factor. That treats adjacent
annuli as decoupled, so the finite-span (aspect-ratio) induced loss — the
drag that a real, load-carrying blade pays for shedding a trailed wake — is
only approximated, never computed.

`--lifting-line` replaces that with a **coupled lifting-line / vortex
model** (`src/lift_line.rs`). It resolves the spanwise coupling directly, so
induced loss emerges from the vortex system rather than a fitted factor.

### The model

The blade is treated as a **radial lifting line** carrying a piecewise-
constant bound circulation, one value `Γ(r)` per station. Because `Γ` changes
across each station boundary, helical trailed vorticity is shed there; after
azimuthal averaging over the `B` blades this is the classic **Trefftz-plane /
Goldstein** vortex system of coaxial rings. The induced velocity the blade
actually feels at radius `r` is then:

- **axial** — the Biot–Savart sum over every trailed ring at radius `a`:
  `u_i(r) ∝ B Σ_j ΔΓ_j K(r, a_j)`, where `K` is the axial velocity in the
  plane of a unit ring (computed by straight-segment Biot–Savart, core-
  regularised near the filament), and the rotor-plane constant converts the
  trailed circulation to the disk inflow (half the fully-developed
  slipstream value — the Glauert/actuator-disk limit);
- **tangential (swirl)** — `v_i(r) = B·Γ(r) / (4 π r)`, the classical swirl
  that opposes the rotor and reduces the relative tangential speed.

These are **coupled**: `u_i` and `v_i` at a station come from *all* stations,
so a heavily-loaded blade (large `Γ`) induces more inflow everywhere, and a
long, lightly-loaded (high-aspect-ratio) blade sheds a weaker, more distributed
wake. `Γ` is made self-consistent with the airfoil polar by iterating
`Γ = ½ c V cl(α)` to convergence (under-relaxed), with the local flow angle
`φ = atan((u_0 + u_i)/(Ωr − v_i))` and `α = θ − φ`.

### Forces

Thrust and torque follow **Kutta–Joukowski + profile drag** per station
(`dT = (ρ Γ u − ½ ρ V c cd u) dr`, `dQ = r (ρ Γ v + ½ ρ V c cd v) dr`),
i.e. lift from `ρ V Γ` and drag anti-parallel to the relative velocity,
resolved along the axial (`u`) and tangential (`v`) directions. Angle of
attack is bounded below the deep-stall region, where the fixed-point
iteration is unreliable.

### Design loop

`Prop::lift_line_design(rpm, thrust, ar)` sizes the blade against this
model:

1. each element runs at its **best-L/D angle** (`argmax cl/cd` on the
   polar — minimum drag for the induced lift it carries) plus a common
   work offset `da`;
2. the local angle of attack is *prescribed* (`alpha = best-L/D + da`), so
   the twist is `phi + alpha` and the solve has no twist↔induction feedback
   to limit-cycle;
3. **smooth chord / per-control optimization**: the chord is a
   **shape-preserving cubic (PCHIP) spline** through `chord_spline_n`
   control points (default 3; `"chord_spline_n": N` in the JSON or
   `--chord-spline-n N`), at radii spread hub→tip.  Each control value is a
   design variable (bounded by the geometrically-allowed chord — the
   `tip_chord·R²/r²` taper capped by blade spacing) and the spline gives a
   kink-free smooth chord at every radius.  The outer level is a
   **NelderMead** over those N control values, seeded at the full chord (so
   it starts at the thrust-capable geometry); for each candidate shape an
   inner **monotone `da` bisection** reliably matches the thrust target
   (`alpha = best-L/D + da`), so the outer only minimises the torque —
   i.e. it finds the most efficient chord *shape* that meets the required
   thrust with no oscillation.  `--ar N` caps the control values (a
   minimum blade aspect ratio), thinning the blade as `N` rises.

Because the wake is coupled, the reported torque *includes* the induced
loss — so the design can make real efficiency trade-offs (induced vs
profile) instead of the empirical-tip-loss shortcut.

### Usage

```sh
cargo run --release -p proply-rs -- --bem --lifting-line [--ar 5] [--chord-spline-n 5] \
    --dir=build/out --param='props/test_prop.json'
```

### Verification

- `biot_savart_ring_matches_closed_form` — straight-segment Biot–Savart on a
  vortex ring reproduces `Γ/(2a)` on-axis.
- `ring_kernel_on_axis_is_one_over_2a` — the trailed-ring kernel recovers
  `1/(2a)` on its own axis.
- `swirl_is_b_gamma_over_4pi_r` — the swirl relation `B·Γ/(4πr)` is exact.
- `lift_line_design_finite_and_thrust_matching` — the coupled solve is
  bounded, produces positive thrust and converges to a stable design.

### Status / limitations

The azimuthally-averaged lifting line couples the radial stations through
the trailed ring system, and the circulation is solved with a **damped
Newton** method (dense analytic Jacobian of the residual
`R_i = Gamma_i − ½ c_i V_i cl(alpha_i)` plus an Armijo backtracking line
search). That converges cleanly to a genuine fixed point (residual
`‖R‖ < 1e-6`, unit-tested) even where a plain Gauss–Seidel fixed point
limit-cycles on the nonlinear stall polars — an earlier iteration of the
solver diverged on the same case and now converges monotonically.

Remaining caveats: the design loop caps the attack-angle offset
(`da ≤ 0.12 rad`) so the blade stays out of deep stall; thrust rises with
`da` only until stall, so a small, lightly-loaded prop at a fixed RPM has a
stall-limited maximum thrust regardless of chord. The chord-scale `s` is
searched up to `s_cap` (the geometric max, or the `--ar`-capped value), and
if the target needs more than the full chord at the stall cap the search
returns the best reachable design. The physics and the solver are validated;
the stall bound is an engineering choice you can relax. The default `--bem`
path remains the shipped design solver; `--lifting-line` is the coupled
vortex alternative described above (see `props/test_prop.json` for a
physically realizable configuration).

## Port notes

- Angle conventions (the design twist of each blade element).  With
  `u_0` the freestream axial velocity, `dv` the axial induction at the
  disc, `a'` the tangential induction, and `omega` the angular velocity,
  the local inflow angle measured from the rotor plane is
  `phi = atan((u_0 + dv) / (omega·r·(1 − a')))` and the angle of attack is
  `alpha = theta − phi`, with `theta` the blade pitch angle from the rotor
  plane.  The design loop solves `theta` per station for the target
  induced velocity at maximum efficiency, bounded to
  `[phi − 8°, phi + 15°]` around the ideal inflow angle (so the design
  angle of attack stays in `[−8°, +15°]`).  Positive `theta` rotates the
  profile around its 0.67·chord point with the leading edge toward −z,
  consistent with the positive-alpha → positive-CL convention of the
  polar model.  Momentum theory is used in its classical second-order
  form: the thrust uses `a(1 + a)` and the torque carries the
  `(1 + a)` axial coupling (`d_t`/`d_m` in `optimize.rs`).
- The Buhl (2005) turbulent-wake `CT(a)` relation (NREL/TP-500-36834,
  Eqs. 1 + 18) is implemented as `ct_buhl`/`a_buhl` in `optimize.rs`
  (golden-tested).  It is *not* wired into the design loop: it is
  published in the decelerating-disk (wind-turbine) convention for
  `a ∈ (0.4, 1)`, whereas the accelerating propeller state
  (`u_disc = u_0 + dv`, `CT = 4a(1 + a)`) has no momentum-theory
  breakdown — the canonical propeller code (XROTOR) likewise solves the
  exact momentum balance without a Glauert-type correction.  The
  functions are provided for brake-state / turbulent-wake momentum-model
  work.
- `optimize.rs` reproduces the Python BEM equations 1:1 (verified against
  numpy reference values in `tests/golden.rs`).  The scipy SLSQP/COBYLA
  optimizer is replaced by a box-constrained Nelder-Mead (with a quadratic
  penalty; all proply constraints are simple variable bounds).  Golden
  tests compare the optimizer result against scipy SLSQP on the same
  objectives.
- `foil_simulator.rs` replicates the polar pipeline: Reynolds number
  snapped to `np.geomspace(3e4, 2e6, 20)`, degree-9 polynomial fit of
  cl/cd over alpha in radians, flat-plate fallback outside the envelope.
  One deliberate deviation from the Python: the airfoil's trailing edge is
  **closed** for the polar solve (the ~0.25 mm design gap has negligible
  aerodynamic effect but destabilizes the boundary-layer convergence below
  Re ≈ 2e5).  The blade geometry keeps the gap.  (Coordinate normalization
  is handled by rust-foil itself — its `set_airfoil` normalizes by default,
  matching canonical XFOIL.)
- The STEP writer (`step_out.rs`) builds the blade as a single watertight
  NURBS solid: cubic B-splines interpolate each station profile (Piegl &
  Tiller A9.1 with averaged knots), and the upper/lower/TE-cap/end-cap
  surfaces are degree-1 ruled lofts.  Edges lie exactly on their faces.
  Verified by round-trip parsing (step-io) and FreeCAD (`Part.read`:
  valid, positive-volume solids).
- Interactive matplotlib plots in the Python are replaced by a console
  design table.

## Verification

```sh
cargo test -p proply-rs                          # unit + golden tests
build/venv/bin/python proply-rs/tests/golden/gen_golden.py # regenerate golden values
freecadcmd build/pyref/check_step.py <file.step> # FreeCAD solid check
build/venv/bin/python build/pyref/run_plate.py props/test_prop.json 30
# Python design-loop reference using the flat-plate polar model
```
