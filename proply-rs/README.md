# proply-rs

A Rust port of the Python [proply](https://github.com/elec-otago/proply)
propeller design package.  Blade element momentum theory sizes each radial
blade element, airfoil polars come from [rust-foil](../rust-foil) (a Rust
port of XFOIL), and the finished propeller is written out as a single
STEP (AP242) file with true NURBS surfaces — no OpenSCAD, no STL.

## Status

First milestone: the NACA path (`--naca --bem`) with STEP output.
Deferred: ARA-D foils (`--arad`), GMSH meshing (`--mesh`).

## Usage

```sh
cargo run --release -p proply-rs -- --naca --bem --n 40 --resolution 30 \
    --dir=build/out --param='props/test_prop.json'
```

This designs the propeller at the motor's maximum-efficiency operating
point, then writes `build/out/<name>.step` — an assembly containing the
hub and `blades` copies of the blade solid.  `--auto` re-runs the design
loop, reducing the thrust target until the torque drops below
`1.5 × Qmax`.

The hub includes the `center_hole` mounting bore as a single closed shell:
outer and bore cylinders plus two annular caps (ring faces with an inner
bound), written so every circular edge is split into two semicircular arcs
shared by exactly two faces.  FreeCAD's OCCT importer reads the solid back
as watertight (`freecadcmd` + `Part.read` reports valid, closed solids).

The first run simulates a polar (an 80-point alpha sweep) for each blade
station; polars are cached in `foil_cache.json` in the working directory,
so reruns are fast.

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
3. **smooth chord**: the chord is a **shape-preserving cubic (PCHIP)
   spline** through `chord_spline_n` control points (default 3; override
   with `"chord_spline_n": N` in the JSON, or `--chord-spline-n N`), at
   radii spread hub→tip.  The control values hold the geometrically-allowed
   chord (the `tip_chord·R²/r²` taper capped by the blade-spacing limit) at
   each reference radius; interpolating them as a cubic spline gives a
   kink-free smooth chord at every radius.  The design sweeps a single
   scale `s` of that smooth spline and, inside each `s`, brackets and
   bisects the single `da` to the thrust target, then keeps the `(s, da)`
   with the lowest torque — pure, independent, exactly-solved evaluations
   with no feedback loop to oscillate in.  `--ar N` caps `s` (a minimum
   blade aspect ratio), thinning the blade as `N` rises.

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
