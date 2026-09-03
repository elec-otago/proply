# Proply: Propeller Design in Rust

This software is a work-in-progress, under heavy development. Please accept my
apologies if things break.

Contact the Author (Tim Molteno) if you have any questions.

`proply-rs` designs propellers automatically: blade element momentum theory
(BEM) sizes each radial blade element, airfoil polars are simulated by
[rust-foil](rust-foil) (a Rust port of XFOIL), and the finished propeller is
written out as a single STEP (AP242) file with true NURBS surfaces — no
OpenSCAD, no STL.

Try the web version (the same pipeline compiled to WebAssembly, running
entirely in your browser tab):
<https://proply-zeta.vercel.app/>

## Repository layout

- `proply-rs/` — the propeller design tool (BEM design loop, STEP output)
- `rust-foil/` — airfoil simulation (Rust port of XFOIL)
- `props/` — JSON propeller parameter files (motor, geometry)
- `images/` — rendered PNG of each prop (see [GALLERY.md](GALLERY.md))
- `legacy/` — the original Python implementation, kept for reference
- [MECHANICAL.md](MECHANICAL.md) — how the mechanical strength calculation
  sizes the airfoil thickness (the `--mech-thickness` design option)

## Building

```sh
cargo build --release
```

## Usage

Specify the prop with a JSON parameter file in `props/` (motor constants,
blade count, radius, chord, ...). See `props/test_prop.json`
for an example.

Design the propeller at the motor's maximum-efficiency operating point:

```sh
cargo run --release -p proply-rs -- --naca --bem --n 40 --element-count 30 \
    --dir=build/out --param='props/test_prop.json'
```

This writes `build/out/<name>.step` — an assembly containing the hub and
`blades` copies of the blade solid — and `build/out/<name>.yml`, a YAML
summary of the design (motor operating point, performance totals, and the
per-station section list).

- Every design converges onto the motor's operating point: the blade
  absorbs the motor's design torque at the design RPM at maximum
  efficiency (the operating point is (torque, RPM); the achieved thrust
  is an output, not an input — `--auto` is implied, and still accepted
  for compatibility).
- Quantity keys in the design JSON may carry unit suffixes as quoted
  strings — lengths in `m`, `cm` or `mm`, pressures in `Pa`, `kPa`, `MPa`
  or `GPa`:

  ```json
      "radius": "6.8cm",
      "tip_chord": "5mm",
      "center_hole": "1.5mm",
      "hub_radius": "6 mm",
      "hub_depth": "6mm",
      "trailing_edge": "0.25mm"
  ```

  A bare number keeps its historical unit — metres for lengths,
  millimetres for `trailing_edge` — so existing files parse exactly as
  before.  `center_hole` may be omitted: it defaults to half the
  `hub_radius`.
- A `motor_torque` + `motor_RPM` pair in the design JSON pins the
  operating point directly (e.g. an engine's rated torque and speed),
  overriding the electric motor model derived from `motor_Kv` and
  `motor_volts`.
- `--step-file <path>` overrides the output STEP file name.
- `--lifting-line [--ar N]` selects the coupled lifting-line / vortex design
  (spanwise-induced losses from the trailed wake instead of the empirical
  tip-loss factor; `--ar` targets a minimum blade aspect ratio). Technique
  described in detail in [proply-rs/README.md](proply-rs/README.md).
- `--mech-thickness` sizes the airfoil thickness mechanically instead of
  with the geometric power law: the blade is treated as a cantilever beam
  anchored at the hub, and the thickness is chosen so the tip deflection
  caused by the blade's own thrust stays within `--deflection-fraction`
  (default 5%) of the radius.  The sizing uses each section's *real*
  bending inertia: the local twist (a pitched section presents the
  z-projection of its chord to the load, so twisted sections need less
  thickness), the actual airfoil shape (a real section is ≈ half as stiff
  as its rectangle), and **camber** (a curved section is stiffer than a
  flat one: 2% camber ≈ 10% thinner sections, 6% ≈ 30%), plus the
  material stiffness via `--modulus` (raise it, e.g. `--modulus=100GPa`
  for carbon, to thin the foils).  How the strength calculation picks the
  section thickness — loads, bending moment, the
  `E·(i_edge·c³t·sin²θ + i_flat·ct³·cos²θ)/12` beam model and its
  parameters — is documented in [MECHANICAL.md](MECHANICAL.md).
- The first run simulates an 80-point alpha sweep per blade station; polars
  are cached in `foil_cache.json` in the working directory, so reruns are
  fast.

`make gallery` designs every prop in `props/` and renders each STEP file to
`images/<name>.png` with FreeCAD (headless, via
[`render-step`](render-step/)); the results are collected in
[GALLERY.md](GALLERY.md), where every prop also links its YAML design
summary (`make summaries` regenerates just those).

The full CLI reference and port notes (deviations from the Python
implementation, deferred features) are in
[proply-rs/README.md](proply-rs/README.md).

## Verification

```sh
cargo test -p proply-rs
```

Golden tests compare the BEM equations and optimizer against numpy/scipy
reference values, and a FreeCAD headless check validates that the STEP
solids are watertight with positive volume.

## TODO

- [ ] **GMSH meshing (`--mesh`)** — not yet ported; the design output is the
      STEP file only.
- [ ] **Version the port-verification tooling** — the reference scripts
      (`build/pyref/`, FreeCAD STEP check) are gitignored and only exist on
      this machine; move them into the repo (e.g. `proply-rs/tests/`) so a
      fresh clone can re-run the checks.
- [ ] **Clean up `images/prop5x3.png`** — leftover asset from the old Python
      README; move it into `legacy/` or delete it.

## Legacy Python version

The original Python `proply` package (scipy optimizer, OpenSCAD/STL output)
lives in [`legacy/`](legacy/). It is not part of the Rust workspace.
