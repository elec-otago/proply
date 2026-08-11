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

The hub is currently a plain cylinder: the Python's `center_hole` mounting
bore is a trivial pocket operation in any CAD package, and step-io's
multi-shell output (rings with inner bounds / `BREP_WITH_VOIDS`) is not
read back correctly by FreeCAD's OCCT importer.

The first run simulates a polar (an 80-point alpha sweep) for each blade
station; polars are cached in `foil_cache.json` in the working directory,
so reruns are fast.

## Port notes

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
