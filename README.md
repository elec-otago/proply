# Proply: Propeller Design in Rust

This software is a work-in-progress, under heavy development. Please accept my
apologies if things break.

Contact the Author (Tim Molteno) if you have any questions.

`proply-rs` designs propellers automatically: blade element momentum theory
(BEM) sizes each radial blade element, airfoil polars are simulated by
[rust-foil](rust-foil) (a Rust port of XFOIL), and the finished propeller is
written out as a single STEP (AP242) file with true NURBS surfaces — no
OpenSCAD, no STL.

## Repository layout

- `proply-rs/` — the propeller design tool (BEM design loop, STEP output)
- `rust-foil/` — airfoil simulation (Rust port of XFOIL)
- `props/` — JSON propeller parameter files (motor, geometry, thrust target)
- `legacy/` — the original Python implementation, kept for reference

## Building

```sh
cargo build --release
```

## Usage

Specify the prop with a JSON parameter file in `props/` (motor constants,
blade count, radius, chord, thrust target, ...). See `props/test_prop.json`
for an example.

Design the propeller at the motor's maximum-efficiency operating point:

```sh
cargo run --release -p proply-rs -- --naca --bem --n 40 --resolution 30 \
    --dir=build/out --param='props/test_prop.json'
```

This writes `build/out/<name>.step` — an assembly containing the hub and
`blades` copies of the blade solid.

- `--auto` re-runs the design loop, reducing the thrust target until the
  torque drops below `1.5 × Qmax`.
- `--step-file <path>` overrides the output STEP file name.
- The first run simulates an 80-point alpha sweep per blade station; polars
  are cached in `foil_cache.json` in the working directory, so reruns are
  fast.

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

## Legacy Python version

The original Python `proply` package (scipy optimizer, OpenSCAD/STL output)
lives in [`legacy/`](legacy/). It is not part of the Rust workspace.
