# rust-foil

A Rust port of the **XFOIL** airfoil analysis program, derived from the
Fortran sources of the `xfoil-python` package.  `rust-foil` is a library that
simulates the inviscid and viscous flow over arbitrary airfoils using the
panel method:

- **Inviscid solution** — linear-vorticity panel method with a Kutta
  condition, solving for surface vortex strengths at any angle of attack.
- **Viscous solution** — integral boundary-layer method with the
  \(e^n\) amplification model for natural transition, marched along both
  surfaces and the wake.
- **Viscous/inviscid coupling** — full Newton iteration on the coupled
  boundary-layer system, with automatic transition-point location and
  forced-transition trips.

The numerical algorithms, closure relations and defaults are a faithful
translation of the Fortran (compiled with `-fdefault-real-8`, i.e. double
precision); the interactive menu and plotting code of the original program
are not ported.

## Usage

```rust
use rust_foil::XFoil;

let mut xf = XFoil::new();

// airfoil coordinates (counterclockwise, starting at the trailing edge)
let x = vec![1.0, 0.95, 0.9, ...];
let y = vec![0.0, 0.05, 0.08, ...];
xf.set_airfoil(&x, &y);

// ... or use a built-in NACA airfoil (generated through a CST fit)
xf.naca(2412);

// ... or define the airfoil directly with CST (Kulfan) parameters —
// the canonical geometry representation (8 weights per side + LEM + TE)
use rust_foil::KulfanParams;
xf.cst(&KulfanParams::default());

// fit CST parameters to an arbitrary coordinate loop, or generate
// coordinates / NACA sections as parameters
let params = KulfanParams::fit_from_coordinates(&x, &y);
let (xs, ys, _nb) = params.coordinates(200);
let (naca_params, _name) = KulfanParams::from_naca(2412, false).unwrap();

xf.set_reynolds(1.0e6);
xf.set_mach(0.0);
xf.set_n_crit(9.0);
xf.set_max_iter(100);

// single operating point at alpha = 10 deg
let (cl, cd, cm, cpmin, converged) = xf.a(10.0);

// or converge to a target lift coefficient
let (alpha, cd, cm, cpmin, converged) = xf.cl(1.0);

// alpha and lift sweeps
let polar = xf.aseq(-20.0, 20.0, 80); // (alpha, CL, CD, CM, Cpmin, conv)
let polar = xf.cseq(-0.5, 0.5, 20);   // (alpha, CL, CD, CM, Cpmin, conv)

// surface pressure distribution
let (xs, cp) = xf.get_cp_distribution();
```

See `examples/lift_curve.rs` for a complete NACA 0012 polar:

```sh
cargo run --release --example lift_curve
```

## Public API

| Method | Description |
| ------ | ----------- |
| `new()` | Create a solver instance |
| `set_airfoil(&[f64], &[f64])` | Set airfoil from coordinates (input points become the panel nodes) |
| `naca(u32)` | Generate and panel a NACA 4- or 5-digit airfoil (via a CST fit) |
| `cst(&KulfanParams)` | Generate and panel an airfoil from CST (Kulfan) parameters |
| `airfoil()` | Get the current buffer airfoil coordinates |
| `repanel(...)` | Re-panel the buffer airfoil with custom bunching parameters |
| `set_reynolds / reynolds` | Set/get the Reynolds number (`0` selects inviscid mode) |
| `set_mach / mach` | Set/get the Mach number |
| `set_n_crit / n_crit` | Set/get the critical amplification ratio |
| `set_xtr / xtr` | Set/get the transition trip x/c (top, bottom) |
| `set_max_iter / max_iter` | Set/get the Newton iteration limit |
| `reset_bls()` | Force BL re-initialization on the next point |
| `a(alpha_deg)` | Solve at an angle of attack → `(CL, CD, CM, Cpmin, conv)` |
| `cl(target)` | Solve at a target CL → `(alpha, CD, CM, Cpmin, conv)` |
| `aseq(a0, a1, n)` | Alpha sweep |
| `cseq(cl0, cl1, n)` | Lift sweep |
| `get_cp_distribution()` | Surface `(x, Cp)` from the current solution |
| `set_show_output(bool)` | Toggle console output |

## Validation

The integration tests in `tests/naca0012.rs` exercise NACA 0012 at
`Re = 1e6`, `M = 0`.  The reference values were originally transcribed from
the `xfoil-python` test suite; with the canonical-XFOIL `hvrat = 0.25`
setting (see below) the rust-foil output matches them to the tolerances
shown, which also cover the small (~1e-4 in CL/CD/CM) shift between the two
settings at low Mach:

| Quantity | Reference | rust-foil (hvrat=0.25) |
| -------- | --------- | ---------------------- |
| `a(10°)` CL | 1.0809 | 1.0809 |
| `a(10°)` CD | 0.0150 | 0.0150 |
| `a(10°)` CM | 0.0053 | 0.0053 |
| `cl(1.0)` alpha | 9.0617° | 9.0617° |
| `aseq(-20°, 20°)` | 80-point polar | matches |
| `cseq(-0.5, 0.5)` | 20-point polar | matches |

Run the tests:

```sh
cargo test
```

## Porting notes

- The Fortran global state (`i_xfoil`, `i_xbl`, ...) is collected in
  `src/state.rs` as a single `Xfoil` struct, so multiple independent solver
  instances can coexist (no global mutable state).
- The `com1`/`com2` pointer aliases of the BL module become a `BlVars` struct
  with 73 named fields in the exact `preptrs` ordering.
- All matrices keep the Fortran column-major layout; the block solver
  `blsolv`, `gauss`, `ludcmp` and `baksub` are translated 1:1 and verified
  against dense solves in `tests/`.
- The interactive menu, plotting, polar file I/O and inverse-design (MDES/
  QDES) paths are not ported.

## Fixes relative to upstream

`rust-foil` targets canonical (Drela) XFOIL rather than the `xfoil-python`
port it was translated from.  The following deviations from a straight 1:1
translation are intentional; all are covered by regression tests under
`tests/`.

- **`hvrat` set to the canonical XFOIL value `0.25`** (`src/state.rs`).  This
  Sutherland-viscosity parameter (`HVRAIN` in `i_blpar.f90`) enters the BL
  Reynolds number `Re_θ` and edge viscosity in `setbl`/`blkin`.  Its `DATA`
  statement was dropped during the SPAG-ification of the Python sources, so
  `xfoil-python` effectively uses `0.0`.  We restore the upstream value so
  results match canonical XFOIL.  At low Mach the impact is small (for NACA
  0012 at `Re=1e6, M=0`, CL/CD/CM shift by <1e-4 between the two settings);
  the shift grows at high Mach where the `herat^1.5 * (1+hvrat)/(herat+hvrat)`
  viscosity correction is more active.  Pinned by the existing
  `tests/naca0012.rs` reference values.
- **`baksub` forward-substitution sentinel** (`src/solve.rs`).  The Numerical
  Recipes `lubksb` pattern uses a "first index where the right-hand side
  became nonzero" variable to suppress the forward-elimination loop until the
  first nonzero appears.  The port used `0` as the "not set" sentinel, which
  collides with the valid index `0`: any system whose first nonzero RHS entry
  lands at row 0 had its forward elimination silently skipped.  The sentinel
  is now `-1` (out of the valid index range).  The bug was latent because
  partial pivoting in `ludcmp` happens to move the affected row out of
  position 0 on the inviscid Kutta system, so the existing NACA 0012 results
  are unchanged; see `tests/baksub.rs` for the failure mode.
- **Wake-gap array indexing** (`src/bl.rs`).  The `wgap` array of wake
  "dead-air" thicknesses is 0-based, but the BL station counter
  `iw = ibl - iblte[is]` is 1-based, so the correct read is `wgap[iw-1]`.  The
  four read sites in `setbl`, `mrchue`, `mrchdu` and `update` all used
  `wgap[iw]`, reading one slot off and, at `iw == nw == IWX`, indexing one
  element past the array end.  The fix is `wgap[iw-1]` at all four sites,
  with the invariant documented on the field in `src/state.rs`; see
  `tests/wake_gap.rs`.  As with `baksub`, the existing sharp-TE NACA 0012
  validation is unaffected because its wake gap is small enough that the
  off-by-one did not trip the array bound.

## License

GPL-3.0-or-later, matching the `xfoil-python` sources this port is derived
from.
