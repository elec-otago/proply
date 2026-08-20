# Design: Kulfan (CST) parametrization for rust-foil

Status: **implemented** (2026-08) — the canonical geometry representation.
`XFoil::naca()` now generates its airfoil through a CST fit; see §3.

## 1. Background — the math

CST (Class/Shape Transformation, Kulfan 2008) represents each surface as a
class function times a Bernstein-weighted shape function:

```
y(x) = C(x) · S(x) + TE terms + LEM term
```

| Term | Formula | Role |
| --- | --- | --- |
| Class function | `C(x) = x^N1 · (1−x)^N2` | N1 = 0.5 gives the √x blunt-LE behavior; N2 = 1.0 gives a sharp TE |
| Shape function | `S(x) = Σᵢ wᵢ · Bᵢ(x)` | Bernstein polynomials `Bᵢ(x) = C(n,i)·xⁱ·(1−x)^(n−i)` weighted by the CST coefficients |
| TE thickness | `y_upper += x·TE/2`, `y_lower −= x·TE/2` | finite TE gap; 0 = sharp TE |
| LEM | `y += w_LE · x·(1−x)^(n+0.5)` | leading-edge camber mode (Kulfan 2020) |

With 8 weights per side + LEM + TE = **18 parameters** — the exact
representation NeuralFoil and AeroSandbox (the "CST with LEM" flavor) use.
The math was verified against the AeroSandbox source
(`get_kulfan_coordinates` / `get_kulfan_parameters` in
`aerosandbox/geometry/airfoil/airfoil_families.py`), including the LEM
exponent `n + 0.5` where `n` = number of weights per side.

## 2. Module `src/cst.rs`

```rust
pub struct KulfanParams {
    pub lower_weights: Vec<f64>,
    pub upper_weights: Vec<f64>,
    pub leading_edge_weight: f64,
    pub te_thickness: f64,
    pub n1: f64,
    pub n2: f64,
}
```

- `Default` → AeroSandbox defaults: `lower = −0.2·ones(8)`,
  `upper = 0.2·ones(8)`, LEM = 0, TE = 0, N1 = 0.5, N2 = 1.0.
- `impl KulfanParams`:
  - `pub fn coordinates(&self, n_points_per_side: usize) -> (Vec<f64>, Vec<f64>, usize)`
    — closed TE→LE→TE loop at cosine-spaced stations (same ordering as
    `naca4`/`naca5`; `nb = 2·n − 1` with a shared LE point).
  - `pub fn upper_y(&self, x: f64) -> f64` and `pub fn lower_y(&self, x: f64) -> f64`
    — single-station evaluation, for design/optimization loops.
  - `pub fn fit_from_coordinates(x: &[f64], y: &[f64]) -> KulfanParams`
    — inverse fit (8 weights/side; see §4).
  - `pub fn fit_from_coordinates_n(..., n_weights: usize)` — as above with an
    explicit weight count.
  - `pub fn from_naca(ides: u32, show_output: bool) -> Option<(KulfanParams, String)>`
    — NACA 4/5-digit → CST parameters by fitting the analytic coordinates
    (`None` for illegal designations).

Private helpers: `cosspace` (`xᵢ = 0.5·(1 − cos(π·i/(n−1)))`), `bernstein`,
`binom`, `class_function`, `lem_term`, `shape_function`.  Re-exported from
the crate root as `rust_foil::KulfanParams` (the module itself is private,
like `naca`).

No new dependencies — the math is a few loops, and the fit reuses the
crate's `solve::gauss`.

## 3. Integration — CST is the default geometry interface

- **`src/xfoil.rs`**: the buffer-setting tail (scalc + segspl×2 + geopar +
  pangen, no `geom::norm` — unit-chord output is already normalized) is
  factored into the private `set_buffer(xf, xb, yb, nb, name)` helper.
  - `pub fn cst(xf: &mut Xfoil, params: &KulfanParams)` — evaluates
    `params.coordinates(IQX/3)` (same `nside` the NACA path used), writes
    the buffer, name `"CST"`.
  - `pub fn naca(xf: &mut Xfoil, ides: i32)` — **generates the airfoil
    through the CST path**: `KulfanParams::from_naca` → `cst`.  The buffer
    geometry is the CST fit of the analytic NACA shape (cosine-spaced
    stations instead of the old TE-bunched stations), with the NACA name
    preserved.
- **`src/lib.rs`**: `pub use cst::KulfanParams;` (the module is private but
  the struct must be nameable by callers of `XFoil::cst`), plus
  `XFoil::cst(&mut self, params: &KulfanParams)` — a one-line delegate to
  `xfoil::cst`, mirroring `XFoil::naca`.

### NACA fidelity

NACA 4/5-digit shapes are *not* exactly CST-representable (the class
function's (1−x) factor fights the analytic TE value; camber lines are
piecewise polynomials), so `naca()` now produces an approximation.
Measured at 8 weights/side:

| Section | max \|Δy/c\| vs analytic | rms |
| --- | --- | --- |
| NACA 0012 | 1.1e-4 | 3.1e-5 |
| NACA 2412 | 1.1e-4 | 3.7e-5 |
| NACA 4412 | 1.2e-4 | 5.0e-5 |
| NACA 23012 | 1.5e-4 | 5.3e-5 |

The fitted TE gap reproduces the analytic open TE (2.57e-3 for 0012), and
the LEM weight absorbs the camber.  The fit is **exact on CST-generated
points** (round-trip to ~1e-15), so design loops in parameter space are
lossless.  The existing `tests/naca0012.rs` reference polars (transcribed
from xfoil-python) still pass unchanged within their original tolerance
bands — the geometry shift is below the panel-method noise floor.

## 4. Inverse fit (coordinates → params)

The fit is **linear in the unknowns** `[upper_weights (n), lower_weights
(n), LEM, TE]`, so it is a linear least-squares problem — no optimizer
needed:

- Normalize input to unit chord (AeroSandbox default).
- Split points at the LE index (`argmin x`): upper = points before it.
- Build `A` (n_coords × (2n+2)): columns are `C(x)·Bᵢ(x)` masked to the
  upper/lower side, the LEM column `x·(1−x)^(n+0.5)`, and the TE column
  `±x/2` (sign by side).
- Solve `A·x = b = y` via normal equations `AᵀA·x = Aᵀb` using the crate's
  existing `gauss` (column-major, NRHS = 1) — reuses the ported `solve.rs`
  machinery (already `pub` and tested in `tests/gauss.rs`).
- If the fitted TE < 0, re-solve with the TE column dropped and TE = 0
  (AeroSandbox behavior).

## 5. Tests `tests/cst.rs`

- **Forward sanity**: default params → LE at (0, 0), TE closure, `2n−1`
  points, upper y > 0 mid-chord, `upper_y`/`lower_y` agree with the loop.
- **Round-trip**: fit arbitrary params' coordinates → regenerate → match
  within ~1e-6; fitted params recover the originals.
- **Symmetry**: NACA 0012 fit gives `lower_weights[i] ≈ −upper_weights[i]`,
  LEM ≈ 0, TE ≈ 2.5e-3.
- **Fidelity**: NACA 0012/2412 CST shapes preserve max thickness (0.12) and
  max camber (0.02 at x/c = 0.4) within 1e-3.
- **Solver**: `XFoil::cst()` on a slightly blunt (te = 0.002) section →
  `set_reynolds(1e6)` → `a(4°)` converges with CL > 0.  (The `Default`
  section has a perfectly sharp TE, on which XFOIL's viscous solve does not
  converge — inviscid and blunt-TE solves are fine.)
- **NACA path**: `XFoil::naca(2412)`'s buffer is bit-identical to
  `XFoil::cst(KulfanParams::from_naca(2412))`.

The existing suite (`naca0012.rs`, `wake_gap.rs`, `normalize.rs`, …) passes
unchanged on the CST-generated NACA geometry.

## 6. Alternatives considered

- **PARSEC** — 12 physics-based params, but nonlinear fit and no
  ecosystem-standard implementation; rejected.
- **B-spline/NURBS** — more params, no canonical form; rejected.
- **Plain AeroSandbox 1:1 port** — the core of this implementation, but
  with a Rust-idiomatic struct + `Default` instead of 6 positional args.
- **NACA kept analytic** (dual path) — considered; rejected in favour of a
  single geometry path (Option A), because the measured polar shift from
  the CST fit is within the existing reference tolerance bands.

## 7. Decisions (formerly open questions)

1. **Naming**: `KulfanParams` (matches AeroSandbox).
2. **Fit normalization**: normalize input to unit chord (AeroSandbox
   default) — always on.
3. **Default weight count**: 8/side (NeuralFoil's canonical 18-param).
4. **NACA path**: generates CST (`from_naca` → `cst`), name preserved;
   buffer is now cosine-spaced CST stations rather than TE-bunched.
