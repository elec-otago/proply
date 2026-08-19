# Design: Kulfan (CST) parametrization for rust-foil

Status: proposal — not yet implemented.

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

## 2. New module `src/cst.rs`

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
    `naca4`/`naca5`).
  - `pub fn upper_y(&self, x: f64) -> f64` and `pub fn lower_y(&self, x: f64) -> f64`
    — single-station evaluation, for design/optimization loops.

Private helpers: `cosspace` (`xᵢ = 0.5·(1 − cos(π·i/(n−1)))`), `bernstein`,
`shape_function` (evaluate S(x), then multiply by C(x)).

No new dependencies — the math is a few loops.

## 3. Integration (mirrors the NACA path exactly)

- **`src/xfoil.rs`**: `pub fn cst(xf: &mut Xfoil, params: &KulfanParams)` —
  evaluate coordinates at `n_points_per_side = IQX/3` (same as NACA's
  `nside = IQX/3`), write `xf.nb/xb/yb`, set `xf.name = "CST"`, then the
  identical `scalc + segspl×2 + geopar + pangen(xf, true)` sequence
  `xfoil::naca` uses. Unit-chord output is already normalized, so it skips
  `geom::norm` (same as `xfoil::naca`).
- **`src/lib.rs`**: `pub fn cst(&mut self, params: &KulfanParams)` on `XFoil`
  — one-line delegate to `xfoil::cst`, mirroring `XFoil::naca(u32)`.
- **`src/lib.rs` module list**: add `mod cst;` (private, like `naca`).

## 4. Inverse fit (coordinates → params)

The fit is **linear in the unknowns** `[lower_weights (n), upper_weights (n),
LEM, TE]`, so it is a linear least-squares problem — no optimizer needed:

- Normalize input to unit chord (AeroSandbox default).
- Split points at the LE index (`argmin x`): upper = points before it.
- Build `A` (n_coords × (2n+2)): columns are `C(x)·Bᵢ(x)` masked to the
  lower/upper side, the LEM row `x·(1−x)^(n+0.5)`, and the TE row `±x/2`
  (sign by side).
- Solve `A·x = b = y` via normal equations `AᵀA·x = Aᵀb` using the crate's
  existing `gauss` (column-major, NRHS = 1) — reuses the ported `solve.rs`
  machinery (already `pub` and tested in `tests/gauss.rs`).
- If the fitted TE < 0, re-solve with the TE column dropped and TE = 0
  (AeroSandbox behavior).

## 5. Tests `tests/cst.rs` (integration-test style, matching the existing suite)

- **Forward sanity**: default params → LE at (0, 0), TE at (1, ±TE/2), 2n−1
  points, closed loop, upper y > 0 mid-chord.
- **Round-trip**: fit NACA 0012 → coordinates → refit → regenerate →
  coordinates match within ~1e-6 (AeroSandbox's own reversibility test).
- **Symmetry**: NACA 0012 fit gives `lower_weights[i] ≈ −upper_weights[i]`,
  LEM ≈ 0, TE ≈ 0.
- **Solver**: `XFoil::cst()` → `set_reynolds(1e6)` → `a(4°)` converges with
  CL > 0 (cambered foil), mirroring `tests/wake_gap.rs`'s `solve_naca` helper.

## 6. Alternatives considered

- **PARSEC** — 12 physics-based params, but nonlinear fit and no
  ecosystem-standard implementation; rejected.
- **B-spline/NURBS** — more params, no canonical form; rejected.
- **Plain AeroSandbox 1:1 port** — the core of this proposal, but with a
  Rust-idiomatic struct + `Default` instead of 6 positional args.

## 7. Open questions

1. **Naming**: `KulfanParams` (matches AeroSandbox) vs `CstParams` — leaning
   `KulfanParams`.
2. **Fit normalization**: normalize input to unit chord before fitting
   (AeroSandbox default) — leaning on.
3. **Default weight count**: 8/side (NeuralFoil's canonical 18-param) —
   leaning 8.
