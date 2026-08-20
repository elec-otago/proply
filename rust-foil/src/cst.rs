// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Kulfan (CST) airfoil parametrization (Kulfan 2008) — the canonical
//! geometry representation for rust-foil.
//!
//! Each surface is written as a class function times a Bernstein-weighted
//! shape function, plus a leading-edge (LEM) and trailing-edge (TE) term:
//!
//! ```text
//! y(x) = C(x) · S(x) + w_LE · x·(1−x)^(n+0.5) ± x·TE/2
//! C(x) = x^N1 · (1−x)^N2
//! S(x) = Σᵢ wᵢ · Bᵢ(x),   Bᵢ(x) = C(n,i)·xⁱ·(1−x)^(n−i)
//! ```
//!
//! where `n + 1` is the number of weights per side and `N1 = 0.5`, `N2 = 1.0`
//! give the classic blunt-LE / sharp-TE behaviour.  The math (including the
//! LEM exponent `n + 0.5`) mirrors `get_kulfan_coordinates` /
//! `get_kulfan_parameters` in
//! `aerosandbox/geometry/airfoil/airfoil_families.py`, the "CST with LEM"
//! flavour NeuralFoil and AeroSandbox use.
//!
//! The inverse fit (coordinates → parameters) is linear in the unknowns
//! `[upper_weights, lower_weights, LEM, TE]`, so it is a plain linear
//! least-squares problem solved via normal equations with the crate's
//! `solve::gauss` (no optimizer).  NACA 4/5-digit sections are obtained by
//! fitting the analytic NACA coordinates ([`KulfanParams::from_naca`]); the
//! resulting CST surface reproduces the NACA shape to ~1e-4 in y/c at the
//! default 8 weights per side.

use crate::naca::{naca4, naca5};
use crate::solve::gauss;
use crate::state::{IQX, PI};

/// CST parameters for one airfoil.
///
/// `lower_weights` and `upper_weights` are the Bernstein shape-function
/// weights for the lower and upper surfaces (one per Bernstein basis
/// function; `n + 1` weights where `n` is the Bernstein degree).  The LEM
/// exponent is `n + 0.5` with `n` = number of weights per side.
#[derive(Clone, Debug, PartialEq)]
pub struct KulfanParams {
    pub lower_weights: Vec<f64>,
    pub upper_weights: Vec<f64>,
    pub leading_edge_weight: f64,
    pub te_thickness: f64,
    pub n1: f64,
    pub n2: f64,
}

impl Default for KulfanParams {
    /// AeroSandbox defaults: 8 weights per side (±0.2), no LEM, sharp TE,
    /// N1 = 0.5, N2 = 1.0.  The result is a symmetric ~14% thick section.
    fn default() -> Self {
        KulfanParams {
            lower_weights: vec![-0.2; 8],
            upper_weights: vec![0.2; 8],
            leading_edge_weight: 0.0,
            te_thickness: 0.0,
            n1: 0.5,
            n2: 1.0,
        }
    }
}

impl KulfanParams {
    /// Evaluates the upper-surface y at normalized chord station `x` in
    /// [0, 1] (single-station evaluation for design/optimization loops).
    pub fn upper_y(&self, x: f64) -> f64 {
        let n = self.upper_weights.len();
        class_function(self.n1, self.n2, x) * shape_function(&self.upper_weights, x)
            + self.leading_edge_weight * lem_term(n, x)
            + self.te_thickness * x / 2.0
    }

    /// Evaluates the lower-surface y at normalized chord station `x` in
    /// [0, 1].
    pub fn lower_y(&self, x: f64) -> f64 {
        let n = self.lower_weights.len();
        class_function(self.n1, self.n2, x) * shape_function(&self.lower_weights, x)
            + self.leading_edge_weight * lem_term(n, x)
            - self.te_thickness * x / 2.0
    }

    /// Generates a closed TE→LE→TE coordinate loop at cosine-spaced
    /// stations (same ordering as `naca4`/`naca5`): the upper surface runs
    /// from the TE to the LE, then the lower surface back to the TE.  The
    /// LE point is shared, so the loop has `2·n_points_per_side − 1` points.
    /// Returns `(x, y, nb)`.
    pub fn coordinates(&self, n_points_per_side: usize) -> (Vec<f64>, Vec<f64>, usize) {
        let n = n_points_per_side;
        assert!(n >= 2, "need at least 2 points per side");
        let x: Vec<f64> = (0..n).map(|i| cosspace(i, n)).collect();

        let mut xb = vec![0.0; 2 * n - 1];
        let mut yb = vec![0.0; 2 * n - 1];
        let mut ib = 0;
        for &xi in x.iter().rev() {
            xb[ib] = xi;
            yb[ib] = self.upper_y(xi);
            ib += 1;
        }
        for &xi in x.iter().skip(1) {
            xb[ib] = xi;
            yb[ib] = self.lower_y(xi);
            ib += 1;
        }
        (xb, yb, ib)
    }

    /// Fits CST parameters (8 weights per side) to a closed TE→LE→TE
    /// coordinate loop, following AeroSandbox's `get_kulfan_parameters`.
    ///
    /// The input is normalized to unit chord; the leading-edge point is
    /// `argmin x`, and points before it are taken as the upper surface.
    /// If the fitted TE thickness is negative, the fit is repeated with the
    /// TE term dropped and `te_thickness = 0` (AeroSandbox behaviour).
    pub fn fit_from_coordinates(x: &[f64], y: &[f64]) -> KulfanParams {
        Self::fit_n(x, y, 8)
    }

    /// As [`fit_from_coordinates`](Self::fit_from_coordinates), but with an
    /// explicit number of weights per side.
    pub fn fit_from_coordinates_n(x: &[f64], y: &[f64], n_weights: usize) -> KulfanParams {
        Self::fit_n(x, y, n_weights)
    }

    /// Builds CST parameters for a NACA 4- or 5-digit airfoil by fitting
    /// the analytic NACA coordinates (sampled at `IQX/3` points per side,
    /// the same resolution the paneling uses).  Returns `None` for illegal
    /// designations (mirroring `naca5`'s error handling).
    pub fn from_naca(ides: u32, show_output: bool) -> Option<(KulfanParams, String)> {
        let nside = IQX / 3;
        let (xb, yb, nb, name) = if ides <= 9999 {
            naca4(ides as i32, nside)
        } else if ides <= 25099 {
            naca5(ides as i32, nside, show_output)
        } else {
            if show_output {
                eprintln!("This designation not implemented.");
            }
            return None;
        };
        if nb == 0 {
            return None; // naca5 illegal designation (message already printed)
        }
        let params = KulfanParams::fit_from_coordinates(&xb, &yb);
        Some((params, name))
    }

    /// The shared fit core: linear least squares via normal equations.
    fn fit_n(x: &[f64], y: &[f64], n_weights: usize) -> KulfanParams {
        assert!(n_weights >= 1, "need at least one weight per side");
        assert_eq!(x.len(), y.len(), "x and y arrays must have equal length");
        assert!(
            x.len() >= 2 * n_weights + 2,
            "need at least 2n+2 coordinates to fit {} weights per side",
            n_weights
        );

        // Normalize to unit chord (AeroSandbox default).
        let x_min = x.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let chord = x_max - x_min;
        assert!(chord > 0.0, "degenerate input: zero chord");

        let n_coords = x.len();
        let mut x0 = vec![0.0; n_coords];
        let mut y0 = vec![0.0; n_coords];
        let mut le = 0;
        let mut x_le = f64::INFINITY;
        for i in 0..n_coords {
            x0[i] = (x[i] - x_min) / chord;
            y0[i] = y[i] / chord;
            if x0[i] < x_le {
                x_le = x0[i];
                le = i;
            }
        }

        let mut params = Self::solve_fit(&x0, &y0, le, n_weights, true);
        if params.te_thickness < 0.0 {
            // Negative TE gap: re-solve with the TE term dropped.
            params = Self::solve_fit(&x0, &y0, le, n_weights, false);
            params.te_thickness = 0.0;
        }
        params
    }

    /// Builds the design matrix and solves the normal equations.  The
    /// unknowns are ordered `[upper_weights (n), lower_weights (n), LEM, TE]`.
    fn solve_fit(
        x0: &[f64],
        y0: &[f64],
        le: usize,
        n_weights: usize,
        with_te: bool,
    ) -> KulfanParams {
        let n_coords = x0.len();
        let n_cols = 2 * n_weights + if with_te { 2 } else { 1 };
        let n = n_weights - 1;

        // Design matrix, column-major (Fortran layout, nsiz = n_cols):
        // element (row i, col c) at c*n_cols + i.
        let mut a = vec![0.0; n_cols * n_coords];
        for i in 0..n_coords {
            let cfun = class_function(0.5, 1.0, x0[i]);
            for w in 0..n_weights {
                let b = bernstein(n, w, x0[i]);
                if i < le {
                    a[w * n_coords + i] = cfun * b;
                } else if i > le {
                    a[(n_weights + w) * n_coords + i] = cfun * b;
                }
            }
            a[(2 * n_weights) * n_coords + i] = lem_term(n_weights, x0[i]);
            if with_te {
                let sign = if i < le { 1.0 } else { -1.0 };
                a[(2 * n_weights + 1) * n_coords + i] = sign * x0[i] / 2.0;
            }
        }

        // Normal equations AᵀA·p = Aᵀy.
        let mut ata = vec![0.0; n_cols * n_cols];
        let mut atb = vec![0.0; n_cols];
        for j in 0..n_cols {
            for i in 0..n_coords {
                atb[j] += a[j * n_coords + i] * y0[i];
            }
            for k in 0..=j {
                let mut s = 0.0;
                for i in 0..n_coords {
                    s += a[j * n_coords + i] * a[k * n_coords + i];
                }
                ata[k * n_cols + j] = s;
                ata[j * n_cols + k] = s;
            }
        }

        gauss(&mut ata, &mut atb, n_cols, n_cols, 1);

        KulfanParams {
            upper_weights: atb[0..n_weights].to_vec(),
            lower_weights: atb[n_weights..2 * n_weights].to_vec(),
            leading_edge_weight: atb[2 * n_weights],
            te_thickness: if with_te { atb[2 * n_weights + 1] } else { 0.0 },
            n1: 0.5,
            n2: 1.0,
        }
    }
}

/// Cosine-spaced station `i` of `n`: `xᵢ = 0.5·(1 − cos(π·i/(n−1)))`.
fn cosspace(i: usize, n: usize) -> f64 {
    0.5 * (1.0 - (PI * i as f64 / (n - 1) as f64).cos())
}

/// Bernstein basis polynomial `Bᵢ(x) = C(n,i)·xⁱ·(1−x)^(n−i)`.
fn bernstein(n: usize, i: usize, x: f64) -> f64 {
    binom(n, i) * x.powi(i as i32) * (1.0 - x).powi((n - i) as i32)
}

/// Binomial coefficient `C(n, k)`.
fn binom(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    let mut c = 1.0;
    for i in 0..k {
        c = c * (n - i) as f64 / (i + 1) as f64;
    }
    c
}

/// The Kulfan class function `C(x) = x^N1·(1−x)^N2`.
fn class_function(n1: f64, n2: f64, x: f64) -> f64 {
    x.powf(n1) * (1.0 - x).powf(n2)
}

/// The leading-edge (camber) mode `x·(1−x)^(n+0.5)`, with `n` = number of
/// weights per side.
fn lem_term(n_weights: usize, x: f64) -> f64 {
    x * (1.0 - x).powf(n_weights as f64 + 0.5)
}

/// The Bernstein-weighted shape function `S(x) = Σᵢ wᵢ·Bᵢ(x)`.
fn shape_function(weights: &[f64], x: f64) -> f64 {
    let n = weights.len() - 1;
    let mut s = 0.0;
    for (i, w) in weights.iter().enumerate() {
        s += w * bernstein(n, i, x);
    }
    s
}
