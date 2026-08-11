//! Least-squares polynomial fitting (replaces `numpy.polyfit`) and Horner
//! evaluation (replaces `numpy.poly1d`).
//!
//! The fit solves the Vandermonde least-squares problem with a column-scaled
//! Householder QR, which is robust for the degree-9 polar fits used by the
//! simulator (x in radians, |x| <= 0.35).

/// Return the coefficients of a degree-`deg` polynomial least-squares fit to
/// `(x, y)`, ordered highest power first (as `np.polyfit` does).
pub fn polyfit(x: &[f64], y: &[f64], deg: usize) -> Vec<f64> {
    assert_eq!(x.len(), y.len(), "x and y arrays must have equal length");
    let m = x.len();
    assert!(m > deg, "at least deg+1 data points required");

    // Vandermonde: A[i][j] = x_i^(deg - j)  (column j is power deg-j)
    let mut a: Vec<f64> = vec![0.0; m * (deg + 1)];
    for i in 0..m {
        let mut p = 1.0;
        for j in (0..=deg).rev() {
            a[i * (deg + 1) + j] = p;
            p *= x[i];
        }
    }
    // a[i][k] now holds x_i^k for k = 0..deg

    // Column scaling for conditioning.
    let n = deg + 1;
    let mut scale = vec![0.0; n];
    for j in 0..n {
        let mut mx: f64 = 0.0;
        for i in 0..m {
            mx = mx.max(a[i * n + j].abs());
        }
        scale[j] = if mx > 0.0 { 1.0 / mx } else { 1.0 };
    }
    for i in 0..m {
        for j in 0..n {
            a[i * n + j] *= scale[j];
        }
    }

    // Householder QR of the scaled A.
    let mut r = a.clone();
    let mut qty: Vec<f64> = y.to_vec();
    for k in 0..n {
        // norm of column k rows k..m
        let mut norm = 0.0;
        for i in k..m {
            norm += r[i * n + k] * r[i * n + k];
        }
        norm = norm.sqrt();
        if norm.abs() < 1e-300 {
            continue;
        }
        let alpha = if r[k * n + k] >= 0.0 { -norm } else { norm };
        // v = r[k..m, k] - alpha * e1
        let mut v: Vec<f64> = (k..m).map(|i| r[i * n + k]).collect();
        v[0] -= alpha;
        let vnorm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if vnorm < 1e-300 {
            continue;
        }
        let v: Vec<f64> = v.iter().map(|x| x / vnorm).collect();

        // Apply H to the remaining columns of R
        for j in k..n {
            let mut dot = 0.0;
            for (i, vi) in v.iter().enumerate() {
                dot += vi * r[(k + i) * n + j];
            }
            for (i, vi) in v.iter().enumerate() {
                r[(k + i) * n + j] -= 2.0 * vi * dot;
            }
        }
        // Apply H to the right-hand side
        let mut dot = 0.0;
        for (i, vi) in v.iter().enumerate() {
            dot += vi * qty[k + i];
        }
        for (i, vi) in v.iter().enumerate() {
            qty[k + i] -= 2.0 * vi * dot;
        }
    }

    // Back substitution on the upper-triangular R (first n rows).
    let mut z = vec![0.0; n];
    for j in (0..n).rev() {
        let mut acc = qty[j];
        for k in (j + 1)..n {
            acc -= r[j * n + k] * z[k];
        }
        z[j] = acc / r[j * n + j];
    }

    // Unscale.
    z.iter().zip(scale.iter()).map(|(zj, s)| zj * s).collect()
}

/// Evaluate a polynomial with coefficients highest-power-first at `x`.
pub fn polyval(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    for c in coeffs {
        acc = acc * x + c;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_line_fit() {
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| 3.0 * v + 2.0).collect();
        let c = polyfit(&x, &y, 1);
        assert!((c[0] - 3.0).abs() < 1e-9, "slope {}", c[0]);
        assert!((c[1] - 2.0).abs() < 1e-9, "intercept {}", c[1]);
    }

    #[test]
    fn exact_parabola_fit() {
        let x: Vec<f64> = (0..20).map(|i| (i as f64 - 10.0) * 0.1).collect();
        let y: Vec<f64> = x.iter().map(|v| 2.0 * v * v - v + 0.5).collect();
        let c = polyfit(&x, &y, 2);
        assert!((c[0] - 2.0).abs() < 1e-8, "c0 {}", c[0]);
        assert!((c[1] + 1.0).abs() < 1e-8, "c1 {}", c[1]);
        assert!((c[2] - 0.5).abs() < 1e-8, "c2 {}", c[2]);
    }

    #[test]
    fn noisy_fit_is_reasonable() {
        // Fit a line through points with small noise; residual should be small.
        let x: Vec<f64> = (0..50).map(|i| i as f64 / 50.0).collect();
        let y: Vec<f64> = x
            .iter()
            .enumerate()
            .map(|(i, v)| 1.0 - 2.0 * v + 0.001 * (i as f64 % 7.0 - 3.0))
            .collect();
        let c = polyfit(&x, &y, 1);
        assert!((c[0] + 2.0).abs() < 1e-3, "slope {}", c[0]);
        assert!((c[1] - 1.0).abs() < 1e-3, "intercept {}", c[1]);
    }

    #[test]
    fn horner_matches_powers() {
        let c = vec![1.0, -2.0, 0.5, 3.0];
        let x: f64 = 0.7;
        let expected = 1.0 * x.powi(3) - 2.0 * x.powi(2) + 0.5 * x + 3.0;
        assert!((polyval(&c, x) - expected).abs() < 1e-12);
    }
}
