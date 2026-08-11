//! NURBS curve and surface construction for the STEP blade loft.
//!
//! The blade side surfaces are lofts through the station profiles: each
//! profile is a cubic B-spline interpolating the foil points (global
//! interpolation, Piegl & Tiller Algo A9.1), and the loft is a B-spline
//! surface of degree 1 in the radial direction with all profiles sharing a
//! common knot vector (the skinning approach of P&T §9.2.4).
//!
//! All curves/surfaces here are non-rational (weights all 1.0), which
//! step-io writes as plain `B_SPLINE_CURVE_WITH_KNOTS` /
//! `B_SPLINE_SURFACE_WITH_KNOTS`.

/// A rational B-spline curve (weights all 1.0).
#[derive(Debug, Clone, PartialEq)]
pub struct NurbsCurve {
    pub degree: usize,
    pub control_points: Vec<[f64; 3]>,
    /// Expanded knot vector (each knot repeated by multiplicity).
    pub knots: Vec<f64>,
}

/// A rational B-spline surface, control grid indexed [u][v].
#[derive(Debug, Clone, PartialEq)]
pub struct NurbsSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    pub control_points: Vec<Vec<[f64; 3]>>,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
}

impl NurbsCurve {
    pub fn weights(&self) -> Vec<f64> {
        vec![1.0; self.control_points.len()]
    }
}

impl NurbsSurface {
    pub fn weights(&self) -> Vec<Vec<f64>> {
        vec![vec![1.0; self.control_points[0].len()]; self.control_points.len()]
    }
}

/// The i-th degree-p B-spline basis function at parameter u (left-continuous
/// except at the last knot, where the final span is closed).
fn basis(i: usize, p: usize, u: f64, knots: &[f64]) -> f64 {
    if p == 0 {
        let last = knots.len() - 1;
        if knots[i] <= u && (u < knots[i + 1] || (u == knots[i + 1] && i + 1 == last)) {
            return 1.0;
        }
        return 0.0;
    }
    let mut n = 0.0;
    let d1 = knots[i + p] - knots[i];
    if d1 > 1e-300 {
        n += (u - knots[i]) / d1 * basis(i, p - 1, u, knots);
    }
    let d2 = knots[i + p + 1] - knots[i + 1];
    if d2 > 1e-300 {
        n += (knots[i + p + 1] - u) / d2 * basis(i + 1, p - 1, u, knots);
    }
    n
}

/// Chord-length parameters for `m` points, normalized to [0, 1].
fn chord_params(points: &[[f64; 3]]) -> Vec<f64> {
    let m = points.len();
    let mut u = vec![0.0; m];
    let mut total = 0.0;
    for i in 1..m {
        let d = (points[i][0] - points[i - 1][0]).hypot(
            (points[i][1] - points[i - 1][1]).hypot(points[i][2] - points[i - 1][2]),
        );
        total += d;
        u[i] = total;
    }
    for v in u.iter_mut() {
        *v /= total.max(1e-300);
    }
    u
}

/// The clamped degree-p knot vector built from interior parameters
/// (averaged per P&T eq. 9.8, which keeps the collocation matrix
/// well-conditioned).
fn clamped_knots(m: usize, degree: usize, u: &[f64]) -> Vec<f64> {
    let mut knots = Vec::with_capacity(m + degree + 1);
    for _ in 0..=degree {
        knots.push(0.0);
    }
    // U_{j+p} = (u_j + u_{j+1} + u_{j+2}) / p, j = 1 .. m-p-1
    let p = degree as f64;
    for j in 1..=(m - degree - 1) {
        let s: f64 = (0..degree).map(|k| u[j + k]).sum();
        knots.push(s / p);
    }
    for _ in 0..=degree {
        knots.push(1.0);
    }
    knots
}

/// Cubic interpolation of `points` (m >= 5): returns the curve with
/// chord-length parametrization.
pub fn interpolate(points: &[[f64; 3]]) -> NurbsCurve {
    let m = points.len();
    assert!(m >= 5, "need at least 5 points for cubic interpolation");
    let u = chord_params(points);
    let knots = clamped_knots(m, 3, &u);
    let control = interpolate_with_knots(points, &u, &knots);
    NurbsCurve {
        degree: 3,
        control_points: control,
        knots,
    }
}

/// Solve the global cubic interpolation system with a fixed collocation
/// parametrization `u` (length = points) and knot vector `knots`.
/// Returns the control points.
pub fn interpolate_with_knots(points: &[[f64; 3]], u: &[f64], knots: &[f64]) -> Vec<[f64; 3]> {
    let m = points.len();
    assert_eq!(u.len(), m);
    let p = 3;
    let n_inner = m - 2; // interior control points (P_1 .. P_{m-2})
    let mut a = vec![0.0; n_inner * n_inner];
    let mut rhs_x = vec![0.0; n_inner];
    let mut rhs_y = vec![0.0; n_inner];
    let mut rhs_z = vec![0.0; n_inner];

    for k in 1..m - 1 {
        let row = k - 1;
        for i in 1..m - 1 {
            a[row * n_inner + (i - 1)] = basis(i, p, u[k], knots);
        }
        // Move the endpoint contributions to the right-hand side.
        let n0 = basis(0, p, u[k], knots);
        let nn = basis(m - 1, p, u[k], knots);
        rhs_x[row] = points[k][0] - n0 * points[0][0] - nn * points[m - 1][0];
        rhs_y[row] = points[k][1] - n0 * points[0][1] - nn * points[m - 1][1];
        rhs_z[row] = points[k][2] - n0 * points[0][2] - nn * points[m - 1][2];
    }

    let inv = gauss_solve(&a, n_inner);
    let px = mat_vec(&inv, n_inner, &rhs_x);
    let py = mat_vec(&inv, n_inner, &rhs_y);
    let pz = mat_vec(&inv, n_inner, &rhs_z);

    let mut control = Vec::with_capacity(m);
    control.push(points[0]);
    for i in 0..n_inner {
        control.push([px[i], py[i], pz[i]]);
    }
    control.push(points[m - 1]);
    control
}

/// Gaussian elimination with partial pivoting; returns the inverse.
fn gauss_solve(a: &[f64], n: usize) -> Vec<f64> {
    // Augment with the identity (row stride 2n).
    let mut m = vec![0.0; n * n * 2];
    for r in 0..n {
        for c in 0..n {
            m[r * 2 * n + c] = a[r * n + c];
        }
        m[r * 2 * n + n + r] = 1.0;
    }
    for col in 0..n {
        // Pivot
        let mut piv = col;
        let mut piv_abs = m[col * 2 * n + col].abs();
        for r in col + 1..n {
            let v = m[r * 2 * n + col].abs();
            if v > piv_abs {
                piv_abs = v;
                piv = r;
            }
        }
        assert!(piv_abs > 1e-12, "singular interpolation matrix");
        if piv != col {
            for c in 0..2 * n {
                m.swap(col * 2 * n + c, piv * 2 * n + c);
            }
        }
        let d = m[col * 2 * n + col];
        for c in 0..2 * n {
            m[col * 2 * n + c] /= d;
        }
        for r in 0..n {
            if r != col {
                let f = m[r * 2 * n + col];
                if f != 0.0 {
                    for c in 0..2 * n {
                        m[r * 2 * n + c] -= f * m[col * 2 * n + c];
                    }
                }
            }
        }
    }
    // Extract the inverse (right half).
    let mut inv = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            inv[r * n + c] = m[r * 2 * n + n + c];
        }
    }
    inv
}

fn mat_vec(a: &[f64], n: usize, x: &[f64]) -> Vec<f64> {
    (0..n).map(|r| (0..n).map(|c| a[r * n + c] * x[c]).sum()).collect()
}

/// A degree-1 polyline curve through `points` (control points = points).
pub fn polyline(points: &[[f64; 3]]) -> NurbsCurve {
    let m = points.len();
    assert!(m >= 2);
    let mut knots = vec![0.0, 0.0];
    for i in 1..m {
        knots.push(i as f64);
    }
    knots.push((m - 1) as f64);
    NurbsCurve {
        degree: 1,
        control_points: points.to_vec(),
        knots,
    }
}

/// A straight line between two points.
pub fn line(a: [f64; 3], b: [f64; 3]) -> NurbsCurve {
    polyline(&[a, b])
}

/// Average the per-profile chord-length parameters to a common
/// parametrization (P&T §9.2.4), for profiles with equal point counts.
pub fn common_params(profiles: &[Vec<[f64; 3]>]) -> Vec<f64> {
    let m = profiles[0].len();
    let mut acc = vec![0.0; m];
    for prof in profiles {
        assert_eq!(prof.len(), m, "profiles must have equal point counts");
        let u = chord_params(prof);
        for (a, v) in acc.iter_mut().zip(u.iter()) {
            *a += v;
        }
    }
    let n = profiles.len() as f64;
    acc.iter().map(|v| v / n).collect()
}

/// The clamped cubic knot vector for the common parametrization.
pub fn common_knots(m: usize, u: &[f64]) -> Vec<f64> {
    clamped_knots(m, 3, u)
}

/// A ruled loft through the station profiles (each already interpolated with
/// the common knots): degree 1 in v, one knot span per station.
pub fn loft(profiles: &[NurbsCurve]) -> NurbsSurface {
    let n_v = profiles.len();
    assert!(n_v >= 2);
    let n_u = profiles[0].control_points.len();
    let knots_u = profiles[0].knots.clone();
    let mut control = vec![vec![[0.0; 3]; n_v]; n_u];
    for (j, prof) in profiles.iter().enumerate() {
        assert_eq!(prof.knots, knots_u, "profiles must share the knot vector");
        assert_eq!(prof.control_points.len(), n_u);
        for (i, pt) in prof.control_points.iter().enumerate() {
            control[i][j] = *pt;
        }
    }
    let mut knots_v = vec![0.0, 0.0];
    for j in 1..n_v {
        knots_v.push(j as f64);
    }
    knots_v.push((n_v - 1) as f64);
    NurbsSurface {
        degree_u: 3,
        degree_v: 1,
        control_points: control,
        knots_u,
        knots_v,
    }
}

/// A ruled surface between two curves (v = 0..1).
pub fn ruled(a: &NurbsCurve, b: &NurbsCurve) -> NurbsSurface {
    assert_eq!(a.knots, b.knots, "ruled curves must share the knot vector");
    assert_eq!(a.control_points.len(), b.control_points.len());
    let n_u = a.control_points.len();
    let mut control = vec![vec![[0.0; 3]; 2]; n_u];
    for i in 0..n_u {
        control[i][0] = a.control_points[i];
        control[i][1] = b.control_points[i];
    }
    NurbsSurface {
        degree_u: a.degree,
        degree_v: 1,
        control_points: control,
        knots_u: a.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval_curve(c: &NurbsCurve, u: f64) -> [f64; 3] {
        // At the last knot a clamped B-spline equals its last control point.
        let last = *c.knots.last().unwrap();
        if u >= last {
            return c.control_points[c.control_points.len() - 1];
        }
        let n = c.control_points.len();
        let mut acc = [0.0; 3];
        for i in 0..n {
            let b = basis(i, c.degree, u, &c.knots);
            acc[0] += b * c.control_points[i][0];
            acc[1] += b * c.control_points[i][1];
            acc[2] += b * c.control_points[i][2];
        }
        acc
    }

    fn eval_surface(s: &NurbsSurface, u: f64, v: f64) -> [f64; 3] {
        let nu = s.control_points.len();
        let nv = s.control_points[0].len();
        let last_u = *s.knots_u.last().unwrap();
        let last_v = *s.knots_v.last().unwrap();
        let mut acc = [0.0; 3];
        for i in 0..nu {
            for j in 0..nv {
                let bu = if u >= last_u && i == nu - 1 {
                    1.0
                } else {
                    basis(i, s.degree_u, u, &s.knots_u)
                };
                let bv = if v >= last_v && j == nv - 1 {
                    1.0
                } else {
                    basis(j, s.degree_v, v, &s.knots_v)
                };
                let b = bu * bv;
                acc[0] += b * s.control_points[i][j][0];
                acc[1] += b * s.control_points[i][j][1];
                acc[2] += b * s.control_points[i][j][2];
            }
        }
        acc
    }

    fn point(x: f64, y: f64) -> [f64; 3] {
        [x, y, 0.0]
    }

    #[test]
    fn cubic_interpolation_passes_through_points() {
        // A smooth arc of points; the interpolant must reproduce them.
        let pts: Vec<[f64; 3]> = (0..15)
            .map(|i| {
                let t = i as f64 / 14.0;
                point(t, 0.3 * (std::f64::consts::PI * t).sin())
            })
            .collect();
        let c = interpolate(&pts);
        let u = chord_params(&pts);
        for (k, p) in pts.iter().enumerate() {
            let ck = eval_curve(&c, u[k]);
            assert!(
                (ck[0] - p[0]).abs() < 1e-9 && (ck[1] - p[1]).abs() < 1e-9,
                "point {}: {:?} vs {:?}",
                k,
                ck,
                p
            );
        }
    }

    #[test]
    fn polyline_interpolates_exactly() {
        let pts = [point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.5)];
        let c = polyline(&pts);
        assert_eq!(eval_curve(&c, 0.0), pts[0]);
        assert_eq!(eval_curve(&c, 1.0), pts[1]);
        assert_eq!(eval_curve(&c, 2.0), pts[2]);
        // midpoint of the first span
        let m = eval_curve(&c, 0.5);
        assert!((m[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn loft_interpolates_each_profile() {
        // Two profiles; the loft surface must reproduce them at v=0 and v=1.
        let prof_a: Vec<[f64; 3]> = (0..12).map(|i| point(i as f64, (i as f64).sin())).collect();
        let prof_b: Vec<[f64; 3]> = (0..12).map(|i| point(i as f64, (i as f64).cos())).collect();
        let u = common_params(&[prof_a.clone(), prof_b.clone()]);
        let knots = common_knots(12, &u);
        let ca = interpolate_with_knots(&prof_a, &u, &knots);
        let cb = interpolate_with_knots(&prof_b, &u, &knots);
        let a = NurbsCurve { degree: 3, control_points: ca, knots: knots.clone() };
        let b = NurbsCurve { degree: 3, control_points: cb, knots };
        let s = ruled(&a, &b);
        for (k, p) in prof_a.iter().enumerate() {
            let on = eval_surface(&s, u[k], 0.0);
            assert!((on[1] - p[1]).abs() < 1e-9, "v=0 point {}: {:?} vs {:?}", k, on, p);
        }
        for (k, p) in prof_b.iter().enumerate() {
            let on = eval_surface(&s, u[k], 1.0);
            assert!((on[1] - p[1]).abs() < 1e-9, "v=1 point {}: {:?} vs {:?}", k, on, p);
        }
        // Linear in v between the two profiles.
        let mid = eval_surface(&s, u[5], 0.5);
        let want = 0.5 * (prof_a[5][1] + prof_b[5][1]);
        assert!((mid[1] - want).abs() < 1e-9, "mid {}", mid[1]);
    }

    #[test]
    fn loft_across_stations_is_exact_at_span_junctions() {
        let n_stations = 4;
        let profiles: Vec<Vec<[f64; 3]>> = (0..n_stations)
            .map(|j| {
                (0..10)
                    .map(|i| point(i as f64, 0.1 * i as f64 + j as f64))
                    .collect()
            })
            .collect();
        let u = common_params(&profiles);
        let knots = common_knots(10, &u);
        let curves: Vec<NurbsCurve> = profiles
            .iter()
            .map(|p| {
                let cp = interpolate_with_knots(p, &u, &knots);
                NurbsCurve { degree: 3, control_points: cp, knots: knots.clone() }
            })
            .collect();
        let s = loft(&curves);
        for (j, prof) in profiles.iter().enumerate() {
            for (k, p) in prof.iter().enumerate() {
                let on = eval_surface(&s, u[k], j as f64);
                assert!(
                    (on[1] - p[1]).abs() < 1e-9,
                    "station {} point {}: {:?} vs {:?}",
                    j,
                    k,
                    on,
                    p
                );
            }
        }
    }
}
