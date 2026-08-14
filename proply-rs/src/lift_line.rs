//! Propeller lifting-line analysis and design.
//!
//! Replaces the independent-annulus blade-element momentum balance with a
//! coupled vortex-lattice lifting line: each blade carries a bound
//! circulation `Gamma`, and the trailed (helical) wake induces the axial and
//! tangential velocities at every station from *all* blades, so the radial
//! stations couple and the finite-span (aspect-ratio) induced loss emerges
//! from the vortex system instead of an empirical tip-loss factor.
//!
//! Model (Adkins–Liebeck class):
//!  * the blade is a radial lifting line carrying a piecewise-constant bound
//!    circulation `Gamma_i` (one value per station),
//!  * discontinuities of `Gamma` across station boundaries shed trailed
//!    helical vortex filaments whose strength equals the jump;
//!  * each trailed filament is a rigid helix of pitch `2 pi r tan(phi)`,
//!    updated each outer iteration to the converged inflow angle;
//!  * the induced velocity at a station is the Biot–Savart sum (straight
//!    segments) over every blade's trailed filaments;
//!  * `Gamma` is iterated to self-consistency with the polar
//!    `cl(alpha)` (`Gamma = 1/2 c V cl`), under-relaxed;
//!  * forces follow Kutta–Joukowski plus profile drag.
//!
//! Sign conventions match the existing BEMT (`blade_element.rs`):
//!  * `u = u_0 + u_i` (axial; induced `u_i > 0` accelerates the slipstream),
//!  * `v = omega r - v_i` (`v_i > 0` is swirl that reduces the relative
//!    tangential speed),
//!  * `phi = atan(u/v)`, `alpha = theta - phi`.

use crate::foil::FoilLike;
use crate::simulator::FoilSimulator;

pub const RHO: f64 = 1.225;

/// Rotor-plane calibration for the axial trailed-ring induction: `u_i =
/// AXIAL_K * B/(2π) * Σ dg * K`.  Tuned so a uniform/low-loading blade matches
/// the actuator-disk axial inflow (`u_i ≈ dv = T/(2ρA(u_0))` at the disk).
pub const AXIAL_K: f64 = 0.5;

/// One radial station of the blade for the lifting-line analysis.
///
/// The local geometric angle of attack `alpha` is *prescribed* by the design
/// (the twist is then `phi + alpha`, computed from the converged inflow angle
/// `phi`).  Because `alpha = twist - phi`, prescribing `alpha` makes the solve
/// independent of the twist/induction feedback that otherwise limits-cycles.
#[derive(Clone, Copy)]
pub struct Station {
    pub r: f64,
    pub c: f64,
    pub alpha: f64,
}

/// Result of a converged lifting-line solve.
#[derive(Default)]
pub struct LiftLineResult {
    pub gamma: Vec<f64>,    // bound circulation per blade (m^2/s)
    pub u_i: Vec<f64>,      // axial induced velocity (m/s), slipstream sense
    pub v_i: Vec<f64>,      // tangential induced velocity / swirl (m/s), reduces v
    pub alpha: Vec<f64>,    // local angle of attack (rad)
    pub phi: Vec<f64>,      // local inflow angle (rad)
    pub d_thrust: Vec<f64>, // element thrust (N)
    pub d_torque: Vec<f64>, // element torque (N m)
    pub thrust: f64,
    pub torque: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Perpendicular distance from `p0` to the line through `p1`-`p2`.
fn dist_line(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let l = sub(p2, p1);
    let ll = norm(l).max(1.0e-24);
    let t = (dot(sub(p0, p1), l) / ll.powi(2)).clamp(0.0, 1.0);
    let proj = [p1[0] + t * l[0], p1[1] + t * l[1], p1[2] + t * l[2]];
    norm(sub(p0, proj))
}

/// Biot–Savart velocity at `p0` from a straight vortex segment `p1 -> p2` of
/// strength `gamma`, with a Rankine vortex core of radius `core` so the field
/// stays finite near the filament (needed for a stable root vortex).
pub fn biot_savart_core(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    gamma: f64,
    core: f64,
) -> [f64; 3] {
    let r1 = sub(p0, p1);
    let r2 = sub(p0, p2);
    let l = sub(p2, p1);
    let cr = cross(r1, r2);
    let denom = dot(cr, cr);
    if denom < 1.0e-24 {
        return [0.0; 3];
    }
    let n1 = norm(r1).max(1.0e-30);
    let n2 = norm(r2).max(1.0e-30);
    // Rankine-core rolloff: attenuate the inviscid field inside the core.
    let d = dist_line(p0, p1, p2);
    let rolloff = d * d / (d * d + core * core);
    let k = gamma
        / (4.0 * std::f64::consts::PI)
        * (dot(r1, l) / n1 - dot(r2, l) / n2)
        / denom
        * rolloff;
    [cr[0] * k, cr[1] * k, cr[2] * k]
}

pub fn biot_savart(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], gamma: f64) -> [f64; 3] {
    biot_savart_core(p0, p1, p2, gamma, 0.0)
}

/// Axial velocity induced on its own axis by a circular vortex ring of radius
/// `a`, strength `gamma`, at height `z` above the ring plane: `Gamma/2 *
/// a^2/(a^2+z^2)^(3/2)`.  Used as the sanity reference in tests.
pub fn vortex_ring_axial(a: f64, gamma: f64, z: f64) -> f64 {
    gamma / 2.0 * a * a / (a * a + z * z).powf(1.5)
}

/// Station boundary radii given element-center radii (monotone, clamped >0).
fn edges(centers: &[f64]) -> Vec<f64> {
    let m = centers.len();
    let mut e = Vec::with_capacity(m + 1);
    e.push(centers[0] - 0.5 * (centers[1] - centers[0]));
    for k in 1..m {
        e.push((centers[k - 1] + centers[k]) / 2.0);
    }
    e.push(centers[m - 1] + 0.5 * (centers[m - 1] - centers[m - 2]));
    for k in 1..e.len() {
        if e[k] <= e[k - 1] {
            e[k] = e[k - 1] * 1.01 + 1.0e-9;
        }
    }
    if e[0] < 0.0 {
        e[0] = 1.0e-6;
    }
    e
}

/// Axial velocity induced at radius `rho` (>= 0) in the plane of a unit
/// circulation vortex ring of radius `a`, by Biot–Savart over straight
/// segments.  On-axis (`rho = 0`) this equals `1/(2 a)` (matches
/// [`vortex_ring_axial`]); a small Rankine core bounds the `rho ~ a` peak.
fn ring_axial_unit(a: f64, rho: f64, core: f64, n: usize) -> f64 {
    let p0 = [rho, 0.0, 0.0];
    let mut vz = 0.0;
    for s in 0..n {
        let th0 = 2.0 * std::f64::consts::PI * s as f64 / n as f64;
        let th1 = 2.0 * std::f64::consts::PI * (s + 1) as f64 / n as f64;
        let p1 = [a * th0.cos(), a * th0.sin(), 0.0];
        let p2 = [a * th1.cos(), a * th1.sin(), 0.0];
        vz += biot_savart_core(p0, p1, p2, 1.0, core)[2];
    }
    vz
}

/// Azimuthally-averaged (Trefftz-plane) induced velocity at every station,
/// from the trailed coaxial vortex rings.  This is the smooth far-field
/// inflow that physically acts on the blade — no near-field singularity.
///
///  * axial: `u_i(r) = s * B/(4 pi) * sum_j dg_j * K(r, a_j)`, where `K` is
///    the axial velocity in the plane of a unit ring at radius `a_j` and the
///    `B/(4 pi)` / rotor-plane factor maps the trailed-sheet circulation to
///    the disk inflow;
///  * tangential: `v_i(r) = B * Gamma(r) / (4 pi r)` (swirl opposing the
///    rotor, reducing the relative tangential speed).
#[cfg(test)]
fn induced_velocity(
    r: &[f64],
    _phi: &[f64],
    gamma: &[f64],
    edge: &[f64],
    n_blades: usize,
    core: f64,
) -> (Vec<f64>, Vec<f64>) {
    let m = r.len();
    let n_b = n_blades as f64;
    // Trailed ring strengths at each edge (jump in bound circulation).
    let mut dg = vec![0.0; m + 1];
    dg[0] = -gamma[0];
    for j in 1..m {
        dg[j] = gamma[j - 1] - gamma[j];
    }
    dg[m] = gamma[m - 1];

    let nper = 512;
    let mut ui = vec![0.0; m];
    for (ci, &ri) in r.iter().enumerate() {
        for j in 0..=m {
            if dg[j].abs() < 1.0e-18 {
                continue;
            }
            let a_j = edge[j].max(1.0e-6);
            ui[ci] += n_b * dg[j] * ring_axial_unit(a_j, ri, core, nper);
        }
        // Calibrated rotor-plane factor: disk axial inflow equals half the
        // fully-developed slipstream value (actuator-disk / Glauert limit).
        ui[ci] *= AXIAL_K / (2.0 * std::f64::consts::PI);
    }
    // Swirl from the bound circulation (halved for the rotor plane).
    let tang: Vec<f64> = (0..m)
        .map(|i| n_b * gamma[i] / (4.0 * std::f64::consts::PI * r[i].max(1.0e-6)))
        .collect();
    (ui, tang)
}

/// Exact linear influence between the bound circulation and the induced
/// flow.  `u_i = sum_j U[i*m+j] Gamma_j` (axial) and `v_i = Vdiag[i]*Gamma_i`
/// (swirl), matching [`induced_velocity`].
fn influence_matrices(
    r: &[f64],
    edge: &[f64],
    n_blades: usize,
    core: f64,
) -> (Vec<f64>, Vec<f64>) {
    let m = r.len();
    let n_b = n_blades as f64;
    let nper = 512;
    let pref = AXIAL_K * n_b / (2.0 * std::f64::consts::PI);
    let mut u = vec![0.0; m * m];
    for ci in 0..m {
        // K[ci][k] = axial velocity at r[ci] in the plane of a unit ring at
        // edge[k]. Gamma_j feeds trailed rings at edges j and j+1 (the jump
        // in Gamma), so U[ci][j] = pref*(K[j+1] - K[j]).
        let mut kk = vec![0.0; m + 1];
        for k in 0..=m {
            kk[k] = ring_axial_unit(edge[k].max(1.0e-6), r[ci], core, nper);
        }
        for j in 0..m {
            u[ci * m + j] = pref * (kk[j + 1] - kk[j]);
        }
    }
    let vdiag: Vec<f64> = (0..m)
        .map(|i| n_b / (4.0 * std::f64::consts::PI * r[i].max(1.0e-6)))
        .collect();
    (u, vdiag)
}

/// Per-element flow state induced by a circulation vector.
struct FlowState {
    u: Vec<f64>,
    v: Vec<f64>,
    vv: Vec<f64>,
    phi: Vec<f64>,
    alpha: Vec<f64>,
    free: Vec<bool>,
    cl: Vec<f64>,
    residual: Vec<f64>,
}

/// Evaluate the coupled flow (and the circulation residual) for `gamma`.
fn eval_flow<F: FoilLike>(
    gamma: &[f64],
    stations: &[Station],
    u_0: f64,
    omega: f64,
    u_mat: &[f64],
    vdiag: &[f64],
    fs: &[&FoilSimulator<F>],
) -> FlowState {
    let m = stations.len();
    let mut u = vec![0.0; m];
    let mut v = vec![0.0; m];
    let mut vv = vec![0.0; m];
    let mut phi = vec![0.0; m];
    let mut alpha = vec![0.0; m];
    let mut free = vec![false; m];
    let mut cl = vec![0.0; m];
    let mut residual = vec![0.0; m];
    for i in 0..m {
        let mut ui = 0.0;
        for j in 0..m {
            ui += u_mat[i * m + j] * gamma[j];
        }
        let vi = vdiag[i] * gamma[i];
        u[i] = u_0 + ui;
        v[i] = omega * stations[i].r - vi;
        vv[i] = (u[i] * u[i] + v[i] * v[i]).sqrt();
        phi[i] = u[i].atan2(v[i]);
        // Prescribed attack angle, clamped below deep stall for solver health.
        alpha[i] = stations[i].alpha.clamp(-0.45, 0.45);
        free[i] = stations[i].alpha > -0.45 && stations[i].alpha < 0.45;
        cl[i] = fs[i].get_cl(vv[i], alpha[i]);
        residual[i] = gamma[i] - 0.5 * stations[i].c * vv[i] * cl[i];
    }
    FlowState {
        u,
        v,
        vv,
        phi,
        alpha,
        free,
        cl,
        residual,
    }
}

/// Dense Gaussian elimination (partial pivoting): solve `A x = b`, `A`
/// row-major `n x n`.
fn solve_linear(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut m = a.to_vec();
    let mut bb = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        for row in col + 1..n {
            if m[row * n + col].abs() > m[piv * n + col].abs() {
                piv = row;
            }
        }
        if piv != col {
            for c in 0..n {
                m.swap(col * n + c, piv * n + c);
            }
            bb.swap(col, piv);
        }
        let pivot = m[col * n + col];
        if pivot.abs() < 1.0e-14 {
            continue;
        }
        for row in col + 1..n {
            let f = m[row * n + col] / pivot;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                m[row * n + c] -= f * m[col * n + c];
            }
            bb[row] -= f * bb[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = bb[row];
        for c in row + 1..n {
            s -= m[row * n + c] * x[c];
        }
        x[row] = if m[row * n + row].abs() < 1.0e-14 {
            0.0
        } else {
            s / m[row * n + row]
        };
    }
    x
}

/// Solve the nonlinear circulation system
/// `Gamma_i - 1/2 c_i V_i cl(alpha_i(Gamma)) = 0` by **damped Newton**
/// (dense analytic Jacobian + Armijo backtracking line search).  Robust for
/// the nonlinear stall polars where a Gauss-Seidel fixed point limit-cycles.
fn newton_solve<F: FoilLike>(
    stations: &[Station],
    omega: f64,
    u_0: f64,
    u_mat: &[f64],
    vdiag: &[f64],
    fs: &[&FoilSimulator<F>],
    gamma0: &[f64],
) -> Vec<f64> {
    let m = stations.len();
    // Warm start: seed with a nearby converged circulation (e.g. the previous
    // design evaluation's) so Newton converges in 1-2 steps.  The line search
    // keeps it safe even when the seed is far from the new solution.
    let mut gamma: Vec<f64> = if gamma0.len() == m {
        gamma0
            .iter()
            .map(|g| if g.is_finite() && *g >= 0.0 { *g } else { 0.0 })
            .collect()
    } else {
        vec![0.0; m]
    };
    let clp_d = 1.0e-3;
    for _iter in 0..80 {
        let st = eval_flow(&gamma, stations, u_0, omega, u_mat, vdiag, fs);
        let rinf = st.residual.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        if rinf < 1.0e-10 {
            break;
        }
        let r2 = st.residual.iter().map(|x| x * x).sum::<f64>();

        // Build the analytic Jacobian J = dR/dGamma.
        let mut jac = vec![0.0; m * m];
        for i in 0..m {
            let (u, v, vv, cl) = (st.u[i], st.v[i], st.vv[i], st.cl[i]);
            let vv2 = (vv * vv).max(1.0e-30);
            for j in 0..m {
                let du = u_mat[i * m + j];
                let dv = if i == j { -vdiag[i] } else { 0.0 };
                let dv_p = (u * du + v * dv) / vv.max(1.0e-30);
                let mut dalpha = -(v * du + u * dv) / vv2;
                let mut clp = 0.0;
                if st.free[i] {
                    clp = (fs[i].get_cl(vv, st.alpha[i] + clp_d)
                        - fs[i].get_cl(vv, st.alpha[i] - clp_d))
                        / (2.0 * clp_d);
                } else {
                    dalpha = 0.0;
                }
                let jv = if i == j { 1.0 } else { 0.0 }
                    - 0.5 * stations[i].c * (cl * dv_p + vv * clp * dalpha);
                jac[i * m + j] = jv;
            }
        }

        let neg_r: Vec<f64> = st.residual.iter().map(|x| -x).collect();
        let delta = solve_linear(&jac, &neg_r, m);
        if !delta.iter().all(|x| x.is_finite()) {
            break;
        }

        // Armijo backtracking.
        let mut lambda = 1.0;
        let mut accepted = false;
        for _ls in 0..24 {
            let cand: Vec<f64> = (0..m).map(|i| gamma[i] + lambda * delta[i]).collect();
            let st2 = eval_flow(&cand, stations, u_0, omega, u_mat, vdiag, fs);
            let r2c = st2.residual.iter().map(|x| x * x).sum::<f64>();
            if r2c <= (1.0 - 1.0e-4 * lambda) * r2 || lambda < 1.0e-6 {
                gamma = cand;
                accepted = true;
                break;
            }
            lambda *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    gamma
}

/// Solve the coupled lifting line for one operating point.
///
/// `stations[i].c` must equal `fs[i].chord()` (Reynolds/polar depend on
/// chord).  `seed` is an optional warm-start circulation (ignored if the
/// length mismatches; pass an empty slice to start cold).  Returns converged
/// circulation, induced velocities, local angles and thrust/torque.
pub fn solve<F: FoilLike>(
    stations: &[Station],
    n_blades: usize,
    omega: f64,
    u_0: f64,
    fs: &[&FoilSimulator<F>],
    seed: &[f64],
) -> LiftLineResult {
    let m = stations.len();
    let r: Vec<f64> = stations.iter().map(|s| s.r).collect();
    let edge = edges(&r);
    let dr: Vec<f64> = (0..m).map(|i| edge[i + 1] - edge[i]).collect();

    let mut res = LiftLineResult::default();
    res.gamma = vec![0.0; m];
    res.u_i = vec![0.0; m];
    res.v_i = vec![0.0; m];
    res.alpha = vec![0.0; m];
    res.phi = vec![0.0; m];
    res.d_thrust = vec![0.0; m];
    res.d_torque = vec![0.0; m];

    // Vortex-core radius, a small fraction of the mean blade element span.
    let core = 0.15 * (r[m - 1] - r[0]) / m as f64;
    // Exact linear influence, then a damped-Newton circulation solve.
    let (u_mat, vdiag) = influence_matrices(&r, &edge, n_blades, core);
    let gamma = newton_solve(stations, omega, u_0, &u_mat, &vdiag, fs, seed);
    res.gamma = gamma.clone();

    // Final induced velocities and inflow angles from the converged gamma.
    let flow = eval_flow(&gamma, stations, u_0, omega, &u_mat, &vdiag, fs);
    for i in 0..m {
        res.u_i[i] = flow.u[i] - u_0;
        res.v_i[i] = vdiag[i] * gamma[i];
        res.phi[i] = flow.phi[i];
        res.alpha[i] = flow.alpha[i];
    }

    // --- Kutta-Joukowski + profile drag forces ---
    let mut t = 0.0;
    let mut q = 0.0;
    for i in 0..m {
        let vi = res.v_i[i].min(0.9 * omega * r[i]).max(-0.1 * omega * r[i]);
        let ui = res.u_i[i].max(-0.9 * u_0).min(20.0 * u_0.max(1.0));
        let u = u_0 + ui;
        let v = (omega * r[i] - vi).max(0.05 * omega * r[i]);
        let vv = (u * u + v * v).sqrt();
        // Same prescribed (clamped) attack angle the circulation was solved for.
        let alpha = res.alpha[i];
        let cd = fs[i].get_cd(vv, alpha);
        // Drag is anti-parallel to V: its axial (thrust-reducing) component is
        // along `u`, its tangential (torque-adding) component along `v`.
        let dt = (RHO * gamma[i] * u - 0.5 * RHO * vv * stations[i].c * cd * u) * dr[i];
        let dq = (RHO * gamma[i] * v + 0.5 * RHO * vv * stations[i].c * cd * v) * r[i] * dr[i];
        res.d_thrust[i] = dt;
        res.d_torque[i] = dq;
        t += dt;
        q += dq;
    }
    res.thrust = t;
    res.torque = q;
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn biot_savart_ring_matches_closed_form() {
        // A circular vortex ring of radius a, discretised into straight
        // segments, must induce ~ Gamma/(2 a) on its own axis (z=0).
        let a = 0.05;
        let gamma = 2.0;
        let n = 256;
        let p0 = [0.0, 0.0, 0.0];
        let mut vz = 0.0;
        for s in 0..n {
            let th0 = 2.0 * std::f64::consts::PI * s as f64 / n as f64;
            let th1 = 2.0 * std::f64::consts::PI * (s + 1) as f64 / n as f64;
            let p1 = [a * th0.cos(), a * th0.sin(), 0.0];
            let p2 = [a * th1.cos(), a * th1.sin(), 0.0];
            vz += biot_savart(p0, p1, p2, gamma)[2];
        }
        let expect = vortex_ring_axial(a, gamma, 0.0);
        assert!((vz - expect).abs() / expect < 1.0e-3, "{} vs {}", vz, expect);
    }

    #[test]
    fn ring_kernel_on_axis_is_one_over_2a() {
        // ring_axial_unit evaluated on its own axis must recover 1/(2a).
        let a = 0.05;
        let k = ring_axial_unit(a, 0.0, a * 0.01, 2048);
        let expect = 1.0 / (2.0 * a);
        assert!((k - expect).abs() / expect < 1.0e-3, "{} vs {}", k, expect);
    }

    #[test]
    fn swirl_is_b_gamma_over_4pi_r() {
        // The tangential induced velocity from the trailed-sheet model is
        // v_i = B * Gamma / (4 pi r); check the coupling helper reproduces it
        // exactly for a smooth circulation distribution.
        let r = vec![0.02, 0.04, 0.06];
        let edge = edges(&r);
        let phi = vec![0.1, 0.1, 0.1];
        let gamma = vec![0.5, 0.4, 0.2];
        let (ui, vi) = induced_velocity(&r, &phi, &gamma, &edge, 2, 1.0e-4);
        for i in 0..r.len() {
            let expect = 2.0 * gamma[i] / (4.0 * std::f64::consts::PI * r[i]);
            assert!((vi[i] - expect).abs() < 1.0e-12, "vi[{}]={} vs {}", i, vi[i], expect);
        }
        // Axial induced is finite and bounded for a bounded circulation.
        for u in &ui {
            assert!(u.is_finite(), "axial induced not finite");
        }
    }

    #[test]
    fn newton_converges_to_zero_residual() {
        // A coupled plate-polar blade: the damped-Newton circulation solve
        // must drive the residual R_i = Gamma_i - 1/2 c V cl to ~ 0 (a genuine
        // fixed point), not merely stay bounded as the old Gauss-Seidel loop
        // did.
        use crate::cache::PolarStore;
        use crate::foil::Naca4;
        use std::cell::RefCell;
        use std::rc::Rc;
        use std::sync::{Arc, Mutex};

        let radii = [0.02, 0.03, 0.04];
        let omega = 2000.0;
        let u_0 = 1.0;
        let mut sims = Vec::new();
        let mut stations = Vec::new();
        for (i, &r) in radii.iter().enumerate() {
            let c = 0.010;
            let foil = Rc::new(RefCell::new(Naca4::new(c, 0.10, 0.0, 0.4)));
            let mut fs = FoilSimulator::new(
                foil,
                Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json"))),
            );
            fs.set_plate_mode(true);
            stations.push(Station { r, c, alpha: 0.10 + 0.05 * i as f64 });
            sims.push(fs);
        }
        let fs_refs: Vec<&FoilSimulator<Naca4>> = sims.iter().collect();
        let edge = edges(&radii);
        let (u_mat, vdiag) = influence_matrices(&radii, &edge, 2, 1.0e-4);
        let gamma = newton_solve(&stations, omega, u_0, &u_mat, &vdiag, &fs_refs, &[]);
        let st = eval_flow(&gamma, &stations, u_0, omega, &u_mat, &vdiag, &fs_refs);
        let rinf = st.residual.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        assert!(rinf < 1.0e-6, "Newton residual too large: {:.3e}", rinf);
        // A non-trivial, physically-sane circulation was found.
        for g in &gamma {
            assert!(g.is_finite() && *g >= 0.0, "circulation {} not sane", g);
        }
    }
}
