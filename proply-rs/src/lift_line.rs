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
//!  * the trailed filaments are azimuthally averaged into coaxial vortex
//!    rings at the station edges; the rotor-plane calibration of the
//!    helix pitch is folded into a single constant ([`AXIAL_K`]), and the
//!    ring field is evaluated a small axial distance downstream
//!    ([`WAKE_OFFSET_K`]) so the blade does not sit in the singular plane
//!    of its own wake;
//!  * the induced velocity at a station is the Biot–Savart sum (straight
//!    segments) over every blade's trailed rings;
//!  * `Gamma` is solved to self-consistency with the polar `cl(alpha)`
//!    (`Gamma = 1/2 c V cl`) by damped Newton, with a hub-loss factor
//!    tapering the circulation to zero at the root (see [`hub_loss`]);
//!  * forces follow Kutta–Joukowski plus profile drag: lift is perpendicular
//!    to the relative velocity, so thrust goes with the tangential speed `v`
//!    and torque with the axial speed `u`, summed over all `B` blades.
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
/// AXIAL_K * B/(2π) * Σ dg * K`.  Calibrated (against the corrected
/// Kutta–Joukowski forces) so a lightly loaded blade's disk-averaged `u_i`
/// matches the actuator-disk value `w/2`, `w = -u_0 + sqrt(u_0² + 2T/(ρA))`
/// — see the `axial_induction_matches_actuator_disk` test.
///
/// Known limitation: a rigorous trailed-helix model carries a per-edge
/// pitch factor `1/tan(phi)` (a semi-infinite helix of pitch
/// `2 pi a tan(phi)` induces `B dg / (2 p)` at the rotor plane), which no
/// constant can reproduce across operating points — this constant is
/// effectively `1/tan(phi)` at the calibration point and under-induces far
/// from it (e.g. in hover).  Iterating the true per-edge pitch inside the
/// design loop made the frozen-pitch Newton subproblems ill-posed (the
/// induction gain exceeds the circulation gain at tight pitch), so the
/// constant is kept until that coupling is solved simultaneously.
pub const AXIAL_K: f64 = 2.22;

/// Magnitude bound on the prescribed attack angle (rad).  20 deg = the range
/// the degree-9 polar fits are fitted over; beyond it the fit is an
/// unreliable extrapolation (`get_cl` only falls back to the flat-plate model
/// past 30 deg).
pub const ALPHA_MAX: f64 = 0.35;

/// Axial offset, as a fraction of the local radius, at which the trailed-ring
/// field is evaluated: `z = WAKE_OFFSET_K * r`.  The trailed helix convects
/// downstream, so the blade does not sit in the plane of its own wake; this
/// offset is the constant stand-in for that (the exact per-edge offset
/// `r*tan(phi)` is pitch-dependent, which the [`AXIAL_K`] doc rules out for
/// the frozen-pitch design loop).
///
/// Without it (z = 0) the in-plane ring kernel has a core-capped O(1/core)
/// spike at `rho ~ a`, so the straddling ring pair at every station
/// ([`Influence::matrices`]) samples the *second difference* of the bound
/// circulation with a huge gain — any station-to-station step in `Gamma`
/// (e.g. a polar Reynolds-bucket jump in `cl`) becomes a multi-m/s spike in
/// `u_i` and the inflow angle oscillates station-to-station.  The offset
/// smooths the kernel at scale `z`, leaving the far-field induction
/// structure intact.
pub const WAKE_OFFSET_K: f64 = 0.1;

/// Fraction of the blade span over which the root loading ramps up from
/// zero (the hub-loss ramp, see [`hub_loss`]).  Deliberately wide: a narrow
/// ramp concentrates the circulation curvature into the hub region, where
/// the trailed-ring kernels are least regularised, and the ramp's own
/// second difference then oscillates the innermost inflow angles.
pub const HUB_RAMP_FRAC: f64 = 0.3;

/// Hub-loss factor: the bound circulation must vanish at the root, where the
/// blade meets the hub/spinner, and ramps smoothly (smoothstep) back to full
/// loading over the first [`HUB_RAMP_FRAC`] of the span — the root analogue
/// of a tip-loss factor.  Without it the model trails a full-strength root
/// vortex ring just inside the first station, whose near-field dominates the
/// innermost stations' inflow angles (a steep, near-singular `phi` at the
/// hub that no amount of wake regularisation removes).
fn hub_loss(r: f64, r_in: f64, r_out: f64) -> f64 {
    let span = (r_out - r_in).max(1.0e-9);
    let x = ((r - r_in) / (HUB_RAMP_FRAC * span)).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

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

/// Axial velocity induced at radius `rho` (>= 0), a height `z` above the
/// plane of a unit-circulation vortex ring of radius `a`, by Biot–Savart over
/// straight segments.  On-axis this equals the closed form
/// `a^2/(2 (a^2+z^2)^(3/2))` (matches [`vortex_ring_axial`]); a small Rankine
/// core bounds the `rho ~ a, z ~ 0` peak.  `z > 0` (see [`WAKE_OFFSET_K`])
/// keeps the evaluation point out of the ring plane, where the kernel is
/// singular.
fn ring_axial_unit(a: f64, rho: f64, z: f64, core: f64, n: usize) -> f64 {
    let p0 = [rho, 0.0, z];
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
///  * axial: `u_i(r) = AXIAL_K * B/(2 pi) * sum_j dg_j * K(r, a_j)`, where
///    `K` is the axial velocity in the plane of a unit ring at radius `a_j`
///    and the prefactor maps the trailed-sheet circulation to the disk
///    inflow (calibrated against actuator-disk momentum; see [`AXIAL_K`]);
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

    // Helix-pitch calibration folded into the rotor-plane prefactor (see
    // AXIAL_K).

    let nper = 512;
    let pref = AXIAL_K * n_b / (2.0 * std::f64::consts::PI);
    let mut ui = vec![0.0; m];
    for (ci, &ri) in r.iter().enumerate() {
        for j in 0..=m {
            if dg[j].abs() < 1.0e-18 {
                continue;
            }
            let a_j = edge[j].max(1.0e-6);
            ui[ci] += pref * dg[j] * ring_axial_unit(a_j, ri, WAKE_OFFSET_K * ri, core, nper);
        }
    }
    // Swirl from the bound circulation (halved for the rotor plane).
    let tang: Vec<f64> = (0..m)
        .map(|i| n_b * gamma[i] / (4.0 * std::f64::consts::PI * r[i].max(1.0e-6)))
        .collect();
    (ui, tang)
}

/// Raw trailed-ring kernels: `kmat[ci*(m+1)+k]` is the axial velocity at
/// station `ci` in the plane of a unit-circulation ring at `edge[k]`.
/// These are pitch-independent (and the expensive part); the helix-pitch
/// factors are applied cheaply in [`Influence::matrices`].
fn ring_kernels(r: &[f64], edge: &[f64], core: f64) -> Vec<f64> {
    let m = r.len();
    let nper = 512;
    let mut kmat = vec![0.0; m * (m + 1)];
    for ci in 0..m {
        let z = WAKE_OFFSET_K * r[ci];
        for k in 0..=m {
            kmat[ci * (m + 1) + k] = ring_axial_unit(edge[k].max(1.0e-6), r[ci], z, core, nper);
        }
    }
    kmat
}

/// Per-element flow state induced by a circulation vector.
struct FlowState {
    u: Vec<f64>,
    v: Vec<f64>,
    vv: Vec<f64>,
    phi: Vec<f64>,
    alpha: Vec<f64>,
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
    let (r_in, r_out) = (stations[0].r, stations[m - 1].r);
    let mut u = vec![0.0; m];
    let mut v = vec![0.0; m];
    let mut vv = vec![0.0; m];
    let mut phi = vec![0.0; m];
    let mut alpha = vec![0.0; m];
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
        // Prescribed attack angle, clamped to the polar fit range (see
        // ALPHA_MAX).
        alpha[i] = stations[i].alpha.clamp(-ALPHA_MAX, ALPHA_MAX);
        cl[i] = fs[i].get_cl(vv[i], alpha[i]);
        // The hub-loss factor scales the circulation target (the blade
        // carries no circulation at the root).
        let f = hub_loss(stations[i].r, r_in, r_out);
        residual[i] = gamma[i] - f * 0.5 * stations[i].c * vv[i] * cl[i];
    }
    FlowState {
        u,
        v,
        vv,
        phi,
        alpha,
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
    for _iter in 0..80 {
        let st = eval_flow(&gamma, stations, u_0, omega, u_mat, vdiag, fs);
        let rinf = st.residual.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        if rinf < 1.0e-10 {
            break;
        }
        // Runaway detector: if a subproblem's induction gain exceeds the
        // circulation gain there is no finite fixed point and the solve
        // limit-cycles instead (an over-loaded candidate).  Such candidates
        // only need a deterministic, monotone "how bad is it" answer for the
        // outer optimizer — return the no-induction circulation instead of a
        // diverging iterate.
        if _iter > 20 && rinf > 1.0e-2 {
            let mut g0 = vec![0.0; m];
            for i in 0..m {
                let v0 =
                    (u_0 * u_0 + (omega * stations[i].r).powi(2)).sqrt();
                let a0 = stations[i].alpha.clamp(-ALPHA_MAX, ALPHA_MAX);
                let f = hub_loss(stations[i].r, stations[0].r, stations[m - 1].r);
                g0[i] = f * 0.5 * stations[i].c * v0 * fs[i].get_cl(v0, a0);
            }
            return g0;
        }
        let r2 = st.residual.iter().map(|x| x * x).sum::<f64>();

        // Build the analytic Jacobian J = dR/dGamma.  alpha is prescribed
        // (independent of Gamma), so cl depends on Gamma only through V; the
        // weak Reynolds dependence of the polar is neglected.  Thus
        // dR_i/dGamma_j = delta_ij - F_i * 1/2 c cl dV/dGamma_j with F the
        // hub-loss factor.
        let mut jac = vec![0.0; m * m];
        for i in 0..m {
            let (u, v, vv, cl) = (st.u[i], st.v[i], st.vv[i], st.cl[i]);
            let f = hub_loss(stations[i].r, stations[0].r, stations[m - 1].r);
            for j in 0..m {
                let du = u_mat[i * m + j];
                let dv = if i == j { -vdiag[i] } else { 0.0 };
                let dv_p = (u * du + v * dv) / vv.max(1.0e-30);
                let jv = if i == j { 1.0 } else { 0.0 }
                    - f * 0.5 * stations[i].c * cl * dv_p;
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

/// Precomputed trailed-wake influence for a fixed set of station radii and
/// blade count: the linear map from the bound circulation to the induced
/// velocities.  Expensive to build (O(m^2) ring kernels); reuse it across
/// solves whose station radii do not change (e.g. a design loop that only
/// varies chord and attack angle).
pub struct Influence {
    /// Unit-ring axial kernels, `kmat[ci*(m+1)+k]` at edge `k` (see
    /// [`ring_kernels`]).
    kmat: Vec<f64>,
    pub vdiag: Vec<f64>,
    n_blades: usize,
    dr: Vec<f64>,
}

impl Influence {
    /// The axial influence matrix `U` (`u_i[i] = sum_j U[i*m+j] Gamma_j`),
    /// matching the trailed-ring model of [`induced_velocity`]: `Gamma_j`
    /// feeds the trailed rings at edges `j` and `j+1` (the jump in bound
    /// circulation) with opposite signs, each scaled by the rotor-plane
    /// calibration [`AXIAL_K`].
    pub fn matrices(&self) -> Vec<f64> {
        let m = self.vdiag.len();
        let pref = AXIAL_K * self.n_blades as f64 / (2.0 * std::f64::consts::PI);
        let mut u = vec![0.0; m * m];
        for ci in 0..m {
            let row = ci * (m + 1);
            for j in 0..m {
                u[ci * m + j] = pref * (self.kmat[row + j + 1] - self.kmat[row + j]);
            }
        }
        u
    }
}

/// Build the [`Influence`] for station radii `r` (monotone increasing).
pub fn influence(r: &[f64], n_blades: usize) -> Influence {
    let m = r.len();
    // edges() and the trailing-vortex model need at least two stations.
    assert!(m >= 2, "lifting line needs >= 2 stations");
    let edge = edges(r);
    let dr: Vec<f64> = (0..m).map(|i| edge[i + 1] - edge[i]).collect();
    // Vortex-core radius, a small fraction of the mean blade element span.
    let core = 0.15 * (r[m - 1] - r[0]) / m as f64;
    let kmat = ring_kernels(r, &edge, core);
    let vdiag: Vec<f64> = (0..m)
        .map(|i| n_blades as f64 / (4.0 * std::f64::consts::PI * r[i].max(1.0e-6)))
        .collect();
    Influence {
        kmat,
        vdiag,
        n_blades,
        dr,
    }
}

/// Solve the coupled lifting line for one operating point, reusing a
/// precomputed [`Influence`] (see [`influence`]).  `stations[i].c` must equal
/// `fs[i].chord()` and `stations[i].r` must match the radii `infl` was built
/// with.  `seed` is an optional warm-start circulation (ignored if the length
/// mismatches; pass an empty slice to start cold).
pub fn solve_with_influence<F: FoilLike>(
    stations: &[Station],
    n_blades: usize,
    omega: f64,
    u_0: f64,
    fs: &[&FoilSimulator<F>],
    seed: &[f64],
    infl: &Influence,
) -> LiftLineResult {
    let m = stations.len();
    let r: Vec<f64> = stations.iter().map(|s| s.r).collect();
    let dr = &infl.dr;

    let mut res = LiftLineResult::default();
    res.gamma = vec![0.0; m];
    res.u_i = vec![0.0; m];
    res.v_i = vec![0.0; m];
    res.alpha = vec![0.0; m];
    res.phi = vec![0.0; m];
    res.d_thrust = vec![0.0; m];
    res.d_torque = vec![0.0; m];

    // Exact linear influence, then a damped-Newton circulation solve.
    let u_mat = infl.matrices();
    let gamma = newton_solve(stations, omega, u_0, &u_mat, &infl.vdiag, fs, seed);
    res.gamma = gamma.clone();

    // Final induced velocities and inflow angles from the converged gamma.
    let flow = eval_flow(&gamma, stations, u_0, omega, &u_mat, &infl.vdiag, fs);
    for i in 0..m {
        res.u_i[i] = flow.u[i] - u_0;
        res.v_i[i] = infl.vdiag[i] * gamma[i];
        res.phi[i] = flow.phi[i];
        res.alpha[i] = flow.alpha[i];
    }

    // --- Kutta-Joukowski + profile drag forces ---
    // Lift is perpendicular to the relative velocity V = (u, v), so its axial
    // (thrust) component goes with `v` and its tangential (torque) component
    // with `u`; drag is anti-parallel to V (axial with `u`, tangential with
    // `v`).  `gamma` is the circulation of ONE blade: scale by the blade count.
    let nb = n_blades as f64;
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
        let dt = nb * (RHO * gamma[i] * v - 0.5 * RHO * vv * stations[i].c * cd * u) * dr[i];
        let dq = nb * (RHO * gamma[i] * u + 0.5 * RHO * vv * stations[i].c * cd * v) * r[i] * dr[i];
        res.d_thrust[i] = dt;
        res.d_torque[i] = dq;
        t += dt;
        q += dq;
    }
    res.thrust = t;
    res.torque = q;
    res
}

/// Solve the coupled lifting line for one operating point (convenience
/// wrapper that builds the [`Influence`] first; see [`solve_with_influence`]).
pub fn solve<F: FoilLike>(
    stations: &[Station],
    n_blades: usize,
    omega: f64,
    u_0: f64,
    fs: &[&FoilSimulator<F>],
    seed: &[f64],
) -> LiftLineResult {
    let r: Vec<f64> = stations.iter().map(|s| s.r).collect();
    let infl = influence(&r, n_blades);
    solve_with_influence(stations, n_blades, omega, u_0, fs, seed, &infl)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A plate-polar blade: `alpha` is uniform and below the clamp, chord
    /// uniform, so the converged solution is easy to reason about.
    fn plate_setup(radii: &[f64], chord: f64, alpha: f64) -> (Vec<Station>, Vec<crate::simulator::FoilSimulator<crate::foil::Naca4>>) {
        use crate::cache::PolarStore;
        use crate::foil::Naca4;
        use std::cell::RefCell;
        use std::rc::Rc;
        use std::sync::{Arc, Mutex};

        let mut sims = Vec::new();
        let mut stations = Vec::new();
        for &r in radii {
            let foil = Rc::new(RefCell::new(Naca4::new(chord, 0.10, 0.0, 0.4)));
            let mut fs = FoilSimulator::new(
                foil,
                Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json"))),
            );
            fs.set_plate_mode(true);
            stations.push(Station { r, c: chord, alpha });
            sims.push(fs);
        }
        (stations, sims)
    }

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
        // ring_axial_unit evaluated on its own axis must recover the closed
        // form vortex_ring_axial(a, 1, z) = a^2 / (2 (a^2+z^2)^(3/2)), both
        // near the ring plane and at the wake offset used by the design.
        let a = 0.05;
        for z in [a * 0.01, WAKE_OFFSET_K * a] {
            let k = ring_axial_unit(a, 0.0, z, a * 0.01, 2048);
            let expect = vortex_ring_axial(a, 1.0, z);
            assert!(
                (k - expect).abs() / expect < 1.0e-3,
                "z={}: {} vs {}",
                z,
                k,
                expect
            );
        }
    }

    #[test]
    fn hub_loss_tapers_root_circulation() {
        // The bound circulation must vanish at the root (hub-loss factor):
        // a full-strength root vortex ring just inside the first station
        // otherwise dominates the innermost stations' inflow angles with a
        // steep, near-singular phi.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.010, 0.10);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let res = solve(
            &stations,
            2,
            crate::optimize::rpm2omega(9000.0),
            5.0,
            &fs_refs,
            &[],
        );
        assert!(res.gamma[0].abs() < 1.0e-9, "root gamma {}", res.gamma[0]);
        assert!(
            res.gamma[1] > 0.0 && res.gamma[1] < res.gamma[3],
            "ramp {} vs {}",
            res.gamma[1],
            res.gamma[3]
        );
        // The factor itself is a smoothstep: mid-ramp is exactly one half.
        let r_mid = 0.02 + 0.5 * HUB_RAMP_FRAC * (0.07 - 0.02);
        let f_mid = hub_loss(r_mid, 0.02, 0.07);
        assert!((f_mid - 0.5).abs() < 1.0e-12, "mid-ramp {}", f_mid);
        assert!(hub_loss(0.02, 0.02, 0.07).abs() < 1.0e-12);
        assert!((hub_loss(0.07, 0.02, 0.07) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn kernel_does_not_amplify_circulation_steps() {
        // Regression: the in-plane (z=0) ring kernel has a core-capped
        // O(1/core) spike at rho ~ a, so the straddling ring pair sampled the
        // *second difference* of Gamma with a huge gain — a 10% step in Gamma
        // at one station (e.g. a polar Reynolds-bucket jump in cl) produced a
        // ~9 m/s spike in u_i and the inflow angle oscillated
        // station-to-station.  With the wake offset the step must induce a
        // smooth, bounded response.
        let m = 43;
        let radii: Vec<f64> = (0..m).map(|i| 0.006 + 0.069 * i as f64 / 42.0).collect();
        let infl = influence(&radii, 3);
        let u_mat = infl.matrices();

        // Smooth (elliptic-ish) circulation with a 10% step at mid-span.
        let gamma: Vec<f64> = radii
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                let x = (r - 0.006) / 0.069;
                let g = 0.12 * (1.0 - x * x);
                g * if i == 21 { 0.9 } else { 1.0 }
            })
            .collect();

        let ui: Vec<f64> = (0..m)
            .map(|i| (0..m).map(|j| u_mat[i * m + j] * gamma[j]).sum())
            .collect();
        // Skip the root-vortex hub stations (their steep u_i gradient is
        // physical, not an oscillation); the step at station 21 is mid-span.
        for w in 5..m {
            let step = (ui[w] - ui[w - 1]).abs();
            assert!(step < 1.0, "u_i step {} m/s at station {}", step, w);
        }
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
        let infl = influence(&radii, 2);
        let u_mat = infl.matrices();
        let vdiag = infl.vdiag.clone();
        let gamma = newton_solve(&stations, omega, u_0, &u_mat, &vdiag, &fs_refs, &[]);
        let st = eval_flow(&gamma, &stations, u_0, omega, &u_mat, &vdiag, &fs_refs);
        let rinf = st.residual.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        assert!(rinf < 1.0e-6, "Newton residual too large: {:.3e}", rinf);
        // A non-trivial, physically-sane circulation was found.
        for g in &gamma {
            assert!(g.is_finite() && *g >= 0.0, "circulation {} not sane", g);
        }
    }

    #[test]
    fn forces_match_blade_element_projection() {
        // The Kutta-Joukowski force on each element must equal the textbook
        // blade-element projection with cos(phi)=v/V, sin(phi)=u/V:
        //   dT = B 1/2 rho V^2 c (cl cos(phi) - cd sin(phi)) dr
        //   dQ = B r 1/2 rho V^2 c (cl sin(phi) + cd cos(phi)) dr
        // (i.e. thrust goes with the *tangential* speed, torque with the
        // *axial* speed).  This pins the projection independently of the
        // implementation in `solve`.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.010, 0.10);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(9000.0);
        let u_0 = 5.0;
        let n_blades = 2;
        let res = solve(&stations, n_blades, omega, u_0, &fs_refs, &[]);

        let edge = edges(&radii);
        let mut t_ref = 0.0;
        let mut q_ref = 0.0;
        for (i, s) in stations.iter().enumerate() {
            let u = u_0 + res.u_i[i];
            let v = omega * s.r - res.v_i[i];
            let vv = (u * u + v * v).sqrt();
            let (cosphi, sinphi) = (v / vv, u / vv);
            // The hub-loss factor scales the circulation (lift) target only;
            // profile drag still acts at zero circulation.
            let f = hub_loss(s.r, radii[0], radii[radii.len() - 1]);
            let cl = f * 2.0 * std::f64::consts::PI * res.alpha[i];
            let cd = 1.28 * res.alpha[i].sin();
            let dr = edge[i + 1] - edge[i];
            let dt = n_blades as f64 * 0.5 * RHO * vv * vv * s.c * (cl * cosphi - cd * sinphi) * dr;
            let dq = n_blades as f64
                * s.r
                * 0.5
                * RHO
                * vv
                * vv
                * s.c
                * (cl * sinphi + cd * cosphi)
                * dr;
            assert!(
                (res.d_thrust[i] - dt).abs() < 1.0e-9 * dt.abs().max(1.0e-9),
                "dT[{}] {} vs {}",
                i,
                res.d_thrust[i],
                dt
            );
            assert!(
                (res.d_torque[i] - dq).abs() < 1.0e-9 * dq.abs().max(1.0e-9),
                "dQ[{}] {} vs {}",
                i,
                res.d_torque[i],
                dq
            );
            t_ref += dt;
            q_ref += dq;
        }
        assert!((res.thrust - t_ref).abs() < 1.0e-9, "{} vs {}", res.thrust, t_ref);
        assert!((res.torque - q_ref).abs() < 1.0e-9, "{} vs {}", res.torque, q_ref);
    }

    #[test]
    fn thrust_survives_the_hover_limit() {
        // u_0 -> 0: the relative wind is almost purely tangential, so the
        // Kutta-Joukowski thrust (rho Gamma v) stays large.  A projection
        // that sends thrust with `u` instead collapses to ~rho Gamma u_i and
        // under-predicts the thrust by ~v/u (order 10-40 here).
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.010, 0.10);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(9000.0);
        let res = solve(&stations, 2, omega, 0.001, &fs_refs, &[]);
        assert!(res.thrust > 0.3, "hover thrust {} collapsed", res.thrust);
    }

    #[test]
    fn doubling_blades_doubles_thrust() {
        // gamma is per-blade: at light loading (large u_0) the converged
        // circulation barely changes with blade count, so the total forces
        // must scale with B.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.006, 0.06);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(6000.0);
        let r1 = solve(&stations, 1, omega, 15.0, &fs_refs, &[]);
        let r2 = solve(&stations, 2, omega, 15.0, &fs_refs, &[]);
        let ratio = r2.thrust / r1.thrust;
        assert!(
            (ratio - 2.0).abs() < 0.3,
            "thrust ratio {} for B=2 vs B=1",
            ratio
        );
        let qratio = r2.torque / r1.torque;
        assert!(
            (qratio - 2.0).abs() < 0.3,
            "torque ratio {} for B=2 vs B=1",
            qratio
        );
    }

    #[test]
    fn axial_induction_matches_actuator_disk() {
        // Calibration check on AXIAL_K: the disk-averaged axial induced
        // velocity of a lightly loaded blade must match the actuator-disk
        // value u_disk = w/2 with w = -u_0 + sqrt(u_0^2 + 2T/(rho A)) to
        // ~25% (the trailed-ring model is a rotor-plane approximation).
        // A constant prefactor cannot hold this at every operating point
        // (see the AXIAL_K doc), but the check trips a gross
        // mis-calibration such as the pre-correction 0.35x under-induction
        // or a missing blade-count factor.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.006, 0.06);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(6000.0);
        let u_0 = 15.0;
        let res = solve(&stations, 2, omega, u_0, &fs_refs, &[]);

        // The hub-loss ramp unloads the root annulus; compare the actuator
        // disk over the *loaded* area and average u_i over the loaded
        // stations only (scaling the whole model to a mean that includes an
        // unloaded root over-induces everything else).
        let r_out = radii[radii.len() - 1];
        let i0 = radii
            .iter()
            .position(|&r| hub_loss(r, radii[0], r_out) >= 1.0)
            .unwrap_or(0);
        let area = std::f64::consts::PI * (r_out.powi(2) - radii[i0].powi(2));
        let w = -u_0 + (u_0 * u_0 + 2.0 * res.thrust / (RHO * area)).sqrt();
        let u_disk = 0.5 * w;
        let loaded: Vec<f64> = (i0..radii.len()).map(|i| res.u_i[i]).collect();
        let u_avg = loaded.iter().sum::<f64>() / loaded.len() as f64;
        assert!(
            u_avg > 0.0 && u_avg.is_finite(),
            "mean axial induction {} not positive/finite",
            u_avg
        );
        assert!(
            (u_avg / u_disk - 1.0).abs() < 0.25,
            "mean axial induction {} vs actuator-disk {}",
            u_avg,
            u_disk
        );
    }
}
