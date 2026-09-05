// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
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
//!    rings at the station edges.  A *tight* wake (hover: helix pitch
//!    small vs the tip radius) azimuth-averages into a stack of ring
//!    planes one pitch apart, and a *loose* wake (cruise) into the
//!    single-plane legacy ring; the two are blended on the wake-pitch
//!    ratio (see [`wake_tightness`], [`Influence::axial_matrix`]);
//!  * the wake pitch is *prescribed* from the momentum slipstream of the
//!    converged load — never from pointwise induced velocities inside the
//!    solve (per-station `ui -> pitch` feedback is unstable at tight
//!    pitch; only the disk-integral thrust may set the pitch);
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

/// Legacy loose-wake prefactor for the axial trailed-ring induction:
/// `u_i = AXIAL_K * B/(2π) * Σ dg * K`.  Calibrated (against the corrected
/// Kutta–Joukowski forces) so a lightly loaded blade's disk-averaged `u_i`
/// matches the actuator-disk value `w/2`, `w = -u_0 + sqrt(u_0² + 2T/(ρA))`
/// — see the `axial_induction_matches_actuator_disk` test.
///
/// With the wake-pitch model (see [`Influence::axial_matrix`]) this
/// prefactor scales the *loose-wake* branch of the regime blend: the
/// legacy single rotor-plane ring at the [`WAKE_OFFSET_K`] offset, used
/// when the wake helix pitch is large compared with the tip radius.  A
/// constant cannot serve the tight/hover branch too: a semi-infinite
/// helix of pitch `p` induces `B dg / (2 p)` at the rotor plane — a
/// `1/tan(phi)`-shaped gain, an order of magnitude larger at hover pitch —
/// so the tight branch carries its own calibration ([`STACK_SCALE`]).
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

/// The trailed-wake axial induction is the azimuth average of the frozen
/// helical wake: each circulation jump trails a semi-infinite helix of
/// pitch `p_j = 2 pi a_j tan(phi_j)` (`tan(phi_j) = u_wake_j / (omega
/// a_j)`, the wake's axial convection).  Azimuth-averaging a tight helix
/// gives a stack of ring planes one pitch apart; a *single* ring in the
/// rotor plane (the pre-pitch model) instead gives the interior of the
/// ring the wrong sign — upwash inboard of the loaded band — which made
/// hover inflow reverse at the hub (`ui/(omega r) ~ -1`).  The wake
/// convection is the momentum far-wake value `u_wake_j = u_0 + 2 u_i`
/// (disk average of the far-wake axial velocity), which makes the pitch
/// a fixed point of the induced flow (see [`solve_with_influence`]).
/// Ring planes are stacked from `z = (n + 1/2) p_j` (the blade sits
/// between turns) out to [`WAKE_REACH`] times the tip radius; beyond that
/// the field has decayed below the calibration noise.
///
/// Invariant (what makes the fixed point single-basin): the wake pitch is
/// prescribed ONLY from the disk-integral momentum of the converged load
/// (one scalar `u_p` per solve round).  Pointwise `ui -> pitch` feedback
/// inside a solve is unstable — it has two attractors (a negative root
/// state and a ~14x oscillation) at tight pitch, where the induction gain
/// exceeds the circulation gain.  Any future refinement (radially varying
/// slipstream, per-edge pitch) must likewise be prescribed from the
/// previous round's converged state, never iterated pointwise within a
/// round.
const WAKE_REACH: f64 = 25.0;
/// Minimum number of ring planes in the stack per trailed edge (also the
/// loose-pitch guard: the far planes of a stretched helix still carry the
/// smooth part of the induction).
const WAKE_PLANES_MIN: usize = 24;
/// Ring-plane z grid samples for the precomputed kernel table (log
/// spaced); per-solve stacks interpolate on it, so the pitch-dependent
/// influence costs O(m^2) per outer iteration, not fresh Biot-Savart.
const N_ZGRID: usize = 72;
/// Clamp on the wake helix angle tangent (the fixed point can wander into
/// absurd territory on non-physical branches).
const TANP_MIN: f64 = 0.02;
const TANP_MAX: f64 = 1.5;

/// Tightness transition of the wake model: the turn-stack (tight-helix
/// vortex-cylinder) is used when the wake pitch `p = 2 pi u_wake/omega`
/// is small compared with the tip radius, the legacy rotor-plane ring
/// (empirically calibrated against the actuator disk in cruise) when the
/// wake is loose, with a smoothstep blend in between.  `p/R` equals
/// `2 pi tan(phi_tip)`; hover wakes run ~0.5-0.9, the cruise calibration
/// point ~2.2.
const WAKE_LOOSE_LO: f64 = 0.9;
const WAKE_LOOSE_HI: f64 = 1.6;

/// Smoothstep tightness weight: 1 (full turn-stack) for `x <= lo`, 0
/// (legacy rotor-plane ring only) for `x >= hi`.
fn wake_tightness(x: f64) -> f64 {
    let t = ((x - WAKE_LOOSE_LO) / (WAKE_LOOSE_HI - WAKE_LOOSE_LO)).clamp(0.0, 1.0);
    1.0 - t * t * (3.0 - 2.0 * t)
}

/// Calibration scale of the turn-stack induction (the tight/hover
/// regime): the legacy prefactor AXIAL_K was calibrated for the loose
/// single-plane ring in cruise and under-induces the multi-plane stack
/// ~5-6x at the hover actuator-disk point; the stack part of the kernel
/// carries this extra scale.
const STACK_SCALE: f64 = 4.3;

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

/// Index of the first station at which the hub-loss ramp has fully
/// recovered (`hub_loss == 1`); stations before it carry little or no
/// circulation.
fn first_loaded_index(r: &[f64]) -> usize {
    r.iter()
        .position(|&rr| hub_loss(rr, r[0], r[r.len() - 1]) >= 1.0)
        .unwrap_or(0)
}

/// Disk area used for the momentum slipstream accounting: the *loaded*
/// annulus from the first fully-loaded station to the tip.  The hub-loss
/// ramp (inner [`HUB_RAMP_FRAC`] of the span) carries negligible load, so
/// counting its area as if it pushed air would dilute the momentum
/// slipstream ~10% in hover.  Shared by the wake-pitch prescription in
/// [`solve_with_influence`] and by the actuator-disk calibration tests, so
/// the two cannot drift apart.
fn momentum_area(r: &[f64]) -> f64 {
    let r_out = r[r.len() - 1];
    std::f64::consts::PI * (r_out.powi(2) - r[first_loaded_index(r)].powi(2))
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
    let k = gamma / (4.0 * std::f64::consts::PI) * (dot(r1, l) / n1 - dot(r2, l) / n2) / denom
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
                let v0 = (u_0 * u_0 + (omega * stations[i].r).powi(2)).sqrt();
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
                let jv = if i == j { 1.0 } else { 0.0 } - f * 0.5 * stations[i].c * cl * dv_p;
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
    /// Station radii (hub -> tip).
    r: Vec<f64>,
    /// Trailed-edge radii (m+1 of them, from [`edges`]).
    edge: Vec<f64>,
    /// Unit-ring axial kernels at the legacy rotor-plane offset
    /// (`z = WAKE_OFFSET_K * r_ci`), `kmat[ci*(m+1)+k]` at edge `k` — the
    /// seed/initial guess of the wake-pitch fixed point (see
    /// [`Influence::matrices`]).
    kmat: Vec<f64>,
    /// Ring-plane z grid (log spaced) of the pitch-dependent kernel table.
    zgrid: Vec<f64>,
    /// Unit-ring axial kernels on the z grid:
    /// `kmat_z[(ci*(m+1)+k)*N_ZGRID + nz]` = ring at edge `k`, station
    /// `ci`, plane height `zgrid[nz]`.
    kmat_z: Vec<f64>,
    pub vdiag: Vec<f64>,
    n_blades: usize,
    dr: Vec<f64>,
}

impl Influence {
    /// The axial influence matrix of the legacy single-plane rotor-plane
    /// model (see [`AXIAL_K`]): the ring-pair combination used to seed the
    /// wake-pitch fixed point, and what the cfg(test) twin
    /// [`induced_velocity`] models.
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

    /// The axial influence matrix of the pitch-resolved wake model: each
    /// trailed edge `j` sheds a semi-infinite helix whose azimuth average
    /// is a stack of ring planes one pitch apart, `z = (n + 1/2) p_j`
    /// with `p_j = 2 pi a_j tan(phi_j)` and `tanp[j] = tan(phi_j)` (see
    /// the module docs).  The table is precomputed on the z grid, so one
    /// matrix build is O(m^2 * planes) — cheap enough for an outer
    /// iteration of the wake-pitch fixed point.
    pub fn axial_matrix(&self, tanp: &[f64]) -> Vec<f64> {
        let m = self.vdiag.len();
        let pref = AXIAL_K * self.n_blades as f64 / (2.0 * std::f64::consts::PI);
        let reach = WAKE_REACH * self.r[self.r.len() - 1];
        // Wake tightness from the tip edge's pitch ratio p/R_tip =
        // 2 pi tan(phi_tip): tight wakes use the turn-stack (vortex
        // cylinder), loose wakes the legacy rotor-plane ring.
        let ratio = 2.0
            * std::f64::consts::PI
            * tanp[m].clamp(TANP_MIN, TANP_MAX)
            * self.edge[m]
            / self.r[m - 1];
        let tight = wake_tightness(ratio);
        // Stacked kernel per (station, edge): sum over the edge helix's
        // ring planes, interpolated on the z table.
        let mut s = vec![0.0; m * (m + 1)];
        for ci in 0..m {
            let row = ci * (m + 1);
            for j in 0..=m {
                let p = 2.0
                    * std::f64::consts::PI
                    * self.edge[j]
                    * tanp[j].clamp(TANP_MIN, TANP_MAX);
                let n_max = (reach / p.max(1.0e-12)).ceil() as usize;
                let n_max = n_max.max(WAKE_PLANES_MIN);
                let mut sum = 0.0;
                for n in 0..n_max {
                    let z = (n as f64 + 0.5) * p;
                    if z > reach {
                        break;
                    }
                    sum += self.z_interp(row + j, z);
                }
                // Blend with the legacy rotor-plane ring (loose wakes).
                let legacy = self.kmat[row + j];
                s[row + j] = tight * STACK_SCALE * sum + (1.0 - tight) * legacy;
            }
        }
        let mut u = vec![0.0; m * m];
        for ci in 0..m {
            let row = ci * (m + 1);
            for j in 0..m {
                u[ci * m + j] = pref * (s[row + j + 1] - s[row + j]);
            }
        }
        u
    }

    /// Linear interpolation (in log z) of the precomputed ring-kernel
    /// table for one (station, edge) pair at plane height `z`.
    fn z_interp(&self, base: usize, z: f64) -> f64 {
        let z0 = self.zgrid[0];
        let z1 = self.zgrid[N_ZGRID - 1];
        if z <= z0 {
            return self.kmat_z[base * N_ZGRID];
        }
        if z >= z1 {
            return self.kmat_z[base * N_ZGRID + N_ZGRID - 1];
        }
        let t = (z.ln() - z0.ln()) / (z1.ln() - z0.ln()) * (N_ZGRID - 1) as f64;
        let i = (t.floor() as usize).min(N_ZGRID - 2);
        let f = t - i as f64;
        let lo = self.kmat_z[base * N_ZGRID + i];
        let hi = self.kmat_z[base * N_ZGRID + i + 1];
        lo + f * (hi - lo)
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
    // Ring-plane z table for the pitch-resolved wake model: log spaced
    // from a small fraction of the mean element width (the tightest
    // realistic hover wake's half-pitch) to WAKE_REACH tip radii.
    let zmin = 0.05 * dr.iter().cloned().fold(f64::INFINITY, f64::min);
    let zmax = WAKE_REACH * r[m - 1];
    let zgrid: Vec<f64> = (0..N_ZGRID)
        .map(|n| {
            let t = n as f64 / (N_ZGRID - 1) as f64;
            (zmin * (zmax / zmin).powf(t)).max(1.0e-9)
        })
        .collect();
    let mut kmat_z = vec![0.0; m * (m + 1) * N_ZGRID];
    for (ci, &rc) in r.iter().enumerate() {
        for (j, &ej) in edge.iter().enumerate() {
            let base = (ci * (m + 1) + j) * N_ZGRID;
            for (nz, &z) in zgrid.iter().enumerate() {
                kmat_z[base + nz] = ring_axial_unit(ej.max(1.0e-6), rc, z, core, 256);
            }
        }
    }
    let vdiag: Vec<f64> = (0..m)
        .map(|i| n_blades as f64 / (4.0 * std::f64::consts::PI * r[i].max(1.0e-6)))
        .collect();
    Influence {
        r: r.to_vec(),
        edge,
        kmat,
        zgrid,
        kmat_z,
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

    let mut res = LiftLineResult {
        gamma: vec![0.0; m],
        u_i: vec![0.0; m],
        v_i: vec![0.0; m],
        alpha: vec![0.0; m],
        phi: vec![0.0; m],
        d_thrust: vec![0.0; m],
        d_torque: vec![0.0; m],
        ..Default::default()
    };

    // Wake-pitch model (see the module docs): the axial influence of the
    // trailed wake is the azimuth average of the frozen helices, whose
    // pitch follows the wake convection.  The convection is the momentum
    // slipstream of the *converged load* (uniform across the disk), so the
    // pitch is prescribed from the previous round's thrust:
    //   w   = -u_0 + sqrt(u_0^2 + 2 T / (rho A))      (far-wake velocity)
    //   u_p = u_0 + w/2                                 (disk convection)
    //   tan(phi_j) = u_p / (omega a_j)
    // A few re-solve rounds converge on T (the circulation is nearly
    // induction-independent: alpha is prescribed).  The regime blend (see
    // [`wake_tightness`]) weights the turn-stack (tight/hover wakes)
    // against the legacy rotor-plane ring (loose wakes) on the wake pitch
    // ratio x = 2 pi u_p / (omega R_tip).
    let area = momentum_area(&r);
    let mut u_mat = infl.matrices();
    let mut gamma = newton_solve(stations, omega, u_0, &u_mat, &infl.vdiag, fs, seed);
    let mut seed_g = gamma.clone();
    let axial_force = |gamma: &[f64]| -> f64 {
        let nb = n_blades as f64;
        (0..m)
            .map(|i| {
                let v = omega * r[i] - infl.vdiag[i] * gamma[i];
                nb * RHO * gamma[i] * v.max(0.0) * infl.dr[i]
            })
            .sum()
    };
    let mut thrust_est = axial_force(&gamma);
    for _ in 0..4 {
        let w = (-u_0 + (u_0 * u_0 + 2.0 * thrust_est.max(0.0) / (RHO * area)).sqrt()).max(0.0);
        // The 0.3 m/s floor on the disk convection only guards the
        // zero-thrust edge of the design space (e.g. feathered descent)
        // against a zero-pitch singularity; it is far below every loaded
        // operating point.
        let u_p = (u_0 + 0.5 * w).max(0.3);
        let tanp: Vec<f64> = infl
            .edge
            .iter()
            .map(|&e| (u_p / (omega * e.max(1.0e-9))).clamp(TANP_MIN, TANP_MAX))
            .collect();
        let um = infl.axial_matrix(&tanp);
        let gm = newton_solve(stations, omega, u_0, &um, &infl.vdiag, fs, &seed_g);
        let t_new = axial_force(&gm);
        u_mat = um;
        gamma = gm.clone();
        seed_g = gm;
        if (t_new - thrust_est).abs() <= 1.0e-3 * thrust_est.abs().max(1.0e-6) {
            break;
        }
        thrust_est = t_new;
    }
    res.gamma = gamma.clone();

    // Final induced velocities and inflow angles from the converged gamma.
    let flow = eval_flow(&gamma, stations, u_0, omega, &u_mat, &infl.vdiag, fs);
    for (i, g) in gamma.iter().enumerate() {
        res.u_i[i] = flow.u[i] - u_0;
        res.v_i[i] = infl.vdiag[i] * g;
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
    fn plate_setup(
        radii: &[f64],
        chord: f64,
        alpha: f64,
    ) -> (
        Vec<Station>,
        Vec<crate::simulator::FoilSimulator<crate::foil::Naca4>>,
    ) {
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
        assert!(
            (vz - expect).abs() / expect < 1.0e-3,
            "{} vs {}",
            vz,
            expect
        );
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
            assert!(
                (vi[i] - expect).abs() < 1.0e-12,
                "vi[{}]={} vs {}",
                i,
                vi[i],
                expect
            );
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
            stations.push(Station {
                r,
                c,
                alpha: 0.10 + 0.05 * i as f64,
            });
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
        assert!(
            (res.thrust - t_ref).abs() < 1.0e-9,
            "{} vs {}",
            res.thrust,
            t_ref
        );
        assert!(
            (res.torque - q_ref).abs() < 1.0e-9,
            "{} vs {}",
            res.torque,
            q_ref
        );
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
    fn hover_axial_induction_matches_actuator_disk() {
        // The hover twin of the cruise calibration: the same plate blade
        // at u0 = 0.  The wake is a tight helix here, so the disk-averaged
        // ui of the converged state must match the actuator-disk value
        // u_disk = w/2, w = sqrt(2T/(rho A)), over the loaded stations.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.006, 0.06);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(6000.0);
        let res = solve(&stations, 2, omega, 0.0, &fs_refs, &[]);
        let i0 = first_loaded_index(&radii);
        let w = (2.0 * res.thrust / (RHO * momentum_area(&radii))).sqrt();
        let u_disk = 0.5 * w;
        let loaded: Vec<f64> = (i0..radii.len()).map(|i| res.u_i[i]).collect();
        let u_avg = loaded.iter().sum::<f64>() / loaded.len() as f64;
        assert!(
            u_avg > 0.0 && u_avg.is_finite(),
            "hover mean axial induction {} not positive/finite",
            u_avg
        );
        assert!(
            (u_avg / u_disk - 1.0).abs() < 0.25,
            "hover mean axial induction {} vs actuator-disk {} (thrust {})",
            u_avg,
            u_disk,
            res.thrust
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
        let i0 = first_loaded_index(&radii);
        let w = -u_0 + (u_0 * u_0 + 2.0 * res.thrust / (RHO * momentum_area(&radii))).sqrt();
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

    #[test]
    fn axial_induction_tracks_momentum_across_wake_regimes() {
        // Sweep u_0 from hover to cruise on the fixed plate blade: the
        // wake-pitch ratio x = 2 pi u_p / (omega R_tip) runs from the tight
        // regime (hover) through the WAKE_LOOSE_LO..HI blend band to the
        // loose regime (cruise).  The two actuator-disk tests pin the band
        // ends; this one probes the middle, where the smoothstep mixes the
        // turn-stack and the legacy ring — the operating territory of
        // mid-size props (honda_gx35 sits near x ~ 1.05).  The disk-averaged
        // ui must stay within a factor 0.5-1.5 of the momentum slipstream
        // w/2 at EVERY step; a blend that sagged or overshot in the band
        // would show up here as a ratio excursion.
        let radii: Vec<f64> = (0..6).map(|i| 0.02 + 0.01 * i as f64).collect();
        let (stations, sims) = plate_setup(&radii, 0.006, 0.06);
        let fs_refs: Vec<&FoilSimulator<crate::foil::Naca4>> = sims.iter().collect();
        let omega = crate::optimize::rpm2omega(6000.0);
        let i0 = first_loaded_index(&radii);
        let area = momentum_area(&radii);
        for u_0 in [0.0, 1.0, 2.5, 4.0, 5.5, 7.0, 8.5, 10.0, 12.0, 15.0] {
            let res = solve(&stations, 2, omega, u_0, &fs_refs, &[]);
            let w = (-u_0 + (u_0 * u_0 + 2.0 * res.thrust / (RHO * area)).sqrt()).max(0.0);
            let u_disk = 0.5 * w;
            let loaded: Vec<f64> = (i0..radii.len()).map(|i| res.u_i[i]).collect();
            let u_avg = loaded.iter().sum::<f64>() / loaded.len() as f64;
            let x = 2.0 * std::f64::consts::PI * (u_0 + u_disk) / (omega * radii[radii.len() - 1]);
            let ratio = u_avg / u_disk.max(1.0e-9);
            assert!(
                (0.5..1.5).contains(&ratio),
                "u_0={}: ui_avg {} vs w/2 {} (pitch ratio x={:.2}, thrust {}) \
                 — mid-regime blend outside 0.5-1.5",
                u_0,
                u_avg,
                u_disk,
                x,
                res.thrust
            );
        }
    }
}
