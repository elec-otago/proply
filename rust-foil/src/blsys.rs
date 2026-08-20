// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Boundary-layer Newton system assembly (port of `m_xblsys.f90`).
//!
//! Secondary variables at each BL station live in `xf.com1` (upstream "1")
//! and `xf.com2` (current "2") as [`crate::state::BlVars`] structs; the
//! primary variables are passed by value.  The Newton-system coefficient
//! arrays `vs1`, `vs2`, `vsm`, `vsr`, `vsrez`, `vsx` are state fields.

use crate::state::Xfoil;

/// Checks whether transition occurs in the current interval; implements the
/// 2nd-order amplification equation.
pub fn trchek(xf: &mut Xfoil) -> bool {
    trchek2(xf)
}

/// Sets average amplification AX over the interval 1..2.
#[allow(clippy::too_many_arguments)]
pub fn axset(
    hk1: f64,
    t1: f64,
    rt1: f64,
    a1: f64,
    hk2: f64,
    t2: f64,
    rt2: f64,
    a2: f64,
    acrit: f64,
    idampv: i32,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    // ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2
    let (ax1, ax1_hk1, ax1_t1, ax1_rt1);
    let (ax2, ax2_hk2, ax2_t2, ax2_rt2);
    if idampv == 0 {
        let r = dampl(hk1, t1, rt1);
        ax1 = r.0;
        ax1_hk1 = r.1;
        ax1_t1 = r.2;
        ax1_rt1 = r.3;
        let r = dampl(hk2, t2, rt2);
        ax2 = r.0;
        ax2_hk2 = r.1;
        ax2_t2 = r.2;
        ax2_rt2 = r.3;
    } else {
        let r = dampl2(hk1, t1, rt1);
        ax1 = r.0;
        ax1_hk1 = r.1;
        ax1_t1 = r.2;
        ax1_rt1 = r.3;
        let r = dampl2(hk2, t2, rt2);
        ax2 = r.0;
        ax2_hk2 = r.1;
        ax2_t2 = r.2;
        ax2_rt2 = r.3;
    }

    // rms-average version
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let (axa, axa_ax1, axa_ax2);
    if axsq <= 0.0 {
        axa = 0.0;
        axa_ax1 = 0.0;
        axa_ax2 = 0.0;
    } else {
        axa = axsq.sqrt();
        axa_ax1 = 0.5 * ax1 / axa;
        axa_ax2 = 0.5 * ax2 / axa;
    }

    // small additional term to ensure dN/dx > 0 near N = Ncrit
    let arg = (20.0 * (acrit - 0.5 * (a1 + a2))).min(20.0);
    let (exn, exn_a1, exn_a2);
    if arg <= 0.0 {
        exn = 1.0;
        exn_a1 = 0.0;
        exn_a2 = 0.0;
    } else {
        exn = (-arg).exp();
        exn_a1 = 20.0 * 0.5 * exn;
        exn_a2 = 20.0 * 0.5 * exn;
    }

    let dax = exn * 0.002 / (t1 + t2);
    let dax_a1 = exn_a1 * 0.002 / (t1 + t2);
    let dax_a2 = exn_a2 * 0.002 / (t1 + t2);
    let dax_t1 = -dax / (t1 + t2);
    let dax_t2 = -dax / (t1 + t2);

    let ax = axa + dax;

    let ax_hk1 = axa_ax1 * ax1_hk1;
    let ax_t1 = axa_ax1 * ax1_t1 + dax_t1;
    let ax_rt1 = axa_ax1 * ax1_rt1;
    let ax_a1 = dax_a1;

    let ax_hk2 = axa_ax2 * ax2_hk2;
    let ax_t2 = axa_ax2 * ax2_t2 + dax_t2;
    let ax_rt2 = axa_ax2 * ax2_rt2;
    let ax_a2 = dax_a2;

    (ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2)
}

/// Second-order transition check.
///
/// Solves the implicit amplification equation for N2.  Sets transition
/// location XT and its sensitivities if transition occurs in the interval;
/// otherwise sets the amplification AMPL2.
pub fn trchek2(xf: &mut Xfoil) -> bool {
    const DAEPS: f64 = 5.0E-5;

    xf.c2sav = xf.com2;

    // calculate average amplification rate AX over X1..X2 interval
    let c1 = xf.com1;
    let (mut ax, mut ax_hk1, mut ax_t1, mut ax_rt1, mut ax_a1, _ax_hk2, _ax_t2, _ax_rt2, _ax_a2) = axset(
        c1.hk,
        c1.t,
        c1.rt,
        c1.ampl,
        xf.com2.hk,
        xf.com2.t,
        xf.com2.rt,
        xf.com2.ampl,
        xf.amcrit,
        xf.idampv,
    );

    // set initial guess for iterate N2 (AMPL2) at X2
    xf.com2.ampl = c1.ampl + ax * (xf.com2.x - c1.x);

    let x1 = c1.x;
    let t1 = c1.t;
    let d1 = c1.d;
    let u1 = c1.u;
    let ampl1 = c1.ampl;

    // variables from the Newton loop that are needed after it (Fortran keeps
    // them in subroutine scope across iterations)
    let mut wf2 = 0.0;
    let mut wf2_a1 = 0.0;
    let mut wf2_a2 = 0.0;
    let mut wf2_x1 = 0.0;
    let mut wf2_x2 = 0.0;
    let mut wf2_xf = 0.0;
    let mut wf1_a1 = 0.0;
    let mut wf1_a2 = 0.0;
    let mut wf1_x1 = 0.0;
    let mut wf1_x2 = 0.0;
    let mut wf1_xf = 0.0;
    let mut xt = 0.0;
    let mut tt = 0.0;
    let mut dt = 0.0;
    let mut ut = 0.0;
    let mut xt_a2 = 0.0;
    let mut tt_a1 = 0.0;
    let mut tt_a2 = 0.0;
    let mut dt_a1 = 0.0;
    let mut dt_a2 = 0.0;
    let mut ut_a1 = 0.0;
    let mut ut_a2 = 0.0;
    let mut tt_t1 = 0.0;
    let mut dt_d1 = 0.0;
    let mut ut_u1 = 0.0;
    let mut tt_t2 = 0.0;
    let mut dt_d2 = 0.0;
    let mut ut_u2 = 0.0;
    let mut tt_x1 = 0.0;
    let mut dt_x1 = 0.0;
    let mut ut_x1 = 0.0;
    let mut tt_x2 = 0.0;
    let mut dt_x2 = 0.0;
    let mut ut_x2 = 0.0;
    let mut tt_xf = 0.0;
    let mut dt_xf = 0.0;
    let mut ut_xf = 0.0;
    let mut amplt = 0.0;
    let mut amplt_a2 = 0.0;
    let mut hkt = 0.0;
    let mut hkt_tt = 0.0;
    let mut hkt_dt = 0.0;
    let mut hkt_ut = 0.0;
    let mut hkt_ms = 0.0;
    let mut rtt = 0.0;
    let mut rtt_tt = 0.0;
    let mut rtt_ut = 0.0;
    let mut rtt_ms = 0.0;
    let mut rtt_re = 0.0;
    let mut ax_hkt = 0.0;
    let mut ax_tt = 0.0;
    let mut ax_rtt = 0.0;
    let mut ax_at = 0.0;
    let mut da2 = 0.0;

    for _itam in 1..=30 {
        // define weighting factors WF1, WF2 for defining "T" quantities from 1,2
        if xf.com2.ampl <= xf.amcrit {
            amplt = xf.com2.ampl;
            amplt_a2 = 1.0;
            let sfa = 1.0;
            let sfa_a1 = 0.0;
            let sfa_a2 = 0.0;
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
            wf2_xf = 0.0;
        } else {
            amplt = xf.amcrit;
            amplt_a2 = 0.0;
            let sfa = (amplt - ampl1) / (xf.com2.ampl - ampl1);
            let sfa_a1 = (sfa - 1.0) / (xf.com2.ampl - ampl1);
            let sfa_a2 = -sfa / (xf.com2.ampl - ampl1);
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
            wf2_xf = 0.0;
        }

        let x2 = xf.com2.x;
        let (sfx, sfx_x1, sfx_x2, sfx_xf);
        if xf.xiforc < x2 {
            sfx = (xf.xiforc - x1) / (x2 - x1);
            sfx_x1 = (sfx - 1.0) / (x2 - x1);
            sfx_x2 = -sfx / (x2 - x1);
            sfx_xf = 1.0 / (x2 - x1);
        } else {
            sfx = 1.0;
            sfx_x1 = 0.0;
            sfx_x2 = 0.0;
            sfx_xf = 0.0;
        }

        // set weighting factor from free or forced transition
        if wf2 < sfx {
            // free (wf2 already set from sfa)
        } else {
            wf2 = sfx;
            wf2_a1 = 0.0;
            wf2_a2 = 0.0;
            wf2_x1 = sfx_x1;
            wf2_x2 = sfx_x2;
            wf2_xf = sfx_xf;
        }

        let wf1 = 1.0 - wf2;
        wf1_a1 = -wf2_a1;
        wf1_a2 = -wf2_a2;
        wf1_x1 = -wf2_x1;
        wf1_x2 = -wf2_x2;
        wf1_xf = -wf2_xf;

        // interpolate BL variables to XT
        xt = x1 * wf1 + x2 * wf2;
        tt = t1 * wf1 + xf.com2.t * wf2;
        dt = d1 * wf1 + xf.com2.d * wf2;
        ut = u1 * wf1 + xf.com2.u * wf2;

        xt_a2 = x1 * wf1_a2 + x2 * wf2_a2;
        tt_a2 = t1 * wf1_a2 + xf.com2.t * wf2_a2;
        dt_a2 = d1 * wf1_a2 + xf.com2.d * wf2_a2;
        ut_a2 = u1 * wf1_a2 + xf.com2.u * wf2_a2;

        // temporarily set "2" variables from "T" for BLKIN
        xf.com2.x = xt;
        xf.com2.t = tt;
        xf.com2.d = dt;
        xf.com2.u = ut;

        // calculate laminar secondary "T" variables HKT, RTT
        blkin(xf);

        hkt = xf.com2.hk;
        hkt_tt = xf.com2.hk_t;
        hkt_dt = xf.com2.hk_d;
        hkt_ut = xf.com2.hk_u;
        hkt_ms = xf.com2.hk_ms;

        rtt = xf.com2.rt;
        rtt_tt = xf.com2.rt_t;
        rtt_ut = xf.com2.rt_u;
        rtt_ms = xf.com2.rt_ms;
        rtt_re = xf.com2.rt_re;

        // restore clobbered "2" variables, except for AMPL2
        let amsave = xf.com2.ampl;
        xf.com2 = xf.c2sav;
        xf.com2.ampl = amsave;

        // calculate amplification rate AX over current X1-XT interval
        let (ax2, ax_hk1_2, ax_t1_2, ax_rt1_2, ax_a1_2, ax_hkt_2, ax_tt_2, ax_rtt_2, ax_at_2) = axset(
            c1.hk,
            t1,
            c1.rt,
            ampl1,
            hkt,
            tt,
            rtt,
            amplt,
            xf.amcrit,
            xf.idampv,
        );
        ax = ax2;
        ax_hk1 = ax_hk1_2;
        ax_t1 = ax_t1_2;
        ax_rt1 = ax_rt1_2;
        ax_a1 = ax_a1_2;
        ax_hkt = ax_hkt_2;
        ax_tt = ax_tt_2;
        ax_rtt = ax_rtt_2;
        ax_at = ax_at_2;

        // punch out early if there is no amplification here
        if ax <= 0.0 {
            break;
        }

        // set sensitivity of AX(A2)
        let ax_a2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2
            + (ax_hkt * hkt_dt) * dt_a2
            + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2
            + ax_at * amplt_a2;

        // residual for implicit AMPL2 definition (amplification equation)
        let res = xf.com2.ampl - ampl1 - ax * (x2 - x1);
        let res_a2 = 1.0 - ax_a2 * (x2 - x1);

        da2 = -res / res_a2;

        let mut rlx = 1.0;
        let dxt = xt_a2 * da2;

        if rlx * (dxt / (x2 - x1)).abs() > 0.05 {
            rlx = 0.05 * ((x2 - x1) / dxt).abs();
        }
        if rlx * da2.abs() > 1.0 {
            rlx = 1.0 * (1.0 / da2).abs();
        }

        // check if converged
        if da2.abs() < DAEPS {
            break;
        }

        if (xf.com2.ampl > xf.amcrit && xf.com2.ampl + rlx * da2 < xf.amcrit)
            || (xf.com2.ampl < xf.amcrit && xf.com2.ampl + rlx * da2 > xf.amcrit)
        {
            // limited Newton step so AMPL2 doesn't step across AMCRIT either way
            xf.com2.ampl = xf.amcrit;
        } else {
            // regular Newton step
            xf.com2.ampl += rlx * da2;
        }
    }

    let x2 = xf.com2.x;
    xf.xt = xt;

    if xf.show_output
        && (x1.is_nan() || xt.is_nan() || x2.is_nan() || ampl1.is_nan() || amplt.is_nan() || xf.com2.ampl.is_nan() || ax.is_nan() || da2.is_nan())
    {
        eprintln!("TRCHEK2: N2 convergence failed.");
        eprintln!("x: {:9.5} {:9.5} {:9.5}  N: {:7.3} {:7.3} {:7.3}  Nx: {:8.3}   dN: {:e}", x1, xt, x2, ampl1, amplt, xf.com2.ampl, ax, da2);
    }

    // Check if ANY of the printed variables contain NaN's. If they do,
    // convergence will never be reached.
    if (x1.is_nan()
        || xt.is_nan()
        || x2.is_nan()
        || ampl1.is_nan()
        || amplt.is_nan()
        || xf.com2.ampl.is_nan()
        || ax.is_nan()
        || da2.is_nan())
        && xf.abort_on_nan
    {
        return false;
    }

    // test for free or forced transition
    xf.trfree = xf.com2.ampl >= xf.amcrit;
    xf.trforc = xf.xiforc > x1 && xf.xiforc <= x2;

    // set transition interval flag
    xf.tran = xf.trforc || xf.trfree;

    // no transition yet: normal early return (trchek2 stays true)
    if xf.abort_on_nan && !xf.tran {
        return true;
    }

    // resolve if both forced and free transition
    if xf.trfree && xf.trforc {
        xf.trforc = xf.xiforc < xf.xt;
        xf.trfree = xf.xiforc >= xf.xt;
    }

    if xf.trforc {
        // if forced transition, then XT is prescribed; no sense calculating
        // the sensitivities, since we know them...
        xf.xt = xf.xiforc;
        xf.xt_a1 = 0.0;
        xf.xt_x1 = 0.0;
        xf.xt_t1 = 0.0;
        xf.xt_d1 = 0.0;
        xf.xt_u1 = 0.0;
        xf.xt_x2 = 0.0;
        xf.xt_t2 = 0.0;
        xf.xt_d2 = 0.0;
        xf.xt_u2 = 0.0;
        xf.xt_ms = 0.0;
        xf.xt_re = 0.0;
        xf.xt_xf = 1.0;
        return true;
    }

    // free transition ... set sensitivities of XT
    let wf1 = 1.0 - wf2;

    // XT( X1 X2 A1 A2 XF ),  TT( T1 T2 A1 A2 X1 X2 XF),  DT( ... )
    xf.xt_x1 = wf1;
    tt_t1 = wf1;
    dt_d1 = wf1;
    ut_u1 = wf1;

    xf.xt_x2 = wf2;
    tt_t2 = wf2;
    dt_d2 = wf2;
    ut_u2 = wf2;

    xf.xt_a1 = x1 * wf1_a1 + x2 * wf2_a1;
    tt_a1 = t1 * wf1_a1 + xf.com2.t * wf2_a1;
    dt_a1 = d1 * wf1_a1 + xf.com2.d * wf2_a1;
    ut_a1 = u1 * wf1_a1 + xf.com2.u * wf2_a1;

    xf.xt_x1 += x1 * wf1_x1 + x2 * wf2_x1;
    tt_x1 = t1 * wf1_x1 + xf.com2.t * wf2_x1;
    dt_x1 = d1 * wf1_x1 + xf.com2.d * wf2_x1;
    ut_x1 = u1 * wf1_x1 + xf.com2.u * wf2_x1;

    xf.xt_x2 += x1 * wf1_x2 + x2 * wf2_x2;
    tt_x2 = t1 * wf1_x2 + xf.com2.t * wf2_x2;
    dt_x2 = d1 * wf1_x2 + xf.com2.d * wf2_x2;
    ut_x2 = u1 * wf1_x2 + xf.com2.u * wf2_x2;

    xf.xt_xf = x1 * wf1_xf + x2 * wf2_xf;
    tt_xf = t1 * wf1_xf + xf.com2.t * wf2_xf;
    dt_xf = d1 * wf1_xf + xf.com2.d * wf2_xf;
    ut_xf = u1 * wf1_xf + xf.com2.u * wf2_xf;

    // at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT )
    // set sensitivities of AX( T1 D1 U1 A1 T2 D2 U2 A2 MS RE )
    let hk1 = c1.hk;
    let rt1 = c1.rt;
    let hk1_t1 = c1.hk_t;
    let hk1_d1 = c1.hk_d;
    let hk1_u1 = c1.hk_u;
    let hk1_ms = c1.hk_ms;
    let rt1_t1 = c1.rt_t;
    let rt1_u1 = c1.rt_u;
    let rt1_ms = c1.rt_ms;
    let rt1_re = c1.rt_re;

    let hk2 = xf.com2.hk;
    let hk2_t2 = xf.com2.hk_t;
    let hk2_d2 = xf.com2.hk_d;
    let hk2_u2 = xf.com2.hk_u;
    let hk2_ms = xf.com2.hk_ms;
    let rt2 = xf.com2.rt;
    let rt2_t2 = xf.com2.rt_t;
    let rt2_u2 = xf.com2.rt_u;
    let rt2_ms = xf.com2.rt_ms;
    let rt2_re = xf.com2.rt_re;

    let ax_t1 = ax_hk1 * hk1_t1 + ax_t1 + ax_rt1 * rt1_t1 + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t1;
    let ax_d1 = ax_hk1 * hk1_d1 + (ax_hkt * hkt_dt) * dt_d1;
    let ax_u1 = ax_hk1 * hk1_u1 + ax_rt1 * rt1_u1 + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u1;
    let ax_a1 = ax_a1 + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a1 + (ax_hkt * hkt_dt) * dt_a1
        + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a1;
    let ax_x1 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x1 + (ax_hkt * hkt_dt) * dt_x1
        + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x1;

    let ax_t2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_t2;
    let ax_d2 = (ax_hkt * hkt_dt) * dt_d2;
    let ax_u2 = (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_u2;
    let ax_a2 = ax_at * amplt_a2 + (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_a2 + (ax_hkt * hkt_dt) * dt_a2
        + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_a2;
    let ax_x2 = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_x2 + (ax_hkt * hkt_dt) * dt_x2
        + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_x2;

    let ax_xf = (ax_hkt * hkt_tt + ax_tt + ax_rtt * rtt_tt) * tt_xf + (ax_hkt * hkt_dt) * dt_xf
        + (ax_hkt * hkt_ut + ax_rtt * rtt_ut) * ut_xf;

    let ax_ms = ax_hkt * hkt_ms + ax_rtt * rtt_ms + ax_hk1 * hk1_ms + ax_rt1 * rt1_ms;
    let ax_re = ax_rtt * rtt_re + ax_rt1 * rt1_re;

    // set sensitivities of residual RES
    // RES = AMPL2 - AMPL1 - AX*(X2-X1)
    let z_ax = -(x2 - x1);

    let z_a1 = z_ax * ax_a1 - 1.0;
    let z_t1 = z_ax * ax_t1;
    let z_d1 = z_ax * ax_d1;
    let z_u1 = z_ax * ax_u1;
    let z_x1 = z_ax * ax_x1 + ax;

    let z_a2 = z_ax * ax_a2 + 1.0;
    let z_t2 = z_ax * ax_t2;
    let z_d2 = z_ax * ax_d2;
    let z_u2 = z_ax * ax_u2;
    let z_x2 = z_ax * ax_x2 - ax;

    let z_xf = z_ax * ax_xf;
    let z_ms = z_ax * ax_ms;
    let z_re = z_ax * ax_re;

    // set sensitivities of XT, with RES being stationary for A2 constraint
    xf.xt_a1 -= (xt_a2 / z_a2) * z_a1;
    xf.xt_t1 = -(xt_a2 / z_a2) * z_t1;
    xf.xt_d1 = -(xt_a2 / z_a2) * z_d1;
    xf.xt_u1 = -(xt_a2 / z_a2) * z_u1;
    xf.xt_x1 -= (xt_a2 / z_a2) * z_x1;
    xf.xt_t2 = -(xt_a2 / z_a2) * z_t2;
    xf.xt_d2 = -(xt_a2 / z_a2) * z_d2;
    xf.xt_u2 = -(xt_a2 / z_a2) * z_u2;
    xf.xt_x2 -= (xt_a2 / z_a2) * z_x2;
    xf.xt_ms = -(xt_a2 / z_a2) * z_ms;
    xf.xt_re = -(xt_a2 / z_a2) * z_re;
    xf.xt_xf = 0.0;

    let _ = (wf1, hk2, hk2_t2, hk2_d2, hk2_u2, hk2_ms, rt2, rt2_t2, rt2_u2, rt2_ms, rt2_re);

    true
}

/// Amplification rate routine for envelope e^n method (original version).
/// Returns (Ax, Ax_hk, Ax_th, Ax_rt).
pub fn dampl(hk: f64, th: f64, rt: f64) -> (f64, f64, f64, f64) {
    const DGR: f64 = 0.08;

    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // log10(Critical Rth) - H correlation for Falkner-Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    let gr = rt.log10();
    let gr_rt = 1.0 / (std::f64::consts::LN_10 * rt);

    if gr < grcrit - DGR {
        // no amplification for Rtheta < Rcrit
        return (0.0, 0.0, 0.0, 0.0);
    }

    // set steep cubic ramp used to turn on AX smoothly as Rtheta exceeds Rcrit
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);

    let (rfac, rfac_hk, rfac_rt);
    if rnorm >= 1.0 {
        rfac = 1.0;
        rfac_hk = 0.0;
        rfac_rt = 0.0;
    } else {
        rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        rfac_hk = rfac_rn * rn_hk;
        rfac_rt = rfac_rn * rn_rt;
    }

    // amplification envelope slope correlation for Falkner-Skan
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;

    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);

    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;

    // new m(H) correlation
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    let af_hk = af_hmi * hmi_hk;

    let ax = (af * dadr / th) * rfac;
    let ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
    let ax_th = -ax / th;
    let ax_rt = (af * dadr / th) * rfac_rt;

    (ax, ax_hk, ax_th, ax_rt)
}

/// Amplification rate routine for the modified envelope e^n method (2nd order).
pub fn dampl2(hk: f64, th: f64, rt: f64) -> (f64, f64, f64, f64) {
    const DGR: f64 = 0.08;
    const HK1: f64 = 3.5;
    const HK2: f64 = 4.0;

    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // log10(Critical Rth) -- H correlation for Falkner-Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    let grc = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    let gr = rt.log10();
    let gr_rt = 1.0 / (std::f64::consts::LN_10 * rt);

    if gr < grc - DGR {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let rnorm = (gr - (grc - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);

    let (rfac, rfac_hk, rfac_rt);
    if rnorm >= 1.0 {
        rfac = 1.0;
        rfac_hk = 0.0;
        rfac_rt = 0.0;
    } else {
        rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
        let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        rfac_hk = rfac_rn * rn_hk;
        rfac_rt = rfac_rn * rn_rt;
    }

    // envelope amplification rate with respect to Rtheta: DADR = d(N)/d(Rtheta) = f(H)
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;

    let ex = (-arg * arg).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);

    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;

    // conversion factor from d/d(Rtheta) to d/dx: AF = Theta d(Rtheta)/dx = f(H)
    let brg = -20.0 * hmi;
    let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi + 0.1 * brg.exp();
    let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi - 2.0 * brg.exp();
    let af_hk = af_hmi * hmi_hk;

    // amplification rate with respect to x, with RFAC shutting off amplification below Rcrit
    let mut ax = (af * dadr / th) * rfac;
    let mut ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
    let mut ax_th = -ax / th;
    let mut ax_rt = (af * dadr / th) * rfac_rt;

    if hk < HK1 {
        return (ax, ax_hk, ax_th, ax_rt);
    }

    // non-envelope max-amplification correction for separated profiles
    let hnorm = (hk - HK1) / (HK2 - HK1);
    let hn_hk = 1.0 / (HK2 - HK1);

    // blending fraction HFAC = 0..1 over HK1 < HK < HK2
    let (hfac, hf_hk);
    if hnorm >= 1.0 {
        hfac = 1.0;
        hf_hk = 0.0;
    } else {
        hfac = 3.0 * hnorm * hnorm - 2.0 * hnorm * hnorm * hnorm;
        hf_hk = (6.0 * hnorm - 6.0 * hnorm * hnorm) * hn_hk;
    }

    // "normal" envelope amplification rate AX1
    let ax1 = ax;
    let ax1_hk = ax_hk;
    let ax1_th = ax_th;
    let ax1_rt = ax_rt;

    // modified amplification rate AX2
    let gr0 = 0.30 + 0.35 * (-0.15 * (hk - 5.0)).exp();
    let gr0_hk = -0.35 * (-0.15 * (hk - 5.0)).exp() * 0.15;

    let tnr = (1.2 * (gr - gr0)).tanh();
    let tnr_rt = (1.0 - tnr * tnr) * 1.2 * gr_rt;
    let tnr_hk = -(1.0 - tnr * tnr) * 1.2 * gr0_hk;

    let mut ax2 = (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / th;
    let mut ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0).powf(2.5)) / th;
    let mut ax2_rt = (0.086 * tnr_rt) / th;
    let mut ax2_th = -ax2 / th;

    if ax2 < 0.0 {
        ax2 = 0.0;
        ax2_hk = 0.0;
        ax2_rt = 0.0;
        ax2_th = 0.0;
    }

    // blend the two amplification rates
    ax = hfac * ax2 + (1.0 - hfac) * ax1;
    ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1);
    ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt;
    ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th;

    (ax, ax_hk, ax_th, ax_rt)
}

/// Calculates kinematic shape parameter (from Whitfield).
/// Returns (Hk, Hk_h, Hk_msq).
pub fn hkin(h: f64, msq: f64) -> (f64, f64, f64) {
    let hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
    let hk_h = 1.0 / (1.0 + 0.113 * msq);
    let hk_msq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq);
    (hk, hk_h, hk_msq)
}

/// Laminar dissipation function (2 CD/H*) from Falkner-Skan.
/// Returns (Di, Di_hk, Di_rt).
pub fn dil(hk: f64, rt: f64) -> (f64, f64, f64) {
    let (di, di_hk);
    if hk < 4.0 {
        di = (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / rt;
        di_hk = (-0.00205 * 5.5 * (4.0 - hk).powf(4.5)) / rt;
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb * hkb;
        di = (-0.0016 * hkb * hkb / den + 0.207) / rt;
        di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / (den * den))) / rt;
    }
    let di_rt = -di / rt;
    (di, di_hk, di_rt)
}

/// Laminar wake dissipation function (2 CD/H*).
/// Returns (Di, Di_hk, Di_rt).
pub fn dilw(hk: f64, rt: f64) -> (f64, f64, f64) {
    let msq = 0.0;
    let (hs, hs_hk, hs_rt, _hs_msq) = hsl(hk, rt, msq);

    let rcd = 1.10 * (1.0 - 1.0 / hk).powi(2) / hk;
    let rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / (hk * hk * hk) - rcd / hk;

    let di = 2.0 * rcd / (hs * rt);
    let di_hk = 2.0 * rcd_hk / (hs * rt) - (di / hs) * hs_hk;
    let di_rt = -di / rt - (di / hs) * hs_rt;

    (di, di_hk, di_rt)
}

/// Laminar HS correlation.  Returns (Hs, Hs_hk, Hs_rt, Hs_msq).
pub fn hsl(hk: f64, _rt: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (hs, hs_hk);
    if hk < 4.35 {
        let tmp = hk - 4.35;
        hs = 0.0111 * tmp * tmp / (hk + 1.0) - 0.0278 * tmp * tmp * tmp / (hk + 1.0)
            + 1.528
            - 0.0002 * (tmp * hk).powi(2);
        hs_hk = 0.0111 * (2.0 * tmp - tmp * tmp / (hk + 1.0)) / (hk + 1.0)
            - 0.0278 * (3.0 * tmp * tmp - tmp * tmp * tmp / (hk + 1.0)) / (hk + 1.0)
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
    } else {
        let hs2 = 0.015;
        hs = hs2 * (hk - 4.35).powi(2) / hk + 1.528;
        hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * (hk - 4.35).powi(2) / (hk * hk);
    }
    (hs, hs_hk, 0.0, 0.0)
}

/// Laminar skin friction function (Cf) from Falkner-Skan.
/// Returns (Cf, Cf_hk, Cf_rt, Cf_msq).
pub fn cfl(hk: f64, rt: f64, _msq: f64) -> (f64, f64, f64, f64) {
    let (cf, cf_hk);
    if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        cf = (0.0727 * tmp - 0.07) / rt;
        cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        cf = (0.015 * tmp * tmp - 0.07) / rt;
        cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5).powi(2)) / rt;
    }
    let cf_rt = -cf / rt;
    (cf, cf_hk, cf_rt, 0.0)
}

/// Turbulent dissipation function (2 CD/H*).
/// Returns (Di, Di_hs, Di_us, Di_cf, Di_st).
#[allow(clippy::too_many_arguments)]
pub fn dit(hs: f64, us: f64, cf: f64, st: f64) -> (f64, f64, f64, f64, f64) {
    let di = (0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs;
    let di_hs = -(0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / (hs * hs);
    let di_us = (0.5 * cf - st * st) * 2.0 / hs;
    let di_cf = (0.5 * us) * 2.0 / hs;
    let di_st = (2.0 * st * (1.0 - us)) * 2.0 / hs;
    (di, di_hs, di_us, di_cf, di_st)
}

/// Turbulent HS correlation.  Returns (Hs, Hs_hk, Hs_rt, Hs_msq).
pub fn hst(hk: f64, rt: f64, msq: f64) -> (f64, f64, f64, f64) {
    const HSMIN: f64 = 1.500;
    const DHSINF: f64 = 0.015;

    let (ho, ho_rt);
    if rt > 400.0 {
        ho = 3.0 + 400.0 / rt;
        ho_rt = -400.0 / (rt * rt);
    } else {
        ho = 4.0;
        ho_rt = 0.0;
    }

    let (rtz, rtz_rt);
    if rt > 200.0 {
        rtz = rt;
        rtz_rt = 1.0;
    } else {
        rtz = 200.0;
        rtz_rt = 0.0;
    }

    let (hs, hs_hk, hs_rt);
    if hk < ho {
        // attached branch (new correlation, from arctan(y+) + Schlichting profiles)
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        hs = (2.0 - HSMIN - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + HSMIN + 4.0 / rtz;
        hs_hk = -(2.0 - HSMIN - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5).powi(2)
            + (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        hs_rt = (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
            + (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / (rtz * rtz) * rtz_rt;
    } else {
        // separated branch
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / (rtmp * rtmp) + DHSINF / hk;
        let htmp_hk = -0.014 * grt / (rtmp * rtmp * rtmp) - DHSINF / (hk * hk);
        let htmp_rt =
            -0.014 * grt / (rtmp * rtmp * rtmp) * (-ho_rt - 4.0 / (grt * grt) / rtz * rtz_rt) + 0.007 / (rtmp * rtmp) / rtz * rtz_rt;
        hs = hdif * hdif * htmp + HSMIN + 4.0 / rtz;
        hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
        hs_rt = hdif * hdif * htmp_rt - 4.0 / (rtz * rtz) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
    }

    // Whitfield's minor additional compressibility correction
    let fm = 1.0 + 0.014 * msq;
    let hs = (hs + 0.028 * msq) / fm;
    let hs_hk = hs_hk / fm;
    let hs_rt = hs_rt / fm;
    let hs_msq = 0.028 / fm - 0.014 * hs / fm;

    (hs, hs_hk, hs_rt, hs_msq)
}

/// Turbulent skin friction function (Cf) (Coles).
/// Returns (Cf, Cf_hk, Cf_rt, Cf_msq).
pub fn cft(hk: f64, rt: f64, msq: f64, cffac: f64) -> (f64, f64, f64, f64) {
    const GAM: f64 = 1.4;
    let gm1 = GAM - 1.0;
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    let mut grt = (rt / fc).ln();
    grt = grt.max(3.0);

    let gex = -1.74 - 0.31 * hk;

    let arg = (-1.33 * hk).max(-20.0);

    let thk = (4.0 - hk / 0.875).tanh();

    let cfo = cffac * 0.3 * arg.exp() * (grt / std::f64::consts::LN_10).powf(gex);
    let cf = (cfo + 1.1E-4 * (thk - 1.0)) / fc;
    let cf_hk = (-1.33 * cfo - 0.31 * (grt / std::f64::consts::LN_10).ln() * cfo - 1.1E-4 * (1.0 - thk * thk) / 0.875) / fc;
    let cf_rt = gex * cfo / (fc * grt) / rt;
    let cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / (fc * fc)) - 0.25 * gm1 * cf / (fc * fc);

    (cf, cf_hk, cf_rt, cf_msq)
}

/// Density shape parameter (from Whitfield).  Returns (Hc, Hc_hk, Hc_msq).
pub fn hct(hk: f64, msq: f64) -> (f64, f64, f64) {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / (hk - 0.8).powi(2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    (hc, hc_hk, hc_msq)
}

/// Sets the primary "2" BL variables from the parameter list (BLPRV).
#[allow(clippy::too_many_arguments)]
pub fn blprv(xf: &mut Xfoil, xsi: f64, ami: f64, cti: f64, thi: f64, dsi: f64, dswaki: f64, uei: f64) {
    xf.com2.x = xsi;
    xf.com2.ampl = ami;
    xf.com2.s = cti;
    xf.com2.t = thi;
    xf.com2.d = dsi - dswaki;
    xf.com2.dw = dswaki;

    xf.com2.u = uei * (1.0 - xf.tkbl) / (1.0 - xf.tkbl * (uei / xf.qinfbl).powi(2));
    xf.com2.u_uei = (1.0 + xf.tkbl * (2.0 * xf.com2.u * uei / xf.qinfbl.powi(2) - 1.0))
        / (1.0 - xf.tkbl * (uei / xf.qinfbl).powi(2));
    xf.com2.u_ms = (xf.com2.u * (uei / xf.qinfbl).powi(2) - uei) * xf.tkbl_ms
        / (1.0 - xf.tkbl * (uei / xf.qinfbl).powi(2));
}

/// Calculates turbulence-independent secondary "2" variables from the
/// primary "2" variables (BLKIN).
pub fn blkin(xf: &mut Xfoil) {
    // set edge Mach number ** 2
    let u2 = xf.com2.u;
    xf.com2.m = u2 * u2 * xf.hstinv / (xf.gm1bl * (1.0 - 0.5 * u2 * u2 * xf.hstinv));
    let tr2 = 1.0 + 0.5 * xf.gm1bl * xf.com2.m;
    xf.com2.m_u = 2.0 * xf.com2.m * tr2 / u2;
    xf.com2.m_ms = u2 * u2 * tr2 / (xf.gm1bl * (1.0 - 0.5 * u2 * u2 * xf.hstinv)) * xf.hstinv_ms;

    // set edge static density (isentropic relation)
    xf.com2.r = xf.rstbl * tr2.powf(-1.0 / xf.gm1bl);
    xf.com2.r_u = -xf.com2.r / tr2 * 0.5 * xf.com2.m_u;
    xf.com2.r_ms = -xf.com2.r / tr2 * 0.5 * xf.com2.m_ms + xf.rstbl_ms * tr2.powf(-1.0 / xf.gm1bl);

    // set shape parameter
    xf.com2.h = xf.com2.d / xf.com2.t;
    xf.com2.h_d = 1.0 / xf.com2.t;
    xf.com2.h_t = -xf.com2.h / xf.com2.t;

    // set edge static/stagnation enthalpy
    let herat = 1.0 - 0.5 * u2 * u2 * xf.hstinv;
    let he_u2 = -u2 * xf.hstinv;
    let he_ms = -0.5 * u2 * u2 * xf.hstinv_ms;

    // set molecular viscosity
    xf.com2.v = herat.powf(1.5) * (1.0 + xf.hvrat) / (herat + xf.hvrat) / xf.reybl;
    let v2_he = xf.com2.v * (1.5 / herat - 1.0 / (herat + xf.hvrat));

    xf.com2.v_u = v2_he * he_u2;
    xf.com2.v_ms = -xf.com2.v / xf.reybl * xf.reybl_ms + v2_he * he_ms;
    xf.com2.v_re = -xf.com2.v / xf.reybl * xf.reybl_re;

    // set kinematic shape parameter
    let (hk2, hk2_h2, hk2_m2) = hkin(xf.com2.h, xf.com2.m);

    xf.com2.hk = hk2;
    xf.com2.hk_u = hk2_m2 * xf.com2.m_u;
    xf.com2.hk_t = hk2_h2 * xf.com2.h_t;
    xf.com2.hk_d = hk2_h2 * xf.com2.h_d;
    xf.com2.hk_ms = hk2_m2 * xf.com2.m_ms;

    // set momentum thickness Reynolds number
    xf.com2.rt = xf.com2.r * u2 * xf.com2.t / xf.com2.v;
    xf.com2.rt_u = xf.com2.rt * (1.0 / u2 + xf.com2.r_u / xf.com2.r - xf.com2.v_u / xf.com2.v);
    xf.com2.rt_t = xf.com2.rt / xf.com2.t;
    xf.com2.rt_ms = xf.com2.rt * (xf.com2.r_ms / xf.com2.r - xf.com2.v_ms / xf.com2.v);
    xf.com2.rt_re = xf.com2.rt * (-xf.com2.v_re / xf.com2.v);
}

/// Calculates all secondary "2" variables from the primary "2" variables
/// (BLVAR).  `ityp` = 1 laminar, 2 turbulent, 3 turbulent wake.
pub fn blvar(xf: &mut Xfoil, ityp: i32) {
    if ityp == 3 {
        xf.com2.hk = xf.com2.hk.max(1.00005);
    } else {
        xf.com2.hk = xf.com2.hk.max(1.05000);
    }

    // density thickness shape parameter (H**)
    let (hc2, hc2_hk2, hc2_m2) = hct(xf.com2.hk, xf.com2.m);
    xf.com2.hc = hc2;
    xf.com2.hc_u = hc2_hk2 * xf.com2.hk_u + hc2_m2 * xf.com2.m_u;
    xf.com2.hc_t = hc2_hk2 * xf.com2.hk_t;
    xf.com2.hc_d = hc2_hk2 * xf.com2.hk_d;
    xf.com2.hc_ms = hc2_hk2 * xf.com2.hk_ms + hc2_m2 * xf.com2.m_ms;

    // set KE thickness shape parameter from H - H* correlations
    let (hs2, hs2_hk2, hs2_rt2, hs2_m2);
    if ityp == 1 {
        let r = hsl(xf.com2.hk, xf.com2.rt, xf.com2.m);
        hs2 = r.0;
        hs2_hk2 = r.1;
        hs2_rt2 = r.2;
        hs2_m2 = r.3;
    } else {
        let r = hst(xf.com2.hk, xf.com2.rt, xf.com2.m);
        hs2 = r.0;
        hs2_hk2 = r.1;
        hs2_rt2 = r.2;
        hs2_m2 = r.3;
    }

    xf.com2.hs = hs2;
    xf.com2.hs_u = hs2_hk2 * xf.com2.hk_u + hs2_rt2 * xf.com2.rt_u + hs2_m2 * xf.com2.m_u;
    xf.com2.hs_t = hs2_hk2 * xf.com2.hk_t + hs2_rt2 * xf.com2.rt_t;
    xf.com2.hs_d = hs2_hk2 * xf.com2.hk_d;
    xf.com2.hs_ms = hs2_hk2 * xf.com2.hk_ms + hs2_rt2 * xf.com2.rt_ms + hs2_m2 * xf.com2.m_ms;
    xf.com2.hs_re = hs2_rt2 * xf.com2.rt_re;

    // normalized slip velocity Us
    xf.com2.us = 0.5 * xf.com2.hs * (1.0 - (xf.com2.hk - 1.0) / (xf.gbcon * xf.com2.h));
    let us2_hs2 = 0.5 * (1.0 - (xf.com2.hk - 1.0) / (xf.gbcon * xf.com2.h));
    let us2_hk2 = 0.5 * xf.com2.hs * (-1.0 / (xf.gbcon * xf.com2.h));
    let us2_h2 = 0.5 * xf.com2.hs * (xf.com2.hk - 1.0) / (xf.gbcon * xf.com2.h.powi(2));

    xf.com2.us_u = us2_hs2 * xf.com2.hs_u + us2_hk2 * xf.com2.hk_u;
    xf.com2.us_t = us2_hs2 * xf.com2.hs_t + us2_hk2 * xf.com2.hk_t + us2_h2 * xf.com2.h_t;
    xf.com2.us_d = us2_hs2 * xf.com2.hs_d + us2_hk2 * xf.com2.hk_d + us2_h2 * xf.com2.h_d;
    xf.com2.us_ms = us2_hs2 * xf.com2.hs_ms + us2_hk2 * xf.com2.hk_ms;
    xf.com2.us_re = us2_hs2 * xf.com2.hs_re;

    if ityp <= 2 && xf.com2.us > 0.95 {
        xf.com2.us = 0.98;
        xf.com2.us_u = 0.0;
        xf.com2.us_t = 0.0;
        xf.com2.us_d = 0.0;
        xf.com2.us_ms = 0.0;
        xf.com2.us_re = 0.0;
    }

    if ityp == 3 && xf.com2.us > 0.99995 {
        xf.com2.us = 0.99995;
        xf.com2.us_u = 0.0;
        xf.com2.us_t = 0.0;
        xf.com2.us_d = 0.0;
        xf.com2.us_ms = 0.0;
        xf.com2.us_re = 0.0;
    }

    // equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
    let mut gcc = 0.0;
    let mut hkc = xf.com2.hk - 1.0;
    let mut hkc_hk2 = 1.0;
    let mut hkc_rt2 = 0.0;
    if ityp == 2 {
        gcc = xf.gccon;
        hkc = xf.com2.hk - 1.0 - gcc / xf.com2.rt;
        hkc_hk2 = 1.0;
        hkc_rt2 = gcc / xf.com2.rt.powi(2);
        if hkc < 0.01 {
            hkc = 0.01;
            hkc_hk2 = 0.0;
            hkc_rt2 = 0.0;
        }
    }

    let hkb = xf.com2.hk - 1.0;
    let usb = 1.0 - xf.com2.us;
    xf.com2.cq = (xf.ctcon * xf.com2.hs * hkb * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2))).sqrt();
    let cq2_hs2 = xf.ctcon * hkb * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) * 0.5 / xf.com2.cq;
    let cq2_us2 = xf.ctcon * xf.com2.hs * hkb * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) / usb * 0.5 / xf.com2.cq;
    let cq2_hk2 = xf.ctcon * xf.com2.hs * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) * 0.5 / xf.com2.cq
        - xf.ctcon * xf.com2.hs * hkb * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(3)) * 2.0 * 0.5 / xf.com2.cq
        + xf.ctcon * xf.com2.hs * hkb * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) * 2.0 * 0.5 / xf.com2.cq * hkc_hk2;
    let cq2_rt2 = xf.ctcon * xf.com2.hs * hkb * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) * 2.0 * 0.5 / xf.com2.cq * hkc_rt2;
    let cq2_h2 = -xf.ctcon * xf.com2.hs * hkb * hkc * hkc / (usb * xf.com2.h * xf.com2.hk.powi(2)) / xf.com2.h * 0.5 / xf.com2.cq;

    xf.com2.cq_u = cq2_hs2 * xf.com2.hs_u + cq2_us2 * xf.com2.us_u + cq2_hk2 * xf.com2.hk_u;
    xf.com2.cq_t = cq2_hs2 * xf.com2.hs_t + cq2_us2 * xf.com2.us_t + cq2_hk2 * xf.com2.hk_t;
    xf.com2.cq_d = cq2_hs2 * xf.com2.hs_d + cq2_us2 * xf.com2.us_d + cq2_hk2 * xf.com2.hk_d;
    xf.com2.cq_ms = cq2_hs2 * xf.com2.hs_ms + cq2_us2 * xf.com2.us_ms + cq2_hk2 * xf.com2.hk_ms;
    xf.com2.cq_re = cq2_hs2 * xf.com2.hs_re + cq2_us2 * xf.com2.us_re;

    xf.com2.cq_u += cq2_rt2 * xf.com2.rt_u;
    xf.com2.cq_t += cq2_h2 * xf.com2.h_t + cq2_rt2 * xf.com2.rt_t;
    xf.com2.cq_d += cq2_h2 * xf.com2.h_d;
    xf.com2.cq_ms += cq2_rt2 * xf.com2.rt_ms;
    xf.com2.cq_re += cq2_rt2 * xf.com2.rt_re;

    // set skin friction coefficient
    let (mut cf2, mut cf2_hk2, mut cf2_rt2, mut cf2_m2);
    if ityp == 3 {
        cf2 = 0.0;
        cf2_hk2 = 0.0;
        cf2_rt2 = 0.0;
        cf2_m2 = 0.0;
    } else if ityp == 1 {
        let r = cfl(xf.com2.hk, xf.com2.rt, xf.com2.m);
        cf2 = r.0;
        cf2_hk2 = r.1;
        cf2_rt2 = r.2;
        cf2_m2 = r.3;
    } else {
        let r = cft(xf.com2.hk, xf.com2.rt, xf.com2.m, xf.cffac);
        cf2 = r.0;
        cf2_hk2 = r.1;
        cf2_rt2 = r.2;
        cf2_m2 = r.3;
        let rl = cfl(xf.com2.hk, xf.com2.rt, xf.com2.m);
        let cf2l = rl.0;
        if cf2l > cf2 {
            cf2 = cf2l;
            cf2_hk2 = rl.1;
            cf2_rt2 = rl.2;
            cf2_m2 = rl.3;
        }
    }

    xf.com2.cf = cf2;
    xf.com2.cf_u = cf2_hk2 * xf.com2.hk_u + cf2_rt2 * xf.com2.rt_u + cf2_m2 * xf.com2.m_u;
    xf.com2.cf_t = cf2_hk2 * xf.com2.hk_t + cf2_rt2 * xf.com2.rt_t;
    xf.com2.cf_d = cf2_hk2 * xf.com2.hk_d;
    xf.com2.cf_ms = cf2_hk2 * xf.com2.hk_ms + cf2_rt2 * xf.com2.rt_ms + cf2_m2 * xf.com2.m_ms;
    xf.com2.cf_re = cf2_rt2 * xf.com2.rt_re;

    // dissipation function 2 CD / H*
    if ityp == 1 {
        // laminar
        let (di2, di2_hk2, di2_rt2) = dil(xf.com2.hk, xf.com2.rt);
        xf.com2.di = di2;
        xf.com2.di_u = di2_hk2 * xf.com2.hk_u + di2_rt2 * xf.com2.rt_u;
        xf.com2.di_t = di2_hk2 * xf.com2.hk_t + di2_rt2 * xf.com2.rt_t;
        xf.com2.di_d = di2_hk2 * xf.com2.hk_d;
        xf.com2.di_s = 0.0;
        xf.com2.di_ms = di2_hk2 * xf.com2.hk_ms + di2_rt2 * xf.com2.rt_ms;
        xf.com2.di_re = di2_rt2 * xf.com2.rt_re;
    } else if ityp == 2 {
        // turbulent wall contribution
        let r = cft(xf.com2.hk, xf.com2.rt, xf.com2.m, xf.cffac);
        let cf2t = r.0;
        let cf2t_hk2 = r.1;
        let cf2t_rt2 = r.2;
        let cf2t_m2 = r.3;
        let cf2t_u2 = cf2t_hk2 * xf.com2.hk_u + cf2t_rt2 * xf.com2.rt_u + cf2t_m2 * xf.com2.m_u;
        let cf2t_t2 = cf2t_hk2 * xf.com2.hk_t + cf2t_rt2 * xf.com2.rt_t;
        let cf2t_d2 = cf2t_hk2 * xf.com2.hk_d;
        let cf2t_ms = cf2t_hk2 * xf.com2.hk_ms + cf2t_rt2 * xf.com2.rt_ms + cf2t_m2 * xf.com2.m_ms;
        let cf2t_re = cf2t_rt2 * xf.com2.rt_re;

        xf.com2.di = (0.5 * cf2t * xf.com2.us) * 2.0 / xf.com2.hs;
        let di2_hs2 = -(0.5 * cf2t * xf.com2.us) * 2.0 / xf.com2.hs.powi(2);
        let di2_us2 = (0.5 * cf2t) * 2.0 / xf.com2.hs;
        let di2_cf2t = (0.5 * xf.com2.us) * 2.0 / xf.com2.hs;

        xf.com2.di_s = 0.0;
        xf.com2.di_u = di2_hs2 * xf.com2.hs_u + di2_us2 * xf.com2.us_u + di2_cf2t * cf2t_u2;
        xf.com2.di_t = di2_hs2 * xf.com2.hs_t + di2_us2 * xf.com2.us_t + di2_cf2t * cf2t_t2;
        xf.com2.di_d = di2_hs2 * xf.com2.hs_d + di2_us2 * xf.com2.us_d + di2_cf2t * cf2t_d2;
        xf.com2.di_ms = di2_hs2 * xf.com2.hs_ms + di2_us2 * xf.com2.us_ms + di2_cf2t * cf2t_ms;
        xf.com2.di_re = di2_hs2 * xf.com2.hs_re + di2_us2 * xf.com2.us_re + di2_cf2t * cf2t_re;

        // set minimum Hk for wake layer to still exist
        let grt = xf.com2.rt.ln();
        let hmin = 1.0 + 2.1 / grt;
        let hm_rt2 = -(2.1 / grt.powi(2)) / xf.com2.rt;

        // set factor DFAC for correcting wall dissipation for very low Hk
        let fl = (xf.com2.hk - 1.0) / (hmin - 1.0);
        let fl_hk2 = 1.0 / (hmin - 1.0);
        let fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

        let tfl = fl.tanh();
        let dfac = 0.5 + 0.5 * tfl;
        let df_fl = 0.5 * (1.0 - tfl * tfl);

        let df_hk2 = df_fl * fl_hk2;
        let df_rt2 = df_fl * fl_rt2;

        xf.com2.di_s *= dfac;
        xf.com2.di_u = xf.com2.di_u * dfac + xf.com2.di * (df_hk2 * xf.com2.hk_u + df_rt2 * xf.com2.rt_u);
        xf.com2.di_t = xf.com2.di_t * dfac + xf.com2.di * (df_hk2 * xf.com2.hk_t + df_rt2 * xf.com2.rt_t);
        xf.com2.di_d = xf.com2.di_d * dfac + xf.com2.di * (df_hk2 * xf.com2.hk_d);
        xf.com2.di_ms = xf.com2.di_ms * dfac + xf.com2.di * (df_hk2 * xf.com2.hk_ms + df_rt2 * xf.com2.rt_ms);
        xf.com2.di_re = xf.com2.di_re * dfac + xf.com2.di * (df_rt2 * xf.com2.rt_re);
        xf.com2.di *= dfac;
    } else {
        // zero wall contribution for wake
        xf.com2.di = 0.0;
        xf.com2.di_s = 0.0;
        xf.com2.di_u = 0.0;
        xf.com2.di_t = 0.0;
        xf.com2.di_d = 0.0;
        xf.com2.di_ms = 0.0;
        xf.com2.di_re = 0.0;
    }

    // add on turbulent outer layer contribution
    if ityp != 1 {
        let mut dd = xf.com2.s.powi(2) * (0.995 - xf.com2.us) * 2.0 / xf.com2.hs;
        let mut dd_hs2 = -xf.com2.s.powi(2) * (0.995 - xf.com2.us) * 2.0 / xf.com2.hs.powi(2);
        let mut dd_us2 = -xf.com2.s.powi(2) * 2.0 / xf.com2.hs;
        let mut dd_s2 = xf.com2.s * 2.0 * (0.995 - xf.com2.us) * 2.0 / xf.com2.hs;

        xf.com2.di += dd;
        xf.com2.di_s = dd_s2;
        xf.com2.di_u += dd_hs2 * xf.com2.hs_u + dd_us2 * xf.com2.us_u;
        xf.com2.di_t += dd_hs2 * xf.com2.hs_t + dd_us2 * xf.com2.us_t;
        xf.com2.di_d += dd_hs2 * xf.com2.hs_d + dd_us2 * xf.com2.us_d;
        xf.com2.di_ms += dd_hs2 * xf.com2.hs_ms + dd_us2 * xf.com2.us_ms;
        xf.com2.di_re += dd_hs2 * xf.com2.hs_re + dd_us2 * xf.com2.us_re;

        // add laminar stress contribution to outer layer CD
        dd = 0.15 * (0.995 - xf.com2.us).powi(2) / xf.com2.rt * 2.0 / xf.com2.hs;
        dd_us2 = -0.15 * (0.995 - xf.com2.us) * 2.0 / xf.com2.rt * 2.0 / xf.com2.hs;
        dd_hs2 = -dd / xf.com2.hs;
        let dd_rt2 = -dd / xf.com2.rt;

        xf.com2.di += dd;
        xf.com2.di_u += dd_hs2 * xf.com2.hs_u + dd_us2 * xf.com2.us_u + dd_rt2 * xf.com2.rt_u;
        xf.com2.di_t += dd_hs2 * xf.com2.hs_t + dd_us2 * xf.com2.us_t + dd_rt2 * xf.com2.rt_t;
        xf.com2.di_d += dd_hs2 * xf.com2.hs_d + dd_us2 * xf.com2.us_d;
        xf.com2.di_ms += dd_hs2 * xf.com2.hs_ms + dd_us2 * xf.com2.us_ms + dd_rt2 * xf.com2.rt_ms;
        xf.com2.di_re += dd_hs2 * xf.com2.hs_re + dd_us2 * xf.com2.us_re + dd_rt2 * xf.com2.rt_re;
    }

    if ityp == 2 {
        let (di2l, di2l_hk2, di2l_rt2) = dil(xf.com2.hk, xf.com2.rt);

        if di2l > xf.com2.di {
            xf.com2.di = di2l;
            xf.com2.di_s = 0.0;
            xf.com2.di_u = di2l_hk2 * xf.com2.hk_u + di2l_rt2 * xf.com2.rt_u;
            xf.com2.di_t = di2l_hk2 * xf.com2.hk_t + di2l_rt2 * xf.com2.rt_t;
            xf.com2.di_d = di2l_hk2 * xf.com2.hk_d;
            xf.com2.di_ms = di2l_hk2 * xf.com2.hk_ms + di2l_rt2 * xf.com2.rt_ms;
            xf.com2.di_re = di2l_rt2 * xf.com2.rt_re;
        }
    }

    if ityp == 3 {
        // laminar wake CD
        let (di2l, di2l_hk2, di2l_rt2) = dilw(xf.com2.hk, xf.com2.rt);
        if di2l > xf.com2.di {
            xf.com2.di = di2l;
            xf.com2.di_s = 0.0;
            xf.com2.di_u = di2l_hk2 * xf.com2.hk_u + di2l_rt2 * xf.com2.rt_u;
            xf.com2.di_t = di2l_hk2 * xf.com2.hk_t + di2l_rt2 * xf.com2.rt_t;
            xf.com2.di_d = di2l_hk2 * xf.com2.hk_d;
            xf.com2.di_ms = di2l_hk2 * xf.com2.hk_ms + di2l_rt2 * xf.com2.rt_ms;
            xf.com2.di_re = di2l_rt2 * xf.com2.rt_re;
        }

        // double dissipation for the wake (two wake halves)
        xf.com2.di *= 2.0;
        xf.com2.di_s *= 2.0;
        xf.com2.di_u *= 2.0;
        xf.com2.di_t *= 2.0;
        xf.com2.di_d *= 2.0;
        xf.com2.di_ms *= 2.0;
        xf.com2.di_re *= 2.0;
    }

    // BL thickness (Delta) from simplified Green's correlation
    xf.com2.de = (3.15 + 1.72 / (xf.com2.hk - 1.0)) * xf.com2.t + xf.com2.d;
    let de2_hk2 = (-1.72 / (xf.com2.hk - 1.0).powi(2)) * xf.com2.t;

    xf.com2.de_u = de2_hk2 * xf.com2.hk_u;
    xf.com2.de_t = de2_hk2 * xf.com2.hk_t + (3.15 + 1.72 / (xf.com2.hk - 1.0));
    xf.com2.de_d = de2_hk2 * xf.com2.hk_d + 1.0;
    xf.com2.de_ms = de2_hk2 * xf.com2.hk_ms;

    let hdmax = 12.0;
    if xf.com2.de > hdmax * xf.com2.t {
        xf.com2.de = hdmax * xf.com2.t;
        xf.com2.de_u = 0.0;
        xf.com2.de_t = hdmax;
        xf.com2.de_d = 0.0;
        xf.com2.de_ms = 0.0;
    }
}

/// Calculates midpoint skin friction CFM (BLMID).  `ityp` = 1 laminar,
/// 2 turbulent, 3 turbulent wake.
pub fn blmid(xf: &mut Xfoil, ityp: i32) {
    // set similarity variables if not defined
    if xf.simi {
        xf.com1.hk = xf.com2.hk;
        xf.com1.hk_t = xf.com2.hk_t;
        xf.com1.hk_d = xf.com2.hk_d;
        xf.com1.hk_u = xf.com2.hk_u;
        xf.com1.hk_ms = xf.com2.hk_ms;
        xf.com1.rt = xf.com2.rt;
        xf.com1.rt_t = xf.com2.rt_t;
        xf.com1.rt_u = xf.com2.rt_u;
        xf.com1.rt_ms = xf.com2.rt_ms;
        xf.com1.rt_re = xf.com2.rt_re;
        xf.com1.m = xf.com2.m;
        xf.com1.m_u = xf.com2.m_u;
        xf.com1.m_ms = xf.com2.m_ms;
    }

    // define stuff for midpoint CF
    let hka = 0.5 * (xf.com1.hk + xf.com2.hk);
    let rta = 0.5 * (xf.com1.rt + xf.com2.rt);
    let ma = 0.5 * (xf.com1.m + xf.com2.m);

    // midpoint skin friction coefficient (zero in wake)
    let mut cfm;
    let mut cfm_hka;
    let mut cfm_rta;
    let mut cfm_ma;
    if ityp == 3 {
        cfm = 0.0;
        cfm_hka = 0.0;
        cfm_rta = 0.0;
        cfm_ma = 0.0;
    } else if ityp == 1 {
        let r = cfl(hka, rta, ma);
        cfm = r.0;
        cfm_hka = r.1;
        cfm_rta = r.2;
        cfm_ma = r.3;
    } else {
        let r = cft(hka, rta, ma, xf.cffac);
        cfm = r.0;
        cfm_hka = r.1;
        cfm_rta = r.2;
        cfm_ma = r.3;
        let rl = cfl(hka, rta, ma);
        let cfml = rl.0;
        if cfml > cfm {
            cfm = cfml;
            cfm_hka = rl.1;
            cfm_rta = rl.2;
            cfm_ma = rl.3;
        }
    }
    xf.cfm = cfm;
    xf.cfm_u1 = 0.5 * (cfm_hka * xf.com1.hk_u + cfm_ma * xf.com1.m_u + cfm_rta * xf.com1.rt_u);
    xf.cfm_t1 = 0.5 * (cfm_hka * xf.com1.hk_t + cfm_rta * xf.com1.rt_t);
    xf.cfm_d1 = 0.5 * (cfm_hka * xf.com1.hk_d);
    xf.cfm_u2 = 0.5 * (cfm_hka * xf.com2.hk_u + cfm_ma * xf.com2.m_u + cfm_rta * xf.com2.rt_u);
    xf.cfm_t2 = 0.5 * (cfm_hka * xf.com2.hk_t + cfm_rta * xf.com2.rt_t);
    xf.cfm_d2 = 0.5 * (cfm_hka * xf.com2.hk_d);
    xf.cfm_ms = 0.5 * (cfm_hka * xf.com1.hk_ms + cfm_ma * xf.com1.m_ms + cfm_rta * xf.com1.rt_ms
        + cfm_hka * xf.com2.hk_ms + cfm_ma * xf.com2.m_ms + cfm_rta * xf.com2.rt_ms);
    xf.cfm_re = 0.5 * (cfm_rta * xf.com1.rt_re + cfm_rta * xf.com2.rt_re);
}

/// Sets up the "dummy" BL system between the airfoil TE point and the first
/// wake point infinitesimally behind the TE (TESYS).
pub fn tesys(xf: &mut Xfoil, cte: f64, tte: f64, dte: f64) {
    for k in 0..4 {
        xf.vsrez[k] = 0.0;
        xf.vsm[k] = 0.0;
        xf.vsr[k] = 0.0;
        xf.vsx[k] = 0.0;
        for l in 0..5 {
            xf.vs1[k][l] = 0.0;
            xf.vs2[k][l] = 0.0;
        }
    }

    blvar(xf, 3);

    xf.vs1[0][0] = -1.0;
    xf.vs2[0][0] = 1.0;
    xf.vsrez[0] = cte - xf.com2.s;

    xf.vs1[1][1] = -1.0;
    xf.vs2[1][1] = 1.0;
    xf.vsrez[1] = tte - xf.com2.t;

    xf.vs1[2][2] = -1.0;
    xf.vs2[2][2] = 1.0;
    xf.vsrez[2] = dte - xf.com2.d - xf.com2.dw;
}

/// Sets up the BL Newton system governing the current interval (BLSYS).
pub fn blsys(xf: &mut Xfoil) {
    // calculate secondary BL variables and their sensitivities
    if xf.wake {
        blvar(xf, 3);
        blmid(xf, 3);
    } else if xf.turb || xf.tran {
        blvar(xf, 2);
        blmid(xf, 2);
    } else {
        blvar(xf, 1);
        blmid(xf, 1);
    }

    // for the similarity station, "1" and "2" variables are the same
    if xf.simi {
        xf.com1 = xf.com2;
    }

    // set up appropriate finite difference system for current interval
    if xf.tran {
        trdif(xf);
    } else if xf.simi {
        bldif(xf, 0);
    } else if !xf.turb {
        bldif(xf, 1);
    } else if xf.wake {
        bldif(xf, 3);
    } else if xf.turb {
        bldif(xf, 2);
    }

    if xf.simi {
        // at similarity station, "1" variables are really "2" variables
        for k in 0..4 {
            for l in 0..5 {
                xf.vs2[k][l] += xf.vs1[k][l];
                xf.vs1[k][l] = 0.0;
            }
        }
    }

    // change system over into incompressible Uei and Mach
    for k in 0..4 {
        // residual derivatives wrt compressible Uec
        let res_u1 = xf.vs1[k][3];
        let res_u2 = xf.vs2[k][3];
        let res_ms = xf.vsm[k];

        // combine with derivatives of compressible U1,U2 = Uec(Uei M)
        xf.vs1[k][3] = res_u1 * xf.com1.u_uei;
        xf.vs2[k][3] = res_u2 * xf.com2.u_uei;
        xf.vsm[k] = res_u1 * xf.com1.u_ms + res_u2 * xf.com2.u_ms + res_ms;
    }
}

/// Sets up the Newton system coefficients and residuals (BLDIF).
///
/// `ityp` = 0 similarity station, 1 laminar interval, 2 turbulent interval,
/// 3 wake interval.
pub fn bldif(xf: &mut Xfoil, ityp: i32) {
    let (xlog, ulog, tlog, hlog, ddlog);
    if ityp == 0 {
        // similarity logarithmic differences (prescribed)
        xlog = 1.0;
        ulog = xf.bule;
        tlog = 0.5 * (1.0 - xf.bule);
        hlog = 0.0;
        ddlog = 0.0;
    } else {
        // usual logarithmic differences
        xlog = (xf.com2.x / xf.com1.x).ln();
        ulog = (xf.com2.u / xf.com1.u).ln();
        tlog = (xf.com2.t / xf.com1.t).ln();
        hlog = (xf.com2.hs / xf.com1.hs).ln();
        ddlog = 1.0;
    }

    for k in 0..4 {
        xf.vsrez[k] = 0.0;
        xf.vsm[k] = 0.0;
        xf.vsr[k] = 0.0;
        xf.vsx[k] = 0.0;
        for l in 0..5 {
            xf.vs1[k][l] = 0.0;
            xf.vs2[k][l] = 0.0;
        }
    }

    // set triggering constant for local upwinding
    let hupwt = 1.0;

    let mut hdcon = 5.0 * hupwt / xf.com2.hk.powi(2);
    let mut hd_hk1 = 0.0;
    let mut hd_hk2 = -hdcon * 2.0 / xf.com2.hk;

    // use less upwinding in the wake
    if ityp == 3 {
        hdcon = hupwt / xf.com2.hk.powi(2);
        hd_hk1 = 0.0;
        hd_hk2 = -hdcon * 2.0 / xf.com2.hk;
    }

    // local upwinding is based on local change in log(Hk-1)
    let arg = ((xf.com2.hk - 1.0) / (xf.com1.hk - 1.0)).abs();
    let hl = arg.ln();
    let hl_hk1 = -1.0 / (xf.com1.hk - 1.0);
    let hl_hk2 = 1.0 / (xf.com2.hk - 1.0);

    let hlsq = (hl * hl).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;
    let upw_hl = ehh * hl * hdcon;
    let upw_hd = 0.5 * ehh * hlsq;

    let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    let upw_u1 = upw_hk1 * xf.com1.hk_u;
    let upw_t1 = upw_hk1 * xf.com1.hk_t;
    let upw_d1 = upw_hk1 * xf.com1.hk_d;
    let upw_u2 = upw_hk2 * xf.com2.hk_u;
    let upw_t2 = upw_hk2 * xf.com2.hk_t;
    let upw_d2 = upw_hk2 * xf.com2.hk_d;
    let upw_ms = upw_hk1 * xf.com1.hk_ms + upw_hk2 * xf.com2.hk_ms;

    let (mut rezc, mut z_cfa, mut z_hka, mut z_da, mut z_sl, mut z_ul, mut z_dxi, mut z_usa, mut z_cqa, mut z_sa, mut z_dea, mut z_upw,
        mut z_de1, mut z_de2, mut z_us1, mut z_us2, mut z_d1, mut z_d2, mut z_u1, mut z_u2, mut z_x1, mut z_x2, mut z_s1, mut z_s2,
        mut z_cq1, mut z_cq2, mut z_cf1, mut z_cf2, mut z_hk1, mut z_hk2);
    z_cfa = 0.0;
    z_hka = 0.0;
    z_da = 0.0;
    z_sl = 0.0;
    z_ul = 0.0;
    z_dxi = 0.0;
    z_usa = 0.0;
    z_cqa = 0.0;
    z_sa = 0.0;
    z_dea = 0.0;
    z_upw = 0.0;
    z_de1 = 0.0;
    z_de2 = 0.0;
    z_us1 = 0.0;
    z_us2 = 0.0;
    z_d1 = 0.0;
    z_d2 = 0.0;
    z_u1 = 0.0;
    z_u2 = 0.0;
    z_x1 = 0.0;
    z_x2 = 0.0;
    z_s1 = 0.0;
    z_s2 = 0.0;
    z_cq1 = 0.0;
    z_cq2 = 0.0;
    z_cf1 = 0.0;
    z_cf2 = 0.0;
    z_hk1 = 0.0;
    z_hk2 = 0.0;
    rezc = 0.0;

    if ityp == 0 {
        // LE point --> set zero amplification factor
        xf.vs2[0][0] = 1.0;
        xf.vsr[0] = 0.0;
        xf.vsrez[0] = -xf.com2.ampl;
    } else if ityp == 1 {
        // laminar part --> set amplification equation
        let (ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2) = axset(
            xf.com1.hk,
            xf.com1.t,
            xf.com1.rt,
            xf.com1.ampl,
            xf.com2.hk,
            xf.com2.t,
            xf.com2.rt,
            xf.com2.ampl,
            xf.amcrit,
            xf.idampv,
        );

        rezc = xf.com2.ampl - xf.com1.ampl - ax * (xf.com2.x - xf.com1.x);
        let z_ax = -(xf.com2.x - xf.com1.x);

        xf.vs1[0][0] = z_ax * ax_a1 - 1.0;
        xf.vs1[0][1] = z_ax * (ax_hk1 * xf.com1.hk_t + ax_t1 + ax_rt1 * xf.com1.rt_t);
        xf.vs1[0][2] = z_ax * (ax_hk1 * xf.com1.hk_d);
        xf.vs1[0][3] = z_ax * (ax_hk1 * xf.com1.hk_u + ax_rt1 * xf.com1.rt_u);
        xf.vs1[0][4] = ax;
        xf.vs2[0][0] = z_ax * ax_a2 + 1.0;
        xf.vs2[0][1] = z_ax * (ax_hk2 * xf.com2.hk_t + ax_t2 + ax_rt2 * xf.com2.rt_t);
        xf.vs2[0][2] = z_ax * (ax_hk2 * xf.com2.hk_d);
        xf.vs2[0][3] = z_ax * (ax_hk2 * xf.com2.hk_u + ax_rt2 * xf.com2.rt_u);
        xf.vs2[0][4] = -ax;
        xf.vsm[0] = z_ax * (ax_hk1 * xf.com1.hk_ms + ax_rt1 * xf.com1.rt_ms + ax_hk2 * xf.com2.hk_ms + ax_rt2 * xf.com2.rt_ms);
        xf.vsr[0] = z_ax * (ax_rt1 * xf.com1.rt_re + ax_rt2 * xf.com2.rt_re);
        xf.vsx[0] = 0.0;
        xf.vsrez[0] = -rezc;
    } else {
        // turbulent part --> set shear lag equation
        let sa = (1.0 - upw) * xf.com1.s + upw * xf.com2.s;
        let cqa = (1.0 - upw) * xf.com1.cq + upw * xf.com2.cq;
        let cfa = (1.0 - upw) * xf.com1.cf + upw * xf.com2.cf;
        let hka = (1.0 - upw) * xf.com1.hk + upw * xf.com2.hk;

        let usa = 0.5 * (xf.com1.us + xf.com2.us);
        let rta = 0.5 * (xf.com1.rt + xf.com2.rt);
        let dea = 0.5 * (xf.com1.de + xf.com2.de);
        let da = 0.5 * (xf.com1.d + xf.com2.d);

        // increased dissipation length in wake (decrease its reciprocal)
        let ald = if ityp == 3 { xf.dlcon } else { 1.0 };

        // set and linearize equilibrium 1/Ue dUe/dx
        let gcc;
        let mut hkc;
        let mut hkc_hka;
        let mut hkc_rta;
        if ityp == 2 {
            gcc = xf.gccon;
            hkc = hka - 1.0 - gcc / rta;
            hkc_hka = 1.0;
            hkc_rta = gcc / (rta * rta);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hka = 0.0;
                hkc_rta = 0.0;
            }
        } else {
            gcc = 0.0;
            hkc = hka - 1.0;
            hkc_hka = 1.0;
            hkc_rta = 0.0;
        }
        let _ = gcc;

        let hr = hkc / (xf.gacon * ald * hka);
        let hr_hka = hkc_hka / (xf.gacon * ald * hka) - hr / hka;
        let hr_rta = hkc_rta / (xf.gacon * ald * hka);

        let uq = (0.5 * cfa - hr * hr) / (xf.gbcon * da);
        let uq_hka = -2.0 * hr * hr_hka / (xf.gbcon * da);
        let uq_rta = -2.0 * hr * hr_rta / (xf.gbcon * da);
        let uq_cfa = 0.5 / (xf.gbcon * da);
        let uq_da = -uq / da;
        let uq_upw = uq_cfa * (xf.com2.cf - xf.com1.cf) + uq_hka * (xf.com2.hk - xf.com1.hk);

        let mut uq_t1 = (1.0 - upw) * (uq_cfa * xf.com1.cf_t + uq_hka * xf.com1.hk_t) + uq_upw * upw_t1;
        let mut uq_d1 = (1.0 - upw) * (uq_cfa * xf.com1.cf_d + uq_hka * xf.com1.hk_d) + uq_upw * upw_d1;
        let mut uq_u1 = (1.0 - upw) * (uq_cfa * xf.com1.cf_u + uq_hka * xf.com1.hk_u) + uq_upw * upw_u1;
        let mut uq_t2 = upw * (uq_cfa * xf.com2.cf_t + uq_hka * xf.com2.hk_t) + uq_upw * upw_t2;
        let mut uq_d2 = upw * (uq_cfa * xf.com2.cf_d + uq_hka * xf.com2.hk_d) + uq_upw * upw_d2;
        let mut uq_u2 = upw * (uq_cfa * xf.com2.cf_u + uq_hka * xf.com2.hk_u) + uq_upw * upw_u2;
        let mut uq_ms = (1.0 - upw) * (uq_cfa * xf.com1.cf_ms + uq_hka * xf.com1.hk_ms)
            + uq_upw * upw_ms
            + upw * (uq_cfa * xf.com2.cf_ms + uq_hka * xf.com2.hk_ms);
        let mut uq_re = (1.0 - upw) * uq_cfa * xf.com1.cf_re + upw * uq_cfa * xf.com2.cf_re;

        uq_t1 += 0.5 * uq_rta * xf.com1.rt_t;
        uq_d1 += 0.5 * uq_da;
        uq_u1 += 0.5 * uq_rta * xf.com1.rt_u;
        uq_t2 += 0.5 * uq_rta * xf.com2.rt_t;
        uq_d2 += 0.5 * uq_da;
        uq_u2 += 0.5 * uq_rta * xf.com2.rt_u;
        uq_ms += 0.5 * uq_rta * xf.com1.rt_ms + 0.5 * uq_rta * xf.com2.rt_ms;
        uq_re += 0.5 * uq_rta * xf.com1.rt_re + 0.5 * uq_rta * xf.com2.rt_re;

        let scc = xf.sccon * 1.333 / (1.0 + usa);
        let scc_usa = -scc / (1.0 + usa);

        let scc_us1 = scc_usa * 0.5;
        let scc_us2 = scc_usa * 0.5;

        let slog = (xf.com2.s / xf.com1.s).ln();
        let dxi = xf.com2.x - xf.com1.x;

        rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * xf.duxcon;

        z_cfa = dea * 2.0 * uq_cfa * dxi * xf.duxcon;
        z_hka = dea * 2.0 * uq_hka * dxi * xf.duxcon;
        z_da = dea * 2.0 * uq_da * dxi * xf.duxcon;
        z_sl = -dea * 2.0;
        z_ul = -dea * 2.0 * xf.duxcon;
        z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * xf.duxcon;
        z_usa = scc_usa * (cqa - sa * ald) * dxi;
        z_cqa = scc * dxi;
        z_sa = -scc * dxi * ald;
        z_dea = 2.0 * ((uq * dxi - ulog) * xf.duxcon - slog);

        z_upw = z_cqa * (xf.com2.cq - xf.com1.cq) + z_sa * (xf.com2.s - xf.com1.s) + z_cfa * (xf.com2.cf - xf.com1.cf)
            + z_hka * (xf.com2.hk - xf.com1.hk);
        z_de1 = 0.5 * z_dea;
        z_de2 = 0.5 * z_dea;
        z_us1 = 0.5 * z_usa;
        z_us2 = 0.5 * z_usa;
        z_d1 = 0.5 * z_da;
        z_d2 = 0.5 * z_da;
        z_u1 = -z_ul / xf.com1.u;
        z_u2 = z_ul / xf.com2.u;
        z_x1 = -z_dxi;
        z_x2 = z_dxi;
        z_s1 = (1.0 - upw) * z_sa - z_sl / xf.com1.s;
        z_s2 = upw * z_sa + z_sl / xf.com2.s;
        z_cq1 = (1.0 - upw) * z_cqa;
        z_cq2 = upw * z_cqa;
        z_cf1 = (1.0 - upw) * z_cfa;
        z_cf2 = upw * z_cfa;
        z_hk1 = (1.0 - upw) * z_hka;
        z_hk2 = upw * z_hka;

        xf.vs1[0][0] = z_s1;
        xf.vs1[0][1] = z_upw * upw_t1 + z_de1 * xf.com1.de_t + z_us1 * xf.com1.us_t;
        xf.vs1[0][2] = z_d1 + z_upw * upw_d1 + z_de1 * xf.com1.de_d + z_us1 * xf.com1.us_d;
        xf.vs1[0][3] = z_u1 + z_upw * upw_u1 + z_de1 * xf.com1.de_u + z_us1 * xf.com1.us_u;
        xf.vs1[0][4] = z_x1;
        xf.vs2[0][0] = z_s2;
        xf.vs2[0][1] = z_upw * upw_t2 + z_de2 * xf.com2.de_t + z_us2 * xf.com2.us_t;
        xf.vs2[0][2] = z_d2 + z_upw * upw_d2 + z_de2 * xf.com2.de_d + z_us2 * xf.com2.us_d;
        xf.vs2[0][3] = z_u2 + z_upw * upw_u2 + z_de2 * xf.com2.de_u + z_us2 * xf.com2.us_u;
        xf.vs2[0][4] = z_x2;
        xf.vsm[0] = z_upw * upw_ms + z_de1 * xf.com1.de_ms + z_us1 * xf.com1.us_ms + z_de2 * xf.com2.de_ms + z_us2 * xf.com2.us_ms;

        xf.vs1[0][1] += z_cq1 * xf.com1.cq_t + z_cf1 * xf.com1.cf_t + z_hk1 * xf.com1.hk_t;
        xf.vs1[0][2] += z_cq1 * xf.com1.cq_d + z_cf1 * xf.com1.cf_d + z_hk1 * xf.com1.hk_d;
        xf.vs1[0][3] += z_cq1 * xf.com1.cq_u + z_cf1 * xf.com1.cf_u + z_hk1 * xf.com1.hk_u;

        xf.vs2[0][1] += z_cq2 * xf.com2.cq_t + z_cf2 * xf.com2.cf_t + z_hk2 * xf.com2.hk_t;
        xf.vs2[0][2] += z_cq2 * xf.com2.cq_d + z_cf2 * xf.com2.cf_d + z_hk2 * xf.com2.hk_d;
        xf.vs2[0][3] += z_cq2 * xf.com2.cq_u + z_cf2 * xf.com2.cf_u + z_hk2 * xf.com2.hk_u;

        xf.vsm[0] += z_cq1 * xf.com1.cq_ms + z_cf1 * xf.com1.cf_ms + z_hk1 * xf.com1.hk_ms + z_cq2 * xf.com2.cq_ms + z_cf2 * xf.com2.cf_ms
            + z_hk2 * xf.com2.hk_ms;
        xf.vsr[0] = z_cq1 * xf.com1.cq_re + z_cf1 * xf.com1.cf_re + z_cq2 * xf.com2.cq_re + z_cf2 * xf.com2.cf_re;
        xf.vsx[0] = 0.0;
        xf.vsrez[0] = -rezc;
    }

    // **** Set up momentum equation ****
    let ha = 0.5 * (xf.com1.h + xf.com2.h);
    let ma = 0.5 * (xf.com1.m + xf.com2.m);
    let xa = 0.5 * (xf.com1.x + xf.com2.x);
    let ta = 0.5 * (xf.com1.t + xf.com2.t);
    let hwa = 0.5 * (xf.com1.dw / xf.com1.t + xf.com2.dw / xf.com2.t);

    // set Cf term, using central value CFM for better accuracy in drag
    let cfx = 0.50 * xf.cfm * xa / ta + 0.25 * (xf.com1.cf * xf.com1.x / xf.com1.t + xf.com2.cf * xf.com2.x / xf.com2.t);
    let cfx_xa = 0.50 * xf.cfm / ta;
    let cfx_ta = -0.50 * xf.cfm * xa / (ta * ta);

    let cfx_x1 = 0.25 * xf.com1.cf / xf.com1.t + cfx_xa * 0.5;
    let cfx_x2 = 0.25 * xf.com2.cf / xf.com2.t + cfx_xa * 0.5;
    let cfx_t1 = -0.25 * xf.com1.cf * xf.com1.x / xf.com1.t.powi(2) + cfx_ta * 0.5;
    let cfx_t2 = -0.25 * xf.com2.cf * xf.com2.x / xf.com2.t.powi(2) + cfx_ta * 0.5;
    let cfx_cf1 = 0.25 * xf.com1.x / xf.com1.t;
    let cfx_cf2 = 0.25 * xf.com2.x / xf.com2.t;
    let cfx_cfm = 0.50 * xa / ta;

    let btmp = ha + 2.0 - ma + hwa;

    let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
    let z_cfx = -xlog * 0.5;
    let z_ha = ulog;
    let z_hwa = ulog;
    let z_ma = -ulog;
    let z_xl = -ddlog * 0.5 * cfx;
    let z_ul = ddlog * btmp;
    let z_tl = ddlog;

    let z_cfm = z_cfx * cfx_cfm;
    let z_cf1 = z_cfx * cfx_cf1;
    let z_cf2 = z_cfx * cfx_cf2;

    let z_t1 = -z_tl / xf.com1.t + z_cfx * cfx_t1 + z_hwa * 0.5 * (-xf.com1.dw / xf.com1.t.powi(2));
    let z_t2 = z_tl / xf.com2.t + z_cfx * cfx_t2 + z_hwa * 0.5 * (-xf.com2.dw / xf.com2.t.powi(2));
    let z_x1 = -z_xl / xf.com1.x + z_cfx * cfx_x1;
    let z_x2 = z_xl / xf.com2.x + z_cfx * cfx_x2;
    let z_u1 = -z_ul / xf.com1.u;
    let z_u2 = z_ul / xf.com2.u;

    xf.vs1[1][1] = 0.5 * z_ha * xf.com1.h_t + z_cfm * xf.cfm_t1 + z_cf1 * xf.com1.cf_t + z_t1;
    xf.vs1[1][2] = 0.5 * z_ha * xf.com1.h_d + z_cfm * xf.cfm_d1 + z_cf1 * xf.com1.cf_d;
    xf.vs1[1][3] = 0.5 * z_ma * xf.com1.m_u + z_cfm * xf.cfm_u1 + z_cf1 * xf.com1.cf_u + z_u1;
    xf.vs1[1][4] = z_x1;
    xf.vs2[1][1] = 0.5 * z_ha * xf.com2.h_t + z_cfm * xf.cfm_t2 + z_cf2 * xf.com2.cf_t + z_t2;
    xf.vs2[1][2] = 0.5 * z_ha * xf.com2.h_d + z_cfm * xf.cfm_d2 + z_cf2 * xf.com2.cf_d;
    xf.vs2[1][3] = 0.5 * z_ma * xf.com2.m_u + z_cfm * xf.cfm_u2 + z_cf2 * xf.com2.cf_u + z_u2;
    xf.vs2[1][4] = z_x2;

    xf.vsm[1] = 0.5 * z_ma * xf.com1.m_ms + z_cfm * xf.cfm_ms + z_cf1 * xf.com1.cf_ms + 0.5 * z_ma * xf.com2.m_ms
        + z_cf2 * xf.com2.cf_ms;
    xf.vsr[1] = z_cfm * xf.cfm_re + z_cf1 * xf.com1.cf_re + z_cf2 * xf.com2.cf_re;
    xf.vsx[1] = 0.0;
    xf.vsrez[1] = -rezt;

    // **** Set up shape parameter equation ****
    let xot1 = xf.com1.x / xf.com1.t;
    let xot2 = xf.com2.x / xf.com2.t;

    let ha = 0.5 * (xf.com1.h + xf.com2.h);
    let hsa = 0.5 * (xf.com1.hs + xf.com2.hs);
    let hca = 0.5 * (xf.com1.hc + xf.com2.hc);
    let hwa = 0.5 * (xf.com1.dw / xf.com1.t + xf.com2.dw / xf.com2.t);

    let dix = (1.0 - upw) * xf.com1.di * xot1 + upw * xf.com2.di * xot2;
    let cfx = (1.0 - upw) * xf.com1.cf * xot1 + upw * xf.com2.cf * xot2;
    let dix_upw = xf.com2.di * xot2 - xf.com1.di * xot1;
    let cfx_upw = xf.com2.cf * xot2 - xf.com1.cf * xot1;

    let btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

    let rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
    let z_cfx = xlog * 0.5;
    let z_dix = -xlog;
    let z_hca = 2.0 * ulog / hsa;
    let z_ha = -ulog;
    let z_hwa = -ulog;
    let z_xl = ddlog * (0.5 * cfx - dix);
    let z_ul = ddlog * btmp;
    let z_hl = ddlog;

    let z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

    let z_hs1 = -hca * ulog / hsa.powi(2) - z_hl / xf.com1.hs;
    let z_hs2 = -hca * ulog / hsa.powi(2) + z_hl / xf.com2.hs;

    let z_cf1 = (1.0 - upw) * z_cfx * xot1;
    let z_cf2 = upw * z_cfx * xot2;
    let z_di1 = (1.0 - upw) * z_dix * xot1;
    let z_di2 = upw * z_dix * xot2;

    let mut z_t1 = (1.0 - upw) * (z_cfx * xf.com1.cf + z_dix * xf.com1.di) * (-xot1 / xf.com1.t);
    let mut z_t2 = upw * (z_cfx * xf.com2.cf + z_dix * xf.com2.di) * (-xot2 / xf.com2.t);
    let z_x1 = (1.0 - upw) * (z_cfx * xf.com1.cf + z_dix * xf.com1.di) / xf.com1.t - z_xl / xf.com1.x;
    let z_x2 = upw * (z_cfx * xf.com2.cf + z_dix * xf.com2.di) / xf.com2.t + z_xl / xf.com2.x;
    let z_u1 = -z_ul / xf.com1.u;
    let z_u2 = z_ul / xf.com2.u;

    z_t1 += z_hwa * 0.5 * (-xf.com1.dw / xf.com1.t.powi(2));
    z_t2 += z_hwa * 0.5 * (-xf.com2.dw / xf.com2.t.powi(2));

    xf.vs1[2][0] = z_di1 * xf.com1.di_s;
    xf.vs1[2][1] = z_hs1 * xf.com1.hs_t + z_cf1 * xf.com1.cf_t + z_di1 * xf.com1.di_t + z_t1;
    xf.vs1[2][2] = z_hs1 * xf.com1.hs_d + z_cf1 * xf.com1.cf_d + z_di1 * xf.com1.di_d;
    xf.vs1[2][3] = z_hs1 * xf.com1.hs_u + z_cf1 * xf.com1.cf_u + z_di1 * xf.com1.di_u + z_u1;
    xf.vs1[2][4] = z_x1;
    xf.vs2[2][0] = z_di2 * xf.com2.di_s;
    xf.vs2[2][1] = z_hs2 * xf.com2.hs_t + z_cf2 * xf.com2.cf_t + z_di2 * xf.com2.di_t + z_t2;
    xf.vs2[2][2] = z_hs2 * xf.com2.hs_d + z_cf2 * xf.com2.cf_d + z_di2 * xf.com2.di_d;
    xf.vs2[2][3] = z_hs2 * xf.com2.hs_u + z_cf2 * xf.com2.cf_u + z_di2 * xf.com2.di_u + z_u2;
    xf.vs2[2][4] = z_x2;
    xf.vsm[2] = z_hs1 * xf.com1.hs_ms + z_cf1 * xf.com1.cf_ms + z_di1 * xf.com1.di_ms + z_hs2 * xf.com2.hs_ms
        + z_cf2 * xf.com2.cf_ms + z_di2 * xf.com2.di_ms;
    xf.vsr[2] = z_hs1 * xf.com1.hs_re + z_cf1 * xf.com1.cf_re + z_di1 * xf.com1.di_re + z_hs2 * xf.com2.hs_re
        + z_cf2 * xf.com2.cf_re + z_di2 * xf.com2.di_re;

    xf.vs1[2][1] += 0.5 * (z_hca * xf.com1.hc_t + z_ha * xf.com1.h_t) + z_upw * upw_t1;
    xf.vs1[2][2] += 0.5 * (z_hca * xf.com1.hc_d + z_ha * xf.com1.h_d) + z_upw * upw_d1;
    xf.vs1[2][3] += 0.5 * (z_hca * xf.com1.hc_u) + z_upw * upw_u1;
    xf.vs2[2][1] += 0.5 * (z_hca * xf.com2.hc_t + z_ha * xf.com2.h_t) + z_upw * upw_t2;
    xf.vs2[2][2] += 0.5 * (z_hca * xf.com2.hc_d + z_ha * xf.com2.h_d) + z_upw * upw_d2;
    xf.vs2[2][3] += 0.5 * (z_hca * xf.com2.hc_u) + z_upw * upw_u2;

    xf.vsm[2] += 0.5 * (z_hca * xf.com1.hc_ms) + z_upw * upw_ms + 0.5 * (z_hca * xf.com2.hc_ms);

    xf.vsx[2] = 0.0;
    xf.vsrez[2] = -rezh;
}

/// Sets up the Newton system governing the transition interval (TRDIF).
/// The laminar part X1 < xi < XT and the turbulent part XT < xi < X2 are
/// simply summed.
pub fn trdif(xf: &mut Xfoil) {
    let mut bl1 = [[0.0f64; 5]; 4];
    let mut bl2 = [[0.0f64; 5]; 4];
    let mut bt1 = [[0.0f64; 5]; 4];
    let mut bt2 = [[0.0f64; 5]; 4];
    let mut blm = [0.0f64; 4];
    let mut blr = [0.0f64; 4];
    let mut blrez = [0.0f64; 4];
    let mut blx = [0.0f64; 4];
    let mut btm = [0.0f64; 4];
    let mut btr = [0.0f64; 4];
    let mut btrez = [0.0f64; 4];
    let mut btx = [0.0f64; 4];

    // save variables and sensitivities for future restoration
    xf.c1sav = xf.com1;
    xf.c2sav = xf.com2;

    // weighting factors for linear interpolation to transition point
    let wf2 = (xf.xt - xf.com1.x) / (xf.com2.x - xf.com1.x);
    let wf2_xt = 1.0 / (xf.com2.x - xf.com1.x);

    let wf2_a1 = wf2_xt * xf.xt_a1;
    let wf2_x1 = wf2_xt * xf.xt_x1 + (wf2 - 1.0) / (xf.com2.x - xf.com1.x);
    let wf2_x2 = wf2_xt * xf.xt_x2 - wf2 / (xf.com2.x - xf.com1.x);
    let wf2_t1 = wf2_xt * xf.xt_t1;
    let wf2_t2 = wf2_xt * xf.xt_t2;
    let wf2_d1 = wf2_xt * xf.xt_d1;
    let wf2_d2 = wf2_xt * xf.xt_d2;
    let wf2_u1 = wf2_xt * xf.xt_u1;
    let wf2_u2 = wf2_xt * xf.xt_u2;
    let wf2_ms = wf2_xt * xf.xt_ms;
    let wf2_re = wf2_xt * xf.xt_re;
    let wf2_xf = wf2_xt * xf.xt_xf;

    let wf1 = 1.0 - wf2;
    let wf1_a1 = -wf2_a1;
    let wf1_x1 = -wf2_x1;
    let wf1_x2 = -wf2_x2;
    let wf1_t1 = -wf2_t1;
    let wf1_t2 = -wf2_t2;
    let wf1_d1 = -wf2_d1;
    let wf1_d2 = -wf2_d2;
    let wf1_u1 = -wf2_u1;
    let wf1_u2 = -wf2_u2;
    let wf1_ms = -wf2_ms;
    let wf1_re = -wf2_re;
    let wf1_xf = -wf2_xf;

    // **** FIRST, do laminar part between X1 and XT ****
    // interpolate primary variables to transition point
    let tt = xf.com1.t * wf1 + xf.com2.t * wf2;
    let tt_a1 = xf.com1.t * wf1_a1 + xf.com2.t * wf2_a1;
    let tt_x1 = xf.com1.t * wf1_x1 + xf.com2.t * wf2_x1;
    let tt_x2 = xf.com1.t * wf1_x2 + xf.com2.t * wf2_x2;
    let tt_t1 = xf.com1.t * wf1_t1 + xf.com2.t * wf2_t1 + wf1;
    let tt_t2 = xf.com1.t * wf1_t2 + xf.com2.t * wf2_t2 + wf2;
    let tt_d1 = xf.com1.t * wf1_d1 + xf.com2.t * wf2_d1;
    let tt_d2 = xf.com1.t * wf1_d2 + xf.com2.t * wf2_d2;
    let tt_u1 = xf.com1.t * wf1_u1 + xf.com2.t * wf2_u1;
    let tt_u2 = xf.com1.t * wf1_u2 + xf.com2.t * wf2_u2;
    let tt_ms = xf.com1.t * wf1_ms + xf.com2.t * wf2_ms;
    let tt_re = xf.com1.t * wf1_re + xf.com2.t * wf2_re;
    let tt_xf = xf.com1.t * wf1_xf + xf.com2.t * wf2_xf;

    let dt = xf.com1.d * wf1 + xf.com2.d * wf2;
    let dt_a1 = xf.com1.d * wf1_a1 + xf.com2.d * wf2_a1;
    let dt_x1 = xf.com1.d * wf1_x1 + xf.com2.d * wf2_x1;
    let dt_x2 = xf.com1.d * wf1_x2 + xf.com2.d * wf2_x2;
    let dt_t1 = xf.com1.d * wf1_t1 + xf.com2.d * wf2_t1;
    let dt_t2 = xf.com1.d * wf1_t2 + xf.com2.d * wf2_t2;
    let dt_d1 = xf.com1.d * wf1_d1 + xf.com2.d * wf2_d1 + wf1;
    let dt_d2 = xf.com1.d * wf1_d2 + xf.com2.d * wf2_d2 + wf2;
    let dt_u1 = xf.com1.d * wf1_u1 + xf.com2.d * wf2_u1;
    let dt_u2 = xf.com1.d * wf1_u2 + xf.com2.d * wf2_u2;
    let dt_ms = xf.com1.d * wf1_ms + xf.com2.d * wf2_ms;
    let dt_re = xf.com1.d * wf1_re + xf.com2.d * wf2_re;
    let dt_xf = xf.com1.d * wf1_xf + xf.com2.d * wf2_xf;

    let ut = xf.com1.u * wf1 + xf.com2.u * wf2;
    let ut_a1 = xf.com1.u * wf1_a1 + xf.com2.u * wf2_a1;
    let ut_x1 = xf.com1.u * wf1_x1 + xf.com2.u * wf2_x1;
    let ut_x2 = xf.com1.u * wf1_x2 + xf.com2.u * wf2_x2;
    let ut_t1 = xf.com1.u * wf1_t1 + xf.com2.u * wf2_t1;
    let ut_t2 = xf.com1.u * wf1_t2 + xf.com2.u * wf2_t2;
    let ut_d1 = xf.com1.u * wf1_d1 + xf.com2.u * wf2_d1;
    let ut_d2 = xf.com1.u * wf1_d2 + xf.com2.u * wf2_d2;
    let ut_u1 = xf.com1.u * wf1_u1 + xf.com2.u * wf2_u1 + wf1;
    let ut_u2 = xf.com1.u * wf1_u2 + xf.com2.u * wf2_u2 + wf2;
    let ut_ms = xf.com1.u * wf1_ms + xf.com2.u * wf2_ms;
    let ut_re = xf.com1.u * wf1_re + xf.com2.u * wf2_re;
    let ut_xf = xf.com1.u * wf1_xf + xf.com2.u * wf2_xf;

    // set primary "T" variables at XT (really placed into "2" variables)
    xf.com2.x = xf.xt;
    xf.com2.t = tt;
    xf.com2.d = dt;
    xf.com2.u = ut;

    xf.com2.ampl = xf.amcrit;
    xf.com2.s = 0.0;

    // calculate laminar secondary "T" variables
    blkin(xf);
    blvar(xf, 1);

    // calculate X1-XT midpoint CFM value
    blmid(xf, 1);

    // at this point, all "2" variables are really "T" variables at XT

    // set up Newton system for dAm, dTh, dDs, dUe, dXi at X1 and XT
    bldif(xf, 1);

    // The current Newton system is in terms of "1" and "T" variables, so
    // calculate its equivalent in terms of "1" and "2" variables.
    for k in 1..3 {
        blrez[k] = xf.vsrez[k];
        blm[k] = xf.vsm[k] + xf.vs2[k][1] * tt_ms + xf.vs2[k][2] * dt_ms + xf.vs2[k][3] * ut_ms + xf.vs2[k][4] * xf.xt_ms;
        blr[k] = xf.vsr[k] + xf.vs2[k][1] * tt_re + xf.vs2[k][2] * dt_re + xf.vs2[k][3] * ut_re + xf.vs2[k][4] * xf.xt_re;
        blx[k] = xf.vsx[k] + xf.vs2[k][1] * tt_xf + xf.vs2[k][2] * dt_xf + xf.vs2[k][3] * ut_xf + xf.vs2[k][4] * xf.xt_xf;

        bl1[k][0] = xf.vs1[k][0] + xf.vs2[k][1] * tt_a1 + xf.vs2[k][2] * dt_a1 + xf.vs2[k][3] * ut_a1 + xf.vs2[k][4] * xf.xt_a1;
        bl1[k][1] = xf.vs1[k][1] + xf.vs2[k][1] * tt_t1 + xf.vs2[k][2] * dt_t1 + xf.vs2[k][3] * ut_t1 + xf.vs2[k][4] * xf.xt_t1;
        bl1[k][2] = xf.vs1[k][2] + xf.vs2[k][1] * tt_d1 + xf.vs2[k][2] * dt_d1 + xf.vs2[k][3] * ut_d1 + xf.vs2[k][4] * xf.xt_d1;
        bl1[k][3] = xf.vs1[k][3] + xf.vs2[k][1] * tt_u1 + xf.vs2[k][2] * dt_u1 + xf.vs2[k][3] * ut_u1 + xf.vs2[k][4] * xf.xt_u1;
        bl1[k][4] = xf.vs1[k][4] + xf.vs2[k][1] * tt_x1 + xf.vs2[k][2] * dt_x1 + xf.vs2[k][3] * ut_x1 + xf.vs2[k][4] * xf.xt_x1;

        bl2[k][0] = 0.0;
        bl2[k][1] = xf.vs2[k][1] * tt_t2 + xf.vs2[k][2] * dt_t2 + xf.vs2[k][3] * ut_t2 + xf.vs2[k][4] * xf.xt_t2;
        bl2[k][2] = xf.vs2[k][1] * tt_d2 + xf.vs2[k][2] * dt_d2 + xf.vs2[k][3] * ut_d2 + xf.vs2[k][4] * xf.xt_d2;
        bl2[k][3] = xf.vs2[k][1] * tt_u2 + xf.vs2[k][2] * dt_u2 + xf.vs2[k][3] * ut_u2 + xf.vs2[k][4] * xf.xt_u2;
        bl2[k][4] = xf.vs2[k][1] * tt_x2 + xf.vs2[k][2] * dt_x2 + xf.vs2[k][3] * ut_x2 + xf.vs2[k][4] * xf.xt_x2;
    }

    // **** SECOND, set up turbulent part between XT and X2 ****
    // calculate equilibrium shear coefficient CQT at transition point
    blvar(xf, 2);

    // set initial shear coefficient value ST at transition point
    let ctr = xf.ctrcon * (-xf.ctrcex / (xf.com2.hk - 1.0)).exp();
    let ctr_hk2 = ctr * xf.ctrcex / (xf.com2.hk - 1.0).powi(2);

    let st = ctr * xf.com2.cq;
    let st_tt = ctr * xf.com2.cq_t + xf.com2.cq * ctr_hk2 * xf.com2.hk_t;
    let st_dt = ctr * xf.com2.cq_d + xf.com2.cq * ctr_hk2 * xf.com2.hk_d;
    let st_ut = ctr * xf.com2.cq_u + xf.com2.cq * ctr_hk2 * xf.com2.hk_u;
    let st_ms = ctr * xf.com2.cq_ms + xf.com2.cq * ctr_hk2 * xf.com2.hk_ms;
    let st_re = ctr * xf.com2.cq_re;

    // calculate ST sensitivities wrt the actual "1" and "2" variables
    let st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
    let st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
    let st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
    let st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
    let st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
    let st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
    let st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
    let st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
    let st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
    let st_ms = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
    let st_re = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
    let st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

    xf.com2.ampl = 0.0;
    xf.com2.s = st;

    // recalculate turbulent secondary "T" variables using proper CTI
    blvar(xf, 2);

    // set "1" variables to "T" variables and reset "2" variables
    xf.com1 = xf.com2;
    xf.com2 = xf.c2sav;

    // calculate XT-X2 midpoint CFM value
    blmid(xf, 2);

    // set up Newton system for dCt, dTh, dDs, dUe, dXi at XT and X2
    bldif(xf, 2);

    // convert sensitivities wrt "T" variables into sensitivities wrt "1" and "2"
    for k in 0..3 {
        btrez[k] = xf.vsrez[k];
        btm[k] = xf.vsm[k]
            + xf.vs1[k][0] * st_ms
            + xf.vs1[k][1] * tt_ms
            + xf.vs1[k][2] * dt_ms
            + xf.vs1[k][3] * ut_ms
            + xf.vs1[k][4] * xf.xt_ms;
        btr[k] = xf.vsr[k]
            + xf.vs1[k][0] * st_re
            + xf.vs1[k][1] * tt_re
            + xf.vs1[k][2] * dt_re
            + xf.vs1[k][3] * ut_re
            + xf.vs1[k][4] * xf.xt_re;
        btx[k] = xf.vsx[k]
            + xf.vs1[k][0] * st_xf
            + xf.vs1[k][1] * tt_xf
            + xf.vs1[k][2] * dt_xf
            + xf.vs1[k][3] * ut_xf
            + xf.vs1[k][4] * xf.xt_xf;

        bt1[k][0] = xf.vs1[k][0] * st_a1 + xf.vs1[k][1] * tt_a1 + xf.vs1[k][2] * dt_a1 + xf.vs1[k][3] * ut_a1 + xf.vs1[k][4] * xf.xt_a1;
        bt1[k][1] = xf.vs1[k][0] * st_t1 + xf.vs1[k][1] * tt_t1 + xf.vs1[k][2] * dt_t1 + xf.vs1[k][3] * ut_t1 + xf.vs1[k][4] * xf.xt_t1;
        bt1[k][2] = xf.vs1[k][0] * st_d1 + xf.vs1[k][1] * tt_d1 + xf.vs1[k][2] * dt_d1 + xf.vs1[k][3] * ut_d1 + xf.vs1[k][4] * xf.xt_d1;
        bt1[k][3] = xf.vs1[k][0] * st_u1 + xf.vs1[k][1] * tt_u1 + xf.vs1[k][2] * dt_u1 + xf.vs1[k][3] * ut_u1 + xf.vs1[k][4] * xf.xt_u1;
        bt1[k][4] = xf.vs1[k][0] * st_x1 + xf.vs1[k][1] * tt_x1 + xf.vs1[k][2] * dt_x1 + xf.vs1[k][3] * ut_x1 + xf.vs1[k][4] * xf.xt_x1;

        bt2[k][0] = xf.vs2[k][0];
        bt2[k][1] = xf.vs2[k][1]
            + xf.vs1[k][0] * st_t2
            + xf.vs1[k][1] * tt_t2
            + xf.vs1[k][2] * dt_t2
            + xf.vs1[k][3] * ut_t2
            + xf.vs1[k][4] * xf.xt_t2;
        bt2[k][2] = xf.vs2[k][2]
            + xf.vs1[k][0] * st_d2
            + xf.vs1[k][1] * tt_d2
            + xf.vs1[k][2] * dt_d2
            + xf.vs1[k][3] * ut_d2
            + xf.vs1[k][4] * xf.xt_d2;
        bt2[k][3] = xf.vs2[k][3]
            + xf.vs1[k][0] * st_u2
            + xf.vs1[k][1] * tt_u2
            + xf.vs1[k][2] * dt_u2
            + xf.vs1[k][3] * ut_u2
            + xf.vs1[k][4] * xf.xt_u2;
        bt2[k][4] = xf.vs2[k][4]
            + xf.vs1[k][0] * st_x2
            + xf.vs1[k][1] * tt_x2
            + xf.vs1[k][2] * dt_x2
            + xf.vs1[k][3] * ut_x2
            + xf.vs1[k][4] * xf.xt_x2;
    }

    // add up laminar and turbulent parts to get final system
    xf.vsrez[0] = btrez[0];
    xf.vsrez[1] = blrez[1] + btrez[1];
    xf.vsrez[2] = blrez[2] + btrez[2];
    xf.vsm[0] = btm[0];
    xf.vsm[1] = blm[1] + btm[1];
    xf.vsm[2] = blm[2] + btm[2];
    xf.vsr[0] = btr[0];
    xf.vsr[1] = blr[1] + btr[1];
    xf.vsr[2] = blr[2] + btr[2];
    xf.vsx[0] = btx[0];
    xf.vsx[1] = blx[1] + btx[1];
    xf.vsx[2] = blx[2] + btx[2];
    for l in 0..5 {
        xf.vs1[0][l] = bt1[0][l];
        xf.vs2[0][l] = bt2[0][l];
        xf.vs1[1][l] = bl1[1][l] + bt1[1][l];
        xf.vs2[1][l] = bl2[1][l] + bt2[1][l];
        xf.vs1[2][l] = bl1[2][l] + bt1[2][l];
        xf.vs2[2][l] = bl2[2][l] + bt2[2][l];
    }

    // restore "1" quantities which got clobbered in all of the numerical
    // gymnastics above
    xf.com1 = xf.c1sav;
}
