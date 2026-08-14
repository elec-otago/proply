//! Boundary-layer marching and Newton update (port of `m_xbl.f90`).
//!
//! `setbl` assembles the coupled viscous-inviscid Newton system by marching
//! both BLs and the wake; `update` applies the solution deltas with
//! under-relaxation; `mrchue`/`mrchdu` initialize and march the BL in
//! direct/mixed mode.

use crate::blsys::{blkin, blmid, blprv, blsys, blvar, tesys, trchek, hkin};
use crate::panel::ueset;
use crate::solve::gauss;
use crate::spline::{seval, sinvrt, splind};
use crate::state::{Xfoil, IZX, NCOM, QOPI};

/// Returns the BL Newton scratch Vecs (pulled out of `xf` via `mem::take` at
/// the top of `setbl`) back to the state, so the buffers are reused on the
/// next call instead of being re-allocated.
macro_rules! restore_bl_scratch {
    ($xf:expr, $usav:expr, $ule1_m:expr, $ule2_m:expr, $ute1_m:expr, $ute2_m:expr, $u1_m:expr, $d1_m:expr, $u2_m:expr, $d2_m:expr) => {
        $xf.bl_usav = $usav;
        $xf.bl_ule1_m = $ule1_m;
        $xf.bl_ule2_m = $ule2_m;
        $xf.bl_ute1_m = $ute1_m;
        $xf.bl_ute2_m = $ute2_m;
        $xf.bl_u1_m = $u1_m;
        $xf.bl_d1_m = $d1_m;
        $xf.bl_u2_m = $u2_m;
        $xf.bl_d2_m = $d2_m;
    };
}

/// Assembles the BL Newton system coefficients for the current BL variables
/// and edge velocities, incorporating the local BL system coefficients into
/// the global Newton system.  Returns false if the BL initialization failed
/// (and abort_on_nan is set).
pub fn setbl(xf: &mut Xfoil) -> bool {
    // set the CL used to define Mach, Reynolds numbers
    let clmr = if xf.lalfa { xf.cl } else { xf.clspec };

    // set current MINF(CL) and REINF(CL)
    let (ma_clmr, re_clmr) = crate::s_xfoil::mrcl(xf, clmr);
    let msq_clmr = 2.0 * xf.minf * ma_clmr;

    // set compressibility parameter TKLAM and derivative TK_MSQ
    crate::s_xfoil::comset(xf);

    // set gas constant (= Cp/Cv)
    xf.gambl = xf.gamma;
    xf.gm1bl = xf.gamm1;

    // set parameters for compressibility correction
    xf.qinfbl = xf.qinf;
    xf.tkbl = xf.tklam;
    xf.tkbl_ms = xf.tkl_msq;

    // stagnation density and 1/enthalpy
    xf.rstbl = (1.0 + 0.5 * xf.gm1bl * xf.minf.powi(2)).powf(1.0 / xf.gm1bl);
    xf.rstbl_ms = 0.5 * xf.rstbl / (1.0 + 0.5 * xf.gm1bl * xf.minf.powi(2));

    xf.hstinv = xf.gm1bl * (xf.minf / xf.qinfbl).powi(2) / (1.0 + 0.5 * xf.gm1bl * xf.minf.powi(2));
    xf.hstinv_ms = xf.gm1bl * (1.0 / xf.qinfbl).powi(2) / (1.0 + 0.5 * xf.gm1bl * xf.minf.powi(2))
        - 0.5 * xf.gm1bl * xf.hstinv / (1.0 + 0.5 * xf.gm1bl * xf.minf.powi(2));

    // set Reynolds number based on freestream density, velocity, viscosity
    let herat = 1.0 - 0.5 * xf.qinfbl.powi(2) * xf.hstinv;
    let herat_ms = -0.5 * xf.qinfbl.powi(2) * xf.hstinv_ms;

    xf.reybl = xf.reinf * herat.powf(1.5) * (1.0 + xf.hvrat) / (herat + xf.hvrat);
    xf.reybl_re = herat.powf(1.5) * (1.0 + xf.hvrat) / (herat + xf.hvrat);
    xf.reybl_ms = xf.reybl * (1.5 / herat - 1.0 / (herat + xf.hvrat)) * herat_ms;

    xf.idampv = xf.idamp;

    // save TE thickness
    xf.dwte = xf.wgap[0];

    if !xf.lbli_ni {
        // initialize BL by marching with Ue (fudge at separation)
        if xf.show_output {
            eprintln!();
            eprintln!("Initializing BL ...");
        }
        let ok = mrchue(xf);
        if ok {
            xf.lbli_ni = true;
        } else if xf.abort_on_nan {
            return false;
        }
    }

    if xf.show_output {
        eprintln!();
    }

    // march BL with current Ue and Ds to establish transition
    let ok = mrchdu(xf);
    if xf.abort_on_nan && !ok {
        return false;
    }

    // Pull the reusable BL Newton scratch Vecs out of the state so they are
    // locally owned during this call (avoids re-allocation on every iteration
    // while sidestepping the &mut xf aliasing that long-lived field borrows
    // would create).  They are returned to the state at the end of setbl.
    let mut usav = std::mem::take(&mut xf.bl_usav);
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            usav[is][ibl] = xf.uedg[is][ibl];
        }
    }

    ueset(xf);

    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let temp = usav[is][ibl];
            usav[is][ibl] = xf.uedg[is][ibl];
            xf.uedg[is][ibl] = temp;
        }
    }

    let ile1 = xf.ipan[0][2] as usize;
    let ile2 = xf.ipan[1][2] as usize;
    let ite1 = xf.ipan[0][xf.iblte[0] as usize] as usize;
    let ite2 = xf.ipan[1][xf.iblte[1] as usize] as usize;

    let jvte1 = xf.isys[0][xf.iblte[0] as usize] as usize;
    let jvte2 = xf.isys[1][xf.iblte[1] as usize] as usize;

    let dule1 = xf.uedg[0][2] - usav[0][2];
    let dule2 = xf.uedg[1][2] - usav[1][2];

    // set LE and TE Ue sensitivities wrt all m values.  The js/jbl loop below
    // writes every jv in 1..=nsys (isys is a bijection onto that range), so the
    // scratch arrays need not be re-zeroed here.
    let mut ule1_m = std::mem::take(&mut xf.bl_ule1_m);
    let mut ule2_m = std::mem::take(&mut xf.bl_ule2_m);
    let mut ute1_m = std::mem::take(&mut xf.bl_ute1_m);
    let mut ute2_m = std::mem::take(&mut xf.bl_ute2_m);
    for js in 0..2 {
        for jbl in 2..=xf.nbl[js] as usize {
            let j = xf.ipan[js][jbl] as usize;
            let jv = xf.isys[js][jbl] as usize;
            ule1_m[jv] = -xf.vti[0][2] * xf.vti[js][jbl] * xf.dij_t[ile1 * IZX + j];
            ule2_m[jv] = -xf.vti[1][2] * xf.vti[js][jbl] * xf.dij_t[ile2 * IZX + j];
            ute1_m[jv] = -xf.vti[0][xf.iblte[0] as usize] * xf.vti[js][jbl] * xf.dij_t[ite1 * IZX + j];
            ute2_m[jv] = -xf.vti[1][xf.iblte[1] as usize] * xf.vti[js][jbl] * xf.dij_t[ite2 * IZX + j];
        }
    }

    let ule1_a = xf.uinv_a[0][2];
    let ule2_a = xf.uinv_a[1][2];

    xf.tindex[0] = 0.0;
    xf.tindex[1] = 0.0;

    // u1_m/d1_m and u2_m/d2_m are also reused via mem::take; they are returned
    // to the state after the per-side sweep.
    let mut u1_m = std::mem::take(&mut xf.bl_u1_m);
    let mut d1_m = std::mem::take(&mut xf.bl_d1_m);
    let mut u2_m = std::mem::take(&mut xf.bl_u2_m);
    let mut d2_m = std::mem::take(&mut xf.bl_d2_m);

    // **** go over each boundary layer/wake ****
    for is in 0..2 {
        // there is no station "1" at similarity, so zero everything out.
        // u1_m/d1_m persist across the ibl sweep (carried via u1_m=u2_m), so
        // they must be cleared once per side; only 0..=nsys is ever read.
        for v in &mut u1_m[..=xf.nsys] {
            *v = 0.0;
        }
        for v in &mut d1_m[..=xf.nsys] {
            *v = 0.0;
        }
        let mut u1_a = 0.0;
        let mut d1_a = 0.0;

        let mut due1 = 0.0;
        let mut dds1 = 0.0;

        // similarity station pressure gradient parameter x/u du/dx
        xf.bule = 1.0;

        xf.amcrit = xf.acrit[is];

        // set forced transition arc length position
        xifset(xf, is);

        xf.tran = false;
        xf.turb = false;

        // **** sweep downstream setting up BL equation linearizations ****
        for ibl in 2..=xf.nbl[is] as usize {
            let iv = xf.isys[is][ibl] as usize;

            xf.simi = ibl == 2;
            xf.wake = ibl > xf.iblte[is] as usize;
            xf.tran = ibl == xf.itran[is] as usize;
            xf.turb = ibl > xf.itran[is] as usize;

            let i = xf.ipan[is][ibl] as usize;

            // set primary variables for current station
            let xsi = xf.xssi[is][ibl];
            let mut ami = 0.0;
            let mut cti = 0.0;
            if ibl < xf.itran[is] as usize {
                ami = xf.ctau[is][ibl];
            }
            if ibl >= xf.itran[is] as usize {
                cti = xf.ctau[is][ibl];
            }
            let uei = xf.uedg[is][ibl];
            let thi = xf.thet[is][ibl];
            let mdi = xf.mass[is][ibl];

            let mut dsi = mdi / uei;

            let dswaki;
            if xf.wake {
                let iw = ibl - xf.iblte[is] as usize;
                // WAKE GAP indexing: wgap is a 0-based array of length IWX.
                // BL station `ibl == iblte[is] + iw` (iw = 1..nw) must read
                // wgap[iw-1].  Using wgap[iw] here would read one element past
                // the intended slot and, at iw == IWX, past the array end.
                dswaki = xf.wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            // set derivatives of DSI (= D2)
            let d2_m2 = 1.0 / uei;
            let d2_u2 = -dsi / uei;

            // u2_m/d2_m are rebuilt for every station; the js/jbl loop writes
            // every jv in 1..=nsys (isys is a bijection), so no re-zero is
            // needed before reuse.
            for js in 0..2 {
                for jbl in 2..=xf.nbl[js] as usize {
                    let j = xf.ipan[js][jbl] as usize;
                    let jv = xf.isys[js][jbl] as usize;
                    u2_m[jv] = -xf.vti[is][ibl] * xf.vti[js][jbl] * xf.dij_t[i * IZX + j];
                    d2_m[jv] = d2_u2 * u2_m[jv];
                }
            }
            d2_m[iv] += d2_m2;

            let u2_a = xf.uinv_a[is][ibl];
            let d2_a = d2_u2 * u2_a;

            // "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
            let due2 = xf.uedg[is][ibl] - usav[is][ibl];
            let dds2 = d2_u2 * due2;

            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);

            // check for transition and set TRAN, XT, etc. if found
            if xf.tran {
                let ok = trchek(xf);
                if xf.abort_on_nan && !ok {
                    restore_bl_scratch!(xf, usav, ule1_m, ule2_m, ute1_m, ute2_m, u1_m, d1_m, u2_m, d2_m);
                    return false;
                }
                ami = xf.com2.ampl;
            }
            if ibl == xf.itran[is] as usize && !xf.tran && xf.show_output {
                eprintln!("SETBL: Xtr???  n1 n2: {} {}", xf.com1.ampl, xf.com2.ampl);
            }

            // assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
            // at the previous "1" station and the current "2" station
            let (tte_tte1, tte_tte2, dte_mte1, dte_ute1, dte_mte2, dte_ute2, cte_cte1, cte_cte2, cte_tte1, cte_tte2);
            if ibl == xf.iblte[is] as usize + 1 {
                // define quantities at start of wake, adding TE base thickness to Dstar
                let tte = xf.thet[0][xf.iblte[0] as usize] + xf.thet[1][xf.iblte[1] as usize];
                let dte = xf.dstr[0][xf.iblte[0] as usize] + xf.dstr[1][xf.iblte[1] as usize] + xf.ante;
                let cte =
                    (xf.ctau[0][xf.iblte[0] as usize] * xf.thet[0][xf.iblte[0] as usize]
                        + xf.ctau[1][xf.iblte[1] as usize] * xf.thet[1][xf.iblte[1] as usize])
                        / tte;
                tesys(xf, cte, tte, dte);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                dte_mte1 = 1.0 / xf.uedg[0][xf.iblte[0] as usize];
                dte_ute1 = -xf.dstr[0][xf.iblte[0] as usize] / xf.uedg[0][xf.iblte[0] as usize];
                dte_mte2 = 1.0 / xf.uedg[1][xf.iblte[1] as usize];
                dte_ute2 = -xf.dstr[1][xf.iblte[1] as usize] / xf.uedg[1][xf.iblte[1] as usize];
                cte_cte1 = xf.thet[0][xf.iblte[0] as usize] / tte;
                cte_cte2 = xf.thet[1][xf.iblte[1] as usize] / tte;
                cte_tte1 = (xf.ctau[0][xf.iblte[0] as usize] - cte) / tte;
                cte_tte2 = (xf.ctau[1][xf.iblte[1] as usize] - cte) / tte;

                // re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
                for js in 0..2 {
                    for jbl in 2..=xf.nbl[js] as usize {
                        let j = xf.ipan[js][jbl] as usize;
                        let jv = xf.isys[js][jbl] as usize;
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                    }
                }
                d1_m[jvte1] += dte_mte1;
                d1_m[jvte2] += dte_mte2;

                // "forced" changes from UEDG --- USAV=UINV+dij*MASS mismatch
                due1 = 0.0;
                dds1 = dte_ute1 * (xf.uedg[0][xf.iblte[0] as usize] - usav[0][xf.iblte[0] as usize])
                    + dte_ute2 * (xf.uedg[1][xf.iblte[1] as usize] - usav[1][xf.iblte[1] as usize]);
            } else {
                blsys(xf);
                tte_tte1 = 0.0;
                tte_tte2 = 0.0;
                dte_mte1 = 0.0;
                dte_ute1 = 0.0;
                dte_mte2 = 0.0;
                dte_ute2 = 0.0;
                cte_cte1 = 0.0;
                cte_cte2 = 0.0;
                cte_tte1 = 0.0;
                cte_tte2 = 0.0;
            }

            // save wall shear and equil. max shear coefficient for plotting output
            xf.tau[is][ibl] = 0.5 * xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.cf;
            xf.dis[is][ibl] = xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.u * xf.com2.di * xf.com2.hs * 0.5;
            xf.ctq[is][ibl] = xf.com2.cq;
            xf.delt[is][ibl] = xf.com2.de;
            xf.uslp[is][ibl] = 1.60 / (1.0 + xf.com2.us);

            // set XI sensitivities wrt LE Ue changes
            let (xi_ule1, xi_ule2);
            if is == 0 {
                xi_ule1 = xf.sst_go;
                xi_ule2 = -xf.sst_gp;
            } else {
                xi_ule1 = -xf.sst_go;
                xi_ule2 = xf.sst_gp;
            }

            // stuff BL system coefficients into main Jacobian matrix.  The vm
            // layout is (i, j, k) with k innermost, so for a fixed iv the
            // nsys*3 elements form one contiguous row.  The original code wrote
            // the k=1,2,3 components in three separate sweeps that all shared
            // the same cache lines; fusing them into a single jv pass over the
            // contiguous row keeps the same values and FP order (bit-identical)
            // while streaming memory once and hoisting the row-offset bounds
            // work out of the loop.
            {
                let base3 = (iv - 1) * IZX * 3;
                for jv in 1..=xf.nsys {
                    let o = base3 + (jv - 1) * 3;
                    xf.vm[o] = xf.vs1[0][2] * d1_m[jv] + xf.vs1[0][3] * u1_m[jv]
                        + xf.vs2[0][2] * d2_m[jv]
                        + xf.vs2[0][3] * u2_m[jv]
                        + (xf.vs1[0][4] + xf.vs2[0][4] + xf.vsx[0]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                    xf.vm[o + 1] = xf.vs1[1][2] * d1_m[jv] + xf.vs1[1][3] * u1_m[jv]
                        + xf.vs2[1][2] * d2_m[jv]
                        + xf.vs2[1][3] * u2_m[jv]
                        + (xf.vs1[1][4] + xf.vs2[1][4] + xf.vsx[1]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                    xf.vm[o + 2] = xf.vs1[2][2] * d1_m[jv] + xf.vs1[2][3] * u1_m[jv]
                        + xf.vs2[2][2] * d2_m[jv]
                        + xf.vs2[2][3] * u2_m[jv]
                        + (xf.vs1[2][4] + xf.vs2[2][4] + xf.vsx[2]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                }
            }

            xf.vb[Xfoil::v_index(iv, 1, 1)] = xf.vs1[0][0];
            xf.vb[Xfoil::v_index(iv, 2, 1)] = xf.vs1[0][1];

            xf.va[Xfoil::v_index(iv, 1, 1)] = xf.vs2[0][0];
            xf.va[Xfoil::v_index(iv, 2, 1)] = xf.vs2[0][1];

            if xf.lalfa {
                xf.vdel[Xfoil::v_index(iv, 2, 1)] = xf.vsr[0] * re_clmr + xf.vsm[0] * msq_clmr;
            } else {
                xf.vdel[Xfoil::v_index(iv, 2, 1)] = (xf.vs1[0][3] * u1_a + xf.vs1[0][2] * d1_a)
                    + (xf.vs2[0][3] * u2_a + xf.vs2[0][2] * d2_a)
                    + (xf.vs1[0][4] + xf.vs2[0][4] + xf.vsx[0]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            xf.vdel[Xfoil::v_index(iv, 1, 1)] = xf.vsrez[0] + (xf.vs1[0][3] * due1 + xf.vs1[0][2] * dds1)
                + (xf.vs2[0][3] * due2 + xf.vs2[0][2] * dds2)
                + (xf.vs1[0][4] + xf.vs2[0][4] + xf.vsx[0]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            xf.vb[Xfoil::v_index(iv, 1, 2)] = xf.vs1[1][0];
            xf.vb[Xfoil::v_index(iv, 2, 2)] = xf.vs1[1][1];

            xf.va[Xfoil::v_index(iv, 1, 2)] = xf.vs2[1][0];
            xf.va[Xfoil::v_index(iv, 2, 2)] = xf.vs2[1][1];

            if xf.lalfa {
                xf.vdel[Xfoil::v_index(iv, 2, 2)] = xf.vsr[1] * re_clmr + xf.vsm[1] * msq_clmr;
            } else {
                xf.vdel[Xfoil::v_index(iv, 2, 2)] = (xf.vs1[1][3] * u1_a + xf.vs1[1][2] * d1_a)
                    + (xf.vs2[1][3] * u2_a + xf.vs2[1][2] * d2_a)
                    + (xf.vs1[1][4] + xf.vs2[1][4] + xf.vsx[1]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            xf.vdel[Xfoil::v_index(iv, 1, 2)] = xf.vsrez[1] + (xf.vs1[1][3] * due1 + xf.vs1[1][2] * dds1)
                + (xf.vs2[1][3] * due2 + xf.vs2[1][2] * dds2)
                + (xf.vs1[1][4] + xf.vs2[1][4] + xf.vsx[1]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            xf.vb[Xfoil::v_index(iv, 1, 3)] = xf.vs1[2][0];
            xf.vb[Xfoil::v_index(iv, 2, 3)] = xf.vs1[2][1];

            xf.va[Xfoil::v_index(iv, 1, 3)] = xf.vs2[2][0];
            xf.va[Xfoil::v_index(iv, 2, 3)] = xf.vs2[2][1];

            if xf.lalfa {
                xf.vdel[Xfoil::v_index(iv, 2, 3)] = xf.vsr[2] * re_clmr + xf.vsm[2] * msq_clmr;
            } else {
                xf.vdel[Xfoil::v_index(iv, 2, 3)] = (xf.vs1[2][3] * u1_a + xf.vs1[2][2] * d1_a)
                    + (xf.vs2[2][3] * u2_a + xf.vs2[2][2] * d2_a)
                    + (xf.vs1[2][4] + xf.vs2[2][4] + xf.vsx[2]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            xf.vdel[Xfoil::v_index(iv, 1, 3)] = xf.vsrez[2] + (xf.vs1[2][3] * due1 + xf.vs1[2][2] * dds1)
                + (xf.vs2[2][3] * due2 + xf.vs2[2][2] * dds2)
                + (xf.vs1[2][4] + xf.vs2[2][4] + xf.vsx[2]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            if ibl == xf.iblte[is] as usize + 1 {
                // redefine coefficients for TTE, DTE, etc
                xf.vz[Xfoil::v_index(1, 1, 1)] = xf.vs1[0][0] * cte_cte1;
                xf.vz[Xfoil::v_index(1, 2, 1)] = xf.vs1[0][0] * cte_tte1 + xf.vs1[0][1] * tte_tte1;
                xf.vb[Xfoil::v_index(iv, 1, 1)] = xf.vs1[0][0] * cte_cte2;
                xf.vb[Xfoil::v_index(iv, 2, 1)] = xf.vs1[0][0] * cte_tte2 + xf.vs1[0][1] * tte_tte2;

                xf.vz[Xfoil::v_index(1, 1, 2)] = xf.vs1[1][0] * cte_cte1;
                xf.vz[Xfoil::v_index(1, 2, 2)] = xf.vs1[1][0] * cte_tte1 + xf.vs1[1][1] * tte_tte1;
                xf.vb[Xfoil::v_index(iv, 1, 2)] = xf.vs1[1][0] * cte_cte2;
                xf.vb[Xfoil::v_index(iv, 2, 2)] = xf.vs1[1][0] * cte_tte2 + xf.vs1[1][1] * tte_tte2;

                xf.vz[Xfoil::v_index(1, 1, 3)] = xf.vs1[2][0] * cte_cte1;
                xf.vz[Xfoil::v_index(1, 2, 3)] = xf.vs1[2][0] * cte_tte1 + xf.vs1[2][1] * tte_tte1;
                xf.vb[Xfoil::v_index(iv, 1, 3)] = xf.vs1[2][0] * cte_cte2;
                xf.vb[Xfoil::v_index(iv, 2, 3)] = xf.vs1[2][0] * cte_tte2 + xf.vs1[2][1] * tte_tte2;
            }

            // turbulent intervals will follow if currently at transition interval
            if xf.tran {
                xf.turb = true;

                // save transition location
                xf.itran[is] = ibl as i32;
                xf.tforce[is] = xf.trforc;
                xf.xssitr[is] = xf.xt;

                // interpolate airfoil geometry to find transition x/c (for user output)
                let str = if is == 0 { xf.sst - xf.xt } else { xf.sst + xf.xt };
                let chx = xf.xte - xf.xle;
                let chy = xf.yte - xf.yle;
                let chsq = chx * chx + chy * chy;
                let xtr = seval(str, &xf.x[..xf.n], &xf.xp[..xf.n], &xf.s[..xf.n]);
                let ytr = seval(str, &xf.y[..xf.n], &xf.yp[..xf.n], &xf.s[..xf.n]);
                xf.xoctr[is] = ((xtr - xf.xle) * chx + (ytr - xf.yle) * chy) / chsq;
                xf.yoctr[is] = ((ytr - xf.yle) * chx - (xtr - xf.xle) * chy) / chsq;
            }

            xf.tran = false;

            if ibl == xf.iblte[is] as usize {
                // set "2" variables at TE to wake correlations for next station
                xf.turb = true;
                xf.wake = true;
                blvar(xf, 3);
                blmid(xf, 3);
            }

            // Roll u2_m/d2_m (the just-computed "2" station) into the "1"
            // station slot by swapping buffers instead of copying all nsys
            // entries.  u2_m is fully overwritten at the top of the next
            // station (the js/jbl rebuild loop above), so whatever stale
            // values it picks up here are discarded before use; u1_m is either
            // read as-is (vm assembly) or fully redefined at the TE, so the
            // swap always presents the correct previous-station values.
            std::mem::swap(&mut u1_m, &mut u2_m);
            std::mem::swap(&mut d1_m, &mut d2_m);

            u1_a = u2_a;
            d1_a = d2_a;

            due1 = due2;
            dds1 = dds2;

            if ibl == xf.itran[is] as usize && xf.com2.x > xf.com1.x {
                if is == 0 {
                    xf.tindex[is] = (xf.ist as f64 - xf.itran[is] as f64 + 3.0) - (xf.xt - xf.com1.x) / (xf.com2.x - xf.com1.x);
                } else {
                    xf.tindex[is] = (xf.ist as f64 + xf.itran[is] as f64 - 2.0) + (xf.xt - xf.com1.x) / (xf.com2.x - xf.com1.x);
                }
            }

            // set BL variables for next station
            xf.com1 = xf.com2;

            // next streamwise station
        }

        if xf.tforce[is] {
            if xf.show_output {
                eprintln!("Side {} forced transition at x/c = {:7.4} {:5}", is + 1, xf.xoctr[is], xf.itran[is]);
            }
        } else {
            if xf.show_output {
                eprintln!("Side {}  free  transition at x/c = {:7.4} {:5}", is + 1, xf.xoctr[is], xf.itran[is]);
            }
        }
    }

    restore_bl_scratch!(xf, usav, ule1_m, ule2_m, ute1_m, ute2_m, u1_m, d1_m, u2_m, d2_m);
    true
}

/// Solves the 4x4 BL Newton system stored in `vs2` (with the 4th row
/// prescribed by the caller), replacing `vsrez` with the solution.  This is
/// the port of `gauss(4, 4, VS2, VSRez, 1)` from the Fortran, which operates
/// on the first four columns of the (4,5) VS2 array.
fn gauss4(vs2: &mut [[f64; 5]; 4], vsrez: &mut [f64; 4]) {
    let mut z = [0.0f64; 16];
    for k in 0..4 {
        for l in 0..4 {
            z[l * 4 + k] = vs2[k][l];
        }
    }
    gauss(&mut z, vsrez, 4, 4, 1);
    for k in 0..4 {
        for l in 0..4 {
            vs2[k][l] = z[l * 4 + k];
        }
    }
}

/// Marches the BLs and wake in direct mode using the UEDG array.  If
/// separation is encountered, a plausible value of Hk extrapolated from
/// upstream is prescribed instead.  Continuous checking of transition onset
/// is performed.
pub fn mrchue(xf: &mut Xfoil) -> bool {
    // shape parameters for separation criteria
    let hlmax = 3.8;
    let htmax = 2.5;

    for is in 0..2 {
        if xf.show_output {
            eprintln!("   side {} ...", is + 1);
        }

        xf.amcrit = xf.acrit[is];

        // set forced transition arc length position
        xifset(xf, is);

        // initialize similarity station with Thwaites' formula
        let mut ami = 0.0;
        let mut cti = 0.03;
        xf.bule = 1.0;
        let xsi = xf.xssi[is][2];
        let uei = xf.uedg[is][2];
        let ucon = uei / xsi.powf(xf.bule);
        let tsq = 0.45 / (ucon * (5.0 * xf.bule + 1.0) * xf.reybl) * xsi.powf(1.0 - xf.bule);
        let mut thi = tsq.sqrt();
        let mut dsi = 2.2 * thi;

        xf.tran = false;
        xf.turb = false;
        xf.itran[is] = xf.iblte[is];

        // march downstream
        for ibl in 2..=xf.nbl[is] as usize {
            let ibm = ibl - 1;

            xf.simi = ibl == 2;
            xf.wake = ibl > xf.iblte[is] as usize;

            // prescribed quantities
            let xsi = xf.xssi[is][ibl];
            let mut uei = xf.uedg[is][ibl];

            let dswaki;
            if xf.wake {
                let iw = ibl - xf.iblte[is] as usize;
                // wgap is 0-based; BL wake station iw (1-based) reads wgap[iw-1].
                dswaki = xf.wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            let mut direct = true;
            let mut htarg = 0.0;
            let mut dmax = 0.0;
            let mut dsw = 0.0;
            let mut hklim = 0.0;
            let mut cte = 0.0;
            let mut tte = 0.0;
            let mut dte = 0.0;

            // Newton iteration loop for current station
            for _itbl in 1..=25 {
                // assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                // check for transition and set appropriate flags and things
                if !xf.simi && !xf.turb {
                    let ok = trchek(xf);
                    if xf.abort_on_nan && !ok {
                        return false;
                    }
                    ami = xf.com2.ampl;

                    if xf.tran {
                        xf.itran[is] = ibl as i32;
                        if cti <= 0.0 {
                            cti = 0.03;
                            xf.com2.s = cti;
                        }
                    } else {
                        xf.itran[is] = ibl as i32 + 2;
                    }
                }

                if ibl == xf.iblte[is] as usize + 1 {
                    tte = xf.thet[0][xf.iblte[0] as usize] + xf.thet[1][xf.iblte[1] as usize];
                    dte = xf.dstr[0][xf.iblte[0] as usize] + xf.dstr[1][xf.iblte[1] as usize] + xf.ante;
                    cte = (xf.ctau[0][xf.iblte[0] as usize] * xf.thet[0][xf.iblte[0] as usize]
                        + xf.ctau[1][xf.iblte[1] as usize] * xf.thet[1][xf.iblte[1] as usize])
                        / tte;
                    tesys(xf, cte, tte, dte);
                } else {
                    blsys(xf);
                }

                if direct {
                    // try direct mode (set dUe = 0 in currently empty 4th line)
                    xf.vs2[3][0] = 0.0;
                    xf.vs2[3][1] = 0.0;
                    xf.vs2[3][2] = 0.0;
                    xf.vs2[3][3] = 1.0;
                    xf.vsrez[3] = 0.0;

                    // solve Newton system for current "2" station
                    gauss4(&mut xf.vs2, &mut xf.vsrez);

                    // determine max changes and underrelax if necessary
                    dmax = (xf.vsrez[1] / thi).abs().max((xf.vsrez[2] / dsi).abs());
                    if ibl < xf.itran[is] as usize {
                        dmax = dmax.max((xf.vsrez[0] / 10.0).abs());
                    }
                    if ibl >= xf.itran[is] as usize {
                        dmax = dmax.max((xf.vsrez[0] / cti).abs());
                    }

                    let mut rlx = 1.0;
                    if dmax > 0.3 {
                        rlx = 0.3 / dmax;
                    }

                    // see if direct mode is not applicable
                    if ibl != xf.iblte[is] as usize + 1 {
                        // calculate resulting kinematic shape parameter Hk
                        let msq = uei * uei * xf.hstinv / (xf.gm1bl * (1.0 - 0.5 * uei * uei * xf.hstinv));
                        let htest = (dsi + rlx * xf.vsrez[2]) / (thi + rlx * xf.vsrez[1]);
                        let (hktest, _dummy, _dummy2) = hkin(htest, msq);

                        // decide whether to do direct or inverse problem based on Hk
                        let hmax = if ibl < xf.itran[is] as usize { hlmax } else { htmax };
                        direct = hktest < hmax;
                    }

                    if direct {
                        // update as usual
                        if ibl >= xf.itran[is] as usize {
                            cti += rlx * xf.vsrez[0];
                        }
                        thi += rlx * xf.vsrez[1];
                        dsi += rlx * xf.vsrez[2];
                    } else {
                        // set prescribed Hk for inverse calculation at the current station
                        if ibl < xf.itran[is] as usize {
                            // laminar case: relatively slow increase in Hk downstream
                            htarg = xf.com1.hk + 0.03 * (xf.com2.x - xf.com1.x) / xf.com1.t;
                        } else if ibl == xf.itran[is] as usize {
                            // transition interval: weighted laminar and turbulent case
                            htarg = xf.com1.hk + (0.03 * (xf.xt - xf.com1.x) - 0.15 * (xf.com2.x - xf.xt)) / xf.com1.t;
                        } else if xf.wake {
                            // turbulent wake case: asymptotic wake behavior with approximate Backward Euler
                            let const0 = 0.03 * (xf.com2.x - xf.com1.x) / xf.com1.t;
                            let mut hk2 = xf.com1.hk;
                            hk2 -= (hk2 + const0 * (hk2 - 1.0).powi(3) - xf.com1.hk) / (1.0 + 3.0 * const0 * (hk2 - 1.0).powi(2));
                            hk2 -= (hk2 + const0 * (hk2 - 1.0).powi(3) - xf.com1.hk) / (1.0 + 3.0 * const0 * (hk2 - 1.0).powi(2));
                            hk2 -= (hk2 + const0 * (hk2 - 1.0).powi(3) - xf.com1.hk) / (1.0 + 3.0 * const0 * (hk2 - 1.0).powi(2));
                            htarg = hk2;
                        } else {
                            // turbulent case: relatively fast decrease in Hk downstream
                            htarg = xf.com1.hk - 0.15 * (xf.com2.x - xf.com1.x) / xf.com1.t;
                        }

                        // limit specified Hk to something reasonable
                        if xf.wake {
                            htarg = htarg.max(1.01);
                        } else {
                            let hmax = if ibl < xf.itran[is] as usize { hlmax } else { htmax };
                            htarg = htarg.max(hmax);
                        }

                        if xf.show_output {
                            eprintln!(" MRCHUE: Inverse mode at {:4}     Hk = {:8.3}", ibl, htarg);
                        }

                        // try again with prescribed Hk
                        continue;
                    }
                } else {
                    // inverse mode (force Hk to prescribed value HTARG)
                    xf.vs2[3][0] = 0.0;
                    xf.vs2[3][1] = xf.com2.hk_t;
                    xf.vs2[3][2] = xf.com2.hk_d;
                    xf.vs2[3][3] = xf.com2.hk_u;
                    xf.vsrez[3] = htarg - xf.com2.hk;

                    gauss4(&mut xf.vs2, &mut xf.vsrez);

                    // added Ue clamp
                    dmax = (xf.vsrez[1] / thi).abs().max((xf.vsrez[2] / dsi).abs()).max((xf.vsrez[3] / uei).abs());
                    if ibl >= xf.itran[is] as usize {
                        dmax = dmax.max((xf.vsrez[0] / cti).abs());
                    }

                    let mut rlx = 1.0;
                    if dmax > 0.3 {
                        rlx = 0.3 / dmax;
                    }

                    // update variables
                    if ibl >= xf.itran[is] as usize {
                        cti += rlx * xf.vsrez[0];
                    }
                    thi += rlx * xf.vsrez[1];
                    dsi += rlx * xf.vsrez[2];
                    uei += rlx * xf.vsrez[3];
                }

                // eliminate absurd transients
                if ibl >= xf.itran[is] as usize {
                    cti = cti.min(0.30);
                    cti = cti.max(0.0000001);
                }

                if ibl <= xf.iblte[is] as usize {
                    hklim = 1.02;
                } else {
                    hklim = 1.00005;
                }
                let msq = uei * uei * xf.hstinv / (xf.gm1bl * (1.0 - 0.5 * uei * uei * xf.hstinv));
                dsw = dsi - dswaki;
                dslim(&mut dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if dmax <= 1.0E-5 {
                    break;
                }
            }

            if dmax > 1.0E-5 {
                if xf.show_output {
                    eprintln!(" MRCHUE: Convergence failed at {:4}  side{}    Res = {:e}", ibl, is + 1, dmax);
                }

                // the current unconverged solution might still be reasonable...
                if dmax > 0.1 {
                    // the current solution is garbage --> extrapolate values instead
                    if ibl > 3 {
                        if ibl <= xf.iblte[is] as usize {
                            thi = xf.thet[is][ibm] * (xf.xssi[is][ibl] / xf.xssi[is][ibm]).sqrt();
                            dsi = xf.dstr[is][ibm] * (xf.xssi[is][ibl] / xf.xssi[is][ibm]).sqrt();
                        } else if ibl == xf.iblte[is] as usize + 1 {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                        } else {
                            thi = xf.thet[is][ibm];
                            let ratlen = (xf.xssi[is][ibl] - xf.xssi[is][ibm]) / (10.0 * xf.dstr[is][ibm]);
                            dsi = (xf.dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                        }
                        if ibl == xf.itran[is] as usize {
                            cti = 0.05;
                        }
                        if ibl > xf.itran[is] as usize {
                            cti = xf.ctau[is][ibm];
                        }

                        uei = xf.uedg[is][ibl];
                        if ibl > 2 && ibl < xf.nbl[is] as usize {
                            uei = 0.5 * (xf.uedg[is][ibl - 1] + xf.uedg[is][ibl + 1]);
                        }
                    }
                }

                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                // check for transition and set appropriate flags and things
                if !xf.simi && !xf.turb {
                    let ok = trchek(xf);
                    if xf.abort_on_nan && !ok {
                        return false;
                    }
                    ami = xf.com2.ampl;
                    if xf.tran {
                        xf.itran[is] = ibl as i32;
                    }
                    if !xf.tran {
                        xf.itran[is] = ibl as i32 + 2;
                    }
                }

                // set all other extrapolated values for current station
                if ibl < xf.itran[is] as usize {
                    blvar(xf, 1);
                }
                if ibl >= xf.itran[is] as usize {
                    blvar(xf, 2);
                }
                if xf.wake {
                    blvar(xf, 3);
                }

                if ibl < xf.itran[is] as usize {
                    blmid(xf, 1);
                }
                if ibl >= xf.itran[is] as usize {
                    blmid(xf, 2);
                }
                if xf.wake {
                    blmid(xf, 3);
                }
            }

            // store primary variables
            if ibl < xf.itran[is] as usize {
                xf.ctau[is][ibl] = ami;
            }
            if ibl >= xf.itran[is] as usize {
                xf.ctau[is][ibl] = cti;
            }
            xf.thet[is][ibl] = thi;
            xf.dstr[is][ibl] = dsi;
            xf.uedg[is][ibl] = uei;
            xf.mass[is][ibl] = dsi * uei;
            xf.tau[is][ibl] = 0.5 * xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.cf;
            xf.dis[is][ibl] = xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.u * xf.com2.di * xf.com2.hs * 0.5;
            xf.ctq[is][ibl] = xf.com2.cq;
            xf.delt[is][ibl] = xf.com2.de;
            xf.tstr[is][ibl] = xf.com2.hs * xf.com2.t;

            // set "1" variables to "2" variables for next streamwise station
            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);
            xf.com1 = xf.com2;

            // turbulent intervals will follow transition interval or TE
            if xf.tran || ibl == xf.iblte[is] as usize {
                xf.turb = true;

                // save transition location
                xf.tforce[is] = xf.trforc;
                xf.xssitr[is] = xf.xt;
            }

            xf.tran = false;

            if ibl == xf.iblte[is] as usize {
                // (TE wake starting values are set by the following station's
                //  tesys call via the last cte/tte/dte values)
                let _ = (cte, tte, dte);
            }
        }
    }

    true
}

/// Marches the BLs and wake in mixed mode using the current Ue and Hk.
/// The calculated Ue and Hk lie along a line quasi-normal to the natural
/// Ue-Hk characteristic line of the current BL so that the Goldstein or
/// Levy-Lees singularity is never encountered.  Continuous checking of
/// transition onset is performed.
pub fn mrchdu(xf: &mut Xfoil) -> bool {
    const DEPS: f64 = 5.0E-6;

    // constant controlling how far Hk is allowed to deviate from the specified value
    let senswt = 1000.0;

    for is in 0..2 {
        xf.amcrit = xf.acrit[is];

        // set forced transition arc length position
        xifset(xf, is);

        // set leading edge pressure gradient parameter x/u du/dx
        xf.bule = 1.0;

        // old transition station
        let itrold = xf.itran[is];

        xf.tran = false;
        xf.turb = false;
        xf.itran[is] = xf.iblte[is];

        let mut sens = 0.0;
        let mut sennew = 0.0;

        // march downstream
        for ibl in 2..=xf.nbl[is] as usize {
            let ibm = ibl - 1;

            xf.simi = ibl == 2;
            xf.wake = ibl > xf.iblte[is] as usize;

            // initialize current station to existing variables
            let xsi = xf.xssi[is][ibl];
            let mut uei = xf.uedg[is][ibl];
            let mut thi = xf.thet[is][ibl];
            let mut dsi = xf.dstr[is][ibl];

            // fixed BUG MD 7 June 99
            let mut ami;
            let mut cti;
            if ibl < itrold as usize {
                ami = xf.ctau[is][ibl];
                cti = 0.03;
            } else {
                ami = 0.0;
                cti = xf.ctau[is][ibl];
                if cti <= 0.0 {
                    cti = 0.03;
                }
            }

            let dswaki;
            if xf.wake {
                let iw = ibl - xf.iblte[is] as usize;
                // wgap is 0-based; BL wake station iw (1-based) reads wgap[iw-1].
                dswaki = xf.wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            if ibl <= xf.iblte[is] as usize {
                dsi = (dsi - dswaki).max(1.02000 * thi) + dswaki;
            }
            if ibl > xf.iblte[is] as usize {
                dsi = (dsi - dswaki).max(1.00005 * thi) + dswaki;
            }

            let mut dmax = 0.0;
            let mut dsw = 0.0;
            let mut hklim = 0.0;
            let mut ueref = 0.0;
            let mut hkref = 0.0;
            let mut cte = 0.0;
            let mut tte = 0.0;
            let mut dte = 0.0;

            // Newton iteration loop for current station
            for itbl in 1..=25 {
                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                // check for transition and set appropriate flags and things
                if !xf.simi && !xf.turb {
                    let ok = trchek(xf);
                    if xf.abort_on_nan && !ok {
                        return false;
                    }
                    ami = xf.com2.ampl;
                    if xf.tran {
                        xf.itran[is] = ibl as i32;
                    }
                    if !xf.tran {
                        xf.itran[is] = ibl as i32 + 2;
                    }
                }

                if ibl == xf.iblte[is] as usize + 1 {
                    tte = xf.thet[0][xf.iblte[0] as usize] + xf.thet[1][xf.iblte[1] as usize];
                    dte = xf.dstr[0][xf.iblte[0] as usize] + xf.dstr[1][xf.iblte[1] as usize] + xf.ante;
                    cte = (xf.ctau[0][xf.iblte[0] as usize] * xf.thet[0][xf.iblte[0] as usize]
                        + xf.ctau[1][xf.iblte[1] as usize] * xf.thet[1][xf.iblte[1] as usize])
                        / tte;
                    tesys(xf, cte, tte, dte);
                } else {
                    blsys(xf);
                }

                // set stuff at first iteration...
                if itbl == 1 {
                    // set "baseline" Ue and Hk for forming Ue(Hk) relation
                    ueref = xf.com2.u;
                    hkref = xf.com2.hk;

                    // if current point IBL was turbulent and is now laminar, then...
                    if ibl < xf.itran[is] as usize && ibl >= itrold as usize {
                        // extrapolate baseline Hk
                        let uem = xf.uedg[is][ibl - 1];
                        let dsm = xf.dstr[is][ibl - 1];
                        let thm = xf.thet[is][ibl - 1];
                        let msq = uem * uem * xf.hstinv / (xf.gm1bl * (1.0 - 0.5 * uem * uem * xf.hstinv));
                        let (hkr, _dummy, _dummy2) = hkin(dsm / thm, msq);
                        hkref = hkr;
                    }

                    // if current point IBL was laminar, then...
                    if ibl < itrold as usize {
                        // reinitialize or extrapolate Ctau if it's now turbulent
                        if xf.tran {
                            xf.ctau[is][ibl] = 0.03;
                        }
                        if xf.turb {
                            xf.ctau[is][ibl] = xf.ctau[is][ibl - 1];
                        }
                        if xf.tran || xf.turb {
                            cti = xf.ctau[is][ibl];
                            xf.com2.s = cti;
                        }
                    }
                }

                if xf.simi || ibl == xf.iblte[is] as usize + 1 {
                    // for similarity station or first wake point, prescribe Ue
                    xf.vs2[3][0] = 0.0;
                    xf.vs2[3][1] = 0.0;
                    xf.vs2[3][2] = 0.0;
                    xf.vs2[3][3] = xf.com2.u_uei;
                    xf.vsrez[3] = ueref - xf.com2.u;
                } else {
                    // ******** calculate Ue-Hk characteristic slope ********
                    let mut vtmp = [[0.0f64; 5]; 4];
                    let mut vztmp = [0.0f64; 4];
                    for k in 0..4 {
                        vztmp[k] = xf.vsrez[k];
                        for l in 0..5 {
                            vtmp[k][l] = xf.vs2[k][l];
                        }
                    }

                    // set unit dHk
                    vtmp[3][0] = 0.0;
                    vtmp[3][1] = xf.com2.hk_t;
                    vtmp[3][2] = xf.com2.hk_d;
                    vtmp[3][3] = xf.com2.hk_u * xf.com2.u_uei;
                    vztmp[3] = 1.0;

                    // calculate dUe response
                    gauss4(&mut vtmp, &mut vztmp);

                    // set SENSWT * (normalized dUe/dHk)
                    sennew = senswt * vztmp[3] * hkref / ueref;
                    if itbl <= 5 {
                        sens = sennew;
                    } else if itbl <= 15 {
                        sens = 0.5 * (sens + sennew);
                    }

                    // set prescribed Ue-Hk combination
                    xf.vs2[3][0] = 0.0;
                    xf.vs2[3][1] = xf.com2.hk_t * hkref;
                    xf.vs2[3][2] = xf.com2.hk_d * hkref;
                    xf.vs2[3][3] = (xf.com2.hk_u * hkref + sens / ueref) * xf.com2.u_uei;
                    xf.vsrez[3] = -(hkref * hkref) * (xf.com2.hk / hkref - 1.0) - sens * (xf.com2.u / ueref - 1.0);
                }

                // solve Newton system for current "2" station
                gauss4(&mut xf.vs2, &mut xf.vsrez);

                // determine max changes and underrelax if necessary (added Ue clamp)
                dmax = (xf.vsrez[1] / thi).abs().max((xf.vsrez[2] / dsi).abs()).max((xf.vsrez[3] / uei).abs());
                if ibl >= xf.itran[is] as usize {
                    dmax = dmax.max((xf.vsrez[0] / (10.0 * cti)).abs());
                }

                let mut rlx = 1.0;
                if dmax > 0.3 {
                    rlx = 0.3 / dmax;
                }

                // update as usual
                if ibl < xf.itran[is] as usize {
                    ami += rlx * xf.vsrez[0];
                }
                if ibl >= xf.itran[is] as usize {
                    cti += rlx * xf.vsrez[0];
                }
                thi += rlx * xf.vsrez[1];
                dsi += rlx * xf.vsrez[2];
                uei += rlx * xf.vsrez[3];

                // eliminate absurd transients
                if ibl >= xf.itran[is] as usize {
                    cti = cti.min(0.30);
                    cti = cti.max(0.0000001);
                }

                if ibl <= xf.iblte[is] as usize {
                    hklim = 1.02;
                } else {
                    hklim = 1.00005;
                }
                let msq = uei * uei * xf.hstinv / (xf.gm1bl * (1.0 - 0.5 * uei * uei * xf.hstinv));
                dsw = dsi - dswaki;
                dslim(&mut dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if dmax <= DEPS {
                    break;
                }
            }

            if dmax > DEPS {
                if xf.show_output {
                    eprintln!(" MRCHDU: Convergence failed at {:4}  side{}    Res = {:e}", ibl, is + 1, dmax);
                }

                // the current unconverged solution might still be reasonable...
                if dmax > 0.1 {
                    // the current solution is garbage --> extrapolate values instead
                    if ibl > 3 {
                        if ibl <= xf.iblte[is] as usize {
                            thi = xf.thet[is][ibm] * (xf.xssi[is][ibl] / xf.xssi[is][ibm]).sqrt();
                            dsi = xf.dstr[is][ibm] * (xf.xssi[is][ibl] / xf.xssi[is][ibm]).sqrt();
                            uei = xf.uedg[is][ibm];
                        } else if ibl == xf.iblte[is] as usize + 1 {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                            uei = xf.uedg[is][ibm];
                        } else {
                            thi = xf.thet[is][ibm];
                            let ratlen = (xf.xssi[is][ibl] - xf.xssi[is][ibm]) / (10.0 * xf.dstr[is][ibm]);
                            dsi = (xf.dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                            uei = xf.uedg[is][ibm];
                        }
                        if ibl == xf.itran[is] as usize {
                            cti = 0.05;
                        }
                        if ibl > xf.itran[is] as usize {
                            cti = xf.ctau[is][ibm];
                        }
                    }
                }

                blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(xf);

                // check for transition and set appropriate flags and things
                if !xf.simi && !xf.turb {
                    let ok = trchek(xf);
                    if xf.abort_on_nan && !ok {
                        return false;
                    }
                    ami = xf.com2.ampl;
                    if xf.tran {
                        xf.itran[is] = ibl as i32;
                    }
                    if !xf.tran {
                        xf.itran[is] = ibl as i32 + 2;
                    }
                }

                // set all other extrapolated values for current station
                if ibl < xf.itran[is] as usize {
                    blvar(xf, 1);
                }
                if ibl >= xf.itran[is] as usize {
                    blvar(xf, 2);
                }
                if xf.wake {
                    blvar(xf, 3);
                }

                if ibl < xf.itran[is] as usize {
                    blmid(xf, 1);
                }
                if ibl >= xf.itran[is] as usize {
                    blmid(xf, 2);
                }
                if xf.wake {
                    blmid(xf, 3);
                }
            }

            sens = sennew;

            // store primary variables
            if ibl < xf.itran[is] as usize {
                xf.ctau[is][ibl] = ami;
            }
            if ibl >= xf.itran[is] as usize {
                xf.ctau[is][ibl] = cti;
            }
            xf.thet[is][ibl] = thi;
            xf.dstr[is][ibl] = dsi;
            xf.uedg[is][ibl] = uei;
            xf.mass[is][ibl] = dsi * uei;
            xf.tau[is][ibl] = 0.5 * xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.cf;
            xf.dis[is][ibl] = xf.com2.r * xf.com2.u * xf.com2.u * xf.com2.u * xf.com2.di * xf.com2.hs * 0.5;
            xf.ctq[is][ibl] = xf.com2.cq;
            xf.delt[is][ibl] = xf.com2.de;
            xf.tstr[is][ibl] = xf.com2.hs * xf.com2.t;

            // set "1" variables to "2" variables for next streamwise station
            blprv(xf, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(xf);
            xf.com1 = xf.com2;

            // turbulent intervals will follow transition interval or TE
            if xf.tran || ibl == xf.iblte[is] as usize {
                xf.turb = true;

                // save transition location
                xf.tforce[is] = xf.trforc;
                xf.xssitr[is] = xf.xt;
            }

            xf.tran = false;
        }
    }

    true
}

/// Sets forced-transition BL coordinate locations (XIFSET).
pub fn xifset(xf: &mut Xfoil, is: usize) {
    if xf.xstrip[is] >= 1.0 {
        xf.xiforc = xf.xssi[is][xf.iblte[is] as usize];
        return;
    }

    let chx = xf.xte - xf.xle;
    let chy = xf.yte - xf.yle;
    let chsq = chx * chx + chy * chy;

    // calculate chord-based x/c, y/c
    for i in 0..xf.n {
        xf.w1[i] = ((xf.x[i] - xf.xle) * chx + (xf.y[i] - xf.yle) * chy) / chsq;
        xf.w2[i] = ((xf.y[i] - xf.yle) * chx - (xf.x[i] - xf.xle) * chy) / chsq;
    }

    splind(&xf.w1[..xf.n], &mut xf.w3[..xf.n], &xf.s[..xf.n], -999.0, -999.0);
    splind(&xf.w2[..xf.n], &mut xf.w4[..xf.n], &xf.s[..xf.n], -999.0, -999.0);

    if is == 0 {
        // set approximate arc length of forced transition point for SINVRT
        let mut str = xf.sle + (xf.s[0] - xf.sle) * xf.xstrip[is];

        // calculate actual arc length
        sinvrt(&mut str, xf.xstrip[is], &xf.w1[..xf.n], &xf.w3[..xf.n], &xf.s[..xf.n], xf.show_output);

        // set BL coordinate value
        xf.xiforc = (xf.sst - str).min(xf.xssi[is][xf.iblte[is] as usize]);
    } else {
        // same for bottom side
        let mut str = xf.sle + (xf.s[xf.n - 1] - xf.sle) * xf.xstrip[is];
        sinvrt(&mut str, xf.xstrip[is], &xf.w1[..xf.n], &xf.w3[..xf.n], &xf.s[..xf.n], xf.show_output);
        xf.xiforc = (str - xf.sst).min(xf.xssi[is][xf.iblte[is] as usize]);
    }

    if xf.xiforc < 0.0 {
        if xf.show_output {
            eprintln!();
            eprintln!(" ***  Stagnation point is past trip on side {}  ***", is + 1);
        }
        xf.xiforc = xf.xssi[is][xf.iblte[is] as usize];
    }
}

/// Adds on Newton deltas to the boundary layer variables, checking for
/// excessive changes and underrelaxing if necessary, and calculating max
/// and rms changes as well as the change in the global variable "AC"
/// (CL if LALFA, else alpha).
pub fn update(xf: &mut Xfoil) {
    // max allowable alpha changes per iteration
    let dalmax = 0.5 * crate::state::DTOR;
    let dalmin = -0.5 * crate::state::DTOR;

    // max allowable CL change per iteration
    let dclmax = 0.5;
    let mut dclmin = -0.5;
    if xf.matyp != 1 {
        dclmin = (-0.5f64).max(-0.9 * xf.cl);
    }

    let hstinv = xf.gamm1 * (xf.minf / xf.qinf).powi(2) / (1.0 + 0.5 * xf.gamm1 * xf.minf.powi(2));

    // calculate new Ue distribution assuming no under-relaxation, and set the
    // sensitivity of Ue wrt alpha or Re.  unew/u_ac and qnew/q_ac are reused
    // across iterations via mem::take; they are returned to the state below.
    let mut unew = std::mem::take(&mut xf.bl_unew);
    let mut u_ac = std::mem::take(&mut xf.bl_u_ac);
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let i = xf.ipan[is][ibl] as usize;

            let mut dui = 0.0;
            let mut dui_ac = 0.0;
            for js in 0..2 {
                for jbl in 2..=xf.nbl[js] as usize {
                    let j = xf.ipan[js][jbl] as usize;
                    let jv = xf.isys[js][jbl] as usize;
                    let ue_m = -xf.vti[is][ibl] * xf.vti[js][jbl] * xf.dij_t[i * IZX + j];
                    dui += ue_m * (xf.mass[js][jbl] + xf.vdel[Xfoil::v_index(jv, 1, 3)]);
                    dui_ac += ue_m * (-xf.vdel[Xfoil::v_index(jv, 2, 3)]);
                }
            }

            // UINV depends on "AC" only if "AC" is alpha
            let uinv_ac = if xf.lalfa { 0.0 } else { xf.uinv_a[is][ibl] };

            unew[is][ibl] = xf.uinv[is][ibl] + dui;
            u_ac[is][ibl] = uinv_ac + dui_ac;
        }
    }

    // set new Qtan from new Ue with appropriate sign change
    let mut qnew = std::mem::take(&mut xf.bl_qnew);
    let mut q_ac = std::mem::take(&mut xf.bl_q_ac);
    for is in 0..2 {
        for ibl in 2..=xf.iblte[is] as usize {
            let i = xf.ipan[is][ibl] as usize;
            qnew[i] = xf.vti[is][ibl] * unew[is][ibl];
            q_ac[i] = xf.vti[is][ibl] * u_ac[is][ibl];
        }
    }

    // calculate new CL from this new Qtan
    let sa = xf.alfa.sin();
    let ca = xf.alfa.cos();

    let beta = (1.0 - xf.minf.powi(2)).sqrt();
    let beta_msq = -0.5 / beta;

    let bfac = 0.5 * xf.minf.powi(2) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    let mut clnew = 0.0;
    let mut cl_a = 0.0;
    let mut cl_ms = 0.0;
    let mut cl_ac = 0.0;

    let mut cpg1;
    let mut cpg1_ms;
    let mut cpg1_ac;
    {
        let i = 0usize;
        let cginc = 1.0 - (qnew[i] / xf.qinf).powi(2);
        cpg1 = cginc / (beta + bfac * cginc);
        cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        let cpi_q = -2.0 * qnew[i] / xf.qinf.powi(2);
        let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
        cpg1_ac = cpc_cpi * cpi_q * q_ac[i];
    }

    for i in 0..xf.n {
        let ip = if i == xf.n - 1 { 0 } else { i + 1 };

        let cginc = 1.0 - (qnew[ip] / xf.qinf).powi(2);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        let cpi_q = -2.0 * qnew[ip] / xf.qinf.powi(2);
        let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

        let dx = (xf.x[ip] - xf.x[i]) * ca + (xf.y[ip] - xf.y[i]) * sa;
        let dx_a = -(xf.x[ip] - xf.x[i]) * sa + (xf.y[ip] - xf.y[i]) * ca;

        let ag = 0.5 * (cpg2 + cpg1);
        let ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
        let ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

        clnew += dx * ag;
        cl_a += dx_a * ag;
        cl_ms += dx * ag_ms;
        cl_ac += dx * ag_ac;

        cpg1 = cpg2;
        cpg1_ms = cpg2_ms;
        cpg1_ac = cpg2_ac;
    }

    // initialize under-relaxation factor
    let mut rlx = 1.0;

    let dac;
    if xf.lalfa {
        // alpha is prescribed: AC is CL
        // set change in Re to account for CL changing, since Re = Re(CL)
        dac = (clnew - xf.cl) / (1.0 - cl_ac - cl_ms * 2.0 * xf.minf * xf.minf_cl);

        // set under-relaxation factor if Re change is too large
        if rlx * dac > dclmax {
            rlx = dclmax / dac;
        }
        if rlx * dac < dclmin {
            rlx = dclmin / dac;
        }
    } else {
        // CL is prescribed: AC is alpha
        // set change in alpha to drive CL to prescribed value
        dac = (clnew - xf.clspec) / (0.0 - cl_ac - cl_a);

        // set under-relaxation factor if alpha change is too large
        if rlx * dac > dalmax {
            rlx = dalmax / dac;
        }
        if rlx * dac < dalmin {
            rlx = dalmin / dac;
        }
    }

    xf.rmsbl = 0.0;
    xf.rmxbl = 0.0;

    let dhi = 1.5;
    let dlo = -0.5;

    // calculate changes in BL variables and under-relaxation if needed
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let iv = xf.isys[is][ibl] as usize;

            // set changes without underrelaxation
            let dctau = xf.vdel[Xfoil::v_index(iv, 1, 1)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 1)];
            let dthet = xf.vdel[Xfoil::v_index(iv, 1, 2)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 2)];
            let dmass = xf.vdel[Xfoil::v_index(iv, 1, 3)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 3)];
            let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - xf.uedg[is][ibl];
            let ddstr = (dmass - xf.dstr[is][ibl] * duedg) / xf.uedg[is][ibl];

            // normalize changes
            let dn1 = if ibl < xf.itran[is] as usize { dctau / 10.0 } else { dctau / xf.ctau[is][ibl] };
            let dn2 = dthet / xf.thet[is][ibl];
            let dn3 = ddstr / xf.dstr[is][ibl];
            let dn4 = duedg.abs() / 0.25;

            // accumulate for rms change
            xf.rmsbl += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

            // see if Ctau needs underrelaxation
            let rdn1 = rlx * dn1;
            if dn1.abs() > xf.rmxbl.abs() {
                xf.rmxbl = dn1;
                if ibl < xf.itran[is] as usize {
                    xf.vmxbl = 'n';
                }
                if ibl >= xf.itran[is] as usize {
                    xf.vmxbl = 'C';
                }
                xf.imxbl = ibl as i32;
                xf.ismxbl = is as i32 + 1;
            }
            if rdn1 > dhi {
                rlx = dhi / dn1;
            }
            if rdn1 < dlo {
                rlx = dlo / dn1;
            }

            // see if Theta needs underrelaxation
            let rdn2 = rlx * dn2;
            if dn2.abs() > xf.rmxbl.abs() {
                xf.rmxbl = dn2;
                xf.vmxbl = 'T';
                xf.imxbl = ibl as i32;
                xf.ismxbl = is as i32 + 1;
            }
            if rdn2 > dhi {
                rlx = dhi / dn2;
            }
            if rdn2 < dlo {
                rlx = dlo / dn2;
            }

            // see if Dstar needs underrelaxation
            let rdn3 = rlx * dn3;
            if dn3.abs() > xf.rmxbl.abs() {
                xf.rmxbl = dn3;
                xf.vmxbl = 'D';
                xf.imxbl = ibl as i32;
                xf.ismxbl = is as i32 + 1;
            }
            if rdn3 > dhi {
                rlx = dhi / dn3;
            }
            if rdn3 < dlo {
                rlx = dlo / dn3;
            }

            // see if Ue needs underrelaxation
            let rdn4 = rlx * dn4;
            if dn4.abs() > xf.rmxbl.abs() {
                xf.rmxbl = dn4;
                xf.vmxbl = 'U';
                xf.imxbl = ibl as i32;
                xf.ismxbl = is as i32 + 1;
            }
            if rdn4 > dhi {
                rlx = dhi / dn4;
            }
            if rdn4 < dlo {
                rlx = dlo / dn4;
            }
        }
    }

    // set true rms change
    xf.rmsbl = (xf.rmsbl / (4.0 * (xf.nbl[0] + xf.nbl[1]) as f64)).sqrt();

    xf.rlx = rlx;

    if xf.lalfa {
        // set underrelaxed change in Reynolds number from change in lift
        xf.cl += rlx * dac;
    } else {
        // set underrelaxed change in alpha
        xf.alfa += rlx * dac;
        xf.adeg = xf.alfa / crate::state::DTOR;
    }

    // update BL variables with underrelaxed changes
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let iv = xf.isys[is][ibl] as usize;

            let dctau = xf.vdel[Xfoil::v_index(iv, 1, 1)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 1)];
            let dthet = xf.vdel[Xfoil::v_index(iv, 1, 2)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 2)];
            let dmass = xf.vdel[Xfoil::v_index(iv, 1, 3)] - dac * xf.vdel[Xfoil::v_index(iv, 2, 3)];
            let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - xf.uedg[is][ibl];
            let ddstr = (dmass - xf.dstr[is][ibl] * duedg) / xf.uedg[is][ibl];

            xf.ctau[is][ibl] += rlx * dctau;
            xf.thet[is][ibl] += rlx * dthet;
            xf.dstr[is][ibl] += rlx * ddstr;
            xf.uedg[is][ibl] += rlx * duedg;

            let dswaki;
            if ibl > xf.iblte[is] as usize {
                let iw = ibl - xf.iblte[is] as usize;
                // wgap is 0-based; BL wake station iw (1-based) reads wgap[iw-1].
                dswaki = xf.wgap[iw - 1];
            } else {
                dswaki = 0.0;
            }

            // eliminate absurd transients
            if ibl >= xf.itran[is] as usize {
                xf.ctau[is][ibl] = xf.ctau[is][ibl].min(0.25);
            }

            let hklim;
            if ibl <= xf.iblte[is] as usize {
                hklim = 1.02;
            } else {
                hklim = 1.00005;
            }
            let msq = xf.uedg[is][ibl].powi(2) * hstinv / (xf.gamm1 * (1.0 - 0.5 * xf.uedg[is][ibl].powi(2) * hstinv));
            let mut dsw = xf.dstr[is][ibl] - dswaki;
            dslim(&mut dsw, xf.thet[is][ibl], xf.uedg[is][ibl], msq, hklim);
            xf.dstr[is][ibl] = dsw + dswaki;

            // set new mass defect (nonlinear update)
            xf.mass[is][ibl] = xf.dstr[is][ibl] * xf.uedg[is][ibl];
        }

        // make sure there are no "islands" of negative Ue
        for ibl in 3..=xf.iblte[is] as usize {
            if xf.uedg[is][ibl - 1] > 0.0 && xf.uedg[is][ibl] <= 0.0 {
                xf.uedg[is][ibl] = xf.uedg[is][ibl - 1];
                xf.mass[is][ibl] = xf.dstr[is][ibl] * xf.uedg[is][ibl];
            }
        }
    }

    // equate upper wake arrays to lower wake arrays
    for kbl in 1..=(xf.nbl[1] - xf.iblte[1]) as usize {
        xf.ctau[0][xf.iblte[0] as usize + kbl] = xf.ctau[1][xf.iblte[1] as usize + kbl];
        xf.thet[0][xf.iblte[0] as usize + kbl] = xf.thet[1][xf.iblte[1] as usize + kbl];
        xf.dstr[0][xf.iblte[0] as usize + kbl] = xf.dstr[1][xf.iblte[1] as usize + kbl];
        xf.uedg[0][xf.iblte[0] as usize + kbl] = xf.uedg[1][xf.iblte[1] as usize + kbl];
        xf.tau[0][xf.iblte[0] as usize + kbl] = xf.tau[1][xf.iblte[1] as usize + kbl];
        xf.dis[0][xf.iblte[0] as usize + kbl] = xf.dis[1][xf.iblte[1] as usize + kbl];
        xf.ctq[0][xf.iblte[0] as usize + kbl] = xf.ctq[1][xf.iblte[1] as usize + kbl];
        xf.delt[0][xf.iblte[0] as usize + kbl] = xf.delt[1][xf.iblte[1] as usize + kbl];
        xf.tstr[0][xf.iblte[0] as usize + kbl] = xf.tstr[1][xf.iblte[1] as usize + kbl];
    }

    // return the scratch Vecs to the state for reuse on the next iteration
    xf.bl_unew = unew;
    xf.bl_u_ac = u_ac;
    xf.bl_qnew = qnew;
    xf.bl_q_ac = q_ac;
}

/// Limits the displacement thickness to keep Hk >= Hklim.
pub fn dslim(dstr: &mut f64, thet: f64, uedg: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet;
    let (hk, hk_h, _hk_m) = hkin(h, msq);

    let dh = (0.0f64).max(hklim - hk) / hk_h;
    *dstr += dh * thet;
}

/// Sets the boundary-layer parameter calibration constants (BLPINI).
pub fn blpini(xf: &mut Xfoil) {
    xf.sccon = 5.6;
    xf.gacon = 6.70;
    xf.gbcon = 0.75;
    xf.gccon = 18.0;
    xf.dlcon = 0.9;

    xf.ctrcon = 1.8;
    xf.ctrcex = 3.3;

    xf.duxcon = 1.0;

    xf.ctcon = 0.5 / (xf.gacon * xf.gacon * xf.gbcon);

    xf.cffac = 1.0;
}

/// Precomputes the compressibility-free momentum-defect influence
/// coefficient (unused placeholder kept for API parity).
pub fn _qopi() -> f64 {
    QOPI
}

// silence unused import warnings for helpers used by tests
#[allow(unused_imports)]
use NCOM as _NCOM;
