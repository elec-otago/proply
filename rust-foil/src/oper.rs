// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Operating-point routines (port of `m_xoper.f90`): specal, speccl,
//! viscal, fcpmin.  The interactive `.OPER` menu is not ported.

use crate::bl::{setbl, update};
use crate::panel::{
    gamqv, iblpan, qdcalc, qiset, qvfue, qwcalc, stfind, stmove, uicalc, xicalc, xywake,
};
use crate::s_xbl::iblsys;
use crate::s_xfoil::{cdcalc, clcalc, comset, cpcalc, mrcl, tecalc};
use crate::solve::blsolv;
use crate::state::{Xfoil, DTOR};

/// Finds the minimum Cp on the distribution (for cavitation work).
pub fn fcpmin(xf: &mut Xfoil) {
    xf.xcpmni = xf.x[0];
    xf.xcpmnv = xf.x[0];
    xf.cpmni = xf.cpi[0];
    xf.cpmnv = xf.cpv[0];

    for i in 1..xf.n + xf.nw {
        if xf.cpi[i] < xf.cpmni {
            xf.xcpmni = xf.x[i];
            xf.cpmni = xf.cpi[i];
        }
        if xf.cpv[i] < xf.cpmnv {
            xf.xcpmnv = xf.x[i];
            xf.cpmnv = xf.cpv[i];
        }
    }

    if xf.lvisc {
        xf.cpmn = xf.cpmnv;
    } else {
        xf.cpmn = xf.cpmni;
        xf.cpmnv = xf.cpmni;
        xf.xcpmnv = xf.xcpmni;
    }
}

/// Converges to the specified alpha.
pub fn specal(xf: &mut Xfoil) {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if !xf.lgamu || !xf.lqaij {
        ggcalc_wrap(xf);
    }

    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();

    // superimpose suitably weighted alpha = 0, 90 distributions
    for i in 0..xf.n {
        xf.gam[i] = xf.cosa * xf.gamu[0][i] + xf.sina * xf.gamu[1][i];
        xf.gam_a[i] = -xf.sina * xf.gamu[0][i] + xf.cosa * xf.gamu[1][i];
    }
    xf.psio = xf.cosa * xf.gamu[0][xf.n] + xf.sina * xf.gamu[1][xf.n];

    tecalc(xf);
    qiset(xf);

    // set initial guess for the Newton variable CLM
    let mut clm = 1.0;

    // set corresponding M(CLM), Re(CLM)
    let (mut minf_clm, mut reinf_clm) = mrcl(xf, clm);
    comset(xf);

    // set corresponding CL(M)
    clcalc(
        &xf.x[..xf.n],
        &xf.y[..xf.n],
        &xf.gam[..xf.n],
        &xf.gam_a[..xf.n],
        xf.alfa,
        xf.minf,
        xf.qinf,
        xf.xcmref,
        xf.ycmref,
        &mut xf.cl,
        &mut xf.cm,
        &mut xf.cdp,
        &mut xf.cl_alf,
        &mut xf.cl_msq,
    );

    // iterate on CLM
    let mut converged = false;
    for _itcl in 1..=20 {
        let msq_clm = 2.0 * xf.minf * minf_clm;
        let dclm = (xf.cl - clm) / (1.0 - xf.cl_msq * msq_clm);

        let clm1 = clm;
        let mut rlx = 1.0;

        // under-relaxation loop to avoid driving M(CL) above 1
        for _irlx in 1..=12 {
            clm = clm1 + rlx * dclm;

            // set new freestream Mach M(CLM)
            let r = mrcl(xf, clm);
            minf_clm = r.0;
            reinf_clm = r.1;

            // if Mach is OK, go do next Newton iteration
            if xf.matyp == 1 || xf.minf == 0.0 || minf_clm != 0.0 {
                break;
            }

            rlx *= 0.5;
        }

        // set new CL(M)
        comset(xf);
        clcalc(
            &xf.x[..xf.n],
            &xf.y[..xf.n],
            &xf.gam[..xf.n],
            &xf.gam_a[..xf.n],
            xf.alfa,
            xf.minf,
            xf.qinf,
            xf.xcmref,
            xf.ycmref,
            &mut xf.cl,
            &mut xf.cm,
            &mut xf.cdp,
            &mut xf.cl_alf,
            &mut xf.cl_msq,
        );

        if dclm.abs() <= 1.0E-6 {
            converged = true;
            break;
        }
    }
    if !converged && xf.show_output {
        eprintln!("SPECAL:  Minf convergence failed");
    }

    // set final Mach, CL, Cp distributions
    mrcl(xf, xf.cl);
    comset(xf);
    clcalc(
        &xf.x[..xf.n],
        &xf.y[..xf.n],
        &xf.gam[..xf.n],
        &xf.gam_a[..xf.n],
        xf.alfa,
        xf.minf,
        xf.qinf,
        xf.xcmref,
        xf.ycmref,
        &mut xf.cl,
        &mut xf.cm,
        &mut xf.cdp,
        &mut xf.cl_alf,
        &mut xf.cl_msq,
    );

    if xf.lvisc {
        cpcalc(
            &xf.qvis[..xf.n + xf.nw],
            xf.qinf,
            xf.minf,
            &mut xf.cpv[..xf.n + xf.nw],
            xf.show_output,
        );
        cpcalc(
            &xf.qinv[..xf.n + xf.nw],
            xf.qinf,
            xf.minf,
            &mut xf.cpi[..xf.n + xf.nw],
            xf.show_output,
        );
    } else {
        cpcalc(
            &xf.qinv[..xf.n],
            xf.qinf,
            xf.minf,
            &mut xf.cpi[..xf.n],
            xf.show_output,
        );
    }
}

fn ggcalc_wrap(xf: &mut Xfoil) {
    crate::panel::ggcalc(xf);
}

/// Converges to the specified inviscid CL.
pub fn speccl(xf: &mut Xfoil) {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if !xf.lgamu || !xf.lqaij {
        ggcalc_wrap(xf);
    }

    // set freestream Mach from specified CL -- Mach will be held fixed
    mrcl(xf, xf.clspec);
    comset(xf);

    // current alpha is the initial guess for Newton variable ALFA
    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();
    for i in 0..xf.n {
        xf.gam[i] = xf.cosa * xf.gamu[0][i] + xf.sina * xf.gamu[1][i];
        xf.gam_a[i] = -xf.sina * xf.gamu[0][i] + xf.cosa * xf.gamu[1][i];
    }
    xf.psio = xf.cosa * xf.gamu[0][xf.n] + xf.sina * xf.gamu[1][xf.n];

    // get corresponding CL, CL_alpha, CL_Mach
    clcalc(
        &xf.x[..xf.n],
        &xf.y[..xf.n],
        &xf.gam[..xf.n],
        &xf.gam_a[..xf.n],
        xf.alfa,
        xf.minf,
        xf.qinf,
        xf.xcmref,
        xf.ycmref,
        &mut xf.cl,
        &mut xf.cm,
        &mut xf.cdp,
        &mut xf.cl_alf,
        &mut xf.cl_msq,
    );

    // Newton loop for alpha to get specified inviscid CL
    let mut converged = false;
    for _ital in 1..=20 {
        let dalfa = (xf.clspec - xf.cl) / xf.cl_alf;
        let mut rlx = 1.0;

        xf.alfa += rlx * dalfa;

        // set new surface speed distribution
        xf.cosa = xf.alfa.cos();
        xf.sina = xf.alfa.sin();
        for i in 0..xf.n {
            xf.gam[i] = xf.cosa * xf.gamu[0][i] + xf.sina * xf.gamu[1][i];
            xf.gam_a[i] = -xf.sina * xf.gamu[0][i] + xf.cosa * xf.gamu[1][i];
        }
        xf.psio = xf.cosa * xf.gamu[0][xf.n] + xf.sina * xf.gamu[1][xf.n];

        // set new CL(alpha)
        clcalc(
            &xf.x[..xf.n],
            &xf.y[..xf.n],
            &xf.gam[..xf.n],
            &xf.gam_a[..xf.n],
            xf.alfa,
            xf.minf,
            xf.qinf,
            xf.xcmref,
            xf.ycmref,
            &mut xf.cl,
            &mut xf.cm,
            &mut xf.cdp,
            &mut xf.cl_alf,
            &mut xf.cl_msq,
        );

        if dalfa.abs() <= 1.0E-6 {
            converged = true;
            break;
        }
    }
    if !converged && xf.show_output {
        eprintln!("SPECCL:  CL convergence failed");
    }

    // set final surface speed and Cp distributions
    tecalc(xf);
    qiset(xf);

    if xf.lvisc {
        cpcalc(
            &xf.qvis[..xf.n + xf.nw],
            xf.qinf,
            xf.minf,
            &mut xf.cpv[..xf.n + xf.nw],
            xf.show_output,
        );
        cpcalc(
            &xf.qinv[..xf.n + xf.nw],
            xf.qinf,
            xf.minf,
            &mut xf.cpi[..xf.n + xf.nw],
            xf.show_output,
        );
    } else {
        cpcalc(
            &xf.qinv[..xf.n],
            xf.qinf,
            xf.minf,
            &mut xf.cpi[..xf.n],
            xf.show_output,
        );
    }
}

/// Converges the viscous operating point.  Returns true if converged.
pub fn viscal(xf: &mut Xfoil, niter1: i32) -> bool {
    // convergence tolerance
    let eps1 = 1.0E-4;

    let niter = niter1;

    // calculate wake trajectory from current inviscid solution if necessary
    if !xf.lwake {
        xywake(xf);
    }

    // set velocities on wake from airfoil vorticity for alpha=0, 90
    qwcalc(xf);

    // set velocities on airfoil and wake for initial alpha
    qiset(xf);

    if !xf.lipan {
        if xf.lbli_ni {
            gamqv(xf);
        }

        // locate stagnation point arc length position and panel index
        stfind(xf);

        // set BL position -> panel position pointers
        iblpan(xf);

        // calculate surface arc length array for current stagnation point location
        xicalc(xf);

        // set BL position -> system line pointers
        iblsys(xf);
    }

    // set inviscid BL edge velocity UINV from QINV
    uicalc(xf);

    if !xf.lbli_ni {
        // set initial Ue from inviscid Ue
        for ibl in 1..=xf.nbl[0] as usize {
            xf.uedg[0][ibl] = xf.uinv[0][ibl];
        }
        for ibl in 1..=xf.nbl[1] as usize {
            xf.uedg[1][ibl] = xf.uinv[1][ibl];
        }
    }

    if xf.lvconv {
        // set correct CL if converged point exists
        qvfue(xf);
        if xf.lvisc {
            cpcalc(
                &xf.qvis[..xf.n + xf.nw],
                xf.qinf,
                xf.minf,
                &mut xf.cpv[..xf.n + xf.nw],
                xf.show_output,
            );
            cpcalc(
                &xf.qinv[..xf.n + xf.nw],
                xf.qinf,
                xf.minf,
                &mut xf.cpi[..xf.n + xf.nw],
                xf.show_output,
            );
        } else {
            cpcalc(
                &xf.qinv[..xf.n],
                xf.qinf,
                xf.minf,
                &mut xf.cpi[..xf.n],
                xf.show_output,
            );
        }
        gamqv(xf);
        clcalc(
            &xf.x[..xf.n],
            &xf.y[..xf.n],
            &xf.gam[..xf.n],
            &xf.gam_a[..xf.n],
            xf.alfa,
            xf.minf,
            xf.qinf,
            xf.xcmref,
            xf.ycmref,
            &mut xf.cl,
            &mut xf.cm,
            &mut xf.cdp,
            &mut xf.cl_alf,
            &mut xf.cl_msq,
        );
        cdcalc(xf);
    }

    // set up source influence matrix if it doesn't exist
    if !xf.lwdij || !xf.ladij {
        qdcalc(xf);
    }

    // Newton iteration for entire BL solution
    if xf.show_output {
        eprintln!();
        eprintln!("Solving BL system ...");
    }

    let mut converged = false;
    for iter in 1..=niter {
        // fill Newton system for BL variables
        let ok = setbl(xf);
        if xf.abort_on_nan && !ok {
            return false;
        }

        // solve Newton system with custom solver
        blsolv(xf);

        // update BL variables
        update(xf);

        if xf.lalfa {
            // set new freestream Mach, Re from new CL
            mrcl(xf, xf.cl);
            comset(xf);
        } else {
            // set new inviscid speeds QINV and UINV for new alpha
            qiset(xf);
            uicalc(xf);
        }

        // calculate edge velocities QVIS(.) from UEDG(..)
        qvfue(xf);

        // set GAM distribution from QVIS
        gamqv(xf);

        // relocate stagnation point
        stmove(xf);

        // set updated CL,CD
        clcalc(
            &xf.x[..xf.n],
            &xf.y[..xf.n],
            &xf.gam[..xf.n],
            &xf.gam_a[..xf.n],
            xf.alfa,
            xf.minf,
            xf.qinf,
            xf.xcmref,
            xf.ycmref,
            &mut xf.cl,
            &mut xf.cm,
            &mut xf.cdp,
            &mut xf.cl_alf,
            &mut xf.cl_msq,
        );
        cdcalc(xf);

        // display changes and test for convergence
        if xf.show_output {
            let cdpdif = xf.cd - xf.cdf;
            if xf.rlx < 1.0 {
                eprintln!();
                eprintln!(
                    "{:3}   rms: {:e}   max: {:e}  {:} at {:4}{:3}   RLX:{:6.3}",
                    iter, xf.rmsbl, xf.rmxbl, xf.vmxbl, xf.imxbl, xf.ismxbl, xf.rlx
                );
            } else {
                eprintln!();
                eprintln!(
                    "{:3}   rms: {:e}   max: {:e}  {:} at {:4}{:3}",
                    iter, xf.rmsbl, xf.rmxbl, xf.vmxbl, xf.imxbl, xf.ismxbl
                );
            }
            eprintln!(
                "         a = {:7.3}      CL = {:8.4}",
                xf.alfa / DTOR,
                xf.cl
            );
            eprintln!(
                "  Cm = {:8.4}     CD = {:9.5}   =>   CDf = {:9.5}    CDp = {:9.5}",
                xf.cm, xf.cd, xf.cdf, cdpdif
            );
        }

        if xf.rmsbl < eps1 {
            xf.lvconv = true;
            xf.avisc = xf.alfa;
            xf.mvisc = xf.minf;
            converged = true;
            break;
        }
    }
    if !converged && xf.show_output {
        eprintln!("VISCAL:  Convergence failed");
    }

    let _ = &mut xf.cpmn;
    cpcalc(
        &xf.qinv[..xf.n + xf.nw],
        xf.qinf,
        xf.minf,
        &mut xf.cpi[..xf.n + xf.nw],
        xf.show_output,
    );
    cpcalc(
        &xf.qvis[..xf.n + xf.nw],
        xf.qinf,
        xf.minf,
        &mut xf.cpv[..xf.n + xf.nw],
        xf.show_output,
    );

    // BL summary output (upper surface only)
    let is = 0usize;
    let mut hkmax = 0.0;
    let mut hkm = 0.0;
    let mut psep = 0.0;
    let mut patt = 0.0;
    for ibl in 2..=xf.iblte[is] as usize {
        let hki = xf.dstr[is][ibl] / xf.thet[is][ibl];
        hkmax = hki.max(hkmax);
        if hkm < 4.0 && hki >= 4.0 {
            let hfrac = (4.0 - hkm) / (hki - hkm);
            let pdefm = xf.uedg[is][ibl - 1].powi(2) * xf.thet[is][ibl - 1];
            let pdefi = xf.uedg[is][ibl].powi(2) * xf.thet[is][ibl];
            psep = pdefm * (1.0 - hfrac) + pdefi * hfrac;
        }
        if hkm > 4.0 && hki < 4.0 {
            let hfrac = (4.0 - hkm) / (hki - hkm);
            let pdefm = xf.uedg[is][ibl - 1].powi(2) * xf.thet[is][ibl - 1];
            let pdefi = xf.uedg[is][ibl].powi(2) * xf.thet[is][ibl];
            patt = pdefm * (1.0 - hfrac) + pdefi * hfrac;
        }
        hkm = hki;
    }
    let delp = patt - psep;

    if xf.show_output {
        eprintln!();
        eprintln!(
            "{:10.3}{:10.4}{:11.6}{:11.6}{:11.6}{:11.6}{:10.4}     #",
            xf.acrit[is],
            hkmax,
            xf.cd,
            2.0 * psep,
            2.0 * patt,
            2.0 * delp,
            xf.xoctr[is]
        );
    }

    converged
}
