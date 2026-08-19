// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Shared coefficient routines (port of `s_xfoil.f90`).

use crate::state::Xfoil;

/// Sets actual Mach and Reynolds numbers from unit-CL values and specified
/// `cls`, depending on the MATYP/RETYP flags.  Returns (M_cls, R_cls).
pub fn mrcl(xf: &mut Xfoil, cls: f64) -> (f64, f64) {
    let cla = cls.max(0.000001);

    if xf.retyp < 1 || xf.retyp > 3 {
        if xf.show_output {
            eprintln!("MRCL:  Illegal Re(CL) dependence trigger.");
            eprintln!("       Setting fixed Re.");
        }
        xf.retyp = 1;
    }
    if xf.matyp < 1 || xf.matyp > 3 {
        if xf.show_output {
            eprintln!("MRCL:  Illegal Mach(CL) dependence trigger.");
            eprintln!("       Setting fixed Mach.");
        }
        xf.matyp = 1;
    }

    let m_cls;
    if xf.matyp == 1 {
        xf.minf = xf.minf1;
        m_cls = 0.0;
    } else if xf.matyp == 2 {
        xf.minf = xf.minf1 / cla.sqrt();
        m_cls = -0.5 * xf.minf / cla;
    } else {
        xf.minf = xf.minf1;
        m_cls = 0.0;
    }

    let r_cls;
    if xf.retyp == 1 {
        xf.reinf = xf.reinf1;
        r_cls = 0.0;
    } else if xf.retyp == 2 {
        xf.reinf = xf.reinf1 / cla.sqrt();
        r_cls = -0.5 * xf.reinf / cla;
    } else {
        xf.reinf = xf.reinf1 / cla;
        r_cls = -xf.reinf / cla;
    }

    if xf.minf >= 0.99 {
        if xf.show_output {
            eprintln!();
            eprintln!("MRCL: CL too low for chosen Mach(CL) dependence");
            eprintln!("      Aritificially limiting Mach to  0.99");
        }
        xf.minf = 0.99;
    }

    let rrat = if xf.reinf1 > 0.0 { xf.reinf / xf.reinf1 } else { 1.0 };

    if rrat > 100.0 {
        if xf.show_output {
            eprintln!();
            eprintln!("MRCL: CL too low for chosen Re(CL) dependence");
            eprintln!("      Aritificially limiting Re to {:e}", xf.reinf1 * 100.0);
        }
        xf.reinf = xf.reinf1 * 100.0;
    }

    (m_cls, r_cls)
}

/// Sets the Karman-Tsien compressibility parameter TKLAM and the sonic
/// pressure coefficient and speed.
pub fn comset(xf: &mut Xfoil) {
    let beta = (1.0 - xf.minf.powi(2)).sqrt();
    let beta_msq = -0.5 / beta;

    xf.tklam = xf.minf.powi(2) / (1.0 + beta).powi(2);
    xf.tkl_msq = 1.0 / (1.0 + beta).powi(2) - 2.0 * xf.tklam / (1.0 + beta) * beta_msq;

    // set sonic pressure coefficient and speed
    if xf.minf == 0.0 {
        xf.cpstar = -999.0;
        xf.qstar = 999.0;
    } else {
        xf.cpstar = 2.0 / (xf.gamma * xf.minf.powi(2))
            * (((1.0 + 0.5 * xf.gamm1 * xf.minf.powi(2)) / (1.0 + 0.5 * xf.gamm1)).powf(xf.gamma / xf.gamm1) - 1.0);
        xf.qstar = xf.qinf / xf.minf
            * ((1.0 + 0.5 * xf.gamm1 * xf.minf.powi(2)) / (1.0 + 0.5 * xf.gamm1)).sqrt();
    }
}

/// Sets compressible Cp from speed.
pub fn cpcalc(q: &[f64], qinf: f64, minf: f64, cp: &mut [f64], show_output: bool) {
    let beta = (1.0 - minf.powi(2)).sqrt();
    let bfac = 0.5 * minf.powi(2) / (1.0 + beta);

    let mut denneg = false;

    for i in 0..q.len() {
        let cpinc = 1.0 - (q[i] / qinf).powi(2);
        let den = beta + bfac * cpinc;
        cp[i] = cpinc / den;
        if den <= 0.0 {
            denneg = true;
        }
    }

    if denneg {
        if show_output {
            eprintln!();
            eprintln!("CPCALC: Local speed too large. Compressibility corrections invalid.");
        }
    }
}

/// Integrates surface pressures to get CL and CM, integrates skin friction to
/// get CDF, and calculates dCL/dAlpha for prescribed-CL routines.
#[allow(clippy::too_many_arguments)]
pub fn clcalc(
    x: &[f64],
    y: &[f64],
    gam: &[f64],
    gam_a: &[f64],
    alfa: f64,
    minf: f64,
    qinf: f64,
    xref: f64,
    yref: f64,
    cl: &mut f64,
    cm: &mut f64,
    cdp: &mut f64,
    cl_alf: &mut f64,
    cl_msq: &mut f64,
) {
    let n = x.len();

    let sa = alfa.sin();
    let ca = alfa.cos();

    let beta = (1.0 - minf.powi(2)).sqrt();
    let beta_msq = -0.5 / beta;

    let bfac = 0.5 * minf.powi(2) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    *cl = 0.0;
    *cm = 0.0;
    *cdp = 0.0;
    *cl_alf = 0.0;
    *cl_msq = 0.0;

    let mut cpg1;
    let mut cpg1_msq;
    let mut cpg1_alf;
    {
        let i = 0usize;
        let cginc = 1.0 - (gam[i] / qinf).powi(2);
        cpg1 = cginc / (beta + bfac * cginc);
        cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        let cpi_gam = -2.0 * gam[i] / qinf.powi(2);
        let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
        cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];
    }

    for i in 0..n {
        let ip = if i == n - 1 { 0 } else { i + 1 };

        let cginc = 1.0 - (gam[ip] / qinf).powi(2);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        let cpi_gam = -2.0 * gam[ip] / qinf.powi(2);
        let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        let dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        let dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        let dg = cpg2 - cpg1;

        let ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        let ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        let ag = 0.5 * (cpg2 + cpg1);

        let dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
        let ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        let ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

        *cl += dx * ag;
        *cdp -= dy * ag;
        *cm -= dx * (ag * ax + dg * dx / 12.0) + dy * (ag * ay + dg * dy / 12.0);

        *cl_alf += dx * ag_alf + ag * dx_alf;
        *cl_msq += dx * ag_msq;

        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
    }
}

/// Calculates the drag coefficient CD (viscous, via Squire-Young wake
/// extrapolation) and the friction drag coefficient CDF.
pub fn cdcalc(xf: &mut Xfoil) {
    let sa = xf.alfa.sin();
    let ca = xf.alfa.cos();

    if xf.lvisc && xf.lbli_ni {
        // set variables at the end of the wake
        let nbl2 = xf.nbl[1] as usize;
        let thwake = xf.thet[1][nbl2];
        let urat = xf.uedg[1][nbl2] / xf.qinf;
        let uewake = xf.uedg[1][nbl2] * (1.0 - xf.tklam) / (1.0 - xf.tklam * urat * urat);
        let shwake = xf.dstr[1][nbl2] / xf.thet[1][nbl2];

        // extrapolate wake to downstream infinity using Squire-Young relation
        xf.cd = 2.0 * thwake * (uewake / xf.qinf).powf(0.5 * (5.0 + shwake));
    } else {
        xf.cd = 0.0;
    }

    // calculate friction drag coefficient
    xf.cdf = 0.0;
    for is in 0..2 {
        for ibl in 3..=xf.iblte[is] as usize {
            let i = xf.ipan[is][ibl] as usize;
            let im = xf.ipan[is][ibl - 1] as usize;
            let dx = (xf.x[i] - xf.x[im]) * ca + (xf.y[i] - xf.y[im]) * sa;
            xf.cdf += 0.5 * (xf.tau[is][ibl] + xf.tau[is][ibl - 1]) * dx * 2.0 / xf.qinf.powi(2);
        }
    }
}

/// Calculates total and projected TE gap areas and TE panel strengths.
pub fn tecalc(xf: &mut Xfoil) {
    let n = xf.n;

    // set TE base vector and TE bisector components
    let dxte = xf.x[0] - xf.x[n - 1];
    let dyte = xf.y[0] - xf.y[n - 1];
    let dxs = 0.5 * (-xf.xp[0] + xf.xp[n - 1]);
    let dys = 0.5 * (-xf.yp[0] + xf.yp[n - 1]);

    // normal and streamwise projected TE gap areas
    xf.ante = dxs * dyte - dys * dxte;
    xf.aste = dxs * dxte + dys * dyte;

    // total TE gap area
    xf.dste = (dxte * dxte + dyte * dyte).sqrt();

    xf.sharp = xf.dste < 0.0001 * xf.chord;

    let (scs, sds);
    if xf.sharp {
        scs = 1.0;
        sds = 0.0;
    } else {
        scs = xf.ante / xf.dste;
        sds = xf.aste / xf.dste;
    }

    // TE panel source and vorticity strengths
    xf.sigte = 0.5 * (xf.gam[0] - xf.gam[n - 1]) * scs;
    xf.gamte = -0.5 * (xf.gam[0] - xf.gam[n - 1]) * sds;

    xf.sigte_a = 0.5 * (xf.gam_a[0] - xf.gam_a[n - 1]) * scs;
    xf.gamte_a = -0.5 * (xf.gam_a[0] - xf.gam_a[n - 1]) * sds;
}
