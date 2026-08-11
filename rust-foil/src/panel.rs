//! Panel method (port of `m_xpanel.f90`).
//!
//! Computes the inviscid streamfunction/vorticity solution and the
//! influence-coefficient machinery used by the viscous coupling.

use crate::solve::{baksub, ludcmp};
use crate::spline::segspl;
use crate::state::{IQX, IWX, Xfoil};
use crate::utils::atanc;

/// Sets angles of airfoil panels.
pub fn apcalc(xf: &mut Xfoil) {
    for i in 0..xf.n - 1 {
        let sx = xf.x[i + 1] - xf.x[i];
        let sy = xf.y[i + 1] - xf.y[i];
        if sx == 0.0 && sy == 0.0 {
            xf.apanel[i] = (-xf.ny[i]).atan2(-xf.nx[i]);
        } else {
            xf.apanel[i] = sx.atan2(-sy);
        }
    }

    // TE panel
    let i = xf.n - 1;
    let ip = 0;
    if xf.sharp {
        xf.apanel[i] = std::f64::consts::PI;
    } else {
        let sx = xf.x[ip] - xf.x[i];
        let sy = xf.y[ip] - xf.y[i];
        xf.apanel[i] = (-sx).atan2(sy) + std::f64::consts::PI;
    }
}

/// Calculates normal unit vector components at airfoil panel nodes.
pub fn ncalc(x: &[f64], y: &[f64], s: &[f64], xn: &mut [f64], yn: &mut [f64]) {
    let n = x.len();
    if n <= 1 {
        return;
    }

    segspl(x, xn, s);
    segspl(y, yn, s);
    for i in 0..n {
        let sx = yn[i];
        let sy = -xn[i];
        let smod = (sx * sx + sy * sy).sqrt();
        if smod == 0.0 {
            xn[i] = -1.0;
            yn[i] = 0.0;
        } else {
            xn[i] = sx / smod;
            yn[i] = sy / smod;
        }
    }

    // average normal vectors at corner points
    for i in 0..n - 1 {
        if s[i] == s[i + 1] {
            let sx = 0.5 * (xn[i] + xn[i + 1]);
            let sy = 0.5 * (yn[i] + yn[i + 1]);
            let smod = (sx * sx + sy * sy).sqrt();
            if smod == 0.0 {
                xn[i] = -1.0;
                yn[i] = 0.0;
                xn[i + 1] = -1.0;
                yn[i + 1] = 0.0;
            } else {
                xn[i] = sx / smod;
                yn[i] = sy / smod;
                xn[i + 1] = sx / smod;
                yn[i + 1] = sy / smod;
            }
        }
    }
}

/// Calculates the current streamfunction Psi at panel node or wake node I due
/// to freestream and all bound vorticity Gam on the airfoil, plus the
/// sensitivity vectors dPsi/dGam (DZDG), dPsi/dn (DZDN), dQtan/dGam (DQDG)
/// and, if `siglin`, the viscous source distribution effects (DZDM, DQDM).
///
/// `i` is 0-based: airfoil 0..N-1, wake N..N+NW-1.
#[allow(clippy::too_many_arguments)]
pub fn psilin(xf: &mut Xfoil, i: usize, xi: f64, yi: f64, nxi: f64, nyi: f64, psi: &mut f64, psi_ni: &mut f64, geolin: bool, siglin: bool) {
    let n = xf.n;
    let nw = xf.nw;

    // distance tolerance for determining if two points are the same
    let seps = (xf.s[n - 1] - xf.s[0]) * 1.0E-5;

    let io = i;

    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();

    for jo in 0..n {
        xf.dzdg[jo] = 0.0;
        xf.dzdn[jo] = 0.0;
        xf.dqdg[jo] = 0.0;
    }
    for jo in 0..n {
        xf.dzdm[jo] = 0.0;
        xf.dqdm[jo] = 0.0;
    }

    xf.z_qinf = 0.0;
    xf.z_alfa = 0.0;
    xf.z_qdof0 = 0.0;
    xf.z_qdof1 = 0.0;
    xf.z_qdof2 = 0.0;
    xf.z_qdof3 = 0.0;

    *psi = 0.0;
    *psi_ni = 0.0;

    xf.qtan1 = 0.0;
    xf.qtan2 = 0.0;
    let mut qtanm = 0.0;

    let (scs, sds);
    if xf.sharp {
        scs = 1.0;
        sds = 0.0;
    } else {
        scs = xf.ante / xf.dste;
        sds = xf.aste / xf.dste;
    }

    // variables carried across the panel loop (needed for TE panel terms)
    let mut g1 = 0.0;
    let mut g2 = 0.0;
    let mut t1 = 0.0;
    let mut t2 = 0.0;
    let mut apan = 0.0;
    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut yy = 0.0;
    let mut x1i = 0.0;
    let mut x2i = 0.0;
    let mut yyi = 0.0;
    let mut x1o = 0.0;
    let mut x1p = 0.0;
    let mut x2o = 0.0;
    let mut x2p = 0.0;
    let mut yyo = 0.0;
    let mut yyp = 0.0;
    let mut jo_last = 0usize;
    let mut jp_last = 1usize;
    let mut te_skip = false;

    for jo in 0..n {
        let mut jp = jo + 1;
        let jm;
        let jq;

        if jo == 0 {
            jm = 0;
            jq = 2;
        } else if jo == n - 2 {
            jm = jo - 1;
            jq = jp;
        } else if jo == n - 1 {
            jp = 0;
            jm = jo - 1;
            jq = 1;
            if (xf.x[jo] - xf.x[jp]).powi(2) + (xf.y[jo] - xf.y[jp]).powi(2) < seps * seps {
                te_skip = true;
                break;
            }
        } else {
            jm = jo - 1;
            jq = jp + 1;
        }

        let dso = ((xf.x[jo] - xf.x[jp]).powi(2) + (xf.y[jo] - xf.y[jp]).powi(2)).sqrt();

        // skip null panel
        if dso != 0.0 {
            let dsio = 1.0 / dso;

            apan = xf.apanel[jo];

            let rx1 = xi - xf.x[jo];
            let ry1 = yi - xf.y[jo];
            let rx2 = xi - xf.x[jp];
            let ry2 = yi - xf.y[jp];

            let sx = (xf.x[jp] - xf.x[jo]) * dsio;
            let sy = (xf.y[jp] - xf.y[jo]) * dsio;

            x1 = sx * rx1 + sy * ry1;
            x2 = sx * rx2 + sy * ry2;
            yy = sx * ry1 - sy * rx1;

            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;

            // set reflection flag SGN to avoid branch problems with arctan
            let sgn;
            if io < n {
                // no problem on airfoil surface
                sgn = 1.0;
            } else {
                // make sure arctan falls between -/+ Pi/2
                sgn = yy.signum();
            }

            // set log(r^2) and arctan(x/y), correcting for reflection if any
            if io != jo && rs1 > 0.0 {
                g1 = rs1.ln();
                t1 = (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * std::f64::consts::PI;
            } else {
                g1 = 0.0;
                t1 = 0.0;
            }

            if io != jp && rs2 > 0.0 {
                g2 = rs2.ln();
                t2 = (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * std::f64::consts::PI;
            } else {
                g2 = 0.0;
                t2 = 0.0;
            }

            x1i = sx * nxi + sy * nyi;
            x2i = sx * nxi + sy * nyi;
            yyi = sx * nyi - sy * nxi;

            if geolin {
                let nxo = xf.nx[jo];
                let nyo = xf.ny[jo];
                let nxp = xf.nx[jp];
                let nyp = xf.ny[jp];

                x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
                x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
                x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
                x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
                yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
                yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
            }

            jo_last = jo;
            jp_last = jp;

            if jo == n - 1 {
                break;
            }

            if siglin {
                // set up midpoint quantities
                let x0 = 0.5 * (x1 + x2);
                let rs0 = x0 * x0 + yy * yy;
                let g0 = rs0.ln();
                let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * std::f64::consts::PI;

                // calculate source contribution to Psi for 1-0 half-panel
                let dxinv = 1.0 / (x1 - x0);
                let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
                let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

                let psx1 = -(t1 - apan);
                let psx0 = t0 - apan;
                let psyy = 0.5 * (g1 - g0);

                let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
                let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
                let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

                let dsm = ((xf.x[jp] - xf.x[jm]).powi(2) + (xf.y[jp] - xf.y[jm]).powi(2)).sqrt();
                let dsim = 1.0 / dsm;

                let ssum = (xf.sig[jp] - xf.sig[jo]) * dsio + (xf.sig[jp] - xf.sig[jm]) * dsim;
                let sdif = (xf.sig[jp] - xf.sig[jo]) * dsio - (xf.sig[jp] - xf.sig[jm]) * dsim;

                *psi += crate::state::QOPI * (psum * ssum + pdif * sdif);

                // dPsi/dm
                xf.dzdm[jm] += crate::state::QOPI * (-psum * dsim + pdif * dsim);
                xf.dzdm[jo] += crate::state::QOPI * (-psum * dsio - pdif * dsio);
                xf.dzdm[jp] += crate::state::QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

                // dPsi/dni
                let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
                let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
                *psi_ni += crate::state::QOPI * (psni * ssum + pdni * sdif);

                qtanm += crate::state::QOPI * (psni * ssum + pdni * sdif);

                xf.dqdm[jm] += crate::state::QOPI * (-psni * dsim + pdni * dsim);
                xf.dqdm[jo] += crate::state::QOPI * (-psni * dsio - pdni * dsio);
                xf.dqdm[jp] += crate::state::QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

                // calculate source contribution to Psi for 0-2 half-panel
                let dxinv = 1.0 / (x0 - x2);
                let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
                let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

                let psx0 = -(t0 - apan);
                let psx2 = t2 - apan;
                let psyy = 0.5 * (g0 - g2);

                let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
                let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
                let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

                let dsp = ((xf.x[jq] - xf.x[jo]).powi(2) + (xf.y[jq] - xf.y[jo]).powi(2)).sqrt();
                let dsip = 1.0 / dsp;

                let ssum = (xf.sig[jq] - xf.sig[jo]) * dsip + (xf.sig[jp] - xf.sig[jo]) * dsio;
                let sdif = (xf.sig[jq] - xf.sig[jo]) * dsip - (xf.sig[jp] - xf.sig[jo]) * dsio;

                *psi += crate::state::QOPI * (psum * ssum + pdif * sdif);

                // dPsi/dm
                xf.dzdm[jo] += crate::state::QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
                xf.dzdm[jp] += crate::state::QOPI * (psum * dsio - pdif * dsio);
                xf.dzdm[jq] += crate::state::QOPI * (psum * dsip + pdif * dsip);

                // dPsi/dni
                let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
                let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
                *psi_ni += crate::state::QOPI * (psni * ssum + pdni * sdif);

                qtanm += crate::state::QOPI * (psni * ssum + pdni * sdif);

                xf.dqdm[jo] += crate::state::QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
                xf.dqdm[jp] += crate::state::QOPI * (psni * dsio - pdni * dsio);
                xf.dqdm[jq] += crate::state::QOPI * (psni * dsip + pdni * dsip);
            }

            // calculate vortex panel contribution to Psi
            let dxinv = 1.0 / (x1 - x2);
            let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

            let psx1 = 0.5 * g1;
            let psx2 = -0.5 * g2;
            let psyy = t1 - t2;

            let pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
            let pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
            let pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

            let gsum1 = xf.gamu[0][jp] + xf.gamu[0][jo];
            let gsum2 = xf.gamu[1][jp] + xf.gamu[1][jo];
            let gdif1 = xf.gamu[0][jp] - xf.gamu[0][jo];
            let gdif2 = xf.gamu[1][jp] - xf.gamu[1][jo];

            let gsum = xf.gam[jp] + xf.gam[jo];
            let gdif = xf.gam[jp] - xf.gam[jo];

            *psi += crate::state::QOPI * (psis * gsum + psid * gdif);

            // dPsi/dGam
            xf.dzdg[jo] += crate::state::QOPI * (psis - psid);
            xf.dzdg[jp] += crate::state::QOPI * (psis + psid);

            // dPsi/dni
            let psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
            let pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
            *psi_ni += crate::state::QOPI * (gsum * psni + gdif * pdni);

            xf.qtan1 += crate::state::QOPI * (gsum1 * psni + gdif1 * pdni);
            xf.qtan2 += crate::state::QOPI * (gsum2 * psni + gdif2 * pdni);

            xf.dqdg[jo] += crate::state::QOPI * (psni - pdni);
            xf.dqdg[jp] += crate::state::QOPI * (psni + pdni);

            if geolin {
                // dPsi/dn
                xf.dzdn[jo] += crate::state::QOPI * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo)
                    + crate::state::QOPI * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
                xf.dzdn[jp] += crate::state::QOPI * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp)
                    + crate::state::QOPI * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);
                // dPsi/dP
                xf.z_qdof0 += crate::state::QOPI * ((psis - psid) * xf.qf0[jo] + (psis + psid) * xf.qf0[jp]);
                xf.z_qdof1 += crate::state::QOPI * ((psis - psid) * xf.qf1[jo] + (psis + psid) * xf.qf1[jp]);
                xf.z_qdof2 += crate::state::QOPI * ((psis - psid) * xf.qf2[jo] + (psis + psid) * xf.qf2[jp]);
                xf.z_qdof3 += crate::state::QOPI * ((psis - psid) * xf.qf3[jo] + (psis + psid) * xf.qf3[jp]);
            }
        }
    }

    let jo = jo_last;
    let jp = jp_last;

    if !te_skip {
        let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

        let psigx1 = -(t1 - apan);
        let psigx2 = t2 - apan;
        let psigyy = 0.5 * (g1 - g2);
        let pgamx1 = 0.5 * g1;
        let pgamx2 = -0.5 * g2;
        let pgamyy = t1 - t2;

        let psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
        let pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

        // TE panel source and vortex strengths
        let sigte1 = 0.5 * scs * (xf.gamu[0][jp] - xf.gamu[0][jo]);
        let sigte2 = 0.5 * scs * (xf.gamu[1][jp] - xf.gamu[1][jo]);
        let gamte1 = -0.5 * sds * (xf.gamu[0][jp] - xf.gamu[0][jo]);
        let gamte2 = -0.5 * sds * (xf.gamu[1][jp] - xf.gamu[1][jo]);

        xf.sigte = 0.5 * scs * (xf.gam[jp] - xf.gam[jo]);
        xf.gamte = -0.5 * sds * (xf.gam[jp] - xf.gam[jo]);

        // TE panel contribution to Psi
        *psi += crate::state::HOPI * (psig * xf.sigte + pgam * xf.gamte);

        // dPsi/dGam
        xf.dzdg[jo] -= crate::state::HOPI * psig * scs * 0.5;
        xf.dzdg[jp] += crate::state::HOPI * psig * scs * 0.5;

        xf.dzdg[jo] += crate::state::HOPI * pgam * sds * 0.5;
        xf.dzdg[jp] -= crate::state::HOPI * pgam * sds * 0.5;

        // dPsi/dni
        *psi_ni += crate::state::HOPI * (psigni * xf.sigte + pgamni * xf.gamte);

        xf.qtan1 += crate::state::HOPI * (psigni * sigte1 + pgamni * gamte1);
        xf.qtan2 += crate::state::HOPI * (psigni * sigte2 + pgamni * gamte2);

        xf.dqdg[jo] -= crate::state::HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
        xf.dqdg[jp] += crate::state::HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

        if geolin {
            // dPsi/dn
            xf.dzdn[jo] += crate::state::HOPI * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * xf.sigte
                + crate::state::HOPI * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * xf.gamte;
            xf.dzdn[jp] += crate::state::HOPI * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * xf.sigte
                + crate::state::HOPI * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * xf.gamte;

            // dPsi/dP
            xf.z_qdof0 += crate::state::HOPI * psig * 0.5 * (xf.qf0[jp] - xf.qf0[jo]) * scs
                - crate::state::HOPI * pgam * 0.5 * (xf.qf0[jp] - xf.qf0[jo]) * sds;
            xf.z_qdof1 += crate::state::HOPI * psig * 0.5 * (xf.qf1[jp] - xf.qf1[jo]) * scs
                - crate::state::HOPI * pgam * 0.5 * (xf.qf1[jp] - xf.qf1[jo]) * sds;
            xf.z_qdof2 += crate::state::HOPI * psig * 0.5 * (xf.qf2[jp] - xf.qf2[jo]) * scs
                - crate::state::HOPI * pgam * 0.5 * (xf.qf2[jp] - xf.qf2[jo]) * sds;
            xf.z_qdof3 += crate::state::HOPI * psig * 0.5 * (xf.qf3[jp] - xf.qf3[jo]) * scs
                - crate::state::HOPI * pgam * 0.5 * (xf.qf3[jp] - xf.qf3[jo]) * sds;
        }
    }

    // freestream terms
    *psi += xf.qinf * (xf.cosa * yi - xf.sina * xi);

    // dPsi/dn
    *psi_ni += xf.qinf * (xf.cosa * nyi - xf.sina * nxi);

    xf.qtan1 += xf.qinf * nyi;
    xf.qtan2 -= xf.qinf * nxi;

    // dPsi/dQinf
    xf.z_qinf += xf.cosa * yi - xf.sina * xi;

    // dPsi/dalfa
    xf.z_alfa -= xf.qinf * (xf.sina * yi + xf.cosa * xi);

    // NOTE: the image-airfoil (ground effect, LIMage) contribution from the
    // original Fortran is omitted; LIMage is never enabled in this port.
    let _ = (qtanm, nw, x2o, x2p, yyp);
}

/// Calculates current streamfunction Psi and tangential velocity Qtan at
/// panel node or wake node I due to freestream and wake sources Sig, plus
/// sensitivity vectors dPsi/dSig (DZDM) and dQtan/dSig (DQDM).
pub fn pswlin(xf: &mut Xfoil, i: usize, xi: f64, yi: f64, nxi: f64, nyi: f64, psi: &mut f64, psi_ni: &mut f64) {
    let n = xf.n;
    let nw = xf.nw;

    let io = i;

    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();

    for jo in n..n + nw {
        xf.dzdm[jo] = 0.0;
        xf.dqdm[jo] = 0.0;
    }

    *psi = 0.0;
    *psi_ni = 0.0;

    for jo in n..n + nw - 1 {
        let jp = jo + 1;

        let jm;
        let jq;
        if jo == n {
            jm = jo;
            jq = jp + 1;
        } else if jo == n + nw - 2 {
            jm = jo - 1;
            jq = jp;
        } else {
            jm = jo - 1;
            jq = jp + 1;
        }

        let dso = ((xf.x[jo] - xf.x[jp]).powi(2) + (xf.y[jo] - xf.y[jp]).powi(2)).sqrt();
        let dsio = 1.0 / dso;

        let apan = xf.apanel[jo];

        let rx1 = xi - xf.x[jo];
        let ry1 = yi - xf.y[jo];
        let rx2 = xi - xf.x[jp];
        let ry2 = yi - xf.y[jp];

        let sx = (xf.x[jp] - xf.x[jo]) * dsio;
        let sy = (xf.y[jp] - xf.y[jo]) * dsio;

        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn;
        if io >= n && io < n + nw {
            sgn = 1.0;
        } else {
            sgn = yy.signum();
        }

        let (g1, t1);
        if io != jo && rs1 > 0.0 {
            g1 = rs1.ln();
            t1 = (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * std::f64::consts::PI;
        } else {
            g1 = 0.0;
            t1 = 0.0;
        }

        let (g2, t2);
        if io != jp && rs2 > 0.0 {
            g2 = rs2.ln();
            t2 = (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * std::f64::consts::PI;
        } else {
            g2 = 0.0;
            t2 = 0.0;
        }

        let x1i = sx * nxi + sy * nyi;
        let x2i = sx * nxi + sy * nyi;
        let yyi = sx * nyi - sy * nxi;

        // set up midpoint quantities
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = rs0.ln();
        let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * std::f64::consts::PI;

        // calculate source contribution to Psi for 1-0 half-panel
        let dxinv = 1.0 / (x1 - x0);
        let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
        let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

        let psx1 = -(t1 - apan);
        let psx0 = t0 - apan;
        let psyy = 0.5 * (g1 - g0);

        let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
        let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
        let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

        let dsm = ((xf.x[jp] - xf.x[jm]).powi(2) + (xf.y[jp] - xf.y[jm]).powi(2)).sqrt();
        let dsim = 1.0 / dsm;

        let ssum = (xf.sig[jp] - xf.sig[jo]) * dsio + (xf.sig[jp] - xf.sig[jm]) * dsim;
        let sdif = (xf.sig[jp] - xf.sig[jo]) * dsio - (xf.sig[jp] - xf.sig[jm]) * dsim;

        *psi += crate::state::QOPI * (psum * ssum + pdif * sdif);

        // dPsi/dm
        xf.dzdm[jm] += crate::state::QOPI * (-psum * dsim + pdif * dsim);
        xf.dzdm[jo] += crate::state::QOPI * (-psum * dsio - pdif * dsio);
        xf.dzdm[jp] += crate::state::QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

        // dPsi/dni
        let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
        let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
        *psi_ni += crate::state::QOPI * (psni * ssum + pdni * sdif);

        xf.dqdm[jm] += crate::state::QOPI * (-psni * dsim + pdni * dsim);
        xf.dqdm[jo] += crate::state::QOPI * (-psni * dsio - pdni * dsio);
        xf.dqdm[jp] += crate::state::QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

        // calculate source contribution to Psi for 0-2 half-panel
        let dxinv = 1.0 / (x0 - x2);
        let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
        let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

        let psx0 = -(t0 - apan);
        let psx2 = t2 - apan;
        let psyy = 0.5 * (g0 - g2);

        let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
        let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
        let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

        let dsp = ((xf.x[jq] - xf.x[jo]).powi(2) + (xf.y[jq] - xf.y[jo]).powi(2)).sqrt();
        let dsip = 1.0 / dsp;

        let ssum = (xf.sig[jq] - xf.sig[jo]) * dsip + (xf.sig[jp] - xf.sig[jo]) * dsio;
        let sdif = (xf.sig[jq] - xf.sig[jo]) * dsip - (xf.sig[jp] - xf.sig[jo]) * dsio;

        *psi += crate::state::QOPI * (psum * ssum + pdif * sdif);

        // dPsi/dm
        xf.dzdm[jo] += crate::state::QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
        xf.dzdm[jp] += crate::state::QOPI * (psum * dsio - pdif * dsio);
        xf.dzdm[jq] += crate::state::QOPI * (psum * dsip + pdif * dsip);

        // dPsi/dni
        let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
        let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
        *psi_ni += crate::state::QOPI * (psni * ssum + pdni * sdif);

        xf.dqdm[jo] += crate::state::QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
        xf.dqdm[jp] += crate::state::QOPI * (psni * dsio - pdni * dsio);
        xf.dqdm[jq] += crate::state::QOPI * (psni * dsip + pdni * dsip);
    }
}

/// Calculates two surface vorticity (gamma) distributions for alpha = 0, 90
/// degrees, superimposed in SPECAL or SPECCL for specified alpha or CL.
pub fn ggcalc(xf: &mut Xfoil) {
    let n = xf.n;

    // distance of internal control point ahead of sharp TE
    // (fraction of smaller panel length adjacent to TE)
    let bwt = 0.1;

    if xf.show_output {
        eprintln!("Calculating unit vorticity distributions ...");
    }

    for i in 0..n {
        xf.gam[i] = 0.0;
        xf.gamu[0][i] = 0.0;
        xf.gamu[1][i] = 0.0;
    }
    xf.psio = 0.0;

    // Set up matrix system for Psi = Psio on airfoil surface.
    // The unknowns are (dGamma)i and dPsio.
    for i in 0..n {
        // calculate Psi and dPsi/dGamma array for current node
        let mut psi = 0.0;
        let mut psi_n = 0.0;
        psilin(xf, i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], &mut psi, &mut psi_n, false, true);

        // RES1 = PSI( 0) - PSIO
        // RES2 = PSI(90) - PSIO
        let res1 = xf.qinf * xf.y[i];
        let res2 = -xf.qinf * xf.x[i];

        // dRes/dGamma
        for j in 0..n {
            xf.aij[Xfoil::m_index(i, j)] = xf.dzdg[j];
        }

        for j in 0..n {
            xf.bij[Xfoil::b_index(i, j)] = -xf.dzdm[j];
        }

        // dRes/dPsio
        xf.aij[Xfoil::m_index(i, n)] = -1.0;

        xf.gamu[0][i] = -res1;
        xf.gamu[1][i] = -res2;
    }

    // set Kutta condition: RES = GAM(1) + GAM(N)
    for j in 0..n + 1 {
        xf.aij[Xfoil::m_index(n, j)] = 0.0;
    }

    xf.aij[Xfoil::m_index(n, 0)] = 1.0;
    xf.aij[Xfoil::m_index(n, n - 1)] = 1.0;

    xf.gamu[0][n] = 0.0;
    xf.gamu[1][n] = 0.0;

    // set up Kutta condition (no direct source influence)
    for j in 0..n {
        xf.bij[Xfoil::b_index(n, j)] = 0.0;
    }

    if xf.sharp {
        // set zero internal velocity in TE corner
        // set TE bisector angle
        let ag1 = (-xf.yp[0]).atan2(-xf.xp[0]);
        let ag2 = atanc(xf.yp[n - 1], xf.xp[n - 1], ag1);
        let abis = 0.5 * (ag1 + ag2);
        let cbis = abis.cos();
        let sbis = abis.sin();

        // minimum panel length adjacent to TE
        let ds1 = ((xf.x[0] - xf.x[1]).powi(2) + (xf.y[0] - xf.y[1]).powi(2)).sqrt();
        let ds2 = ((xf.x[n - 1] - xf.x[n - 2]).powi(2) + (xf.y[n - 1] - xf.y[n - 2]).powi(2)).sqrt();
        let dsmin = ds1.min(ds2);

        // control point on bisector just ahead of TE point
        let xbis = xf.xte - bwt * dsmin * cbis;
        let ybis = xf.yte - bwt * dsmin * sbis;

        // set velocity component along bisector line
        let mut psi = 0.0;
        let mut qbis = 0.0;
        psilin(xf, 0, xbis, ybis, -sbis, cbis, &mut psi, &mut qbis, false, true);

        // RES = QDg*Gam + QDm*Mass + QINF*(COSA*CBIS + SINA*SBIS)
        let res = qbis;

        // dRes/dGamma
        for j in 0..n {
            xf.aij[Xfoil::m_index(n - 1, j)] = xf.dqdg[j];
        }

        // -dRes/dMass
        for j in 0..n {
            xf.bij[Xfoil::b_index(n - 1, j)] = -xf.dqdm[j];
        }

        // dRes/dPsio
        xf.aij[Xfoil::m_index(n - 1, n)] = 0.0;

        // -dRes/dUinf
        xf.gamu[0][n - 1] = -cbis;

        // -dRes/dVinf
        xf.gamu[1][n - 1] = -sbis;
    }

    // LU-factor coefficient matrix AIJ
    ludcmp(&mut xf.aij, xf.aijpiv.as_mut_slice(), n + 1);
    xf.lqaij = true;

    // solve system for the two vorticity distributions
    baksub(&xf.aij, &xf.aijpiv, &mut xf.gamu[0], n + 1);
    baksub(&xf.aij, &xf.aijpiv, &mut xf.gamu[1], n + 1);

    // set inviscid alpha=0,90 surface speeds for this geometry
    for i in 0..n {
        xf.qinvu[0][i] = xf.gamu[0][i];
        xf.qinvu[1][i] = xf.gamu[1][i];
    }

    xf.lgamu = true;
}

/// Sets inviscid tangential velocity for alpha = 0, 90 on the wake due to
/// freestream and airfoil surface vorticity.
pub fn qwcalc(xf: &mut Xfoil) {
    let n = xf.n;

    // first wake point (same as TE)
    xf.qinvu[0][n] = xf.qinvu[0][n - 1];
    xf.qinvu[1][n] = xf.qinvu[1][n - 1];

    // rest of wake
    for i in n + 1..n + xf.nw {
        let mut psi = 0.0;
        let mut psi_ni = 0.0;
        psilin(xf, i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], &mut psi, &mut psi_ni, false, false);
        xf.qinvu[0][i] = xf.qtan1;
        xf.qinvu[1][i] = xf.qtan2;
    }
}

/// Calculates source panel influence coefficient matrix for current airfoil
/// and wake geometry.
pub fn qdcalc(xf: &mut Xfoil) {
    let n = xf.n;
    let nw = xf.nw;

    if xf.show_output {
        eprintln!("Calculating source influence matrix ...");
    }

    if !xf.ladij {
        // calculate source influence matrix for airfoil surface if it doesn't exist
        for j in 0..n {
            // multiply each dPsi/Sig vector by inverse of factored dPsi/dGam matrix
            baksub(&xf.aij, &xf.aijpiv, &mut xf.bij[j * IQX..], n + 1);

            // store resulting dGam/dSig = dQtan/dSig vector
            for i in 0..n {
                xf.dij[Xfoil::d_index(i, j)] = xf.bij[Xfoil::b_index(i, j)];
            }
        }
        xf.ladij = true;
    }

    // set up coefficient matrix of dPsi/dm on airfoil surface
    for i in 0..n {
        let mut psi = 0.0;
        let mut psi_n = 0.0;
        pswlin(xf, i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], &mut psi, &mut psi_n);
        for j in n..n + nw {
            xf.bij[Xfoil::b_index(i, j)] = -xf.dzdm[j];
        }
    }

    // set up Kutta condition (no direct source influence)
    for j in n..n + nw {
        xf.bij[Xfoil::b_index(n, j)] = 0.0;
    }

    // sharp TE gamma extrapolation also has no source influence
    if xf.sharp {
        for j in n..n + nw {
            xf.bij[Xfoil::b_index(n - 1, j)] = 0.0;
        }
    }

    // multiply by inverse of factored dPsi/dGam matrix
    for j in n..n + nw {
        baksub(&xf.aij, &xf.aijpiv, &mut xf.bij[j * IQX..], n + 1);
    }

    // set the source influence matrix for the wake sources
    for i in 0..n {
        for j in n..n + nw {
            xf.dij[Xfoil::d_index(i, j)] = xf.bij[Xfoil::b_index(i, j)];
        }
    }

    // Now calculate the influence of sources on the wake velocities:
    // dQtan/dGam and dQtan/dSig at the wake points
    for i in n..n + nw {
        let iw = i - n;
        if iw > IWX {
            break;
        }

        // airfoil contribution at wake panel node
        let mut psi = 0.0;
        let mut psi_n = 0.0;
        psilin(xf, i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], &mut psi, &mut psi_n, false, true);

        for j in 0..n {
            xf.cij[j * IWX + iw] = xf.dqdg[j];
        }

        for j in 0..n {
            xf.dij[Xfoil::d_index(i, j)] = xf.dqdm[j];
        }

        // wake contribution
        pswlin(xf, i, xf.x[i], xf.y[i], xf.nx[i], xf.ny[i], &mut psi, &mut psi_n);

        for j in n..n + nw {
            xf.dij[Xfoil::d_index(i, j)] = xf.dqdm[j];
        }
    }

    // add on effect of all sources on airfoil vorticity which effects wake Qtan
    for i in n..n + nw {
        let iw = i - n;

        // airfoil surface source contribution first
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += xf.cij[k * IWX + iw] * xf.dij[Xfoil::d_index(k, j)];
            }
            xf.dij[Xfoil::d_index(i, j)] += sum;
        }

        // wake source contribution next
        for j in n..n + nw {
            let mut sum = 0.0;
            for k in 0..n {
                sum += xf.cij[k * IWX + iw] * xf.bij[Xfoil::b_index(k, j)];
            }
            xf.dij[Xfoil::d_index(i, j)] += sum;
        }
    }

    // make sure first wake point has same velocity as trailing edge
    for j in 0..n + nw {
        xf.dij[Xfoil::d_index(n, j)] = xf.dij[Xfoil::d_index(n - 1, j)];
    }

    xf.lwdij = true;
}

/// Sets wake coordinate array for current surface vorticity and/or mass
/// source distributions.
pub fn xywake(xf: &mut Xfoil) {
    let n = xf.n;

    if xf.show_output {
        eprintln!("Calculating wake trajectory ...");
    }

    // number of wake points
    xf.nw = n / 12 + (10.0 * xf.waklen) as usize;
    if xf.nw > IWX {
        if xf.show_output {
            eprintln!("Array size (IWX) too small.  Last wake point index reduced.");
        }
        xf.nw = IWX;
    }

    let ds1 = 0.5 * (xf.s[1] - xf.s[0] + xf.s[n - 1] - xf.s[n - 2]);
    crate::utils::setexp(&mut xf.snew[n..n + xf.nw], ds1, xf.waklen * xf.chord, xf.show_output);

    xf.xte = 0.5 * (xf.x[0] + xf.x[n - 1]);
    xf.yte = 0.5 * (xf.y[0] + xf.y[n - 1]);

    // set first wake point a tiny distance behind TE
    let i = n;
    let sx = 0.5 * (xf.yp[n - 1] - xf.yp[0]);
    let sy = 0.5 * (xf.xp[0] - xf.xp[n - 1]);
    let smod = (sx * sx + sy * sy).sqrt();
    xf.nx[i] = sx / smod;
    xf.ny[i] = sy / smod;
    xf.x[i] = xf.xte - 0.0001 * xf.ny[i];
    xf.y[i] = xf.yte + 0.0001 * xf.nx[i];
    xf.s[i] = xf.s[n - 1];

    // calculate streamfunction gradient components at first point
    let mut psi = 0.0;
    let mut psi_x = 0.0;
    let mut psi_y = 0.0;
    psilin(xf, i, xf.x[i], xf.y[i], 1.0, 0.0, &mut psi, &mut psi_x, false, false);
    psilin(xf, i, xf.x[i], xf.y[i], 0.0, 1.0, &mut psi, &mut psi_y, false, false);

    // set unit vector normal to wake at first point
    xf.nx[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
    xf.ny[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

    // set angle of wake panel normal
    xf.apanel[i] = psi_y.atan2(psi_x);

    // set rest of wake points
    for i in n + 1..n + xf.nw {
        let ds = xf.snew[i] - xf.snew[i - 1];

        // set new point DS downstream of last point
        xf.x[i] = xf.x[i - 1] - ds * xf.ny[i];
        xf.y[i] = xf.y[i - 1] + ds * xf.nx[i];
        xf.s[i] = xf.s[i - 1] + ds;

        if i != n + xf.nw - 1 {
            // calculate normal vector for next point
            let mut psi = 0.0;
            let mut psi_x = 0.0;
            let mut psi_y = 0.0;
            psilin(xf, i, xf.x[i], xf.y[i], 1.0, 0.0, &mut psi, &mut psi_x, false, false);
            psilin(xf, i, xf.x[i], xf.y[i], 0.0, 1.0, &mut psi, &mut psi_y, false, false);

            xf.nx[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
            xf.ny[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

            // set angle of wake panel normal
            xf.apanel[i] = psi_y.atan2(psi_x);
        }
    }

    // set wake presence flag and corresponding alpha
    xf.lwake = true;
    xf.awake = xf.alfa;

    // old source influence matrix is invalid for the new wake geometry
    xf.lwdij = false;
}

/// Locates stagnation point arc length location SST and panel index IST.
pub fn stfind(xf: &mut Xfoil) {
    let n = xf.n;

    let mut i = n - 1;
    let mut found = false;
    for ii in 0..n - 1 {
        if xf.gam[ii] >= 0.0 && xf.gam[ii + 1] < 0.0 {
            i = ii;
            found = true;
            break;
        }
    }

    if !found {
        if xf.show_output {
            eprintln!("STFIND: Stagnation point not found. Continuing ...");
        }
        i = n / 2;
    }

    xf.ist = i;
    let dgam = xf.gam[i + 1] - xf.gam[i];
    let ds = xf.s[i + 1] - xf.s[i];

    // evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
    if xf.gam[i] < -xf.gam[i + 1] {
        xf.sst = xf.s[i] - ds * (xf.gam[i] / dgam);
    } else {
        xf.sst = xf.s[i + 1] - ds * (xf.gam[i + 1] / dgam);
    }

    // tweak stagnation point if it falls right on a node (very unlikely)
    if xf.sst <= xf.s[i] {
        xf.sst = xf.s[i] + 1.0E-7;
    }
    if xf.sst >= xf.s[i + 1] {
        xf.sst = xf.s[i + 1] - 1.0E-7;
    }

    xf.sst_go = (xf.sst - xf.s[i + 1]) / dgam;
    xf.sst_gp = (xf.s[i] - xf.sst) / dgam;
}

/// Sets BL location -> panel location pointer array IPAN.
pub fn iblpan(xf: &mut Xfoil) {
    // top surface first (Fortran increments ibl before storing, so entry 1
    // is left unused as the "similarity" station)
    let mut ibl = 1usize;
    for i in (0..=xf.ist).rev() {
        ibl += 1;
        xf.ipan[0][ibl] = i as i32;
        xf.vti[0][ibl] = 1.0;
    }
    xf.iblte[0] = ibl as i32;
    xf.nbl[0] = ibl as i32;

    // bottom surface next
    ibl = 1;
    for i in xf.ist + 1..xf.n {
        ibl += 1;
        xf.ipan[1][ibl] = i as i32;
        xf.vti[1][ibl] = -1.0;
    }

    // wake
    xf.iblte[1] = ibl as i32;
    for iw in 1..=xf.nw {
        let i = xf.n + iw - 1;
        ibl = xf.iblte[1] as usize + iw;
        xf.ipan[1][ibl] = i as i32;
        xf.vti[1][ibl] = -1.0;
    }
    xf.nbl[1] = xf.iblte[1] + xf.nw as i32;

    // upper wake pointers (for plotting only)
    for iw in 1..=xf.nw {
        xf.ipan[0][xf.iblte[0] as usize + iw] = xf.ipan[1][xf.iblte[1] as usize + iw];
        xf.vti[0][xf.iblte[0] as usize + iw] = 1.0;
    }

    let iblmax = xf.iblte[0].max(xf.iblte[1]) as usize + xf.nw;
    if iblmax > crate::state::IVX {
        if xf.show_output {
            eprintln!(" ***  BL array overflow.");
            eprintln!(" ***  Increase IVX to at least {}", iblmax);
        }
        panic!("IBLPAN: BL array overflow");
    }

    xf.lipan = true;
}

/// Sets BL arc length array on each airfoil side and wake.
pub fn xicalc(xf: &mut Xfoil) {
    let n = xf.n;

    // minimum xi node arc length near stagnation point
    let xeps = 1.0E-7 * (xf.s[n - 1] - xf.s[0]);

    let mut is = 0usize; // side index (0-based)
    xf.xssi[is][0] = 0.0;
    for ibl in 2..=xf.iblte[is] as usize {
        let i = xf.ipan[is][ibl] as usize;
        xf.xssi[is][ibl] = (xf.sst - xf.s[i]).max(xeps);
    }

    is = 1;
    xf.xssi[is][0] = 0.0;
    for ibl in 2..=xf.iblte[is] as usize {
        let i = xf.ipan[is][ibl] as usize;
        xf.xssi[is][ibl] = (xf.s[i] - xf.sst).max(xeps);
    }

    let is1 = 0usize;
    let is2 = 1usize;

    let ibl1 = xf.iblte[is1] as usize + 1;
    xf.xssi[is1][ibl1] = xf.xssi[is1][ibl1 - 1];

    let ibl2 = xf.iblte[is2] as usize + 1;
    xf.xssi[is2][ibl2] = xf.xssi[is2][ibl2 - 1];

    for ibl in xf.iblte[is] as usize + 2..=xf.nbl[is] as usize {
        let i = xf.ipan[is][ibl] as usize;
        let dxssi = ((xf.x[i] - xf.x[i - 1]).powi(2) + (xf.y[i] - xf.y[i - 1]).powi(2)).sqrt();

        let ibl1 = xf.iblte[is1] as usize + ibl - xf.iblte[is] as usize;
        let ibl2 = xf.iblte[is2] as usize + ibl - xf.iblte[is] as usize;
        xf.xssi[is1][ibl1] = xf.xssi[is1][ibl1 - 1] + dxssi;
        xf.xssi[is2][ibl2] = xf.xssi[is2][ibl2 - 1] + dxssi;
    }

    // trailing edge flap length to TE gap ratio
    let telrat = 2.50;

    // set up parameters for TE flap cubics
    let crosp = (xf.xp[0] * xf.yp[n - 1] - xf.yp[0] * xf.xp[n - 1])
        / ((xf.xp[0].powi(2) + xf.yp[0].powi(2)) * (xf.xp[n - 1].powi(2) + xf.yp[n - 1].powi(2))).sqrt();
    let mut dwdxte = crosp / (1.0 - crosp * crosp).sqrt();

    // limit cubic to avoid absurd TE gap widths
    dwdxte = dwdxte.max(-3.0 / telrat);
    dwdxte = dwdxte.min(3.0 / telrat);

    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;

    if xf.sharp {
        for wgap in xf.wgap.iter_mut() {
            *wgap = 0.0;
        }
    } else {
        // set TE flap (wake gap) array
        let is = 1usize;
        for iw in 1..=xf.nw {
            let ibl = xf.iblte[is] as usize + iw;
            let zn = 1.0 - (xf.xssi[is][ibl] - xf.xssi[is][xf.iblte[is] as usize]) / (telrat * xf.ante);
            xf.wgap[iw - 1] = 0.0;
            if zn >= 0.0 {
                xf.wgap[iw - 1] = xf.ante * (aa + bb * zn) * zn * zn;
            }
        }
    }
}

/// Sets inviscid Ue from panel inviscid tangential velocity.
pub fn uicalc(xf: &mut Xfoil) {
    for is in 0..2 {
        xf.uinv[is][0] = 0.0;
        xf.uinv_a[is][0] = 0.0;
        for ibl in 2..=xf.nbl[is] as usize {
            let i = xf.ipan[is][ibl] as usize;
            xf.uinv[is][ibl] = xf.vti[is][ibl] * xf.qinv[i];
            xf.uinv_a[is][ibl] = xf.vti[is][ibl] * xf.qinv_a[i];
        }
    }
}

/// Sets viscous Ue from panel viscous tangential velocity.
pub fn uecalc(xf: &mut Xfoil) {
    for is in 0..2 {
        xf.uedg[is][0] = 0.0;
        for ibl in 2..=xf.nbl[is] as usize {
            let i = xf.ipan[is][ibl] as usize;
            xf.uedg[is][ibl] = xf.vti[is][ibl] * xf.qvis[i];
        }
    }
}

/// Sets panel viscous tangential velocity from viscous Ue.
pub fn qvfue(xf: &mut Xfoil) {
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let i = xf.ipan[is][ibl] as usize;
            xf.qvis[i] = xf.vti[is][ibl] * xf.uedg[is][ibl];
        }
    }
}

/// Sets inviscid panel tangential velocity for current alpha.
pub fn qiset(xf: &mut Xfoil) {
    xf.cosa = xf.alfa.cos();
    xf.sina = xf.alfa.sin();

    for i in 0..xf.n + xf.nw {
        xf.qinv[i] = xf.cosa * xf.qinvu[0][i] + xf.sina * xf.qinvu[1][i];
        xf.qinv_a[i] = -xf.sina * xf.qinvu[0][i] + xf.cosa * xf.qinvu[1][i];
    }
}

/// Sets GAM from the viscous tangential velocity QVIS.
pub fn gamqv(xf: &mut Xfoil) {
    for i in 0..xf.n {
        xf.gam[i] = xf.qvis[i];
        xf.gam_a[i] = xf.qinv_a[i];
    }
}

/// Moves stagnation point location to new panel.
pub fn stmove(xf: &mut Xfoil) {
    let istold = xf.ist;
    stfind(xf);

    if istold == xf.ist {
        // recalculate new arc length array
        xicalc(xf);
    } else {
        // set new BL position -> panel position pointers
        iblpan(xf);

        // set new inviscid BL edge velocity UINV from QINV
        uicalc(xf);

        // recalculate new arc length array
        xicalc(xf);

        // set BL position -> system line pointers
        crate::s_xbl::iblsys(xf);

        if xf.ist > istold {
            // increase in number of points on top side (IS=1)
            let idif = xf.ist - istold;

            xf.itran[0] += idif as i32;
            xf.itran[1] -= idif as i32;

            // move top side BL variables downstream
            for ibl in (idif + 2..=xf.nbl[0] as usize).rev() {
                xf.ctau[0][ibl] = xf.ctau[0][ibl - idif];
                xf.thet[0][ibl] = xf.thet[0][ibl - idif];
                xf.dstr[0][ibl] = xf.dstr[0][ibl - idif];
                xf.uedg[0][ibl] = xf.uedg[0][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = xf.uedg[0][idif + 2] / xf.xssi[0][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                xf.ctau[0][ibl] = xf.ctau[0][idif + 2];
                xf.thet[0][ibl] = xf.thet[0][idif + 2];
                xf.dstr[0][ibl] = xf.dstr[0][idif + 2];
                xf.uedg[0][ibl] = dudx * xf.xssi[0][ibl];
            }

            // move bottom side BL variables upstream
            for ibl in 2..=xf.nbl[1] as usize {
                xf.ctau[1][ibl] = xf.ctau[1][ibl + idif];
                xf.thet[1][ibl] = xf.thet[1][ibl + idif];
                xf.dstr[1][ibl] = xf.dstr[1][ibl + idif];
                xf.uedg[1][ibl] = xf.uedg[1][ibl + idif];
            }
        } else {
            // increase in number of points on bottom side (IS=2)
            let idif = istold - xf.ist;

            xf.itran[0] -= idif as i32;
            xf.itran[1] += idif as i32;

            // move bottom side BL variables downstream
            for ibl in (idif + 2..=xf.nbl[1] as usize).rev() {
                xf.ctau[1][ibl] = xf.ctau[1][ibl - idif];
                xf.thet[1][ibl] = xf.thet[1][ibl - idif];
                xf.dstr[1][ibl] = xf.dstr[1][ibl - idif];
                xf.uedg[1][ibl] = xf.uedg[1][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = xf.uedg[1][idif + 2] / xf.xssi[1][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                xf.ctau[1][ibl] = xf.ctau[1][idif + 2];
                xf.thet[1][ibl] = xf.thet[1][idif + 2];
                xf.dstr[1][ibl] = xf.dstr[1][idif + 2];
                xf.uedg[1][ibl] = dudx * xf.xssi[1][ibl];
            }

            // move top side BL variables upstream
            for ibl in 2..=xf.nbl[0] as usize {
                xf.ctau[0][ibl] = xf.ctau[0][ibl + idif];
                xf.thet[0][ibl] = xf.thet[0][ibl + idif];
                xf.dstr[0][ibl] = xf.dstr[0][ibl + idif];
                xf.uedg[0][ibl] = xf.uedg[0][ibl + idif];
            }
        }

        // tweak Ue so it's not zero, in case stag. point is right on node
        let ueps = 1.0E-7;
        for is in 0..2 {
            for ibl in 2..=xf.nbl[is] as usize {
                let i = xf.ipan[is][ibl] as usize;
                if xf.uedg[is][ibl] <= ueps {
                    xf.uedg[is][ibl] = ueps;
                    xf.qvis[i] = xf.vti[is][ibl] * ueps;
                    xf.gam[i] = xf.vti[is][ibl] * ueps;
                }
            }
        }
    }

    // set new mass array since Ue has been tweaked
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            xf.mass[is][ibl] = xf.dstr[is][ibl] * xf.uedg[is][ibl];
        }
    }
}

/// Sets Ue from inviscid Ue plus all source influence.
pub fn ueset(xf: &mut Xfoil) {
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            let i = xf.ipan[is][ibl] as usize;

            let mut dui = 0.0;
            for js in 0..2 {
                for jbl in 2..=xf.nbl[js] as usize {
                    let j = xf.ipan[js][jbl] as usize;
                    let ue_m = -xf.vti[is][ibl] * xf.vti[js][jbl] * xf.dij[Xfoil::d_index(i, j)];
                    dui += ue_m * xf.mass[js][jbl];
                }
            }

            xf.uedg[is][ibl] = xf.uinv[is][ibl] + dui;
        }
    }
}

/// Sets displacement thickness DSTR from mass defect MASS.
pub fn dsset(xf: &mut Xfoil) {
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            xf.dstr[is][ibl] = xf.mass[is][ibl] / xf.uedg[is][ibl];
        }
    }
}
