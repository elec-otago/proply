//! Top-level XFOIL driver routines (port of `m_xfoil.f90`): init, naca,
//! pangen.  The interactive top-level menu is not ported.

use crate::bl::blpini;
use crate::geom::{cang, geopar, lefind};
use crate::naca::{naca4, naca5};
use crate::panel::{apcalc, ncalc};
use crate::spline::{curv, deval, scalc, segspl, seval, trisol};
use crate::s_xfoil::{comset, mrcl, tecalc};
use crate::state::{IQX, Xfoil};

/// Variable initialization/default routine.
pub fn init(xf: &mut Xfoil) {
    // set unity freestream speed
    xf.qinf = 1.0;

    // initialize freestream Mach number to zero
    xf.matyp = 1;
    xf.minf1 = 0.0;

    xf.alfa = 0.0;
    xf.cosa = 1.0;
    xf.sina = 0.0;

    for i in 0..IQX {
        xf.gamu[0][i] = 0.0;
        xf.gamu[1][i] = 0.0;
        xf.gam[i] = 0.0;
        xf.gam_a[i] = 0.0;
    }
    xf.psio = 0.0;

    xf.cl = 0.0;
    xf.cm = 0.0;
    xf.cd = 0.0;

    xf.sigte = 0.0;
    xf.gamte = 0.0;
    xf.sigte_a = 0.0;
    xf.gamte_a = 0.0;

    for i in 0..crate::state::IZX {
        xf.sig[i] = 0.0;
    }

    xf.awake = 0.0;
    xf.avisc = 0.0;

    xf.kimage = 1;
    xf.yimage = -10.0;
    xf.limage = false;

    xf.lgamu = false;
    xf.lqinu = false;
    xf.lvisc = false;
    xf.lwake = false;
    xf.lbli_ni = false;
    xf.lipan = false;
    xf.lqaij = false;
    xf.ladij = false;
    xf.lwdij = false;
    xf.lvconv = false;
    xf.lbflap = false;
    xf.lflap = false;
    xf.leiw = false;
    xf.lscini = false;
    xf.lgsame = false;

    // input airfoil will be normalized to unit chord (XFOIL default)
    xf.lnorm = true;

    // buffer and current airfoil flap hinge coordinates
    xf.xbf = 0.0;
    xf.ybf = 0.0;
    xf.xof = 0.0;
    xf.yof = 0.0;

    // circle plane array size (257, or largest 2^n + 1 that will fit array size)
    let ann = ((2 * IQX - 1) as f64).ln() / 2.0f64.ln();
    let nn = (ann + 0.00001) as i32;
    xf.nc1 = ((1usize << nn) + 1).min(257);

    // default paneling parameters
    xf.npan = 160;
    xf.cvpar = 1.0;
    xf.cterat = 0.15;
    xf.ctrrat = 0.2;

    // default paneling refinement zone x/c endpoints
    xf.xsref1 = 1.0;
    xf.xsref2 = 1.0;
    xf.xpref1 = 1.0;
    xf.xpref2 = 1.0;

    // default Cm reference location
    xf.xcmref = 0.25;
    xf.ycmref = 0.0;

    // default viscous parameters
    xf.retyp = 1;
    xf.reinf1 = 0.0;
    xf.acrit = [9.0; crate::state::ISX];
    xf.xstrip = [1.0; crate::state::ISX];

    // initialize BL parameter constants
    blpini(xf);

    // set MINF, REINF, based on current CL-dependence
    mrcl(xf, 1.0);
    comset(xf);

    let _ = &mut xf.cl_alf;
    xf.cl_alf = 0.0;
    xf.cl_msq = 0.0;
}

/// Sets the buffer airfoil to the specified NACA 4- or 5-digit airfoil and
/// panels it.
pub fn naca(xf: &mut Xfoil, ides1: i32) {
    // number of points per side
    let nside = IQX / 3;

    let mut ides = ides1;

    let mut itype = 0;
    if ides <= 25099 {
        itype = 5;
    }
    if ides <= 9999 {
        itype = 4;
    }

    if itype == 0 {
        if xf.show_output {
            eprintln!("This designation not implemented.");
        }
        return;
    }

    let (xb, yb, nb, name) = if itype == 4 {
        naca4(ides, nside)
    } else {
        naca5(ides, nside, xf.show_output)
    };

    if nb == 0 {
        return; // naca5 illegal designation
    }

    let _ = &mut ides;

    xf.name = name;

    xf.lclock = false;

    xf.xbf = 0.0;
    xf.ybf = 0.0;
    xf.lbflap = false;

    xf.nb = nb;
    for i in 0..nb {
        xf.xb[i] = xb[i];
        xf.yb[i] = yb[i];
    }

    scalc(&xf.xb[..nb], &xf.yb[..nb], &mut xf.sb[..nb]);
    segspl(&xf.xb[..nb], &mut xf.xbp[..nb], &xf.sb[..nb]);
    segspl(&xf.yb[..nb], &mut xf.ybp[..nb], &xf.sb[..nb]);

    geopar(
        &xf.xb[..nb],
        &xf.xbp[..nb],
        &xf.yb[..nb],
        &xf.ybp[..nb],
        &xf.sb[..nb],
        &mut xf.w1[..nb],
        &mut xf.sble,
        &mut xf.chordb,
        &mut xf.areab,
        &mut xf.radble,
        &mut xf.angbte,
        &mut xf.ei11ba,
        &mut xf.ei22ba,
        &mut xf.apx1ba,
        &mut xf.apx2ba,
        &mut xf.ei11bt,
        &mut xf.ei22bt,
        &mut xf.apx1bt,
        &mut xf.apx2bt,
        &mut xf.thickb,
        &mut xf.cambrb,
        xf.show_output,
    );

    if xf.show_output {
        eprintln!();
        eprintln!(" Buffer airfoil set using {:4} points", nb);
    }

    // set paneling
    pangen(xf, true);
}

/// Sets the paneling distribution from the buffer airfoil geometry, thus
/// creating the current airfoil.  If `shopar` is true, paneling parameters
/// are printed.
pub fn pangen(xf: &mut Xfoil, shopar: bool) {
    let nb = xf.nb;

    if nb < 2 {
        if xf.show_output {
            eprintln!("PANGEN: Buffer airfoil not available.");
        }
        xf.n = 0;
        return;
    }

    // number of temporary nodes for panel distribution calculation
    let ipfac = 5;

    // number of airfoil panel points
    xf.n = xf.npan;

    // set arc length spline parameter
    scalc(&xf.xb[..nb], &xf.yb[..nb], &mut xf.sb[..nb]);

    // spline raw airfoil coordinates
    segspl(&xf.xb[..nb], &mut xf.xbp[..nb], &xf.sb[..nb]);
    segspl(&xf.yb[..nb], &mut xf.ybp[..nb], &xf.sb[..nb]);

    // normalizing length (~ chord)
    let sbref = 0.5 * (xf.sb[nb - 1] - xf.sb[0]);

    // set up curvature array
    for i in 0..nb {
        xf.w5[i] = curv(xf.sb[i], &xf.xb[..nb], &xf.xbp[..nb], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]).abs() * sbref;
    }

    // locate LE point arc length value and the normalized curvature there
    lefind(
        &mut xf.sble,
        &xf.xb[..nb],
        &xf.xbp[..nb],
        &xf.yb[..nb],
        &xf.ybp[..nb],
        &xf.sb[..nb],
        xf.show_output,
    );
    let cvle = curv(xf.sble, &xf.xb[..nb], &xf.xbp[..nb], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]).abs() * sbref;

    // check for doubled point (sharp corner) at LE
    let mut ible = 0usize; // 1-based index into SB
    for i in 1..nb {
        if xf.sble == xf.sb[i - 1] && xf.sble == xf.sb[i] {
            ible = i;
            if xf.show_output {
                eprintln!();
                eprintln!("Sharp leading edge");
            }
            break;
        }
    }

    // set LE, TE points
    let xble = seval(xf.sble, &xf.xb[..nb], &xf.xbp[..nb], &xf.sb[..nb]);
    let yble = seval(xf.sble, &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]);
    let xbte = 0.5 * (xf.xb[0] + xf.xb[nb - 1]);
    let ybte = 0.5 * (xf.yb[0] + xf.yb[nb - 1]);
    let chbsq = (xbte - xble).powi(2) + (ybte - yble).powi(2);

    // set average curvature over 2*NK+1 points within Rcurv of LE point
    let nk = 3i32;
    let mut cvsum = 0.0;
    for k in -nk..=nk {
        let frac = k as f64 / nk as f64;
        let sbk = xf.sble + frac * sbref / cvle.max(20.0);
        let cvk = curv(sbk, &xf.xb[..nb], &xf.xbp[..nb], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]).abs() * sbref;
        cvsum += cvk;
    }
    let mut cvavg = cvsum / (2 * nk + 1) as f64;

    // dummy curvature for sharp LE
    if ible != 0 {
        cvavg = 10.0;
    }

    // set curvature attraction coefficient actually used
    let cc = 6.0 * xf.cvpar;

    // set artificial curvature at TE to bunch panels there
    let cvte = cvavg * xf.cterat;
    xf.w5[0] = cvte;
    xf.w5[nb - 1] = cvte;

    // set smoothing length = 1 / averaged LE curvature, but no more than 5%
    // of chord and no less than 1/4 average panel spacing
    let smool = (1.0 / cvavg.max(20.0)).max(0.25 / (xf.npan / 2) as f64);

    let smoosq = (smool * sbref).powi(2);

    // set up tri-diagonal system for smoothed curvatures
    xf.w2[0] = 1.0;
    xf.w3[0] = 0.0;
    for i in 1..nb - 1 {
        let dsm = xf.sb[i] - xf.sb[i - 1];
        let dsp = xf.sb[i + 1] - xf.sb[i];
        let dso = 0.5 * (xf.sb[i + 1] - xf.sb[i - 1]);

        if dsm == 0.0 || dsp == 0.0 {
            // leave curvature at corner point unchanged
            xf.w1[i] = 0.0;
            xf.w2[i] = 1.0;
            xf.w3[i] = 0.0;
        } else {
            xf.w1[i] = smoosq * (-1.0 / dsm) / dso;
            xf.w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            xf.w3[i] = smoosq * (-1.0 / dsp) / dso;
        }
    }
    xf.w1[nb - 1] = 0.0;
    xf.w2[nb - 1] = 1.0;

    // fix curvature at LE point by modifying equations adjacent to LE
    for i in 1..nb - 1 {
        if xf.sb[i] == xf.sble || i == ible || i == ible + 1 {
            // if node falls right on LE point, fix curvature there
            xf.w1[i] = 0.0;
            xf.w2[i] = 1.0;
            xf.w3[i] = 0.0;
            xf.w5[i] = cvle;
        } else if xf.sb[i - 1] < xf.sble && xf.sb[i] > xf.sble {
            // modify equation at node just before LE point
            let dsm = xf.sb[i - 1] - xf.sb[i - 2];
            let dsp = xf.sble - xf.sb[i - 1];
            let dso = 0.5 * (xf.sble - xf.sb[i - 2]);

            xf.w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
            xf.w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            xf.w3[i - 1] = 0.0;
            xf.w5[i - 1] += smoosq * cvle / (dsp * dso);

            // modify equation at node just after LE point
            let dsm = xf.sb[i] - xf.sble;
            let dsp = xf.sb[i + 1] - xf.sb[i];
            let dso = 0.5 * (xf.sb[i + 1] - xf.sble);
            xf.w1[i] = 0.0;
            xf.w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            xf.w3[i] = smoosq * (-1.0 / dsp) / dso;
            xf.w5[i] += smoosq * cvle / (dsm * dso);

            break;
        }
    }

    // set artificial curvature at bunching points and fix it there
    for i in 1..nb - 1 {
        // chord-based x/c coordinate
        let xoc = ((xf.xb[i] - xble) * (xbte - xble) + (xf.yb[i] - yble) * (ybte - yble)) / chbsq;

        if xf.sb[i] < xf.sble {
            // check if top side point is in refinement area
            if xoc > xf.xsref1 && xoc < xf.xsref2 {
                xf.w1[i] = 0.0;
                xf.w2[i] = 1.0;
                xf.w3[i] = 0.0;
                xf.w5[i] = cvle * xf.ctrrat;
            }
        } else if xoc > xf.xpref1 && xoc < xf.xpref2 {
            // check if bottom side point is in refinement area
            xf.w1[i] = 0.0;
            xf.w2[i] = 1.0;
            xf.w3[i] = 0.0;
            xf.w5[i] = cvle * xf.ctrrat;
        }
    }

    // solve for smoothed curvature array W5
    if ible == 0 {
        trisol(&mut xf.w2[..nb], &xf.w1[..nb], &mut xf.w3[..nb], &mut xf.w5[..nb]);
    } else {
        trisol(&mut xf.w2[..ible], &xf.w1[..ible], &mut xf.w3[..ible], &mut xf.w5[..ible]);
        trisol(
            &mut xf.w2[ible..nb],
            &xf.w1[ible..nb],
            &mut xf.w3[ible..nb],
            &mut xf.w5[ible..nb],
        );
    }

    // find max curvature
    let mut cvmax = 0.0f64;
    for i in 0..nb {
        cvmax = cvmax.max(xf.w5[i].abs());
    }

    // normalize curvature array
    for i in 0..nb {
        xf.w5[i] /= cvmax;
    }

    // spline curvature array
    segspl(&xf.w5[..nb], &mut xf.w6[..nb], &xf.sb[..nb]);

    // set initial guess for node positions uniform in s.  More nodes than
    // specified (by factor of IPFAC) are temporarily used for more reliable
    // convergence.
    let nn = ipfac * (xf.n - 1) + 1;

    // ratio of lengths of panel at TE to one away from the TE
    let rdste = 0.667;
    let rtf = (rdste - 1.0) * 2.0 + 1.0;

    let mut nn1 = 0;
    if ible == 0 {
        let dsavg = (xf.sb[nb - 1] - xf.sb[0]) / (nn as f64 - 3.0 + 2.0 * rtf);
        xf.snew[0] = xf.sb[0];
        for i in 2..=nn - 1 {
            xf.snew[i - 1] = xf.sb[0] + dsavg * (i as f64 - 2.0 + rtf);
        }
        xf.snew[nn - 1] = xf.sb[nb - 1];
    } else {
        let nfrac1 = (xf.n * ible) / nb;

        nn1 = ipfac * (nfrac1 - 1) + 1;
        let dsavg1 = (xf.sble - xf.sb[0]) / (nn1 as f64 - 2.0 + rtf);
        xf.snew[0] = xf.sb[0];
        for i in 2..=nn1 {
            xf.snew[i - 1] = xf.sb[0] + dsavg1 * (i as f64 - 2.0 + rtf);
        }

        let nn2 = nn - nn1 + 1;
        let dsavg2 = (xf.sb[nb - 1] - xf.sble) / (nn2 as f64 - 2.0 + rtf);
        for i in 2..nn2 - 1 {
            xf.snew[i - 1 + nn1 - 1] = xf.sble + dsavg2 * (i as f64 - 2.0 + rtf);
        }
        xf.snew[nn - 1] = xf.sb[nb - 1];
    }

    // Newton iteration loop for new node positions
    let mut converged = false;
    for _iter in 1..=20 {
        // set up tri-diagonal system for node position deltas
        let mut cv1 = seval(xf.snew[0], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);
        let mut cv2 = seval(xf.snew[1], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);
        let mut cvs1 = deval(xf.snew[0], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);
        let mut cvs2 = deval(xf.snew[1], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);

        let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
        let (mut cavm_s1, mut cavm_s2);
        if cavm == 0.0 {
            cavm_s1 = 0.0;
            cavm_s2 = 0.0;
        } else {
            cavm_s1 = cvs1 * cv1 / cavm;
            cavm_s2 = cvs2 * cv2 / cavm;
        }

        for i in 1..nn - 1 {
            let dsm = xf.snew[i] - xf.snew[i - 1];
            let dsp = xf.snew[i] - xf.snew[i + 1];
            let cv3 = seval(xf.snew[i + 1], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);
            let cvs3 = deval(xf.snew[i + 1], &xf.w5[..nb], &xf.w6[..nb], &xf.sb[..nb]);

            let cavp = (cv3 * cv3 + cv2 * cv2).sqrt();
            let (cavp_s2, cavp_s3);
            if cavp == 0.0 {
                cavp_s2 = 0.0;
                cavp_s3 = 0.0;
            } else {
                cavp_s2 = cvs2 * cv2 / cavp;
                cavp_s3 = cvs3 * cv3 / cavp;
            }

            let fm = cc * cavm + 1.0;
            let fp = cc * cavp + 1.0;

            let rez = dsp * fp + dsm * fm;

            // lower, main, and upper diagonals
            xf.w1[i] = -fm + cc * dsm * cavm_s1;
            xf.w2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            xf.w3[i] = -fp + cc * dsp * cavp_s3;

            // residual, requiring that (1 + C*curv)*deltaS is equal on both sides of node i
            xf.w4[i] = -rez;

            cv1 = cv2;
            cv2 = cv3;
            cvs1 = cvs2;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        // fix endpoints (at TE)
        xf.w2[0] = 1.0;
        xf.w3[0] = 0.0;
        xf.w4[0] = 0.0;
        xf.w1[nn - 1] = 0.0;
        xf.w2[nn - 1] = 1.0;
        xf.w4[nn - 1] = 0.0;

        if rtf != 1.0 {
            // fudge equations adjacent to TE to get TE panel length ratio RTF
            let i = 1; // 0-based (Fortran i=2)
            xf.w4[i] = -((xf.snew[i] - xf.snew[i - 1]) + rtf * (xf.snew[i] - xf.snew[i + 1]));
            xf.w1[i] = -1.0;
            xf.w2[i] = 1.0 + rtf;
            xf.w3[i] = -rtf;

            let i = nn - 2; // 0-based (Fortran i=nn-1)
            xf.w4[i] = -((xf.snew[i] - xf.snew[i + 1]) + rtf * (xf.snew[i] - xf.snew[i - 1]));
            xf.w3[i] = -1.0;
            xf.w2[i] = 1.0 + rtf;
            xf.w1[i] = -rtf;
        }

        // fix sharp LE point
        if ible != 0 {
            let i = nn1 - 1; // 0-based (Fortran i=nn1)
            xf.w1[i] = 0.0;
            xf.w2[i] = 1.0;
            xf.w3[i] = 0.0;
            xf.w4[i] = xf.sble - xf.snew[i];
        }

        // solve for changes W4 in node position arc length values
        trisol(&mut xf.w2[..nn], &xf.w1[..nn], &mut xf.w3[..nn], &mut xf.w4[..nn]);

        // find under-relaxation factor to keep nodes from changing order
        let mut rlx = 1.0;
        let mut dmax = 0.0f64;
        for i in 0..nn - 1 {
            let ds = xf.snew[i + 1] - xf.snew[i];
            let dds = xf.w4[i + 1] - xf.w4[i];
            let dsrat = 1.0 + rlx * dds / ds;
            if dsrat > 4.0 {
                rlx = (4.0 - 1.0) * ds / dds;
            }
            if dsrat < 0.2 {
                rlx = (0.2 - 1.0) * ds / dds;
            }
            dmax = dmax.max(xf.w4[i].abs());
        }

        // update node position
        for i in 1..nn - 1 {
            xf.snew[i] += rlx * xf.w4[i];
        }

        if dmax.abs() < 1.0E-3 {
            converged = true;
            break;
        }
    }
    if !converged && xf.show_output {
        eprintln!("Paneling convergence failed.  Continuing anyway...");
    }

    // set new panel node coordinates
    for i in 0..xf.n {
        let ind = ipfac * i;
        xf.s[i] = xf.snew[ind];
        xf.x[i] = seval(xf.snew[ind], &xf.xb[..nb], &xf.xbp[..nb], &xf.sb[..nb]);
        xf.y[i] = seval(xf.snew[ind], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]);
    }

    // go over buffer airfoil again, checking for corners (double points)
    let mut n = xf.n;
    for ib in 1..nb {
        if xf.sb[ib - 1] == xf.sb[ib] {
            // found one!
            let xbcorn = xf.xb[ib - 1];
            let ybcorn = xf.yb[ib - 1];
            let sbcorn = xf.sb[ib - 1];

            // find current-airfoil panel which contains corner
            for i in 0..n {
                // keep stepping until first node past corner
                if xf.s[i] > sbcorn {
                    // move remainder of panel nodes to make room for additional node
                    for j in (i..n).rev() {
                        xf.x[j + 1] = xf.x[j];
                        xf.y[j + 1] = xf.y[j];
                        xf.s[j + 1] = xf.s[j];
                    }
                    n += 1;

                    if n > IQX - 1 {
                        panic!("PANEL: Too many panels. Increase IQX in XFOIL.INC");
                    }

                    xf.x[i] = xbcorn;
                    xf.y[i] = ybcorn;
                    xf.s[i] = sbcorn;

                    // shift nodes adjacent to corner to keep panel sizes comparable
                    if i >= 2 {
                        xf.s[i - 1] = 0.5 * (xf.s[i] + xf.s[i - 2]);
                        xf.x[i - 1] = seval(xf.s[i - 1], &xf.xb[..nb], &xf.xbp[..nb], &xf.sb[..nb]);
                        xf.y[i - 1] = seval(xf.s[i - 1], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]);
                    }

                    if i + 2 < n {
                        xf.s[i + 1] = 0.5 * (xf.s[i] + xf.s[i + 2]);
                        xf.x[i + 1] = seval(xf.s[i + 1], &xf.xb[..nb], &xf.xbp[..nb], &xf.sb[..nb]);
                        xf.y[i + 1] = seval(xf.s[i + 1], &xf.yb[..nb], &xf.ybp[..nb], &xf.sb[..nb]);
                    }

                    // go on to next input geometry point to check for corner
                    break;
                }
            }
        }
    }
    xf.n = n;

    scalc(&xf.x[..xf.n], &xf.y[..xf.n], &mut xf.s[..xf.n]);
    segspl(&xf.x[..xf.n], &mut xf.xp[..xf.n], &xf.s[..xf.n]);
    segspl(&xf.y[..xf.n], &mut xf.yp[..xf.n], &xf.s[..xf.n]);
    lefind(
        &mut xf.sle,
        &xf.x[..xf.n],
        &xf.xp[..xf.n],
        &xf.y[..xf.n],
        &xf.yp[..xf.n],
        &xf.s[..xf.n],
        xf.show_output,
    );

    xf.xle = seval(xf.sle, &xf.x[..xf.n], &xf.xp[..xf.n], &xf.s[..xf.n]);
    xf.yle = seval(xf.sle, &xf.y[..xf.n], &xf.yp[..xf.n], &xf.s[..xf.n]);
    xf.xte = 0.5 * (xf.x[0] + xf.x[xf.n - 1]);
    xf.yte = 0.5 * (xf.y[0] + xf.y[xf.n - 1]);
    xf.chord = ((xf.xte - xf.xle).powi(2) + (xf.yte - xf.yle).powi(2)).sqrt();

    // set various flags for new airfoil
    xf.lgamu = false;
    xf.lqinu = false;
    xf.lwake = false;
    xf.lqaij = false;
    xf.ladij = false;
    xf.lwdij = false;
    xf.lipan = false;
    xf.lbli_ni = false;
    xf.lvconv = false;
    xf.lscini = false;
    xf.lgsame = false;

    if xf.lbflap {
        xf.xof = xf.xbf;
        xf.yof = xf.ybf;
        xf.lflap = true;
    }

    // determine if TE is blunt or sharp, calculate TE geometry parameters
    tecalc(xf);

    // calculate normal vectors
    ncalc(&xf.x[..xf.n], &xf.y[..xf.n], &xf.s[..xf.n], &mut xf.nx[..xf.n], &mut xf.ny[..xf.n]);

    // calculate panel angles for panel routines
    apcalc(xf);

    if xf.show_output {
        if xf.sharp {
            eprintln!();
            eprintln!("Sharp trailing edge");
        } else {
            let gap = ((xf.x[0] - xf.x[xf.n - 1]).powi(2) + (xf.y[0] - xf.y[xf.n - 1]).powi(2)).sqrt();
            eprintln!();
            eprintln!("Blunt trailing edge.  Gap = {:9.5}", gap);
        }
    }

    if shopar && xf.show_output {
        eprintln!();
        eprintln!(" Paneling parameters used...");
        eprintln!("   Number of panel nodes      {:4}", xf.npan);
        eprintln!("   Panel bunching parameter   {:6.3}", xf.cvpar);
        eprintln!("   TE/LE panel density ratio  {:6.3}", xf.cterat);
        eprintln!("   Refined-area/LE panel density ratio   {:6.3}", xf.ctrrat);
        eprintln!("   Top    side refined area x/c limits {:6.3}{:6.3}", xf.xsref1, xf.xsref2);
        eprintln!("   Bottom side refined area x/c limits {:6.3}{:6.3}", xf.xpref1, xf.xpref2);
    }

    let _ = &mut cang_placeholder(xf);
}

// cang is kept for API completeness (max panel corner angle check)
fn cang_placeholder(xf: &mut Xfoil) -> (f64, usize) {
    let mut amax = 0.0;
    let mut imax = 0;
    cang(&xf.x[..xf.n], &xf.y[..xf.n], 0, &mut amax, &mut imax, xf.show_output);
    (amax, imax)
}
