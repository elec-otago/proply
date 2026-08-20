// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Airfoil geometry routines (port of `m_xgeom.f90`).
//!
//! Only the routines used by the analysis path are ported (lefind, sopps,
//! norm, geopar, aecalc, tccalc, cang).  Design/flap/plot utilities are
//! omitted.

use crate::spline::{curv, d2val, deval, scalc, segspl, seval};
use crate::utils::atanc;

/// Locates leading edge spline-parameter value `sle`.
///
/// The defining condition is that the surface tangent is normal to the chord
/// line connecting X(SLE),Y(SLE) and the TE point.
pub fn lefind(
    sle: &mut f64,
    x: &[f64],
    xp: &[f64],
    y: &[f64],
    yp: &[f64],
    s: &[f64],
    show_output: bool,
) {
    let n = s.len();

    // convergence tolerance
    let dseps = (s[n - 1] - s[0]) * 1.0E-5;

    // set trailing edge point coordinates
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);

    // get first guess for SLE
    let mut i = 2usize; // Fortran: i from 3..N-2 inclusive -> Rust 2..=n-3
    let mut dotp = 0.0;
    loop {
        if i > n - 3 {
            break;
        }
        let dxte = x[i] - xte;
        let dyte = y[i] - yte;
        let dx = x[i + 1] - x[i];
        let dy = y[i + 1] - y[i];
        dotp = dxte * dx + dyte * dy;
        if dotp < 0.0 {
            break;
        }
        i += 1;
    }

    *sle = s[i];

    // check for sharp LE case
    if s[i] == s[i - 1] {
        return;
    }

    // Newton iteration to get exact SLE value
    for _ in 1..=50 {
        let xle = seval(*sle, x, xp, s);
        let yle = seval(*sle, y, yp, s);
        let dxds = deval(*sle, x, xp, s);
        let dyds = deval(*sle, y, yp, s);
        let dxdd = d2val(*sle, x, xp, s);
        let dydd = d2val(*sle, y, yp, s);

        let xchord = xle - xte;
        let ychord = yle - yte;

        // drive dot product between chord line and LE tangent to zero
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        // Newton delta for SLE
        let mut dsle = -res / ress;

        dsle = dsle.max(-0.02 * (xchord + ychord).abs());
        dsle = dsle.min(0.02 * (xchord + ychord).abs());
        *sle += dsle;
        if dsle.abs() < dseps {
            return;
        }
    }
    if show_output {
        eprintln!("LEFIND:  LE point not found.  Continuing...");
    }
    *sle = s[i];
}

/// Calculates arc length `sopp` of the point which is opposite of point `si`,
/// on the other side of the airfoil baseline.
#[allow(clippy::too_many_arguments)] // XFOIL port signature
pub fn sopps(
    sopp: &mut f64,
    si: f64,
    x: &[f64],
    xp: &[f64],
    y: &[f64],
    yp: &[f64],
    s: &[f64],
    sle: f64,
    show_output: bool,
) {
    let n = s.len();

    // reference length for testing convergence
    let slen = s[n - 1] - s[0];

    // set chordline vector
    let xle = seval(sle, x, xp, s);
    let yle = seval(sle, y, yp, s);
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;

    let (inn, inopp);
    if si < sle {
        inn = 0usize;
        inopp = n - 1;
    } else {
        inn = n - 1;
        inopp = 0;
    }
    let sfrac = (si - sle) / (s[inn] - sle);
    *sopp = sle + sfrac * (s[inopp] - sle);

    if sfrac.abs() <= 1.0E-5 {
        *sopp = sle;
        return;
    }

    // XBAR = x coordinate in chord-line axes
    let xi = seval(si, x, xp, s);
    let yi = seval(si, y, yp, s);
    let xle = seval(sle, x, xp, s);
    let yle = seval(sle, y, yp, s);
    let xbar = (xi - xle) * dxc + (yi - yle) * dyc;

    // converge on exact opposite point with same XBAR value
    for _ in 1..=12 {
        let xopp = seval(*sopp, x, xp, s);
        let yopp = seval(*sopp, y, yp, s);
        let xoppd = deval(*sopp, x, xp, s);
        let yoppd = deval(*sopp, y, yp, s);

        let res = (xopp - xle) * dxc + (yopp - yle) * dyc - xbar;
        let resd = xoppd * dxc + yoppd * dyc;

        if (res / slen).abs() < 1.0E-5 {
            return;
        }
        if resd == 0.0 {
            break;
        }

        let dsopp = -res / resd;
        *sopp += dsopp;

        if (dsopp / slen).abs() < 1.0E-5 {
            return;
        }
    }
    if show_output {
        eprintln!("SOPPS: Opposite-point location failed. Continuing...");
    }
    *sopp = sle + sfrac * (s[inopp] - sle);
}

/// Scales coordinates to get unit chord.
pub fn norm(x: &mut [f64], xp: &mut [f64], y: &mut [f64], yp: &mut [f64], s: &mut [f64]) {
    let n = x.len();

    scalc(x, y, s);
    segspl(x, xp, s);
    segspl(y, yp, s);

    let mut sle = 0.0;
    lefind(&mut sle, x, xp, y, yp, s, false);

    let xmax = 0.5 * (x[0] + x[n - 1]);
    let xmin = seval(sle, x, xp, s);
    let ymin = seval(sle, y, yp, s);

    let fudge = 1.0 / (xmax - xmin);
    for i in 0..n {
        x[i] = (x[i] - xmin) * fudge;
        y[i] = (y[i] - ymin) * fudge;
        s[i] *= fudge;
    }
}

/// Sets geometric parameters for airfoil shape (chord, LE radius, TE angle,
/// area, inertias, max thickness/camber).
#[allow(clippy::too_many_arguments)]
pub fn geopar(
    x: &[f64],
    xp: &[f64],
    y: &[f64],
    yp: &[f64],
    s: &[f64],
    t: &mut [f64],
    sle: &mut f64,
    chord: &mut f64,
    area: &mut f64,
    radle: &mut f64,
    angte: &mut f64,
    ei11a: &mut f64,
    ei22a: &mut f64,
    apx1a: &mut f64,
    apx2a: &mut f64,
    ei11t: &mut f64,
    ei22t: &mut f64,
    apx1t: &mut f64,
    apx2t: &mut f64,
    thick: &mut f64,
    cambr: &mut f64,
    show_output: bool,
) {
    let n = s.len();

    lefind(sle, x, xp, y, yp, s, show_output);

    let xle = seval(*sle, x, xp, s);
    let yle = seval(*sle, y, yp, s);
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);

    let chsq = (xte - xle).powi(2) + (yte - yle).powi(2);
    *chord = chsq.sqrt();

    let curvle = curv(*sle, x, xp, y, yp, s);

    *radle = 0.0;
    if curvle.abs() > 0.001 * (s[n - 1] - s[0]) {
        *radle = 1.0 / curvle;
    }

    let ang1 = (-yp[0]).atan2(-xp[0]);
    let ang2 = atanc(yp[n - 1], xp[n - 1], ang1);
    *angte = ang2 - ang1;

    for t in t.iter_mut() {
        *t = 1.0;
    }

    let mut xcena = 0.0;
    let mut ycena = 0.0;
    aecalc(
        x, y, t, 1, area, &mut xcena, &mut ycena, ei11a, ei22a, apx1a, apx2a,
    );

    let mut slen = 0.0;
    let mut xcent = 0.0;
    let mut ycent = 0.0;
    aecalc(
        x, y, t, 2, &mut slen, &mut xcent, &mut ycent, ei11t, ei22t, apx1t, apx2t,
    );

    // Old, approximate thickness,camber routine (on discrete points only)
    tccalc(x, xp, y, yp, s, thick, cambr, show_output);

    if show_output {
        // NOTE: original prints thickness/camber maxima via tccalc outputs;
        // x-locations are not returned, so print only the values.
        eprintln!(" Max thickness = {:12.6}", *thick);
        eprintln!(" Max camber    = {:12.6}", *cambr);
    }
}

/// Calculates geometric properties of shape X,Y (centroid, inertias,
/// principal-axis angles).  `itype = 1` integrates over the whole area dx dy,
/// `itype = 2` integrates over skin area t ds.
#[allow(clippy::too_many_arguments)] // XFOIL port signature
pub fn aecalc(
    x: &[f64],
    y: &[f64],
    t: &[f64],
    itype: i32,
    area: &mut f64,
    xcen: &mut f64,
    ycen: &mut f64,
    ei11: &mut f64,
    ei22: &mut f64,
    apx1: &mut f64,
    apx2: &mut f64,
) {
    let n = x.len();
    const PI: f64 = std::f64::consts::PI;

    let mut sint = 0.0;
    let mut aint = 0.0;
    let mut xint = 0.0;
    let mut yint = 0.0;
    let mut xxint = 0.0;
    let mut xyint = 0.0;
    let mut yyint = 0.0;

    for io in 0..n {
        let ip = if io == n - 1 { 0 } else { io + 1 };

        let dx = x[io] - x[ip];
        let dy = y[io] - y[ip];
        let xa = (x[io] + x[ip]) * 0.50;
        let ya = (y[io] + y[ip]) * 0.50;
        let ta = (t[io] + t[ip]) * 0.50;

        let ds = (dx * dx + dy * dy).sqrt();
        sint += ds;

        let da;
        if itype == 1 {
            // integrate over airfoil cross-section
            da = ya * dx;
            aint += da;
            xint += xa * da;
            yint += ya * da / 2.0;
            xxint += xa * xa * da;
            xyint += xa * ya * da / 2.0;
            yyint += ya * ya * da / 3.0;
        } else {
            // integrate over skin thickness
            da = ta * ds;
            aint += da;
            xint += xa * da;
            yint += ya * da;
            xxint += xa * xa * da;
            xyint += xa * ya * da;
            yyint += ya * ya * da;
        }
    }

    *area = aint;

    if aint == 0.0 {
        *xcen = 0.0;
        *ycen = 0.0;
        *ei11 = 0.0;
        *ei22 = 0.0;
        *apx1 = 0.0;
        *apx2 = 1.0f64.atan2(0.0);
        return;
    }

    // calculate centroid location
    *xcen = xint / aint;
    *ycen = yint / aint;

    // calculate inertias
    let eixx = yyint - *ycen * *ycen * aint;
    let eixy = xyint - *xcen * *ycen * aint;
    let eiyy = xxint - *xcen * *xcen * aint;

    // set principal-axis inertias, EI11 is closest to "up-down" bending inertia
    let eisq = 0.25 * (eixx - eiyy).powi(2) + eixy * eixy;
    let sgn = (eiyy - eixx).signum();
    *ei11 = 0.5 * (eixx + eiyy) - sgn * eisq.sqrt();
    *ei22 = 0.5 * (eixx + eiyy) + sgn * eisq.sqrt();

    if *ei11 == 0.0 || *ei22 == 0.0 {
        // vanishing section stiffness
        *apx1 = 0.0;
        *apx2 = 1.0f64.atan2(0.0);
    } else if eisq / (*ei11 * *ei22) < (0.001 * sint).powi(4) {
        // rotationally-invariant section (circle, square, etc.)
        *apx1 = 0.0;
        *apx2 = 1.0f64.atan2(0.0);
    } else {
        // normal airfoil section
        let c1 = eixy;
        let s1 = eixx - *ei11;

        let c2 = eixy;
        let s2 = eixx - *ei22;

        if s1.abs() > s2.abs() {
            *apx1 = s1.atan2(c1);
            *apx2 = *apx1 + 0.5 * PI;
        } else {
            *apx2 = s2.atan2(c2);
            *apx1 = *apx2 - 0.5 * PI;
        }

        if *apx1 < -0.5 * PI {
            *apx1 += PI;
        }
        if *apx1 > 0.5 * PI {
            *apx1 -= PI;
        }
        if *apx2 < -0.5 * PI {
            *apx2 += PI;
        }
        if *apx2 > 0.5 * PI {
            *apx2 -= PI;
        }
    }
}

/// Calculates max thickness and camber at airfoil points (discrete-point
/// approximation).
#[allow(clippy::too_many_arguments)] // XFOIL port signature
pub fn tccalc(
    x: &[f64],
    xp: &[f64],
    y: &[f64],
    yp: &[f64],
    s: &[f64],
    thick: &mut f64,
    cambr: &mut f64,
    show_output: bool,
) {
    let n = s.len();

    let mut sle = 0.0;
    lefind(&mut sle, x, xp, y, yp, s, show_output);
    let xle = seval(sle, x, xp, s);
    let yle = seval(sle, y, yp, s);
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();

    // set unit chord-line vector
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;

    *thick = 0.0;
    *cambr = 0.0;

    // go over each point, finding the y-thickness and camber
    for i in 0..n {
        let xbar = (x[i] - xle) * dxc + (y[i] - yle) * dyc;
        let ybar = (y[i] - yle) * dxc - (x[i] - xle) * dyc;

        // set point on the opposite side with the same chord x value
        let mut sopp = 0.0;
        sopps(&mut sopp, s[i], x, xp, y, yp, s, sle, show_output);
        let xopp = seval(sopp, x, xp, s);
        let yopp = seval(sopp, y, yp, s);

        let ybarop = (yopp - yle) * dxc - (xopp - xle) * dyc;

        let yc = 0.5 * (ybar + ybarop);
        let yt = (ybar - ybarop).abs();

        if yc.abs() > cambr.abs() {
            *cambr = yc;
        }
        if yt.abs() > thick.abs() {
            *thick = yt;
        }
        let _ = xbar;
    }
}

/// Computes the maximum panel node corner angle (in degrees).
pub fn cang(
    x: &[f64],
    y: &[f64],
    iprint: i32,
    amax: &mut f64,
    imax: &mut usize,
    show_output: bool,
) {
    let n = x.len();

    *amax = 0.0;
    *imax = 0;

    for i in 1..n - 1 {
        let mut dx1 = x[i] - x[i - 1];
        let mut dy1 = y[i] - y[i - 1];
        let mut dx2 = x[i] - x[i + 1];
        let mut dy2 = y[i] - y[i + 1];

        // allow for doubled points
        if dx1 == 0.0 && dy1 == 0.0 {
            dx1 = x[i] - x[i - 2];
            dy1 = y[i] - y[i - 2];
        }
        if dx2 == 0.0 && dy2 == 0.0 {
            dx2 = x[i] - x[i + 2];
            dy2 = y[i] - y[i + 2];
        }

        let crossp =
            (dx2 * dy1 - dy2 * dx1) / ((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2)).sqrt();
        let angl = crossp.asin() * (180.0 / std::f64::consts::PI);
        if iprint == 2 && show_output {
            eprintln!("{:3} {:9.4} {:9.4} {:9.3}", i + 1, x[i], y[i], angl);
        }
        if angl.abs() > amax.abs() {
            *amax = angl;
            *imax = i;
        }
    }

    if iprint >= 1 && show_output {
        eprintln!(
            " Maximum panel corner angle = {:7.3}   at  i,x,y  = {:3} {:9.4} {:9.4}",
            *amax,
            *imax + 1,
            x[*imax],
            y[*imax]
        );
    }
}
