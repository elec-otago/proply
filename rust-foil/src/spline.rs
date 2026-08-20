// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Spline utilities (port of `m_spline.f90`).
//!
//! All routines work on slices; the number of points is the slice length.
//! Fortran 1-based indexing maps to Rust 0-based slices: `S(i)` -> `s[i-1]`.

/// Calculates spline coefficients for X(S) with zero 2nd derivative end
/// conditions. `xs` (dX/dS) is computed.
pub fn spline(x: &[f64], xs: &mut [f64], s: &[f64]) {
    let n = x.len();
    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];

    for i in 1..n - 1 {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }

    // set zero second derivative end conditions
    a[0] = 2.0;
    c[0] = 1.0;
    xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);
    b[n - 1] = 1.0;
    a[n - 1] = 2.0;
    xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);

    trisol(&mut a, &b, &mut c, xs);
}

/// Calculates spline coefficients for X(S) with specified endpoint derivatives.
///
/// `xs1`, `xs2` endpoint derivatives; `999.0` selects the usual zero 2nd
/// derivative condition, `-999.0` selects zero 3rd derivative condition.
pub fn splind(x: &[f64], xs: &mut [f64], s: &[f64], xs1: f64, xs2: f64) {
    let n = x.len();
    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];

    for i in 1..n - 1 {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }

    if xs1 == 999.0 {
        // zero second derivative end condition
        a[0] = 2.0;
        c[0] = 1.0;
        xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);
    } else if xs1 == -999.0 {
        // zero third derivative end condition
        a[0] = 1.0;
        c[0] = 1.0;
        xs[0] = 2.0 * (x[1] - x[0]) / (s[1] - s[0]);
    } else {
        // specified first derivative end condition
        a[0] = 1.0;
        c[0] = 0.0;
        xs[0] = xs1;
    }

    if xs2 == 999.0 {
        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    } else if xs2 == -999.0 {
        b[n - 1] = 1.0;
        a[n - 1] = 1.0;
        xs[n - 1] = 2.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    } else {
        a[n - 1] = 1.0;
        b[n - 1] = 0.0;
        xs[n - 1] = xs2;
    }

    if n == 2 && xs1 == -999.0 && xs2 == -999.0 {
        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    }

    trisol(&mut a, &b, &mut c, xs);
}

/// Spline with simple averaging of adjacent segment slopes (non-oscillatory).
pub fn splina(x: &[f64], xs: &mut [f64], s: &[f64]) {
    let n = x.len();
    let mut xs1 = 0.0;
    let mut xs2 = 0.0;
    let mut lend = true;

    for i in 0..n - 1 {
        let ds = s[i + 1] - s[i];
        if ds == 0.0 {
            xs[i] = xs1;
            lend = true;
        } else {
            let dx = x[i + 1] - x[i];
            xs2 = dx / ds;
            if lend {
                xs[i] = xs2;
                lend = false;
            } else {
                xs[i] = 0.5 * (xs1 + xs2);
            }
        }
        xs1 = xs2;
    }
    xs[n - 1] = xs1;
}

/// Solves a KK long tridiagonal system.  The righthand side `d` is replaced
/// by the solution; `a`, `c` are destroyed.
pub fn trisol(a: &mut [f64], b: &[f64], c: &mut [f64], d: &mut [f64]) {
    let kk = a.len();
    for k in 1..kk {
        let km = k - 1;
        c[km] /= a[km];
        d[km] /= a[km];
        a[k] -= b[k] * c[km];
        d[k] -= b[k] * d[km];
    }

    d[kk - 1] /= a[kk - 1];

    for k in (0..kk - 1).rev() {
        d[k] -= c[k] * d[k + 1];
    }
}

/// Returns the index `i` (0-based) of the bracketing interval of `ss`, i.e.
/// the interval [s[i-1], s[i]] containing ss.
fn locate(ss: f64, s: &[f64]) -> usize {
    let n = s.len();
    let mut ilow = 0usize;
    let mut i = n - 1;
    while i - ilow > 1 {
        let imid = (i + ilow) / 2;
        if ss < s[imid] {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    i
}

/// Calculates X(SS).  `xs` must have been calculated by [`spline`].
pub fn seval(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = locate(ss, s);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
}

/// Calculates dX/dS(SS).  `xs` must have been calculated by [`spline`].
pub fn deval(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = locate(ss, s);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let deval = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    deval / ds
}

/// Calculates d2X/dS2(SS).  `xs` must have been calculated by [`spline`].
pub fn d2val(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = locate(ss, s);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let d2val = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    d2val / (ds * ds)
}

/// Calculates curvature of splined 2-D curve at S = SS.
pub fn curv(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
    let i = locate(ss, s);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;

    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

    let mut sd = (xd * xd + yd * yd).sqrt();
    sd = sd.max(0.001 * ds);

    (xd * ydd - yd * xdd) / (sd * sd * sd)
}

/// Calculates curvature derivative of splined 2-D curve at S = SS.
pub fn curvs(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
    let i = locate(ss, s);
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;

    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    let xddd = 6.0 * cx1 + 6.0 * cx2;

    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
    let yddd = 6.0 * cy1 + 6.0 * cy2;

    let mut sd = (xd * xd + yd * yd).sqrt();
    sd = sd.max(0.001 * ds);

    let bot = sd * sd * sd;
    let dbotdt = 3.0 * sd * (xd * xdd + yd * ydd);

    let top = xd * ydd - yd * xdd;
    let dtopdt = xd * yddd - yd * xddd;

    (dtopdt * bot - dbotdt * top) / (bot * bot)
}

/// Calculates the "inverse" spline function S(X).  `si` must be a
/// sufficiently good initial guess (input/output).
pub fn sinvrt(si: &mut f64, xi: f64, x: &[f64], xs: &[f64], s: &[f64], show_output: bool) {
    let sisav = *si;

    for _ in 1..=10 {
        let res = seval(*si, x, xs, s) - xi;
        let resp = deval(*si, x, xs, s);
        let ds = -res / resp;
        *si += ds;
        if (ds / (s[s.len() - 1] - s[0])).abs() < 1.0E-5 {
            return;
        }
    }
    if show_output {
        eprintln!("SINVRT: spline inversion failed. Input value returned.");
    }
    *si = sisav;
}

/// Calculates the arc length array S for a 2-D array of points (X, Y).
pub fn scalc(x: &[f64], y: &[f64], s: &mut [f64]) {
    let n = x.len();
    s[0] = 0.0;
    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
}

/// Splines 2-D shape X(S), Y(S), along with true arc length parameter S.
pub fn splnxy(x: &[f64], xs: &mut [f64], y: &[f64], ys: &mut [f64], s: &mut [f64], show_output: bool) {
    const KMAX: usize = 32;
    let n = x.len();
    let kk = KMAX;
    let npass = 10;

    // set first estimate of arc length parameter
    scalc(x, y, s);

    // spline X(S) and Y(S)
    segspl(x, xs, s);
    segspl(y, ys, s);

    // re-integrate true arc length
    let mut xt = [0.0f64; KMAX + 1];
    let mut yt = [0.0f64; KMAX + 1];

    for _ipass in 1..=npass {
        let mut serr = 0.0f64;

        let mut ds = s[1] - s[0];
        for i in 1..n {
            let dx = x[i] - x[i - 1];
            let dy = y[i] - y[i - 1];

            let cx1 = ds * xs[i - 1] - dx;
            let cx2 = ds * xs[i] - dx;
            let cy1 = ds * ys[i - 1] - dy;
            let cy2 = ds * ys[i] - dy;

            xt[0] = 0.0;
            yt[0] = 0.0;
            for k in 1..kk {
                let t = k as f64 / kk as f64;
                xt[k] = t * dx + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
                yt[k] = t * dy + (t - t * t) * ((1.0 - t) * cy1 - t * cy2);
            }
            xt[kk] = dx;
            yt[kk] = dy;

            let mut sint1 = 0.0;
            for k in 1..=kk {
                let ddx = xt[k] - xt[k - 1];
                let ddy = yt[k] - yt[k - 1];
                sint1 += (ddx * ddx + ddy * ddy).sqrt();
            }

            let mut sint2 = 0.0;
            let mut k = 2;
            while k <= kk {
                let ddx = xt[k] - xt[k - 2];
                let ddy = yt[k] - yt[k - 2];
                sint2 += (ddx * ddx + ddy * ddy).sqrt();
                k += 2;
            }

            let sint = (4.0 * sint1 - sint2) / 3.0;

            if (sint - ds).abs() > serr.abs() {
                serr = sint - ds;
            }

            if i < n - 1 {
                ds = s[i + 1] - s[i];
            }

            s[i] = s[i - 1] + sint.sqrt();
        }

        serr /= s[n - 1] - s[0];
        if show_output {
            eprintln!("  {}  {:e}", _ipass, serr);
        }

        // re-spline X(S) and Y(S)
        segspl(x, xs, s);
        segspl(y, ys, s);

        if serr.abs() < 1.0E-7 {
            return;
        }
    }
}

/// Splines X(S) like [`spline`], allowing derivative discontinuities at
/// segment joints (defined by identical successive S values).
pub fn segspl(x: &[f64], xs: &mut [f64], s: &[f64]) {
    let n = x.len();
    if s[0] == s[1] {
        panic!("SEGSPL: First input point duplicated");
    }
    if s[n - 1] == s[n - 2] {
        panic!("SEGSPL: Last input point duplicated");
    }

    let mut iseg0 = 0usize;
    for iseg in 1..n - 2 {
        if s[iseg] == s[iseg + 1] {
            let nseg = iseg - iseg0 + 1;
            splind(&x[iseg0..], &mut xs[iseg0..], &s[iseg0..], -999.0, -999.0);
            iseg0 = iseg + 1;
        }
    }

    let nseg = n - iseg0;
    splind(&x[iseg0..], &mut xs[iseg0..], &s[iseg0..], -999.0, -999.0);
    let _ = nseg;
}

/// Like [`segspl`] with specified endpoint derivatives.
pub fn segspld(x: &[f64], xs: &mut [f64], s: &[f64], xs1: f64, xs2: f64) {
    let n = x.len();
    if s[0] == s[1] {
        panic!("SEGSPL: First input point duplicated");
    }
    if s[n - 1] == s[n - 2] {
        panic!("SEGSPL: Last input point duplicated");
    }

    let mut iseg0 = 0usize;
    for iseg in 1..n - 2 {
        if s[iseg] == s[iseg + 1] {
            let nseg = iseg - iseg0 + 1;
            splind(&x[iseg0..], &mut xs[iseg0..], &s[iseg0..], xs1, xs2);
            iseg0 = iseg + 1;
        }
    }

    splind(&x[iseg0..], &mut xs[iseg0..], &s[iseg0..], xs1, xs2);
}
