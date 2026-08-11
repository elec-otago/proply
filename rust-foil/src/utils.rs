//! Miscellaneous utilities (port of `m_xutils.f90`).

/// Sets a geometrically stretched array S:
///
/// ```text
/// S(i+1) - S(i) = r * [S(i) - S(i-1)]
/// ```
///
/// `ds1` is the first increment, `smax` the final value, `nn` the number of
/// points.
pub fn setexp(s: &mut [f64], ds1: f64, smax: f64, show_output: bool) {
    let nn = s.len();
    let sigma = smax / ds1;
    let nex = nn - 1;
    let rnex = nex as f64;
    let rni = 1.0 / rnex;

    // solve quadratic for initial geometric ratio guess
    let aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    let bbb = rnex * (rnex - 1.0) / 2.0;
    let ccc = rnex - sigma;

    let disc = (bbb * bbb - 4.0 * aaa * ccc).max(0.0);

    let mut ratio;
    if nex <= 1 {
        panic!("SETEXP: Cannot fill array.  N too small.");
    } else if nex == 2 {
        ratio = -ccc / bbb + 1.0;
    } else {
        ratio = (-bbb + disc.sqrt()) / (2.0 * aaa) + 1.0;
    }

    if ratio != 1.0 {
        // Newton iteration for actual geometric ratio
        let mut converged = false;
        for _ in 1..=100 {
            let sigman = (ratio.powf(nex as f64) - 1.0) / (ratio - 1.0);
            let res = sigman.powf(rni) - sigma.powf(rni);
            let dresdr =
                rni * sigman.powf(rni) * (rnex * ratio.powf((nex - 1) as f64) - sigman) / (ratio.powf(nex as f64) - 1.0);

            let dratio = -res / dresdr;
            ratio += dratio;

            if dratio.abs() < 1.0E-5 {
                converged = true;
                break;
            }
        }
        if !converged && show_output {
            eprintln!("SETEXP: Convergence failed.  Continuing anyway ...");
        }
    }

    // set up stretched array using converged geometric ratio
    s[0] = 0.0;
    let mut ds = ds1;
    for n in 1..nn {
        s[n] = s[n - 1] + ds;
        ds *= ratio;
    }
}

/// ATAN2 function with branch cut checking.
///
/// Increments position angle of point (x, y) from some previous value
/// `thold`, ensuring that the position change does not cross the ATAN2
/// branch cut (which is in the -x direction).
pub fn atanc(y: f64, x: f64, thold: f64) -> f64 {
    const PI: f64 = 3.1415926535897932384;
    const TPI: f64 = 6.2831853071795864769;

    // set new position angle, ignoring branch cut in ATAN2 for now
    let thnew = y.atan2(x);
    let dthet = thnew - thold;

    // angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi
    let dtcorr = dthet - TPI * ((dthet + PI * dthet.signum()) / TPI) as i64 as f64;

    thold + dtcorr
}
