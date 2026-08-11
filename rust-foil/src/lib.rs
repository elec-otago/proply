//! rust-foil: a Rust port of the XFOIL airfoil analysis program.
//!
//! The library computes inviscid and viscous (integral boundary layer) flow
//! solutions over arbitrary airfoil geometries using the panel method, and
//! mirrors the public API of the Python `xfoil-python` package.
//!
//! Several routines from the original Fortran are ported but not exercised by
//! the current analysis path (e.g. inverse-design helpers); these carry a
//! crate-wide dead-code allowance.  The Fortran sources also compute many
//! sensitivity derivatives that are only stored, never used; those produce
//! benign assignment warnings, silenced here.
#![allow(dead_code)]
#![allow(unused_assignments)]
#![allow(unused_mut)]
#![allow(unused_variables)]

mod bl;
mod blsys;
mod gdes;
mod geom;
mod naca;
mod oper;
mod panel;
mod s_xbl;
mod s_xfoil;
mod solve;
mod spline;
pub mod state;
mod utils;
mod xfoil;

pub use solve::blsolv;
pub use solve::gauss;
// LU factor/back-substitution are internal to the inviscid solve, but exposed
// (hidden from docs) so the integration tests can exercise them directly.
#[doc(hidden)]
pub use solve::baksub;
#[doc(hidden)]
pub use solve::ludcmp;

use rayon::prelude::*;

use state::Xfoil;

/// A thin, safe wrapper around the XFOIL engine state providing the
/// high-level airfoil analysis API.
#[derive(Clone)]
pub struct XFoil {
    state: Xfoil,
}

impl Default for XFoil {
    fn default() -> Self {
        Self::new()
    }
}

impl XFoil {
    /// Creates a new XFOIL instance with default settings.
    pub fn new() -> Self {
        let mut xf = Xfoil::new();
        xfoil::init(&mut xf);
        XFoil { state: xf }
    }

    /// Enables or disables console output.
    pub fn set_show_output(&mut self, setting: bool) {
        self.state.show_output = setting;
    }

    /// Control whether `set_airfoil` normalizes the input coordinates to unit
    /// chord (default on, matching canonical XFOIL).  The panel and viscous
    /// analysis assumes a unit-chord airfoil with the Reynolds number set
    /// separately, so chord-scaled input silently produces wrong polars when
    /// normalization is off.
    pub fn set_normalize(&mut self, on: bool) {
        self.state.lnorm = on;
    }

    /// Returns whether console output is enabled.
    pub fn get_show_output(&self) -> bool {
        self.state.show_output
    }

    /// Sets the airfoil from coordinate arrays (counterclockwise ordering is
    /// enforced; input coordinates become the panel nodes).
    ///
    /// Returns the LE point, chord and TE point of the stored buffer airfoil.
    pub fn set_airfoil(&mut self, x: &[f64], y: &[f64]) -> ([f64; 2], f64, [f64; 2]) {
        assert_eq!(x.len(), y.len(), "x and y arrays must have equal length");
        assert!(x.len() >= 3, "at least 3 airfoil coordinates required");
        assert!(
            x.len() <= state::IBX,
            "too many input coordinates (max {})",
            state::IBX
        );

        let n_in = x.len();
        let xf = &mut self.state;

        for i in 0..n_in {
            xf.xb[i] = x[i];
            xf.yb[i] = y[i];
        }
        xf.nb = n_in;

        // calculate airfoil area assuming counterclockwise ordering
        let mut area = 0.0;
        for i in 0..xf.nb {
            let ip = if i == xf.nb - 1 { 0 } else { i + 1 };
            area += 0.5 * (xf.yb[i] + xf.yb[ip]) * (xf.xb[i] - xf.xb[ip]);
        }

        if area >= 0.0 {
            xf.lclock = false;
            if xf.show_output {
                eprintln!();
                eprintln!(" Number of input coordinate points: {:4}", xf.nb);
                eprintln!(" Counterclockwise ordering");
            }
        } else {
            // if area is negative (clockwise order), reverse coordinate order
            xf.lclock = true;
            if xf.show_output {
                eprintln!();
                eprintln!(" Number of input coordinate points: {:4}", xf.nb);
                eprintln!(" Clockwise ordering");
            }
            for i in 0..xf.nb / 2 {
                let xtmp = xf.xb[xf.nb - i - 1];
                let ytmp = xf.yb[xf.nb - i - 1];
                xf.xb[xf.nb - i - 1] = xf.xb[i];
                xf.yb[xf.nb - i - 1] = xf.yb[i];
                xf.xb[i] = xtmp;
                xf.yb[i] = ytmp;
            }
        }

        if xf.lnorm {
            geom::norm(
                &mut xf.xb[..xf.nb],
                &mut xf.xbp[..xf.nb],
                &mut xf.yb[..xf.nb],
                &mut xf.ybp[..xf.nb],
                &mut xf.sb[..xf.nb],
            );
            if xf.show_output {
                eprintln!();
                eprintln!(" Airfoil has been normalized");
            }
        }

        spline::scalc(&xf.xb[..xf.nb], &xf.yb[..xf.nb], &mut xf.sb[..xf.nb]);
        spline::segspl(&xf.xb[..xf.nb], &mut xf.xbp[..xf.nb], &xf.sb[..xf.nb]);
        spline::segspl(&xf.yb[..xf.nb], &mut xf.ybp[..xf.nb], &xf.sb[..xf.nb]);

        geom::geopar(
            &xf.xb[..xf.nb],
            &xf.xbp[..xf.nb],
            &xf.yb[..xf.nb],
            &xf.ybp[..xf.nb],
            &xf.sb[..xf.nb],
            &mut xf.w1[..xf.nb],
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

        let xble = spline::seval(xf.sble, &xf.xb[..xf.nb], &xf.xbp[..xf.nb], &xf.sb[..xf.nb]);
        let yble = spline::seval(xf.sble, &xf.yb[..xf.nb], &xf.ybp[..xf.nb], &xf.sb[..xf.nb]);
        let xbte = 0.5 * (xf.xb[0] + xf.xb[xf.nb - 1]);
        let ybte = 0.5 * (xf.yb[0] + xf.yb[xf.nb - 1]);

        if xf.show_output {
            eprintln!();
            eprintln!("  LE  x,y  = {:10.5}{:10.5}  |   Chord = {:10.5}", xble, yble, xf.chordb);
            eprintln!("  TE  x,y  = {:10.5}{:10.5}  |", xbte, ybte);
        }

        // wipe out old flap hinge location
        xf.xbf = 0.0;
        xf.ybf = 0.0;
        xf.lbflap = false;

        gdes::abcopy(xf, true);

        // check for excessive panel corner angles
        let mut amax = 0.0;
        let mut imax = 0;
        geom::cang(&xf.x[..xf.n], &xf.y[..xf.n], 0, &mut amax, &mut imax, xf.show_output);
        if amax.abs() > 40.0 && xf.show_output {
            eprintln!();
            eprintln!(" WARNING: Poor input coordinate distribution");
            eprintln!("          Excessive panel angle {:7.1}  at i = {:4}", amax, imax + 1);
        }

        ([xble, yble], xf.chordb, [xbte, ybte])
    }

    /// Generates and panels the specified NACA 4- or 5-digit airfoil.
    pub fn naca(&mut self, spec: u32) {
        if spec == 0 {
            eprintln!("Invalid NACA specifier. Specify a NACA 4 or 5 series airfoil code.");
        } else {
            xfoil::naca(&mut self.state, spec as i32);
        }
    }

    /// Returns the current buffer airfoil coordinates (input points).
    pub fn airfoil(&self) -> (Vec<f64>, Vec<f64>) {
        let nb = self.state.nb;
        (self.state.xb[..nb].to_vec(), self.state.yb[..nb].to_vec())
    }

    /// Sets the Reynolds number.  A value of 0 selects inviscid mode.
    pub fn set_reynolds(&mut self, re: f64) {
        let xf = &mut self.state;
        if re == 0.0 {
            xf.lvisc = false;
        } else {
            xf.lvisc = true;
        }
        xf.reinf1 = re;
        xf.lbli_ni = false;
        xf.lipan = false;
        xf.lvconv = false;
    }

    /// Returns the current Reynolds number.
    pub fn reynolds(&self) -> f64 {
        self.state.reinf1
    }

    /// Sets the freestream Mach number.
    pub fn set_mach(&mut self, m: f64) {
        let xf = &mut self.state;
        if m != xf.minf {
            xf.minf = m;

            s_xfoil::comset(xf);
            if xf.minf > 0.0 && xf.show_output {
                eprintln!();
                eprintln!(" Sonic Cp = {:10.2}      Sonic Q/Qinf = {:10.3}", xf.cpstar, xf.qstar / xf.qinf);
            }

            s_xfoil::cpcalc(&xf.qinv[..xf.n], xf.qinf, xf.minf, &mut xf.cpi[..xf.n], xf.show_output);
            if xf.lvisc {
                s_xfoil::cpcalc(
                    &xf.qvis[..xf.n + xf.nw],
                    xf.qinf,
                    xf.minf,
                    &mut xf.cpv[..xf.n + xf.nw],
                    xf.show_output,
                );
            }
            s_xfoil::clcalc(
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
            s_xfoil::cdcalc(xf);
        }
        xf.lvconv = false;
    }

    /// Returns the current Mach number.
    pub fn mach(&self) -> f64 {
        self.state.minf
    }

    /// Sets the transition trip x/c locations on the top and bottom surfaces.
    pub fn set_xtr(&mut self, xtr_top: f64, xtr_bot: f64) {
        self.state.xstrip[0] = xtr_top;
        self.state.xstrip[1] = xtr_bot;
        self.state.lvconv = false;
    }

    /// Returns the transition trip x/c locations (top, bottom).
    pub fn xtr(&self) -> (f64, f64) {
        (self.state.xstrip[0], self.state.xstrip[1])
    }

    /// Sets the critical amplification ratio (Ncrit) for the e^n transition
    /// model on both surfaces.
    pub fn set_n_crit(&mut self, n_crit: f64) {
        self.state.acrit[0] = n_crit;
        self.state.acrit[1] = n_crit;
        self.state.lvconv = false;
    }

    /// Returns the critical amplification ratio.
    pub fn n_crit(&self) -> f64 {
        self.state.acrit[0]
    }

    /// Sets the maximum number of BL Newton iterations.
    pub fn set_max_iter(&mut self, max_iter: i32) {
        self.state.itmax = max_iter;
    }

    /// Returns the maximum number of BL Newton iterations.
    pub fn max_iter(&self) -> i32 {
        self.state.itmax
    }

    /// Forces BL re-initialization on the next operating point.
    pub fn reset_bls(&mut self) {
        self.state.lbli_ni = false;
        self.state.lipan = false;
        if self.state.show_output {
            eprintln!("BLs will be initialized on next point");
        }
    }

    /// Re-panels the current buffer airfoil with the specified paneling
    /// parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn repanel(
        &mut self,
        n_panel: i32,
        cv_par: f64,
        cte_ratio: f64,
        ctr_ratio: f64,
        xs_ref1: f64,
        xs_ref2: f64,
        xp_ref1: f64,
        xp_ref2: f64,
    ) {
        let xf = &mut self.state;
        xf.npan = (n_panel as usize).min(state::IQX - 6);
        xf.cvpar = cv_par;
        xf.cterat = cte_ratio;
        xf.ctrrat = ctr_ratio;
        xf.xsref1 = xs_ref1;
        xf.xsref2 = xs_ref2;
        xf.xpref1 = xp_ref1;
        xf.xpref2 = xp_ref2;

        xfoil::pangen(xf, true);
        if xf.n > 0 {
            let mut amax = 0.0;
            let mut imax = 0;
            geom::cang(&xf.x[..xf.n], &xf.y[..xf.n], 1, &mut amax, &mut imax, xf.show_output);
        }
    }

    /// Converges the flow at the specified angle of attack (degrees) and
    /// returns (CL, CD, CM, Cpmin, converged).
    pub fn a(&mut self, alpha_deg: f64) -> (f64, f64, f64, f64, bool) {
        let xf = &mut self.state;
        xf.adeg = alpha_deg;
        xf.alfa = alpha_deg * state::DTOR;
        xf.qinf = 1.0;
        xf.lalfa = true;

        oper::specal(xf);

        if (xf.alfa - xf.awake).abs() > 1.0E-5 {
            xf.lwake = false;
        }
        if (xf.alfa - xf.avisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }
        if (xf.minf - xf.mvisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }

        let mut conv;
        if xf.lvisc {
            conv = oper::viscal(xf, xf.itmax);
            conv = xf.lvconv && conv;
        } else {
            conv = true;
        }

        oper::fcpmin(xf);

        (xf.cl, xf.cd, xf.cm, xf.cpmn, conv)
    }

    /// Converges the flow at the specified lift coefficient and returns
    /// (alpha_deg, CD, CM, Cpmin, converged).
    pub fn cl(&mut self, target_cl: f64) -> (f64, f64, f64, f64, bool) {
        let xf = &mut self.state;
        xf.clspec = target_cl;
        xf.alfa = 0.0;
        xf.qinf = 1.0;
        xf.lalfa = false;

        oper::speccl(xf);

        xf.adeg = xf.alfa / state::DTOR;
        if (xf.alfa - xf.awake).abs() > 1.0E-5 {
            xf.lwake = false;
        }
        if (xf.alfa - xf.avisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }
        if (xf.minf - xf.mvisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }

        let mut conv;
        if xf.lvisc {
            conv = oper::viscal(xf, xf.itmax);
            conv = xf.lvconv && conv;
        } else {
            conv = true;
        }

        oper::fcpmin(xf);

        (xf.alfa / state::DTOR, xf.cd, xf.cm, xf.cpmn, conv)
    }

    /// Runs an alpha sweep from `a_start` to `a_end` (degrees) in `n` equal
    /// steps.  Returns (alpha, CL, CD, CM, Cpmin, converged) per point.
    ///
    /// Each point warm-starts from the previous point's converged state (the
    /// canonical XFOIL sweep behaviour).  For sweeps large enough to be
    /// worth splitting, see [`aseq_par`](Self::aseq_par).
    pub fn aseq(&mut self, a_start: f64, a_end: f64, n: usize) -> Vec<(f64, f64, f64, f64, f64, bool)> {
        let xf = &mut self.state;
        let a0 = a_start * state::DTOR;
        let da = (a_end - a_start) / n as f64 * state::DTOR;
        xf.lalfa = true;
        aseq_range(xf, a0, da, 0..n)
    }

    /// Parallel alpha sweep: identical stepping to [`aseq`], but the sweep is
    /// split into contiguous chunks solved on separate rayon threads.
    ///
    /// Each chunk warm-starts point-to-point within itself; only the first
    /// point of each chunk cold-starts from the configured state (the same
    /// initialization the serial sweep's first point does).  At Reynolds
    /// numbers of 1e6 and above the boundary-layer solve is robust enough
    /// that this matches the serial sweep bit-for-bit (verified for NACA
    /// 0012, Re = 1e6).  Below 1e6 the BL convergence is history-sensitive
    /// and the cold-started chunk boundaries can land on different
    /// (non-)converged states, so the sweep falls back to the serial path.
    pub fn aseq_par(&self, a_start: f64, a_end: f64, n: usize) -> Vec<(f64, f64, f64, f64, f64, bool)> {
        let threads = rayon::current_num_threads();
        let a0 = a_start * state::DTOR;
        let da = (a_end - a_start) / n as f64 * state::DTOR;
        if n < 4 * threads || self.state.reinf1 < 1.0e6 {
            // Too small to split, or the viscous convergence is
            // history-sensitive: run the serial path.
            let mut xf = self.state.clone();
            xf.lalfa = true;
            return aseq_range(&mut xf, a0, da, 0..n);
        }
        // Equal-sized chunks, one per rayon thread: each chunk boundary
        // cold-starts, and smaller chunks (more boundaries) cost more than
        // the work-stealing load balance gains.
        let chunk = n.div_ceil(threads);
        (0..n)
            .into_par_iter()
            .chunks(chunk)
            .flat_map(|range| {
                let mut xf = self.state.clone();
                xf.lalfa = true;
                aseq_range(&mut xf, a0, da, range.into_iter())
            })
            .collect()
    }

    /// Runs a CL sweep from `cl_start` to `cl_end` in `n` equal steps.
    /// Returns (alpha, CL, CD, CM, Cpmin, converged) per point.
    pub fn cseq(&mut self, cl_start: f64, cl_end: f64, n: usize) -> Vec<(f64, f64, f64, f64, f64, bool)> {
        let xf = &mut self.state;
        let cl0 = cl_start;
        let dcl = (cl_end - cl_start) / n as f64;
        xf.lalfa = false;

        let mut results = Vec::with_capacity(n);
        for i in 0..n {
            xf.clspec = cl0 + dcl * i as f64;

            oper::speccl(xf);

            if (xf.alfa - xf.awake).abs() > 1.0E-5 {
                xf.lwake = false;
            }
            if (xf.alfa - xf.avisc).abs() > 1.0E-5 {
                xf.lvconv = false;
            }
            if (xf.minf - xf.mvisc).abs() > 1.0E-5 {
                xf.lvconv = false;
            }

            let itmaxs = xf.itmax + 5;
            let mut conv = false;
            if xf.lvisc {
                conv = oper::viscal(xf, itmaxs);
            }
            xf.adeg = xf.alfa / state::DTOR;

            oper::fcpmin(xf);

            if (xf.lvconv && conv) || !xf.lvisc {
                conv = true;
            } else if xf.lvisc && !(xf.lvconv && conv) {
                conv = false;
            }

            results.push((xf.adeg, xf.cl, xf.cd, xf.cm, xf.cpmn, conv));
        }
        results
    }

    /// Returns the current number of panel nodes.
    pub fn get_n_cp(&self) -> usize {
        self.state.n
    }

    /// Returns true if the current airfoil was treated as having a sharp
    /// (closed) trailing edge, false if it has a finite TE gap.
    pub fn has_sharp_te(&self) -> bool {
        self.state.sharp
    }

    /// Returns the trailing-edge wake-gap thickness held by the most recent
    /// converged viscous solution (`dwte = wgap[0]`).  This is zero for a
    /// sharp-TE airfoil and a small positive value for a blunt TE.
    pub fn wake_gap_te(&self) -> f64 {
        self.state.dwte
    }

    /// Returns the current (x, Cp) distribution along the airfoil surface,
    /// computed from the surface vorticity.
    pub fn get_cp_distribution(&mut self) -> (Vec<f64>, Vec<f64>) {
        let xf = &mut self.state;
        s_xfoil::comset(xf);

        let n = xf.n;
        // Reuse cpcalc for the Karman-Tsien compressibility correction (it
        // hoists beta/bfac out of the loop and warns on supersonic points).
        s_xfoil::cpcalc(&xf.gam[..n], xf.qinf, xf.minf, &mut xf.cpi[..n], xf.show_output);

        let mut x_out = vec![0.0; n];
        let mut cp_out = vec![0.0; n];
        x_out.copy_from_slice(&xf.x[..n]);
        cp_out.copy_from_slice(&xf.cpi[..n]);
        (x_out, cp_out)
    }
}

/// The per-point body of `aseq`, factored out so the parallel sweep can run
/// it on per-chunk engine clones.  Points are solved sequentially over
/// `range`, warm-starting from the previous point within the same engine.
fn aseq_range(
    xf: &mut Xfoil,
    a0: f64,
    da: f64,
    range: impl ExactSizeIterator<Item = usize>,
) -> Vec<(f64, f64, f64, f64, f64, bool)> {
    let mut results = Vec::with_capacity(range.size_hint().0);
    for i in range {
        xf.alfa = a0 + da * i as f64;
        if (xf.alfa - xf.awake).abs() > 1.0E-5 {
            xf.lwake = false;
        }
        if (xf.alfa - xf.avisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }
        if (xf.minf - xf.mvisc).abs() > 1.0E-5 {
            xf.lvconv = false;
        }

        oper::specal(xf);

        let itmaxs = xf.itmax + 5;
        let mut conv = false;
        if xf.lvisc {
            conv = oper::viscal(xf, itmaxs);
        }
        xf.adeg = xf.alfa / state::DTOR;

        oper::fcpmin(xf);

        if (xf.lvconv && conv) || !xf.lvisc {
            conv = true;
        } else if xf.lvisc && !(xf.lvconv && conv) {
            conv = false;
        }

        results.push((xf.adeg, xf.cl, xf.cd, xf.cm, xf.cpmn, conv));
    }
    results
}
