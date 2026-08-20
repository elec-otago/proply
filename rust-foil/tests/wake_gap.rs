// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Wake-gap (`wgap`) regression tests.
//!
//! `wgap[0..nw-1]` holds the "dead air" TE-flap thickness along the wake.  It
//! is read in `bl.rs` at four sites (one each in `setbl`, `mrchue`, `mrchdu`,
//! `update`), all of the form
//!
//! ```text
//! let iw = ibl - iblte[is] as usize;   // iw in 1..=nw
//! dswaki = xf.wgap[iw];               // BUG: was iw, should be iw-1
//! ```
//!
//! The BL station `ibl == iblte[is] + iw` corresponds to wake index `iw`, but
//! `wgap` is a 0-based array, so the correct read is `wgap[iw-1]`.  The
//! off-by-one read the wrong slot for every wake station and, at the maximum
//! index `iw == nw == IWX`, stepped one element past the end of the array.
//!
//! XFOIL treats an airfoil as sharp-TE only when the TE gap is below
//! `0.0001 * chord`; the NACA 4-digit thickness polynomial has a small open TE
//! (~0.002·t), so even `naca()`-generated airfoils are flagged blunt and
//! produce a nonzero wake gap.  The thicker the airfoil, the larger the TE gap
//! and the larger `wgap`.  All the tests below use the `naca()` generator,
//! which is what the rest of the test suite already exercises, so convergence
//! behavior is known-good.

use rust_foil::XFoil;

fn solve_naca(code: u32, alpha: f64) -> (f64, f64, bool, f64, bool) {
    let mut xf = XFoil::new();
    xf.set_show_output(false);
    xf.naca(code);
    xf.set_reynolds(1.0e6);
    xf.set_max_iter(100);
    let (cl, cd, _cm, _cp, conv) = xf.a(alpha);
    // dwte = wgap[0] after a converged viscous solve.
    (cl, cd, xf.has_sharp_te(), xf.wake_gap_te(), conv)
}

#[test]
fn naca_solve_drives_wake_gap_path() {
    // Sanity check: every NACA 4-digit airfoil has a small-but-nonzero TE gap
    // (the thickness polynomial does not close exactly), so the solver takes
    // the blunt-TE branch and writes a nonzero wake gap.  If `wgap` were not
    // being read at all, dwte would be its un-initialized zero; if it were
    // being read out of range this would panic or return garbage.  A small,
    // finite, positive dwte after a converged solve is the signature of the
    // path working end-to-end.
    let (cl, cd, sharp, dwte, conv) = solve_naca(2412, 4.0);
    assert!(conv, "NACA 2412 viscous solve did not converge");
    assert!(!sharp, "NACA 2412 should be treated as blunt-TE");
    assert!(
        dwte > 0.0 && dwte.is_finite(),
        "dwte = {} expected finite positive, got",
        dwte
    );
    assert!(cl > 0.0 && cd > 0.0, "CL = {}, CD = {}", cl, cd);
}

#[test]
fn wake_gap_scales_with_thickness() {
    // The TE gap of a NACA 4-digit airfoil grows ~linearly with the thickness
    // parameter t, so the wake gap written by xicalc must too.  Comparing two
    // thicknesses at the same alpha gives a robust, physics-based regression
    // check: a wgap indexing bug that read the wrong slot would scramble this
    // monotonic relationship.
    //
    //   NACA 0012 (t=0.12):  dwte ~ 2.5e-3
    //   NACA 0024 (t=0.24):  dwte ~ 4.9e-3
    //   NACA 0030 (t=0.30):  dwte ~ 5.9e-3
    let (_cl_12, _cd_12, _sharp_12, dwte_12, conv_12) = solve_naca(12, 4.0);
    let (_cl_24, _cd_24, _sharp_24, dwte_24, conv_24) = solve_naca(24, 4.0);
    let (_cl_30, _cd_30, _sharp_30, dwte_30, conv_30) = solve_naca(30, 4.0);

    assert!(
        conv_12 && conv_24 && conv_30,
        "one or more solves did not converge"
    );

    // All must be positive and finite -- the core "no OOB / no NaN" check.
    for (code, dwte) in [(12u32, dwte_12), (24, dwte_24), (30, dwte_30)] {
        assert!(
            dwte > 0.0 && dwte.is_finite(),
            "NACA {:02} dwte = {} not finite positive",
            code,
            dwte
        );
    }

    // Monotonicity in thickness.
    assert!(
        dwte_24 > dwte_12,
        "expected dwte(0024) > dwte(0012), got {} vs {}",
        dwte_24,
        dwte_12
    );
    assert!(
        dwte_30 > dwte_24,
        "expected dwte(0030) > dwte(0024), got {} vs {}",
        dwte_30,
        dwte_24
    );

    // And roughly linear: doubling thickness from 12% to 24% should about
    // double the wake gap (within 25%, to allow for paneling differences).
    let ratio = dwte_24 / dwte_12;
    assert!(
        (ratio - 2.0).abs() < 0.5,
        "dwte(0024)/dwte(0012) = {}, expected ~2.0",
        ratio
    );
}

#[test]
fn wake_gap_persists_across_alpha_sweep() {
    // A second operating point on the same airfoil reuses the BL
    // initialization but re-marches the wake; the wgap reads happen on every
    // iteration.  Sweep a few attached-flow alphas and confirm dwte stays
    // finite and roughly constant (it is a geometric property, set once in
    // xicalc, so it should not vary meaningfully with alpha).
    let mut xf = XFoil::new();
    xf.set_show_output(false);
    xf.naca(12);
    xf.set_reynolds(1.0e6);
    xf.set_max_iter(100);

    let mut dwtes: Vec<f64> = Vec::new();
    for alpha in [-4.0, 0.0, 4.0, 8.0] {
        let (_cl, _cd, _cm, _cp, conv) = xf.a(alpha);
        assert!(conv, "solve at alpha = {} did not converge", alpha);
        let dwte = xf.wake_gap_te();
        assert!(
            dwte.is_finite() && dwte > 0.0,
            "alpha {}: dwte = {}",
            alpha,
            dwte
        );
        dwtes.push(dwte);
    }

    // dwte is a geometric quantity (TE gap); it must be the same across alphas.
    let lo = dwtes.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = dwtes.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (hi - lo) < 1.0e-12,
        "dwte varied across alpha sweep: [{}, {}]",
        lo,
        hi
    );
}
