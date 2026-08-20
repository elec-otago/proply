// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Kulfan (CST) parametrization tests: forward geometry, the inverse fit
//! (coordinates → params), the NACA-to-CST path, and end-to-end solver
//! convergence on a CST-generated airfoil.

use rust_foil::{KulfanParams, XFoil};

fn close(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() < tol,
        "{}: got {} expected {} (diff {})",
        what,
        a,
        b,
        (a - b).abs()
    );
}

#[test]
fn default_params_forward_sanity() {
    // Default params: upper = +0.2·ones(8), lower = −0.2·ones(8), LEM = 0,
    // TE = 0 — a symmetric ~14%-thick section.
    let p = KulfanParams::default();
    assert_eq!(p.upper_weights.len(), 8);
    assert_eq!(p.lower_weights.len(), 8);
    close(p.n1, 0.5, 1e-12, "N1");
    close(p.n2, 1.0, 1e-12, "N2");

    let n = 213;
    let (xb, yb, nb) = p.coordinates(n);
    assert_eq!(nb, 2 * n - 1, "closed loop with shared LE");
    assert_eq!(xb.len(), nb);
    assert_eq!(yb.len(), nb);

    // Closed loop: first and last point share the TE.
    close(xb[0], 1.0, 1e-12, "TE x (start)");
    close(xb[nb - 1], 1.0, 1e-12, "TE x (end)");
    close(yb[0], yb[nb - 1], 1e-12, "TE y closure");

    // LE at (0, 0) and positive upper thickness mid-chord.
    let le = nb / 2;
    close(xb[le], 0.0, 1e-12, "LE x");
    close(yb[le], 0.0, 1e-12, "LE y");
    assert!(yb[le / 2] > 0.0, "upper surface y > 0 mid-chord");

    // Single-station evaluation agrees with the loop points.
    close(p.upper_y(0.5), yb[le / 2], 1e-12, "upper_y vs loop");
    close(
        p.lower_y(0.5),
        yb[nb - 1 - le / 2],
        1e-12,
        "lower_y vs loop",
    );
}

#[test]
fn fit_round_trip_is_exact() {
    // Generate coordinates from random-ish params, fit them back, and
    // regenerate: the fit is exact on CST-generated points, so the
    // coordinates must match to solver precision.
    let p = KulfanParams {
        lower_weights: vec![-0.30, -0.12, 0.05, -0.08, 0.12, -0.05, 0.02, 0.10],
        upper_weights: vec![0.28, 0.10, -0.04, 0.15, -0.10, 0.03, -0.02, 0.18],
        leading_edge_weight: 0.07,
        te_thickness: 0.013,
        n1: 0.5,
        n2: 1.0,
    };
    let n = 213;
    let (x, y, nb) = p.coordinates(n);

    let fitted = KulfanParams::fit_from_coordinates(&x, &y);
    let (x2, y2, nb2) = fitted.coordinates(n);
    assert_eq!(nb2, nb);

    let mut max_err = 0.0f64;
    for i in 0..nb {
        max_err = max_err.max((x2[i] - x[i]).abs().max((y2[i] - y[i]).abs()));
    }
    assert!(
        max_err < 1e-6,
        "round-trip coordinate error too large: {}",
        max_err
    );
}

#[test]
fn fit_recovers_known_params() {
    let p = KulfanParams {
        lower_weights: vec![-0.30, -0.12, 0.05, -0.08, 0.12, -0.05, 0.02, 0.10],
        upper_weights: vec![0.28, 0.10, -0.04, 0.15, -0.10, 0.03, -0.02, 0.18],
        leading_edge_weight: 0.07,
        te_thickness: 0.013,
        n1: 0.5,
        n2: 1.0,
    };
    let (x, y, _nb) = p.coordinates(213);
    let f = KulfanParams::fit_from_coordinates(&x, &y);
    for i in 0..8 {
        close(
            f.upper_weights[i],
            p.upper_weights[i],
            1e-6,
            &format!("upper[{}]", i),
        );
        close(
            f.lower_weights[i],
            p.lower_weights[i],
            1e-6,
            &format!("lower[{}]", i),
        );
    }
    close(f.leading_edge_weight, p.leading_edge_weight, 1e-6, "LEM");
    close(f.te_thickness, p.te_thickness, 1e-6, "TE");
}

#[test]
fn naca_0012_fit_is_symmetric() {
    // The NACA 0012 shape is symmetric: the fit must give lower = −upper,
    // LEM ≈ 0, and a small positive TE gap (the analytic NACA thickness
    // polynomial has a small open TE, which the fit captures as the TE
    // thickness term).
    let (p, name) = KulfanParams::from_naca(12, false).expect("NACA 0012 is valid");
    assert_eq!(name, "NACA0012");

    for i in 0..8 {
        close(
            p.upper_weights[i],
            -p.lower_weights[i],
            1e-6,
            &format!("antisymmetry[{}]", i),
        );
    }
    close(p.leading_edge_weight, 0.0, 1e-6, "LEM");
    assert!(
        p.te_thickness > 0.0 && p.te_thickness < 5e-3,
        "TE gap = {} (NACA 0012 ~ 2.5e-3)",
        p.te_thickness
    );
}

#[test]
fn naca_fit_fidelity_is_bounded() {
    // NACA 4/5-digit shapes are not exactly CST-representable; the fit must
    // reproduce the analytic geometry to ~1e-3 in y/c or better (measured
    // max error at 8 weights/side is ~1e-4), so known NACA geometry
    // quantities must survive the CST path.
    let (p0012, _) = KulfanParams::from_naca(12, false).unwrap();
    let (p2412, _) = KulfanParams::from_naca(2412, false).unwrap();

    // Max thickness: NACA 0012/2412 are 12% thick.
    for p in [&p0012, &p2412] {
        let mut tmax = 0.0f64;
        for i in 0..1000 {
            let x = i as f64 / 999.0;
            tmax = tmax.max(p.upper_y(x) - p.lower_y(x));
        }
        close(tmax, 0.12, 1e-3, "max thickness");
    }

    // Max camber: NACA 2412 has 2% camber at x/c = 0.4.
    let mut cmax = 0.0f64;
    let mut x_cmax = 0.0;
    for i in 0..1000 {
        let x = i as f64 / 999.0;
        let cam = 0.5 * (p2412.upper_y(x) + p2412.lower_y(x));
        if cam > cmax {
            cmax = cam;
            x_cmax = x;
        }
    }
    close(cmax, 0.02, 1e-3, "max camber");
    close(x_cmax, 0.4, 0.05, "max camber x/c");
}

#[test]
fn naca5_illegal_designation() {
    // The 5-digit camber is only defined for first-three-digit codes 210–250.
    assert!(KulfanParams::from_naca(10012, false).is_none());
}

#[test]
fn cst_solver_converges() {
    // A CST-generated airfoil at Re = 1e6 must converge with positive CL at
    // 4° incidence.  The default params (te = 0) give a perfectly sharp TE,
    // on which XFOIL's viscous solve does not converge, so the test uses a
    // small blunt TE gap — like every airfoil the NACA path produces.
    let params = KulfanParams {
        te_thickness: 0.002,
        ..KulfanParams::default()
    };

    let mut xf = XFoil::new();
    xf.set_show_output(false);
    xf.cst(&params);
    xf.set_reynolds(1.0e6);
    xf.set_max_iter(100);
    let (cl, cd, _cm, _cp, conv) = xf.a(4.0);
    assert!(conv, "CST viscous solve did not converge");
    assert!(cl > 0.0, "CL = {} expected positive", cl);
    assert!(cd > 0.0, "CD = {} expected positive", cd);
}

#[test]
fn naca_through_cst_matches_cst_path() {
    // XFoil::naca now generates geometry through the CST fit: the buffer
    // must be identical to what XFoil::cst produces from the fitted params.
    let (params, _name) = KulfanParams::from_naca(2412, false).unwrap();

    let mut xf_naca = XFoil::new();
    xf_naca.set_show_output(false);
    xf_naca.naca(2412);
    let (xn, yn) = xf_naca.airfoil();

    let mut xf_cst = XFoil::new();
    xf_cst.set_show_output(false);
    xf_cst.cst(&params);
    let (xc, yc) = xf_cst.airfoil();

    assert_eq!(xn.len(), xc.len());
    for i in 0..xn.len() {
        close(xn[i], xc[i], 1e-12, &format!("x[{}]", i));
        close(yn[i], yc[i], 1e-12, &format!("y[{}]", i));
    }
}
