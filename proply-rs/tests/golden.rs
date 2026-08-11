//! Golden-value tests: compare the Rust port against reference values
//! generated from the Python implementation (numpy/scipy) by
//! `build/golden/gen_golden.py`.

use proply_rs::foil::{FoilLike, Naca4};
use proply_rs::motor::Motor;
use proply_rs::optimize::{bem_iterate, optimize_all, PlateSim};
use proply_rs::pchip::Pchip;
use proply_rs::polyfit::{polyfit, polyval};
use proply_rs::smooth::smooth;
use serde_json::Value;

fn golden(name: &str) -> Value {
    let path = format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name);
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read golden file {}: {}", path, e));
    serde_json::from_str(&data).expect("invalid golden JSON")
}

fn near(a: f64, b: f64, tol: f64, ctx: &str) {
    assert!((a - b).abs() <= tol, "{}: got {}, want {}", ctx, a, b);
}

fn near_array(got: &[f64], want: &[f64], tol: f64, ctx: &str) {
    assert_eq!(got.len(), want.len(), "{}: length mismatch", ctx);
    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        near(*g, *w, tol, &format!("{}[{}]", ctx, i));
    }
}

fn arr(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect()
}

#[test]
fn naca4_golden() {
    let g = golden("naca4.json");

    // Symmetric NACA 0012, unit chord.
    let f = Naca4::new(1.0, 0.12, 0.0, 0.4);
    let (xl, yl, xu, yu) = f.get_shape_points(42);
    let s = &g["symmetric_12"];
    near_array(&xl, &arr(&s["xl"]), 1e-12, "xl");
    near_array(&yl, &arr(&s["yl"]), 1e-12, "yl");
    near_array(&xu, &arr(&s["xu"]), 1e-12, "xu");
    near_array(&yu, &arr(&s["yu"]), 1e-12, "yu");

    // Cambered NACA 6415 with a trailing-edge gap, chord 0.1.
    let mut f = Naca4::new(0.1, 0.15, 0.06, 0.4);
    f.set_trailing_edge(0.01);
    let (xl, yl, xu, yu) = f.get_shape_points(42);
    let s = &g["cambered_15"];
    near_array(&xl, &arr(&s["xl"]), 1e-12, "cambered xl");
    near_array(&yl, &arr(&s["yl"]), 1e-12, "cambered yl");
    near_array(&xu, &arr(&s["xu"]), 1e-12, "cambered xu");
    near_array(&yu, &arr(&s["yu"]), 1e-12, "cambered yu");

    // Bounding box of a rotated foil (chord 0.02, twist 0.3 rad).
    let f = Naca4::new(0.02, 0.15, 0.06, 0.4);
    let (x0, x1, y0, y1) = f.get_bounding_box(0.3);
    let b = &g["bounding_box"];
    near(x0, b["x0"].as_f64().unwrap(), 1e-12, "bb x0");
    near(x1, b["x1"].as_f64().unwrap(), 1e-12, "bb x1");
    near(y0, b["y0"].as_f64().unwrap(), 1e-12, "bb y0");
    near(y1, b["y1"].as_f64().unwrap(), 1e-12, "bb y1");
}

#[test]
fn pchip_golden() {
    let g = golden("pchip.json");
    let radius = 0.0625;
    let hub_r = 0.005;

    // Scimitar offset dataset.
    let x = [0.0, hub_r, radius * 0.8, radius];
    let y = [0.0, 0.0, radius * (-5.0 / 100.0), 0.0];
    let p = Pchip::new(&x, &y);
    let q = arr(&g["scimitar"]["x"]);
    let want = arr(&g["scimitar"]["y"]);
    let got: Vec<f64> = q.iter().map(|t| p.eval(*t)).collect();
    near_array(&got, &want, 1e-9, "scimitar");

    // Chord smoothing: smooth() output then PCHIP evaluation.
    let extra = arr(&g["chord_smooth"]["smooth_input"]);
    let sm = smooth(&extra, 11, "hanning");
    let want_sm = arr(&g["chord_smooth"]["smooth_output"]);
    near_array(&sm, &want_sm, 1e-12, "smooth");

    let mut c_pts = vec![0.0, hub_r / 2.0, 0.9 * hub_r];
    c_pts.extend((0..40).map(|i| hub_r + (radius - hub_r) * i as f64 / 39.0));
    let pc = Pchip::new(&c_pts, &sm);
    let rq = arr(&g["chord_smooth"]["r"]);
    let want_c = arr(&g["chord_smooth"]["chord"]);
    let got_c: Vec<f64> = rq.iter().map(|t| pc.eval(*t)).collect();
    near_array(&got_c, &want_c, 1e-9, "chord pchip");
}

#[test]
fn polyfit_golden() {
    let g = golden("polyfit.json");
    let alpha: Vec<f64> = (-40..=40).map(|i| (i as f64 * 0.5).to_radians()).collect();
    let cl: Vec<f64> = alpha
        .iter()
        .map(|a| 0.62 * (2.0 * a).sin() + 0.05 * a + 0.4 * a * a - 3.0 * a.powi(3))
        .collect();
    let cd: Vec<f64> = alpha
        .iter()
        .map(|a| 0.008 + 0.2 * a * a + 1.5 * a.powi(4))
        .collect();

    let cl9 = polyfit(&alpha, &cl, 9);
    let cd9 = polyfit(&alpha, &cd, 9);

    let ev = arr(&g["eval_rad"]);
    let want_cl = arr(&g["eval_cl"]);
    let want_cd = arr(&g["eval_cd"]);
    let got_cl: Vec<f64> = ev.iter().map(|a| polyval(&cl9, *a)).collect();
    let got_cd: Vec<f64> = ev.iter().map(|a| polyval(&cd9, *a)).collect();
    near_array(&got_cl, &want_cl, 1e-6, "eval cl");
    near_array(&got_cd, &want_cd, 1e-8, "eval cd");

    // Degree-4 twist fit.
    let r: Vec<f64> = (0..40).map(|i| 0.0625 - (0.0625 - 0.005) * i as f64 / 39.0).collect();
    let twist: Vec<f64> = r
        .iter()
        .map(|ri| (25.0 * (ri / 0.0625).powf(0.6) + 5.0 * (ri * 40.0).sin()).to_radians())
        .collect();
    let rrev: Vec<f64> = r.iter().rev().copied().collect();
    let c4 = polyfit(&rrev, &twist, 4);
    let want_tw = arr(&g["twist_eval"]);
    let got_tw: Vec<f64> = rrev.iter().map(|ri| polyval(&c4, *ri)).collect();
    near_array(&got_tw, &want_tw, 1e-7, "twist eval");
}

#[test]
fn motor_golden() {
    let g = golden("motor.json");
    let m = Motor::new(1900.0, 0.5, 0.405);
    near(m.get_imax(11.0), g["Imax"].as_f64().unwrap(), 1e-9, "Imax");
    let (q, rpm) = m.get_qmax(11.0);
    near(q, g["Qmax"].as_f64().unwrap(), 1e-9, "Qmax");
    near(rpm, g["RPMmax"].as_f64().unwrap(), 1e-9, "RPMmax");
    near(m.get_pmax(11.0), g["Pmax"].as_f64().unwrap(), 1e-9, "Pmax");
    near(m.get_torque(3.0), g["torque_at_3A"].as_f64().unwrap(), 1e-12, "torque@3A");
    near(m.get_rpm(0.01), g["rpm_at_0_01"].as_f64().unwrap(), 1e-9, "rpm@0.01");
}

#[test]
fn bem_equations_golden() {
    let g = golden("bem.json");
    let fs = PlateSim { chord: 0.008 };
    let (dv, ap, theta, rpm, r, dr, u_0, b) = (5.0, 0.05, 28.0_f64.to_radians(), 12000.0, 0.03, 0.002, 1.0, 3.0);
    let omega = proply_rs::optimize::rpm2omega(rpm);
    let p = &g["precalc"];
    near(omega, p["omega"].as_f64().unwrap(), 1e-9, "omega");
    let (cl, cd, phi) = proply_rs::optimize::precalc(&fs, dv, ap, theta, omega, r, dr, u_0, b as f64);
    near(cl, p["CL"].as_f64().unwrap(), 1e-9, "CL");
    near(cd, p["CD"].as_f64().unwrap(), 1e-9, "CD");
    near(phi, p["phi"].as_f64().unwrap(), 1e-9, "phi");
    let (dv_new, ap_new) = proply_rs::optimize::iterate(&fs, fs.chord, dv, ap, theta, omega, r, dr, u_0, b as f64);
    let it = &g["iterate"];
    near(dv_new, it["dv_new"].as_f64().unwrap(), 1e-9, "dv_new");
    near(ap_new, it["a_prime_new"].as_f64().unwrap(), 1e-9, "a_prime_new");
    let f = &g["forces"];
    near(proply_rs::optimize::d_t(dv, r, dr, u_0), f["dT"].as_f64().unwrap(), 1e-9, "dT");
    near(proply_rs::optimize::d_m(dv, ap, r, dr, omega, u_0), f["dM"].as_f64().unwrap(), 1e-9, "dM");
    near(
        proply_rs::optimize::dv_from_thrust(0.3, 0.05, 1.0),
        f["dv_from_thrust"].as_f64().unwrap(),
        1e-9,
        "dv_from_thrust",
    );
}

#[test]
fn optimizer_matches_slsqp_reference() {
    // The Python version solves these with scipy SLSQP; the Nelder-Mead port
    // should land on the same interior optimum.  Tolerances reflect that
    // Nelder-Mead is a different (derivative-free) algorithm.
    let g = golden("bem.json");
    let fs = PlateSim { chord: 0.008 };
    let (dv_goal, theta, rpm, r, dr, u_0, b) = (5.0, 28.0_f64.to_radians(), 12000.0, 0.03, 0.002, 1.0, 3.0);
    let (dv, ap, err) = bem_iterate(&fs, dv_goal, theta, rpm, r, dr, u_0, b as f64);
    assert!(err < 1e-6, "bem err {}", err);
    let x = &g["slsqp_bem"]["x"];
    near(dv, x[0].as_f64().unwrap(), 0.05, "bem dv vs SLSQP");
    near(ap, x[1].as_f64().unwrap(), 0.005, "bem a_prime vs SLSQP");

    let (x4, fun) = optimize_all(&fs, dv_goal, rpm, r, dr, u_0, b as f64, 0.012);
    let xa = &g["slsqp_all"]["x"];
    // The objective is flat in theta and chord near the optimum, so those
    // may differ while the objective value matches (NM actually finds a
    // slightly better point here); dv and a_prime track SLSQP.
    near(x4[0], xa[0].as_f64().unwrap(), 0.08, "all theta vs SLSQP");
    near(x4[1], xa[1].as_f64().unwrap(), 0.3, "all dv vs SLSQP");
    near(x4[2], xa[2].as_f64().unwrap(), 0.02, "all a_prime vs SLSQP");
    near(x4[3], xa[3].as_f64().unwrap(), 0.005, "all chord vs SLSQP");
    // ...and never materially worse than the SLSQP optimum.
    assert!(
        fun <= g["slsqp_all"]["fun"].as_f64().unwrap() + 0.02,
        "optimize_all fun {} vs SLSQP {}",
        fun,
        g["slsqp_all"]["fun"].as_f64().unwrap()
    );
}

#[test]
fn buhl_golden() {
    // Buhl (2005) NREL/TP-500-36834 turbulent-wake CT(a) relation, Eqs. 1
    // and 18 (decelerating-disk convention), verified against numpy in
    // build/golden/gen_golden.py.
    let g = golden("buhl.json");

    // Forward relation on a grid, for F = 1 and F = 0.8.
    for (name, f) in [("F1", 1.0), ("F08", 0.8)] {
        let s = &g[name];
        let a = arr(&s["a"]);
        let want = arr(&s["ct"]);
        let got: Vec<f64> = a.iter().map(|ai| proply_rs::optimize::ct_buhl(*ai, f)).collect();
        near_array(&got, &want, 1e-12, &format!("ct_buhl {}", name));
    }

    // Inverse a_buhl(CT, F) across both branches.
    let inv = &g["invert"];
    let ct = arr(&inv["ct"]);
    for (name, f) in [("F1", 1.0), ("F08", 0.8)] {
        let want = arr(&inv[name]);
        let got: Vec<f64> = ct.iter().map(|c| proply_rs::optimize::a_buhl(*c, f)).collect();
        near_array(&got, &want, 1e-12, &format!("a_buhl {}", name));
    }

    // Round-trip: ct_buhl(a_buhl(CT)) == CT on both branches.
    for f in [1.0, 0.8] {
        for c in [0.2, 0.5, 0.9, 0.96, 1.0, 1.2, 1.5, 1.9] {
            let a = proply_rs::optimize::a_buhl(c, f);
            near(proply_rs::optimize::ct_buhl(a, f), c, 1e-9, &format!("round-trip F={} CT={}", f, c));
        }
    }
}
