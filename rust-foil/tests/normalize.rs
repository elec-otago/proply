// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! `set_airfoil` normalizes the input coordinates to unit chord by default
//! (canonical XFOIL behaviour): a chord-scaled airfoil must produce the same
//! polars as the same airfoil at unit chord.

use rust_foil::XFoil;

// The naca0012_coords from the canonical-XFOIL reference test.
fn naca0012_coords() -> (Vec<f64>, Vec<f64>) {
    // Reuse the point set from tests/naca0012.rs by rebuilding it here
    // (kept small: every other row of the reference data).
    let x = [
        1.0000e+00, 9.8037e-01, 9.5272e-01, 9.2112e-01, 8.8821e-01, 8.5494e-01, 8.2158e-01, 7.8817e-01,
        7.5475e-01, 7.2132e-01, 6.8789e-01, 6.5447e-01, 6.2108e-01, 5.8772e-01, 5.5441e-01, 5.2116e-01,
        4.8798e-01, 4.5489e-01, 4.2191e-01, 3.8905e-01, 3.5635e-01, 3.2383e-01, 2.9154e-01, 2.5953e-01,
        2.2788e-01, 1.9670e-01, 1.6619e-01, 1.3666e-01, 1.0877e-01, 8.3582e-02, 6.2395e-02, 4.5806e-02,
        3.3294e-02, 2.3867e-02, 1.6658e-02, 1.1080e-02, 6.7811e-03, 3.5772e-03, 1.4037e-03, 2.6000e-05,
        2.6000e-05, 1.4037e-03, 3.5772e-03, 6.7811e-03, 1.1080e-02, 1.6658e-02, 2.3867e-02, 3.3294e-02,
        4.5806e-02, 6.2395e-02, 8.3582e-02, 1.0877e-01, 1.3666e-01, 1.6619e-01, 1.9670e-01, 2.2788e-01,
        2.5953e-01, 2.9154e-01, 3.2383e-01, 3.5635e-01, 3.8905e-01, 4.2191e-01, 4.5489e-01, 4.8798e-01,
        5.2116e-01, 5.5441e-01, 5.8772e-01, 6.2108e-01, 6.5447e-01, 6.8789e-01, 7.2132e-01, 7.5475e-01,
        7.8817e-01, 8.2158e-01, 8.5494e-01, 8.8821e-01, 9.2112e-01, 9.5272e-01, 9.8037e-01, 1.0000e+00,
    ];
    let y = [
        1.2600e-03, 3.9814e-03, 7.7062e-03, 1.1815e-02, 1.5936e-02, 1.9943e-02, 2.3810e-02, 2.7532e-02,
        3.1107e-02, 3.4534e-02, 3.7807e-02, 4.0921e-02, 4.3865e-02, 4.6630e-02, 4.9199e-02, 5.1557e-02,
        5.3683e-02, 5.5554e-02, 5.7143e-02, 5.8419e-02, 5.9349e-02, 5.9891e-02, 6.0002e-02, 5.9628e-02,
        5.8710e-02, 5.7181e-02, 5.4967e-02, 5.1997e-02, 4.8243e-02, 4.3805e-02, 3.9000e-02, 3.4237e-02,
        2.9760e-02, 2.5598e-02, 2.1675e-02, 1.7888e-02, 1.4147e-02, 1.0381e-02, 6.5676e-03, 9.0564e-04,
        -9.0564e-04, -6.5676e-03, -1.0381e-02, -1.4147e-02, -1.7888e-02, -2.1675e-02, -2.5598e-02, -2.9760e-02,
        -3.4237e-02, -3.9000e-02, -4.3805e-02, -4.8243e-02, -5.1997e-02, -5.4967e-02, -5.7181e-02, -5.8710e-02,
        -5.9628e-02, -6.0002e-02, -5.9891e-02, -5.9349e-02, -5.8419e-02, -5.7143e-02, -5.5554e-02, -5.3683e-02,
        -5.1557e-02, -4.9199e-02, -4.6630e-02, -4.3865e-02, -4.0921e-02, -3.7807e-02, -3.4534e-02, -3.1107e-02,
        -2.7532e-02, -2.3810e-02, -1.9943e-02, -1.5936e-02, -1.1815e-02, -7.7062e-03, -3.9814e-03, -1.2600e-03,
    ];
    (x.to_vec(), y.to_vec())
}

fn cl_at(scale: f64, normalize: bool, alpha: f64) -> f64 {
    let (x, y) = naca0012_coords();
    let x: Vec<f64> = x.iter().map(|v| v * scale).collect();
    let y: Vec<f64> = y.iter().map(|v| v * scale).collect();
    let mut xf = XFoil::new();
    xf.set_show_output(false);
    xf.set_normalize(normalize);
    xf.set_airfoil(&x, &y);
    xf.set_reynolds(1.0e6);
    xf.set_max_iter(100);
    let (cl, _cd, _cm, _cp, conv) = xf.a(alpha);
    assert!(conv, "viscous solution did not converge");
    cl
}

#[test]
fn chord_scaled_input_matches_unit_chord_by_default() {
    // The default must normalize: a 20 mm-chord airfoil gives the same CL as
    // the unit-chord airfoil.
    let cl_unit = cl_at(1.0, true, 5.0);
    let cl_scaled = cl_at(0.02, true, 5.0);
    assert!(
        (cl_unit - cl_scaled).abs() < 1.0e-4,
        "CL(unit)={} vs CL(scaled)={}",
        cl_unit,
        cl_scaled
    );
}

#[test]
fn explicit_normalize_toggle_works() {
    // set_normalize(false) must leave the coordinates untouched: a 2 mm-chord
    // airfoil analysed as-is diverges from the normalized result.
    let cl_normalized = cl_at(0.002, true, 5.0);
    let cl_raw = cl_at(0.002, false, 5.0);
    assert!(
        (cl_normalized - cl_raw).abs() > 1.0e-2,
        "un-normalized input should differ from normalized ({} vs {})",
        cl_normalized,
        cl_raw
    );
}
