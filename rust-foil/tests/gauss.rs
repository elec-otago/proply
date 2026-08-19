// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
use rust_foil::gauss;

#[test]
fn gauss_small_system() {
    // 2x2 system embedded in a 4x4 matrix (like the BL similarity solve):
    // -118145*Th + 131440*Ds = 0.776
    //  156771*Th -  95057*Ds = -0.323
    // dAm = 0, dUe = 0
    let mut z = [0.0f64; 16];
    let mut r = [0.0f64; 4];
    z[0] = 1.0; // row0 col0
    z[1 * 4 + 1] = -118145.0; // row1 col1
    z[2 * 4 + 1] = 131440.0; // row1 col2
    z[1 * 4 + 2] = 156771.0; // row2 col1
    z[2 * 4 + 2] = -95057.0; // row2 col2
    z[3 * 4 + 3] = 1.0; // row3 col3
    r[1] = 0.776;
    r[2] = -0.323;

    gauss(&mut z, &mut r, 4, 4, 1);

    // expected: dTh ~ 3.34e-6, dDs ~ 8.9e-6
    assert!(
        (r[1] - 3.34e-6).abs() < 1e-8,
        "dTh = {} (expected ~3.34e-6)",
        r[1]
    );
    assert!(
        (r[2] - 8.9e-6).abs() < 1e-8,
        "dDs = {} (expected ~8.9e-6)",
        r[2]
    );
    assert!(r[0].abs() < 1e-12);
    assert!(r[3].abs() < 1e-12);
}

#[test]
fn gauss_identity() {
    // Solve I x = b
    let mut z = [0.0f64; 16];
    let mut r = [1.0f64, 2.0, 3.0, 4.0];
    for i in 0..4 {
        z[i * 4 + i] = 1.0;
    }
    gauss(&mut z, &mut r, 4, 4, 1);
    assert!((r[0] - 1.0).abs() < 1e-12);
    assert!((r[1] - 2.0).abs() < 1e-12);
    assert!((r[2] - 3.0).abs() < 1e-12);
    assert!((r[3] - 4.0).abs() < 1e-12);
}

#[test]
fn gauss_three_by_three() {
    // x + 2y + 3z = 14; 2x + y + z = 7; 3x + y + 2z = 11 -> x=1,y=2,z=3
    let mut z = [0.0f64; 16];
    let mut r = [14.0f64, 7.0, 11.0, 0.0];
    z[0] = 1.0;
    z[1] = 2.0;
    z[2] = 3.0;
    z[4] = 2.0;
    z[5] = 1.0;
    z[6] = 1.0;
    z[8] = 3.0;
    z[9] = 1.0;
    z[10] = 2.0;
    gauss(&mut z, &mut r, 4, 3, 1);
    assert!((r[0] - 1.0).abs() < 1e-8, "x = {}", r[0]);
    assert!((r[1] - 2.0).abs() < 1e-8, "y = {}", r[1]);
    assert!((r[2] - 3.0).abs() < 1e-8, "z = {}", r[2]);
}

#[test]
fn gauss_with_pivoting() {
    // System that requires a row swap:
    //  0*x + 1*y = 2 ;  1*x + 0*y = 1   -> x=1, y=2
    let mut z = [0.0f64; 16];
    let mut r = [2.0f64, 1.0, 0.0, 0.0];
    z[0 * 4 + 0] = 0.0;
    z[0 * 4 + 1] = 1.0; // row1 col0
    z[1 * 4 + 0] = 1.0; // row0 col1
    z[1 * 4 + 1] = 0.0;
    gauss(&mut z, &mut r, 4, 2, 1);
    assert!((r[0] - 1.0).abs() < 1e-8, "x = {}", r[0]);
    assert!((r[1] - 2.0).abs() < 1e-8, "y = {}", r[1]);
}

#[test]
fn gauss_random_system() {
    // 3x3 with a known solution x = [2, -1, 0.5]
    // A = [[4, 1, 2], [1, 3, 1], [2, 1, 5]]
    // b = [8 + (-1) + 1, 2 + (-3) + 0.5, 4 + (-1) + 2.5] = [8, -0.5, 5.5]
    let mut z = [0.0f64; 16];
    let mut r = [8.0f64, -0.5, 5.5, 0.0];
    let a = [[4.0f64, 1.0, 2.0], [1.0, 3.0, 1.0], [2.0, 1.0, 5.0]];
    for col in 0..3 {
        for row in 0..3 {
            z[col * 4 + row] = a[row][col];
        }
    }
    gauss(&mut z, &mut r, 4, 3, 1);
    assert!((r[0] - 2.0).abs() < 1e-8, "x = {}", r[0]);
    assert!((r[1] + 1.0).abs() < 1e-8, "y = {}", r[1]);
    assert!((r[2] - 0.5).abs() < 1e-8, "z = {}", r[2]);
}
