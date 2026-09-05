// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Verifies the custom block solver (blsolv) against a dense gauss solve on
//! a hand-built 2-station Newton system with the same block structure that
//! setbl produces.

use rust_foil::gauss;
use rust_foil::state::Xfoil;

/// The blsolv matrix for NSYs=2 has columns
///   [dC1, dT1, dD1, dC2, dT2, dD2]
/// and rows
///   row 1: [ VA1(3x2)  | VM11(3x1) |     0     | VM21(3x1) ]
///   row 2: [ VB2(3x2)  | VM12(3x1) | VA2(3x2)  | VM22(3x1) ]
/// where V* blocks are stored column-first (VA1[c][r]).
fn build_system() -> (Xfoil, [[f64; 6]; 6], [f64; 6]) {
    let mut xf = Xfoil::new();
    xf.vaccel = 0.01;
    xf.n = 3;
    xf.s = vec![0.0, 0.5, 1.0];

    xf.nsys = 2;
    xf.isys[0][2] = 1;
    xf.isys[1][3] = 2;
    xf.isys[1][4] = 2; // wake station pointer used by the VZ elimination
    xf.iblte[0] = 2;
    xf.iblte[1] = 3;

    // blocks stored as [column][row]
    let va1 = [[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]];
    let va2 = [[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]];
    let vb2 = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0]];

    let m11 = [0.1, 0.2, 1.0]; // VM(., 1, 1): row 1, mass col of st1
    let m12 = [0.3, 0.4, 1.0]; // VM(., 1, 2): row 2, mass col of st1
    let m21 = [0.05, 0.06, 0.5]; // VM(., 2, 1): row 1, mass col of st2
    let m22 = [0.2, 0.3, 1.5]; // VM(., 2, 2): row 2, mass col of st2

    let b1 = [1.0, 2.0, 3.0];
    let b2 = [4.0, 5.0, 6.0];

    for k in 0..3 {
        xf.va[Xfoil::v_index(1, 1, k + 1)] = va1[0][k];
        xf.va[Xfoil::v_index(1, 2, k + 1)] = va1[1][k];
        xf.va[Xfoil::v_index(2, 1, k + 1)] = va2[0][k];
        xf.va[Xfoil::v_index(2, 2, k + 1)] = va2[1][k];
        xf.vb[Xfoil::v_index(2, 1, k + 1)] = vb2[0][k];
        xf.vb[Xfoil::v_index(2, 2, k + 1)] = vb2[1][k];
        xf.vm[Xfoil::vm_index(1, 1, k + 1, 2)] = m11[k];
        xf.vm[Xfoil::vm_index(2, 1, k + 1, 2)] = m12[k];
        xf.vm[Xfoil::vm_index(1, 2, k + 1, 2)] = m21[k];
        xf.vm[Xfoil::vm_index(2, 2, k + 1, 2)] = m22[k];
        xf.vdel[Xfoil::v_index(1, 1, k + 1)] = b1[k];
        xf.vdel[Xfoil::v_index(2, 1, k + 1)] = b2[k];
    }

    // dense assembly (row-major): A[row][col]
    let mut dense = [[0.0f64; 6]; 6];
    for r in 0..3 {
        // row 1 equations (block row 1)
        dense[r][0] = va1[0][r];
        dense[r][1] = va1[1][r];
        dense[r][2] = m11[r];
        dense[r][5] = m21[r];
        // row 2 equations (block row 2)
        dense[3 + r][0] = vb2[0][r];
        dense[3 + r][1] = vb2[1][r];
        dense[3 + r][2] = m12[r];
        dense[3 + r][3] = va2[0][r];
        dense[3 + r][4] = va2[1][r];
        dense[3 + r][5] = m22[r];
    }

    let rhs = [b1[0], b1[1], b1[2], b2[0], b2[1], b2[2]];

    (xf, dense, rhs)
}

#[test]
fn blsolv_matches_dense_solve() {
    let (mut xf, dense, rhs) = build_system();

    // reference solution via dense gauss
    let mut z = [0.0f64; 36];
    for c in 0..6 {
        for r in 0..6 {
            z[c * 6 + r] = dense[r][c];
        }
    }
    let mut r = rhs;
    gauss(&mut z, &mut r, 6, 6, 1);
    let expected = r;

    // blsolv solution (into VDEl(1..3, 1, 1..2))
    rust_foil::blsolv(&mut xf);

    let mut got = [0.0f64; 6];
    for k in 0..3 {
        got[k] = xf.vdel[Xfoil::v_index(1, 1, k + 1)];
        got[3 + k] = xf.vdel[Xfoil::v_index(2, 1, k + 1)];
    }

    for i in 0..6 {
        assert!(
            (got[i] - expected[i]).abs() < 1e-8,
            "mismatch at {}: blsolv={} dense={}",
            i,
            got[i],
            expected[i]
        );
    }
}

/// 3-station version exercising the "eliminate lower VM column" path
/// (kv = iv+2..NSYs) as well as the VB block elimination at two stations.
#[test]
fn blsolv_three_stations() {
    let nst = 3usize;
    let mut xf = Xfoil::new();
    xf.vaccel = 0.01;
    xf.n = 3;
    xf.s = vec![0.0, 0.5, 1.0];
    xf.nsys = nst;
    for (is, ibl) in [(0usize, 2usize), (0, 3), (1, 2), (1, 3), (1, 4)] {
        xf.isys[is][ibl] = 0;
    }
    xf.isys[0][2] = 1;
    xf.isys[0][3] = 2;
    xf.isys[1][2] = 3;
    xf.isys[1][5] = 3; // wake station pointer used by VZ elimination
    xf.iblte[0] = 3;
    xf.iblte[1] = 4;

    // blocks: VA[station][col][row], VB[station][col][row], M[station][masscol][row]
    let mut va = vec![[[0.0f64; 3]; 2]; nst + 1];
    let mut vb = vec![[[0.0f64; 3]; 2]; nst + 1];
    let mut m = vec![[[0.0f64; 3]; 3]; nst + 1];
    for st in 1..=nst {
        for c in 0..2 {
            for r in 0..3 {
                va[st][c][r] = 1.0 + st as f64 * 10.0 + (c * 3 + r) as f64;
                if st > 1 {
                    vb[st][c][r] = 0.1 + st as f64 + (c * 3 + r) as f64 * 0.1;
                }
            }
        }
        for mc in 1..=nst {
            for (r, mval) in m[st][mc - 1].iter_mut().enumerate() {
                *mval = 0.5 + st as f64 * 0.7 + mc as f64 * 0.3 + r as f64 * 0.05;
            }
        }
    }

    // dense: columns [dC1 dT1 dD1 dC2 dT2 dD2 dC3 dT3 dD3]
    let mut dense = [[0.0f64; 9]; 9];
    let mut rhs = [0.0f64; 9];
    for st in 1..=nst {
        // row block st: [VB(st) | VM(st,1) ... VA(st) ... ]
        for r in 0..3 {
            let row = (st - 1) * 3 + r;
            // VB columns: st-1's dC,dT
            if st > 1 {
                dense[row][(st - 2) * 3] = vb[st][0][r];
                dense[row][(st - 2) * 3 + 1] = vb[st][1][r];
            }
            // mass columns for all stations (upper triangular coupling)
            for mc in 1..=nst {
                dense[row][(mc - 1) * 3 + 2] = m[st][mc - 1][r];
            }
            // VA columns: st's dC,dT
            dense[row][(st - 1) * 3] = va[st][0][r];
            dense[row][(st - 1) * 3 + 1] = va[st][1][r];
            rhs[row] = (st as f64) * 10.0 + r as f64;
        }
    }

    // fill block arrays
    for st in 1..=nst {
        for c in 0..2 {
            for r in 0..3 {
                xf.va[Xfoil::v_index(st, c + 1, r + 1)] = va[st][c][r];
                if st > 1 {
                    xf.vb[Xfoil::v_index(st, c + 1, r + 1)] = vb[st][c][r];
                }
            }
        }
        for mc in 1..=nst {
            for (r, mval) in m[st][mc - 1].iter().enumerate() {
                xf.vm[Xfoil::vm_index(st, mc, r + 1, nst)] = *mval;
            }
        }
        for r in 0..3 {
            xf.vdel[Xfoil::v_index(st, 1, r + 1)] = rhs[(st - 1) * 3 + r];
        }
    }

    // reference via dense gauss
    let mut z = [0.0f64; 81];
    for c in 0..9 {
        for r in 0..9 {
            z[c * 9 + r] = dense[r][c];
        }
    }
    let mut r = rhs;
    gauss(&mut z, &mut r, 9, 9, 1);
    let expected = r;

    rust_foil::blsolv(&mut xf);

    let mut ok = true;
    for st in 1..=nst {
        for r in 0..3 {
            let got = xf.vdel[Xfoil::v_index(st, 1, r + 1)];
            let want = expected[(st - 1) * 3 + r];
            if (got - want).abs() > 1e-8 {
                eprintln!("mismatch st={} r={}: blsolv={} dense={}", st, r, got, want);
                ok = false;
            }
        }
    }
    assert!(ok, "3-station blsolv mismatch");
}
