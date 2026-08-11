//! Regression tests for the LU factor/back-substitution pair
//! (`ludcmp`/`baksub`).
//!
//! These tests exist primarily to lock down the forward-substitution
//! "first nonzero rhs" logic in `baksub`, which previously used `0` as the
//! "no nonzero yet" sentinel.  With that buggy sentinel the forward-elimination
//! inner loop was skipped whenever the first nonzero entry of `b` landed at
//! index 0, silently producing wrong solutions.
//!
//! Each test builds a small system A x = b, LU-factorizes A with `ludcmp`, and
//! solves with `baksub`.  The expected solution `x` is computed by an
//! independent dense Gaussian elimination over row-major storage, so the test
//! is not just checking `baksub` against itself.

use rust_foil::{baksub, ludcmp};

/// Reference solver: standard partial-pivot Gaussian elimination on a dense
/// row-major `n x n` matrix.  Solves A x = b in place, returning x.
fn dense_solve(a_in: [[f64; 4]; 4], b_in: [f64; 4], n: usize) -> [f64; 4] {
    let mut a = a_in;
    let mut b = b_in;

    for k in 0..n {
        // partial pivot
        let mut p = k;
        for r in (k + 1)..n {
            if a[r][k].abs() > a[p][k].abs() {
                p = r;
            }
        }
        if p != k {
            a.swap(p, k);
            b.swap(p, k);
        }
        // eliminate
        for r in (k + 1)..n {
            let m = a[r][k] / a[k][k];
            for c in k..n {
                a[r][c] -= m * a[k][c];
            }
            b[r] -= m * b[k];
        }
    }

    let mut x = [0.0f64; 4];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    x
}

/// Packs a row-major `n x n` matrix into the Fortran column-major layout used
/// by the solver (`element (i,j) at j*nsiz + i`), with `nsiz = 4`.
fn pack_col_major(a: [[f64; 4]; 4], n: usize) -> Vec<f64> {
    let mut z = vec![0.0f64; 16];
    for i in 0..n {
        for j in 0..n {
            z[j * 4 + i] = a[i][j];
        }
    }
    z
}

fn check(a: [[f64; 4]; 4], b: [f64; 4], n: usize) {
    let expected = dense_solve(a, b, n);

    let mut z = pack_col_major(a, n);
    let mut indx = vec![0i32; n];
    ludcmp(&mut z, &mut indx, n);

    let mut x = vec![0.0f64; 4];
    for i in 0..n {
        x[i] = b[i];
    }
    baksub(&z, &indx, &mut x, n);

    for i in 0..n {
        assert!(
            (x[i] - expected[i]).abs() < 1.0e-9,
            "solution component {}: baksub = {}, dense = {}",
            i,
            x[i],
            expected[i]
        );
    }
}

#[test]
fn baksub_rhs_nonzero_at_zero() {
    // This is the case the old sentinel bug got wrong: the first nonzero entry
    // of b is at index 0.  With a non-trivial pivot row order the buggy
    // implementation skipped forward elimination for every subsequent row.
    //
    //   [ 2 1 0 ] [x]   [ 5 ]
    //   [ 1 3 1 ] [y] = [10 ]
    //   [ 0 1 2 ] [z]   [ 5 ]
    //
    // b[0] = 5 != 0 -> the trigger.  Expected solution: [1.25, 2.5, 1.25].
    let a = [
        [2.0, 1.0, 0.0, 0.0],
        [1.0, 3.0, 1.0, 0.0],
        [0.0, 1.0, 2.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    let b = [5.0, 10.0, 5.0, 0.0];
    check(a, b, 3);

    // Confirm the actual solution is the expected one, not just "matches an
    // independent buggy solver" -- both implementations here are correct, but
    // pinning the numeric answer catches future regressions directly.
    let mut z = pack_col_major(a, 3);
    let mut indx = vec![0i32; 3];
    ludcmp(&mut z, &mut indx, 3);
    let mut x = vec![5.0_f64, 10.0, 5.0];
    baksub(&z, &indx, &mut x, 3);
    assert!((x[0] - 1.25).abs() < 1.0e-9, "x[0] = {}", x[0]);
    assert!((x[1] - 2.50).abs() < 1.0e-9, "x[1] = {}", x[1]);
    assert!((x[2] - 1.25).abs() < 1.0e-9, "x[2] = {}", x[2]);
}

#[test]
fn baksub_first_rhs_zero() {
    // Complementary case: b[0] = 0 so the first nonzero entry appears later.
    // This case worked even with the old sentinel; included so the suite covers
    // both branches.
    let a = [
        [4.0, 1.0, 2.0, 0.0],
        [1.0, 3.0, 1.0, 0.0],
        [2.0, 1.0, 5.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    let b = [0.0, -2.5, 5.5, 0.0];
    check(a, b, 3);
}

#[test]
fn baksub_requires_pivot_from_zero() {
    // System that forces a row swap to bring a nonzero pivot into row 0:
    //   [ 0 1 ] [x]   [2]
    //   [ 1 0 ] [y] = [1]
    // Expected: x = 1, y = 2.
    let a = [
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ];
    let b = [2.0, 1.0, 0.0, 0.0];
    check(a, b, 2);

    let mut z = pack_col_major(a, 2);
    let mut indx = vec![0i32; 2];
    ludcmp(&mut z, &mut indx, 2);
    let mut x = vec![2.0_f64, 1.0];
    baksub(&z, &indx, &mut x, 2);
    assert!((x[0] - 1.0).abs() < 1.0e-9, "x = {}", x[0]);
    assert!((x[1] - 2.0).abs() < 1.0e-9, "y = {}", x[1]);
}

#[test]
fn baksub_identity() {
    let a = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let b = [7.0, -3.0, 4.5, 2.0];
    check(a, b, 4);
}

#[test]
fn baksub_4x4_all_rhs_nonzero() {
    // Larger dense system where every b[i] is nonzero from the start -- the
    // strongest stress test for the forward-elimination loop.
    let a = [
        [3.0, 1.0, -1.0, 2.0],
        [1.0, 4.0, 1.0, -1.0],
        [-1.0, 1.0, 5.0, 1.0],
        [2.0, -1.0, 1.0, 6.0],
    ];
    let b = [1.0, 2.0, 3.0, 4.0];
    check(a, b, 4);
}
