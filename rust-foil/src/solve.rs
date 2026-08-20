// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Linear solvers (port of `m_xsolve.f90`).
//!
//! Matrices use the Fortran column-major layout: element (i, j) [1-based] of
//! an `Nsiz x Nsiz` matrix lives at offset `(j-1)*Nsiz + (i-1)`.

use crate::state::Xfoil;

/// Solves general NxN system in NN unknowns with an arbitrary number (NRHS)
/// of righthand sides.  Assumes the system is invertible.  `z` is the
/// coefficient matrix (destroyed during solution), `r` the righthand side(s)
/// (replaced by the solution vector(s)).  Matrices use the Fortran
/// column-major layout: element (row i, col j) [1-based] at `(j-1)*nsiz + (i-1)`.
pub fn gauss(z: &mut [f64], r: &mut [f64], nsiz: usize, nn: usize, nrhs: usize) {
    for np in 0..nn - 1 {
        let np1 = np + 1;

        // find max pivot index NX
        let mut nx = np;
        for n in np1..nn {
            if z[np * nsiz + n].abs() > z[np * nsiz + nx].abs() {
                nx = n;
            }
        }

        let pivot = 1.0 / z[np * nsiz + nx];

        // switch pivots
        z[np * nsiz + nx] = z[np * nsiz + np];

        // switch rows & normalize pivot row
        for l in np1..nn {
            let temp = z[l * nsiz + nx] * pivot;
            z[l * nsiz + nx] = z[l * nsiz + np];
            z[l * nsiz + np] = temp;
        }

        for l in 0..nrhs {
            let temp = r[l * nsiz + nx] * pivot;
            r[l * nsiz + nx] = r[l * nsiz + np];
            r[l * nsiz + np] = temp;
        }

        // forward eliminate everything
        for k in np1..nn {
            let ztmp = z[np * nsiz + k];

            for l in np1..nn {
                z[l * nsiz + k] -= ztmp * z[l * nsiz + np];
            }
            for l in 0..nrhs {
                r[l * nsiz + k] -= ztmp * r[l * nsiz + np];
            }
        }
    }

    // solve for last row
    for l in 0..nrhs {
        r[l * nsiz + nn - 1] /= z[(nn - 1) * nsiz + nn - 1];
    }

    // back substitute everything
    for np in (0..nn - 1).rev() {
        let np1 = np + 1;
        for l in 0..nrhs {
            for k in np1..nn {
                r[l * nsiz + np] -= z[k * nsiz + np] * r[l * nsiz + k];
            }
        }
    }
}

/// Factors a full NxN matrix into an LU form.  Subr. `baksub` can
/// back-substitute it with some RHS.  `a` is replaced with its LU factors.
pub fn ludcmp(a: &mut [f64], indx: &mut [i32], n: usize) {
    let nsiz = (a.len() as f64).sqrt() as usize;
    let mut vv = vec![0.0f64; n];

    for i in 0..n {
        let mut aamax = 0.0;
        for j in 0..n {
            aamax = a[j * nsiz + i].abs().max(aamax);
        }
        vv[i] = 1.0 / aamax;
    }

    let mut imax = 0usize;
    for j in 0..n {
        for i in 0..j {
            let mut sum = a[j * nsiz + i];
            for k in 0..i {
                sum -= a[k * nsiz + i] * a[j * nsiz + k];
            }
            a[j * nsiz + i] = sum;
        }

        let mut aamax = 0.0;
        for i in j..n {
            let mut sum = a[j * nsiz + i];
            for k in 0..j {
                sum -= a[k * nsiz + i] * a[j * nsiz + k];
            }
            a[j * nsiz + i] = sum;

            let dum = vv[i] * sum.abs();
            if dum >= aamax {
                imax = i;
                aamax = dum;
            }
        }

        if j != imax {
            for k in 0..n {
                a.swap(k * nsiz + imax, k * nsiz + j);
            }
            vv[imax] = vv[j];
        }

        indx[j] = imax as i32;
        if j != n - 1 {
            let dum = 1.0 / a[j * nsiz + j];
            for i in j + 1..n {
                a[j * nsiz + i] *= dum;
            }
        }
    }
}

/// Back-substitutes a previously LU-factored matrix (see [`ludcmp`]).
pub fn baksub(a: &[f64], indx: &[i32], b: &mut [f64], n: usize) {
    let nsiz = (a.len() as f64).sqrt() as usize;

    // First index at which b became nonzero, or -1 if none has yet.  Using a
    // sentinel outside the valid index range (rather than 0) is required for
    // correctness: if the first nonzero entry appears at i == 0, the loop
    // suppression must still be enabled for all later rows.  The original
    // Numerical-Recipes implementation used ii = 0 as "not set", which is wrong.
    let mut ii: isize = -1;
    for i in 0..n {
        let ll = indx[i] as usize;
        let mut sum = b[ll];
        b[ll] = b[i];
        if ii >= 0 {
            for j in ii as usize..i {
                sum -= a[j * nsiz + i] * b[j];
            }
        } else if sum != 0.0 {
            ii = i as isize;
        }
        b[i] = sum;
    }

    for i in (0..n).rev() {
        let mut sum = b[i];
        if i < n - 1 {
            for j in i + 1..n {
                sum -= a[j * nsiz + i] * b[j];
            }
        }
        b[i] = sum / a[i * nsiz + i];
    }
}

/// Custom solver for the coupled viscous-inviscid Newton system:
///
/// ```text
///   A  |  |  .  |  |  .  |    d       R       S
///   B  A  |  .  |  |  .  |    d       R       S
///   |  B  A  .  |  |  .  |    d       R       S
///   .  .  .  .  |  |  .  |    d   =   R - dRe S
///   |  |  |  B  A  |  .  |    d       R       S
///   |  Z  |  |  B  A  .  |    d       R       S
///   .  .  .  .  .  .  .  |    d       R       S
///   |  |  |  |  |  |  B  A    d       R       S
/// ```
///
/// A, B, Z are 3x3 blocks; | are 3x1 mass-defect influence vectors; d, R, S
/// are 3x1 unknown, residual and Re-influence vectors.  All block indices in
/// this routine are 1-based (matching the Fortran), converted by the state
/// index helpers.
pub fn blsolv(xf: &mut Xfoil) {
    let ivte1 = xf.isys[0][xf.iblte[0] as usize] as usize;

    let vacc1 = xf.vaccel;
    let vacc2 = xf.vaccel * 2.0 / (xf.s[xf.n - 1] - xf.s[0]);
    let vacc3 = xf.vaccel * 2.0 / (xf.s[xf.n - 1] - xf.s[0]);

    for iv in 1..=xf.nsys {
        let ivp = iv + 1;

        // ====== Invert VA(IV) block ======

        // normalize first row
        let mut pivot = 1.0 / xf.va[Xfoil::v_index(iv, 1, 1)];
        xf.va[Xfoil::v_index(iv, 2, 1)] *= pivot;
        for l in iv..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 1)] *= pivot;
        }
        xf.vdel[Xfoil::v_index(iv, 1, 1)] *= pivot;
        xf.vdel[Xfoil::v_index(iv, 2, 1)] *= pivot;

        // eliminate lower first column in VA block
        for k in 2..=3 {
            let vtmp = xf.va[Xfoil::v_index(iv, 1, k)];
            xf.va[Xfoil::v_index(iv, 2, k)] -= vtmp * xf.va[Xfoil::v_index(iv, 2, 1)];
            for l in iv..=xf.nsys {
                xf.vm[Xfoil::vm_index(iv, l, k)] -= vtmp * xf.vm[Xfoil::vm_index(iv, l, 1)];
            }
            xf.vdel[Xfoil::v_index(iv, 1, k)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 1, 1)];
            xf.vdel[Xfoil::v_index(iv, 2, k)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 2, 1)];
        }

        // normalize second row
        pivot = 1.0 / xf.va[Xfoil::v_index(iv, 2, 2)];
        for l in iv..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 2)] *= pivot;
        }
        xf.vdel[Xfoil::v_index(iv, 1, 2)] *= pivot;
        xf.vdel[Xfoil::v_index(iv, 2, 2)] *= pivot;

        // eliminate lower second column in VA block
        let vtmp = xf.va[Xfoil::v_index(iv, 2, 3)];
        for l in iv..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 3)] -= vtmp * xf.vm[Xfoil::vm_index(iv, l, 2)];
        }
        xf.vdel[Xfoil::v_index(iv, 1, 3)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 1, 2)];
        xf.vdel[Xfoil::v_index(iv, 2, 3)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 2, 2)];

        // normalize third row
        pivot = 1.0 / xf.vm[Xfoil::vm_index(iv, iv, 3)];
        for l in ivp..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 3)] *= pivot;
        }
        xf.vdel[Xfoil::v_index(iv, 1, 3)] *= pivot;
        xf.vdel[Xfoil::v_index(iv, 2, 3)] *= pivot;

        // eliminate upper third column in VA block
        let vtmp1 = xf.vm[Xfoil::vm_index(iv, iv, 1)];
        let vtmp2 = xf.vm[Xfoil::vm_index(iv, iv, 2)];
        for l in ivp..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 1)] -= vtmp1 * xf.vm[Xfoil::vm_index(iv, l, 3)];
            xf.vm[Xfoil::vm_index(iv, l, 2)] -= vtmp2 * xf.vm[Xfoil::vm_index(iv, l, 3)];
        }
        xf.vdel[Xfoil::v_index(iv, 1, 1)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
        xf.vdel[Xfoil::v_index(iv, 1, 2)] -= vtmp2 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
        xf.vdel[Xfoil::v_index(iv, 2, 1)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 2, 3)];
        xf.vdel[Xfoil::v_index(iv, 2, 2)] -= vtmp2 * xf.vdel[Xfoil::v_index(iv, 2, 3)];

        // eliminate upper second column in VA block
        let vtmp = xf.va[Xfoil::v_index(iv, 2, 1)];
        for l in ivp..=xf.nsys {
            xf.vm[Xfoil::vm_index(iv, l, 1)] -= vtmp * xf.vm[Xfoil::vm_index(iv, l, 2)];
        }
        xf.vdel[Xfoil::v_index(iv, 1, 1)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 1, 2)];
        xf.vdel[Xfoil::v_index(iv, 2, 1)] -= vtmp * xf.vdel[Xfoil::v_index(iv, 2, 2)];

        if iv != xf.nsys {
            // ====== Eliminate VB(IV+1) block, rows 1 -> 3 ======
            for k in 1..=3 {
                let vtmp1 = xf.vb[Xfoil::v_index(ivp, 1, k)];
                let vtmp2 = xf.vb[Xfoil::v_index(ivp, 2, k)];
                let vtmp3 = xf.vm[Xfoil::vm_index(ivp, iv, k)];
                for l in ivp..=xf.nsys {
                    xf.vm[Xfoil::vm_index(ivp, l, k)] -= vtmp1 * xf.vm[Xfoil::vm_index(iv, l, 1)]
                        + vtmp2 * xf.vm[Xfoil::vm_index(iv, l, 2)]
                        + vtmp3 * xf.vm[Xfoil::vm_index(iv, l, 3)];
                }
                xf.vdel[Xfoil::v_index(ivp, 1, k)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 1, 1)]
                    + vtmp2 * xf.vdel[Xfoil::v_index(iv, 1, 2)]
                    + vtmp3 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
                xf.vdel[Xfoil::v_index(ivp, 2, k)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 2, 1)]
                    + vtmp2 * xf.vdel[Xfoil::v_index(iv, 2, 2)]
                    + vtmp3 * xf.vdel[Xfoil::v_index(iv, 2, 3)];
            }

            if iv == ivte1 {
                // eliminate VZ block
                let ivz = xf.isys[1][xf.iblte[1] as usize + 1] as usize;

                for k in 1..=3 {
                    let vtmp1 = xf.vz[Xfoil::v_index(1, 1, k)];
                    let vtmp2 = xf.vz[Xfoil::v_index(1, 2, k)];
                    for l in ivp..=xf.nsys {
                        xf.vm[Xfoil::vm_index(ivz, l, k)] -= vtmp1
                            * xf.vm[Xfoil::vm_index(iv, l, 1)]
                            + vtmp2 * xf.vm[Xfoil::vm_index(iv, l, 2)];
                    }
                    xf.vdel[Xfoil::v_index(ivz, 1, k)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 1, 1)]
                        + vtmp2 * xf.vdel[Xfoil::v_index(iv, 1, 2)];
                    xf.vdel[Xfoil::v_index(ivz, 2, k)] -= vtmp1 * xf.vdel[Xfoil::v_index(iv, 2, 1)]
                        + vtmp2 * xf.vdel[Xfoil::v_index(iv, 2, 2)];
                }
            }

            if ivp != xf.nsys {
                // ====== Eliminate lower VM column ======
                for kv in iv + 2..=xf.nsys {
                    let vtmp1 = xf.vm[Xfoil::vm_index(kv, iv, 1)];
                    let vtmp2 = xf.vm[Xfoil::vm_index(kv, iv, 2)];
                    let vtmp3 = xf.vm[Xfoil::vm_index(kv, iv, 3)];

                    if vtmp1.abs() > vacc1 {
                        for l in ivp..=xf.nsys {
                            xf.vm[Xfoil::vm_index(kv, l, 1)] -=
                                vtmp1 * xf.vm[Xfoil::vm_index(iv, l, 3)];
                        }
                        xf.vdel[Xfoil::v_index(kv, 1, 1)] -=
                            vtmp1 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
                        xf.vdel[Xfoil::v_index(kv, 2, 1)] -=
                            vtmp1 * xf.vdel[Xfoil::v_index(iv, 2, 3)];
                    }

                    if vtmp2.abs() > vacc2 {
                        for l in ivp..=xf.nsys {
                            xf.vm[Xfoil::vm_index(kv, l, 2)] -=
                                vtmp2 * xf.vm[Xfoil::vm_index(iv, l, 3)];
                        }
                        xf.vdel[Xfoil::v_index(kv, 1, 2)] -=
                            vtmp2 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
                        xf.vdel[Xfoil::v_index(kv, 2, 2)] -=
                            vtmp2 * xf.vdel[Xfoil::v_index(iv, 2, 3)];
                    }

                    if vtmp3.abs() > vacc3 {
                        for l in ivp..=xf.nsys {
                            xf.vm[Xfoil::vm_index(kv, l, 3)] -=
                                vtmp3 * xf.vm[Xfoil::vm_index(iv, l, 3)];
                        }
                        xf.vdel[Xfoil::v_index(kv, 1, 3)] -=
                            vtmp3 * xf.vdel[Xfoil::v_index(iv, 1, 3)];
                        xf.vdel[Xfoil::v_index(kv, 2, 3)] -=
                            vtmp3 * xf.vdel[Xfoil::v_index(iv, 2, 3)];
                    }
                }
            }
        }
    }

    for iv in (2..=xf.nsys).rev() {
        // eliminate upper VM columns
        let vtmp = xf.vdel[Xfoil::v_index(iv, 1, 3)];
        for kv in (1..iv).rev() {
            xf.vdel[Xfoil::v_index(kv, 1, 1)] -= xf.vm[Xfoil::vm_index(kv, iv, 1)] * vtmp;
            xf.vdel[Xfoil::v_index(kv, 1, 2)] -= xf.vm[Xfoil::vm_index(kv, iv, 2)] * vtmp;
            xf.vdel[Xfoil::v_index(kv, 1, 3)] -= xf.vm[Xfoil::vm_index(kv, iv, 3)] * vtmp;
        }

        let vtmp = xf.vdel[Xfoil::v_index(iv, 2, 3)];
        for kv in (1..iv).rev() {
            xf.vdel[Xfoil::v_index(kv, 2, 1)] -= xf.vm[Xfoil::vm_index(kv, iv, 1)] * vtmp;
            xf.vdel[Xfoil::v_index(kv, 2, 2)] -= xf.vm[Xfoil::vm_index(kv, iv, 2)] * vtmp;
            xf.vdel[Xfoil::v_index(kv, 2, 3)] -= xf.vm[Xfoil::vm_index(kv, iv, 3)] * vtmp;
        }
    }
}
