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

/// Flat index into the 3x2xIZX BL block arrays (va/vb/vz/vdel).  All indices
/// are 1-based, matching the Fortran `va(3, 2, IZX)` / `vdel(3, 2, IZX)`
/// layout.
#[inline(always)]
fn vi(iv: usize, l: usize, k: usize) -> usize {
    ((iv - 1) * 2 + (l - 1)) * 3 + (k - 1)
}

/// Returns mutable slices for two distinct rows of one VM plane over the
/// same column range `a..b` (0-based offsets within the plane).  `ra`/`rb`
/// are plane-relative row base offsets; the first returned slice is row `ra`,
/// the second row `rb`.
#[inline(always)]
fn two_rows(p: &mut [f64], ra: usize, rb: usize, a: usize, b: usize) -> (&mut [f64], &mut [f64]) {
    debug_assert!(ra != rb);
    if ra < rb {
        let (h, t) = p.split_at_mut(rb);
        (&mut h[ra + a..ra + b], &mut t[a..b])
    } else {
        let (h, t) = p.split_at_mut(ra);
        (&mut t[a..b], &mut h[rb + a..rb + b])
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
/// this routine are 1-based (matching the Fortran).
///
/// The `vm` array is stored as three row-major nsys x nsys planes packed by
/// the current system size (one per component k; see [`Xfoil::vm_index`]),
/// so every row sweep below is a
/// unit-stride slice axpy/scale that auto-vectorizes.  The elimination and
/// back-substitution perform exactly the same floating-point operations in
/// the same order per element as the scalar Fortran translation, so results
/// are bit-identical.
pub fn blsolv(xf: &mut Xfoil) {
    let nsys = xf.nsys;
    let ivte1 = xf.isys[0][xf.iblte[0] as usize] as usize;

    let vacc1 = xf.vaccel;
    let vacc2 = xf.vaccel * 2.0 / (xf.s[xf.n - 1] - xf.s[0]);
    let vacc3 = xf.vaccel * 2.0 / (xf.s[xf.n - 1] - xf.s[0]);

    // vm is packed by nsys: row stride nsys, plane size nsys*nsys
    let stride = nsys;
    let plane_sz: usize = stride * stride;
    let (m1, rest) = xf.vm.split_at_mut(plane_sz);
    let (m2, m3) = rest.split_at_mut(plane_sz);
    let va = &mut xf.va;
    let vb = &mut xf.vb;
    let vz = &mut xf.vz;
    let vdel = &mut xf.vdel;

    for iv in 1..=nsys {
        let ivp = iv + 1;
        let r = (iv - 1) * stride; // row-iv base offset within a plane
        let rp = iv * stride; // row-ivp base offset

        // ====== Invert VA(IV) block ======

        // normalize first row
        let mut pivot = 1.0 / va[vi(iv, 1, 1)];
        va[vi(iv, 2, 1)] *= pivot;
        for v in &mut m1[r + iv - 1..r + nsys] {
            *v *= pivot;
        }
        vdel[vi(iv, 1, 1)] *= pivot;
        vdel[vi(iv, 2, 1)] *= pivot;

        // eliminate lower first column in VA block, k = 2
        let vtmp = va[vi(iv, 1, 2)];
        va[vi(iv, 2, 2)] -= vtmp * va[vi(iv, 2, 1)];
        for (d, s) in m2[r + iv - 1..r + nsys]
            .iter_mut()
            .zip(&m1[r + iv - 1..r + nsys])
        {
            *d -= vtmp * s;
        }
        vdel[vi(iv, 1, 2)] -= vtmp * vdel[vi(iv, 1, 1)];
        vdel[vi(iv, 2, 2)] -= vtmp * vdel[vi(iv, 2, 1)];

        // eliminate lower first column in VA block, k = 3
        let vtmp = va[vi(iv, 1, 3)];
        va[vi(iv, 2, 3)] -= vtmp * va[vi(iv, 2, 1)];
        for (d, s) in m3[r + iv - 1..r + nsys]
            .iter_mut()
            .zip(&m1[r + iv - 1..r + nsys])
        {
            *d -= vtmp * s;
        }
        vdel[vi(iv, 1, 3)] -= vtmp * vdel[vi(iv, 1, 1)];
        vdel[vi(iv, 2, 3)] -= vtmp * vdel[vi(iv, 2, 1)];

        // normalize second row
        pivot = 1.0 / va[vi(iv, 2, 2)];
        for v in &mut m2[r + iv - 1..r + nsys] {
            *v *= pivot;
        }
        vdel[vi(iv, 1, 2)] *= pivot;
        vdel[vi(iv, 2, 2)] *= pivot;

        // eliminate lower second column in VA block
        let vtmp = va[vi(iv, 2, 3)];
        for (d, s) in m3[r + iv - 1..r + nsys]
            .iter_mut()
            .zip(&m2[r + iv - 1..r + nsys])
        {
            *d -= vtmp * s;
        }
        vdel[vi(iv, 1, 3)] -= vtmp * vdel[vi(iv, 1, 2)];
        vdel[vi(iv, 2, 3)] -= vtmp * vdel[vi(iv, 2, 2)];

        // normalize third row
        pivot = 1.0 / m3[r + iv - 1];
        for v in &mut m3[r + iv..r + nsys] {
            *v *= pivot;
        }
        vdel[vi(iv, 1, 3)] *= pivot;
        vdel[vi(iv, 2, 3)] *= pivot;

        // eliminate upper third column in VA block
        let vtmp1 = m1[r + iv - 1];
        let vtmp2 = m2[r + iv - 1];
        {
            let src = &m3[r + iv..r + nsys];
            for ((d1, d2), s) in m1[r + iv..r + nsys]
                .iter_mut()
                .zip(m2[r + iv..r + nsys].iter_mut())
                .zip(src)
            {
                *d1 -= vtmp1 * s;
                *d2 -= vtmp2 * s;
            }
        }
        vdel[vi(iv, 1, 1)] -= vtmp1 * vdel[vi(iv, 1, 3)];
        vdel[vi(iv, 1, 2)] -= vtmp2 * vdel[vi(iv, 1, 3)];
        vdel[vi(iv, 2, 1)] -= vtmp1 * vdel[vi(iv, 2, 3)];
        vdel[vi(iv, 2, 2)] -= vtmp2 * vdel[vi(iv, 2, 3)];

        // eliminate upper second column in VA block
        let vtmp = va[vi(iv, 2, 1)];
        for (d, s) in m1[r + iv..r + nsys]
            .iter_mut()
            .zip(&m2[r + iv..r + nsys])
        {
            *d -= vtmp * s;
        }
        vdel[vi(iv, 1, 1)] -= vtmp * vdel[vi(iv, 1, 2)];
        vdel[vi(iv, 2, 1)] -= vtmp * vdel[vi(iv, 2, 2)];

        if iv != nsys {
            // ====== Eliminate VB(IV+1) block, rows 1 -> 3 ======
            // k = 1 (destination row ivp of plane 1)
            {
                let vtmp1 = vb[vi(ivp, 1, 1)];
                let vtmp2 = vb[vi(ivp, 2, 1)];
                let vtmp3 = m1[rp + iv - 1];
                let (sa, dst) = two_rows(&mut *m1, r, rp, iv, nsys);
                for (((d, a), b), c) in dst
                    .iter_mut()
                    .zip(&*sa)
                    .zip(&m2[r + iv..r + nsys])
                    .zip(&m3[r + iv..r + nsys])
                {
                    *d -= vtmp1 * a + vtmp2 * b + vtmp3 * c;
                }
                vdel[vi(ivp, 1, 1)] -= vtmp1 * vdel[vi(iv, 1, 1)]
                    + vtmp2 * vdel[vi(iv, 1, 2)]
                    + vtmp3 * vdel[vi(iv, 1, 3)];
                vdel[vi(ivp, 2, 1)] -= vtmp1 * vdel[vi(iv, 2, 1)]
                    + vtmp2 * vdel[vi(iv, 2, 2)]
                    + vtmp3 * vdel[vi(iv, 2, 3)];
            }
            // k = 2
            {
                let vtmp1 = vb[vi(ivp, 1, 2)];
                let vtmp2 = vb[vi(ivp, 2, 2)];
                let vtmp3 = m2[rp + iv - 1];
                let (sb, dst) = two_rows(&mut *m2, r, rp, iv, nsys);
                for (((d, a), b), c) in dst
                    .iter_mut()
                    .zip(&m1[r + iv..r + nsys])
                    .zip(&*sb)
                    .zip(&m3[r + iv..r + nsys])
                {
                    *d -= vtmp1 * a + vtmp2 * b + vtmp3 * c;
                }
                vdel[vi(ivp, 1, 2)] -= vtmp1 * vdel[vi(iv, 1, 1)]
                    + vtmp2 * vdel[vi(iv, 1, 2)]
                    + vtmp3 * vdel[vi(iv, 1, 3)];
                vdel[vi(ivp, 2, 2)] -= vtmp1 * vdel[vi(iv, 2, 1)]
                    + vtmp2 * vdel[vi(iv, 2, 2)]
                    + vtmp3 * vdel[vi(iv, 2, 3)];
            }
            // k = 3
            {
                let vtmp1 = vb[vi(ivp, 1, 3)];
                let vtmp2 = vb[vi(ivp, 2, 3)];
                let vtmp3 = m3[rp + iv - 1];
                let (sc, dst) = two_rows(&mut *m3, r, rp, iv, nsys);
                for (((d, a), b), c) in dst
                    .iter_mut()
                    .zip(&m1[r + iv..r + nsys])
                    .zip(&m2[r + iv..r + nsys])
                    .zip(&*sc)
                {
                    *d -= vtmp1 * a + vtmp2 * b + vtmp3 * c;
                }
                vdel[vi(ivp, 1, 3)] -= vtmp1 * vdel[vi(iv, 1, 1)]
                    + vtmp2 * vdel[vi(iv, 1, 2)]
                    + vtmp3 * vdel[vi(iv, 1, 3)];
                vdel[vi(ivp, 2, 3)] -= vtmp1 * vdel[vi(iv, 2, 1)]
                    + vtmp2 * vdel[vi(iv, 2, 2)]
                    + vtmp3 * vdel[vi(iv, 2, 3)];
            }

            if iv == ivte1 {
                // eliminate VZ block
                let ivz = xf.isys[1][xf.iblte[1] as usize + 1] as usize;
                let rz = (ivz - 1) * stride;

                // k = 1
                {
                    let vtmp1 = vz[vi(1, 1, 1)];
                    let vtmp2 = vz[vi(1, 2, 1)];
                    let (sa, dst) = two_rows(&mut *m1, r, rz, iv, nsys);
                    for ((d, a), b) in dst.iter_mut().zip(&*sa).zip(&m2[r + iv..r + nsys]) {
                        *d -= vtmp1 * a + vtmp2 * b;
                    }
                    vdel[vi(ivz, 1, 1)] -=
                        vtmp1 * vdel[vi(iv, 1, 1)] + vtmp2 * vdel[vi(iv, 1, 2)];
                    vdel[vi(ivz, 2, 1)] -=
                        vtmp1 * vdel[vi(iv, 2, 1)] + vtmp2 * vdel[vi(iv, 2, 2)];
                }
                // k = 2
                {
                    let vtmp1 = vz[vi(1, 1, 2)];
                    let vtmp2 = vz[vi(1, 2, 2)];
                    let (sb, dst) = two_rows(&mut *m2, r, rz, iv, nsys);
                    for ((d, a), b) in dst.iter_mut().zip(&m1[r + iv..r + nsys]).zip(&*sb) {
                        *d -= vtmp1 * a + vtmp2 * b;
                    }
                    vdel[vi(ivz, 1, 2)] -=
                        vtmp1 * vdel[vi(iv, 1, 1)] + vtmp2 * vdel[vi(iv, 1, 2)];
                    vdel[vi(ivz, 2, 2)] -=
                        vtmp1 * vdel[vi(iv, 2, 1)] + vtmp2 * vdel[vi(iv, 2, 2)];
                }
                // k = 3
                {
                    let vtmp1 = vz[vi(1, 1, 3)];
                    let vtmp2 = vz[vi(1, 2, 3)];
                    let (sc, dst) = two_rows(&mut *m3, r, rz, iv, nsys);
                    for ((d, a), b) in dst.iter_mut().zip(&m1[r + iv..r + nsys]).zip(&*sc) {
                        *d -= vtmp1 * a + vtmp2 * b;
                    }
                    vdel[vi(ivz, 1, 3)] -=
                        vtmp1 * vdel[vi(iv, 1, 1)] + vtmp2 * vdel[vi(iv, 1, 2)];
                    vdel[vi(ivz, 2, 3)] -=
                        vtmp1 * vdel[vi(iv, 2, 1)] + vtmp2 * vdel[vi(iv, 2, 2)];
                }
            }

            if ivp != nsys {
                // ====== Eliminate lower VM column ======
                // The vdel[iv] source terms are loop-invariant here (only
                // vdel[kv] entries are written), so hoist them.
                let d13 = vdel[vi(iv, 1, 3)];
                let d23 = vdel[vi(iv, 2, 3)];

                for kv in iv + 2..=nsys {
                    let c = (kv - 1) * stride + iv - 1;
                    let vtmp1 = m1[c];
                    let vtmp2 = m2[c];
                    let vtmp3 = m3[c];

                    let h1 = vtmp1.abs() > vacc1;
                    let h2 = vtmp2.abs() > vacc2;
                    let h3 = vtmp3.abs() > vacc3;
                    if !(h1 || h2 || h3) {
                        continue;
                    }


                    let rk = (kv - 1) * stride;
                    // all three components eliminate against row iv of plane 3
                    let (src3, dst3) = two_rows(&mut *m3, r, rk, iv, nsys);
                    let src3: &[f64] = src3;

                    if h1 && h2 && h3 {
                        // common case: all three gates hit — fuse into one
                        // pass so the shared source row is loaded once
                        for ((d1, d2), (d3, s)) in m1[rk + iv..rk + nsys]
                            .iter_mut()
                            .zip(m2[rk + iv..rk + nsys].iter_mut())
                            .zip(dst3.iter_mut().zip(src3))
                        {
                            *d1 -= vtmp1 * s;
                            *d2 -= vtmp2 * s;
                            *d3 -= vtmp3 * s;
                        }
                        vdel[vi(kv, 1, 1)] -= vtmp1 * d13;
                        vdel[vi(kv, 2, 1)] -= vtmp1 * d23;
                        vdel[vi(kv, 1, 2)] -= vtmp2 * d13;
                        vdel[vi(kv, 2, 2)] -= vtmp2 * d23;
                        vdel[vi(kv, 1, 3)] -= vtmp3 * d13;
                        vdel[vi(kv, 2, 3)] -= vtmp3 * d23;
                        continue;
                    }

                    if h1 {
                        for (d, s) in m1[rk + iv..rk + nsys].iter_mut().zip(src3) {
                            *d -= vtmp1 * s;
                        }
                        vdel[vi(kv, 1, 1)] -= vtmp1 * d13;
                        vdel[vi(kv, 2, 1)] -= vtmp1 * d23;
                    }

                    if h2 {
                        for (d, s) in m2[rk + iv..rk + nsys].iter_mut().zip(src3) {
                            *d -= vtmp2 * s;
                        }
                        vdel[vi(kv, 1, 2)] -= vtmp2 * d13;
                        vdel[vi(kv, 2, 2)] -= vtmp2 * d23;
                    }

                    if h3 {
                        for (d, s) in dst3.iter_mut().zip(src3) {
                            *d -= vtmp3 * s;
                        }
                        vdel[vi(kv, 1, 3)] -= vtmp3 * d13;
                        vdel[vi(kv, 2, 3)] -= vtmp3 * d23;
                    }
                }
            }
        }
    }

    // ====== Back-substitution: eliminate upper VM columns ======
    // Loops are interchanged vs the Fortran (kv outer, iv inner descending);
    // each vdel[kv] element still accumulates its updates in the original
    // iv-descending order, so results are bit-identical, while the vm reads
    // become unit-stride row sweeps.
    for kv in (1..nsys).rev() {
        let rk = (kv - 1) * stride;
        let s1 = &m1[rk + kv..rk + nsys];
        let s2 = &m2[rk + kv..rk + nsys];
        let s3 = &m3[rk + kv..rk + nsys];
        let b = (kv - 1) * 6;
        for j in (0..nsys - kv).rev() {
            let t1 = vdel[(kv + j) * 6 + 2];
            let t2 = vdel[(kv + j) * 6 + 5];
            vdel[b] -= s1[j] * t1;
            vdel[b + 1] -= s2[j] * t1;
            vdel[b + 2] -= s3[j] * t1;
            vdel[b + 3] -= s1[j] * t2;
            vdel[b + 4] -= s2[j] * t2;
            vdel[b + 5] -= s3[j] * t2;
        }
    }
}
