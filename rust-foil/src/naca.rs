// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! NACA airfoil generation (port of `m_naca.f90`).
//!
//! Generates a 4- or 5-digit NACA airfoil into the buffer arrays (XB, YB).

/// TE point bunching parameter
const AN: f64 = 1.5;

/// Generates a 4-digit NACA airfoil (e.g. 2412 -> m=0.02, p=0.4, t=0.12).
/// Returns the airfoil name.
pub fn naca4(ides: i32, nside: usize) -> (Vec<f64>, Vec<f64>, usize, String) {
    let n4 = ides / 1000;
    let n3 = (ides - n4 * 1000) / 100;
    let n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;

    let m = n4 as f64 / 100.0;
    let p = n3 as f64 / 10.0;
    let t = (n2 * 10 + n1) as f64 / 100.0;

    let anp = AN + 1.0;
    let mut xx = vec![0.0f64; nside];
    let mut yt = vec![0.0f64; nside];
    let mut yc = vec![0.0f64; nside];
    for i in 0..nside {
        let frac = i as f64 / (nside - 1) as f64;
        if i == nside - 1 {
            xx[i] = 1.0;
        } else {
            xx[i] = 1.0 - anp * frac * (1.0 - frac).powf(AN) - (1.0 - frac).powf(anp);
        }
        yt[i] = (0.29690 * xx[i].sqrt() - 0.12600 * xx[i] - 0.35160 * xx[i].powi(2)
            + 0.28430 * xx[i].powi(3)
            - 0.10150 * xx[i].powi(4))
            * t
            / 0.20;
        if xx[i] < p {
            yc[i] = m / p.powi(2) * (2.0 * p * xx[i] - xx[i].powi(2));
        } else {
            yc[i] = m / (1.0 - p).powi(2) * ((1.0 - 2.0 * p) + 2.0 * p * xx[i] - xx[i].powi(2));
        }
    }

    let mut xb = vec![0.0f64; 2 * nside];
    let mut yb = vec![0.0f64; 2 * nside];
    let mut ib = 0;
    for i in (0..nside).rev() {
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
        ib += 1;
    }
    for i in 1..nside {
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
        ib += 1;
    }
    let nb = ib;

    let name = format!("NACA{}{}{}{}", n4, n3, n2, n1);
    (xb, yb, nb, name)
}

/// Generates a 5-digit NACA airfoil (e.g. 23012).  Returns the airfoil name,
/// or an empty name if the designation is illegal.
pub fn naca5(ides: i32, nside: usize, show_output: bool) -> (Vec<f64>, Vec<f64>, usize, String) {
    let n5 = ides / 10000;
    let n4 = (ides - n5 * 10000) / 1000;
    let n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    let n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;

    let n543 = 100 * n5 + 10 * n4 + n3;

    let (m, c);
    if n543 == 210 {
        m = 0.0580;
        c = 361.4;
    } else if n543 == 220 {
        m = 0.1260;
        c = 51.64;
    } else if n543 == 230 {
        m = 0.2025;
        c = 15.957;
    } else if n543 == 240 {
        m = 0.2900;
        c = 6.643;
    } else if n543 == 250 {
        m = 0.3910;
        c = 3.230;
    } else {
        if show_output {
            eprintln!("Illegal 5-digit designation");
            eprintln!("First three digits must be 210, 220, ... 250");
        }
        return (vec![0.0; 2 * nside], vec![0.0; 2 * nside], 0, String::new());
    }

    let t = (n2 * 10 + n1) as f64 / 100.0;

    let anp = AN + 1.0;
    let mut xx = vec![0.0f64; nside];
    let mut yt = vec![0.0f64; nside];
    let mut yc = vec![0.0f64; nside];
    for i in 0..nside {
        let frac = i as f64 / (nside - 1) as f64;
        if i == nside - 1 {
            xx[i] = 1.0;
        } else {
            xx[i] = 1.0 - anp * frac * (1.0 - frac).powf(AN) - (1.0 - frac).powf(anp);
        }
        yt[i] = (0.29690 * xx[i].sqrt() - 0.12600 * xx[i] - 0.35160 * xx[i].powi(2)
            + 0.28430 * xx[i].powi(3)
            - 0.10150 * xx[i].powi(4))
            * t
            / 0.20;
        if xx[i] < m {
            yc[i] =
                (c / 6.0) * (xx[i].powi(3) - 3.0 * m * xx[i].powi(2) + m * m * (3.0 - m) * xx[i]);
        } else {
            yc[i] = (c / 6.0) * m.powi(3) * (1.0 - xx[i]);
        }
    }

    let mut xb = vec![0.0f64; 2 * nside];
    let mut yb = vec![0.0f64; 2 * nside];
    let mut ib = 0;
    for i in (0..nside).rev() {
        xb[ib] = xx[i];
        yb[ib] = yc[i] + yt[i];
        ib += 1;
    }
    for i in 1..nside {
        xb[ib] = xx[i];
        yb[ib] = yc[i] - yt[i];
        ib += 1;
    }
    let nb = ib;

    let name = format!("NACA{}{}{}{}{}", n5, n4, n3, n2, n1);
    (xb, yb, nb, name)
}
