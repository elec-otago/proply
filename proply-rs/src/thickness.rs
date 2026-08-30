// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Mechanical blade-thickness sizing: choose the airfoil thickness so the
//! blade — treated as a cantilever beam anchored at the hub — does not
//! deform too much under its own thrust.
//!
//! The design's station loads are the *annular* element thrusts of the
//! converged BEM/lifting-line solve (newtons over the whole disk annulus).
//! Dividing by the panel width and the blade count gives the load per unit
//! span on one blade, `q(r)`, acting in the z direction (normal to the
//! rotor disk, the direction thrust actually tugs the blade).  The bending
//! moment follows by integrating that load from the tip inwards:
//!
//! ```text
//! M(r) = ∫_r^R q(s) (s - r) ds .
//! ```
//!
//! The blade bends as an Euler–Bernoulli beam with a rectangular-section
//! approximation — chord `c`, thickness `t`, second moment of area
//! `I = c t³ / 12` about the flapping axis.  The chord enters the
//! stiffness (`I ∝ c t³`), so a wider blade is stiffer and needs less
//! thickness; the twist enters through the design's chord distribution
//! (the geometric chord cap and the lifting-line chord optimisation both
//! depend on the local twist), i.e. through `c(r)` itself:
//!
//! ```text
//! w''(r) = M(r) / (E I(r)),    w(r_hub) = 0,  w'(r_hub) = 0 .
//! ```
//!
//! The thickness is laid out for a constant-curvature bend: with `w'' = κ`
//! and the tip deflection closed at the allowed value `δ` over the beam
//! length `L = R - r_hub` (`κ = 2δ/L²`), the section thickness satisfying
//! both is closed form:
//!
//! ```text
//! t(r) = (6 M(r) L² / (E c(r) δ))^(1/3) .
//! ```
//!
//! A uniform-curvature bend has no curvature concentration — the blade
//! bends as a smooth arc — and the cubic keeps the sizing stable against
//! the load (a blade with double the thrust needs only 2^(1/3) more
//! thickness).  The result is floored at a minimum fraction of the local
//! chord (the tip, where `M → 0`, would otherwise taper to a knife edge;
//! the default floor is the thinnest ARA-D table, 6%).
//!
//! The hub thickness (`hub_depth`) is deliberately not part of this: it
//! describes the hub mounting, and the blade's airfoil thickness is a
//! separate structural question answered by the load.  The law is returned
//! as a radius → *thickness/chord* curve — the foils carry `t/c`, and the
//! chord is the design's own output, so the sized absolute thickness
//! scales with the final chord exactly.

use crate::pchip::Pchip;

/// The sized thickness law: station radii (hub → tip, m) and the
/// thickness-to-chord fraction at each, plus the beam-model bookkeeping
/// (the predicted and allowed tip deflections it was closed on, m).
pub struct ThicknessLaw {
    pub rr: Vec<f64>,
    pub t_over_c: Vec<f64>,
    /// Predicted tip deflection of the sized blade (m), from the
    /// trapezoid double-integration of `M/(E I)`.
    pub tip_deflection: f64,
    /// Allowed tip deflection (m): `deflection_fraction * R`.
    pub deflection_limit: f64,
}

/// Size the blade thickness from the thrust load.
///
/// * `rr` — station radii, hub → tip, strictly increasing (m).
/// * `chords` — the converged design's chord at each station (m).
/// * `thrust` — the *annular* element thrust at each station (N, all
///   blades; it is divided by `n_blades` to get one blade's load).
/// * `n_blades` — blade count.
/// * `modulus` — blade material Young's modulus `E` (Pa).
/// * `deflection_fraction` — allowed tip deflection as a fraction of the
///   prop radius (the beam length is still `R - r_hub`).
/// * `thickness_floor` — minimum thickness as a fraction of the local
///   chord.
///
/// Returns `None` when there is nothing to size (fewer than two stations,
/// ragged input, a non-positive modulus/deflection, or zero/negative total
/// thrust — a design that produced no load cannot be sized).
pub fn size_mechanical_thickness(
    rr: &[f64],
    chords: &[f64],
    thrust: &[f64],
    n_blades: usize,
    modulus: f64,
    deflection_fraction: f64,
    thickness_floor: f64,
) -> Option<ThicknessLaw> {
    let n = rr.len();
    if n < 2
        || rr.len() != chords.len()
        || rr.len() != thrust.len()
        || n_blades == 0
        || !modulus.is_finite()
        || modulus <= 0.0
        || !deflection_fraction.is_finite()
        || deflection_fraction <= 0.0
        || !thickness_floor.is_finite()
        || !(0.0..1.0).contains(&thickness_floor)
    {
        return None;
    }
    // Station panel widths (central differences; endpoints one-sided).
    let mut dr = vec![0.0; n];
    dr[0] = rr[1] - rr[0];
    for i in 1..n - 1 {
        dr[i] = 0.5 * (rr[i + 1] - rr[i - 1]);
    }
    dr[n - 1] = rr[n - 1] - rr[n - 2];
    if dr.iter().any(|d| !d.is_finite() || *d <= 0.0) {
        return None;
    }

    // Per-blade load per unit span, and the total annular thrust sanity
    // check (a zero-load design — e.g. a failed solve — cannot be sized).
    let mut q = vec![0.0; n];
    let mut total = 0.0;
    for i in 0..n {
        q[i] = thrust[i] / dr[i] / n_blades as f64;
        total += thrust[i] * dr[i];
    }
    if !total.is_finite() || total <= 0.0 {
        return None;
    }

    // Bending moment at each station: trapezoid integration of the load
    // from the tip inboards — `M(r_i) = ∫ q(s)(s - r_i) ds`.
    let mut m = vec![0.0; n];
    for i in 0..n {
        let mut acc = 0.0;
        for j in i..n - 1 {
            let h = rr[j + 1] - rr[j];
            acc += 0.5 * (q[j] * (rr[j] - rr[i]) + q[j + 1] * (rr[j + 1] - rr[i])) * h;
        }
        m[i] = acc;
    }

    // Constant-curvature closed form, floored at a minimum fraction of the
    // local chord (the tip moment vanishes, so without the floor the tip
    // would taper to a knife edge).
    let l = rr[n - 1] - rr[0];
    let delta = deflection_fraction * rr[n - 1];
    let mut t = vec![0.0; n];
    for i in 0..n {
        let c = chords[i].abs().max(1.0e-9);
        let floor = thickness_floor * c;
        let sized = if m[i] > 0.0 {
            (6.0 * m[i] * l * l / (modulus * c * delta)).cbrt()
        } else {
            0.0
        };
        t[i] = sized.max(floor);
    }

    // The achieved tip deflection (trapezoid double integration of the
    // curvature), with the floor applied — with the floor the blade is
    // stiffer than the closed form everywhere it binds, so the achieved
    // deflection is at most the allowed one.
    let mut wpp = vec![0.0; n];
    for i in 0..n {
        let i_sec = chords[i].abs().max(1.0e-9) * t[i].powi(3) / 12.0;
        wpp[i] = m[i] / (modulus * i_sec.max(1.0e-30));
    }
    let mut wp = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 1..n {
        let h = rr[i] - rr[i - 1];
        wp[i] = wp[i - 1] + 0.5 * (wpp[i] + wpp[i - 1]) * h;
        w[i] = w[i - 1] + 0.5 * (wp[i] + wp[i - 1]) * h;
    }

    Some(ThicknessLaw {
        rr: rr.to_vec(),
        t_over_c: t
            .iter()
            .zip(chords.iter())
            .map(|(ti, ci)| (ti / ci.abs().max(1.0e-9)).clamp(0.0, 1.0))
            .collect(),
        tip_deflection: w[n - 1],
        deflection_limit: delta,
    })
}

/// Build the interpolant the design loop evaluates: the thickness/chord
/// fractions of a sized law as a function of radius.
pub fn law_interpolant(law: &ThicknessLaw) -> Pchip {
    Pchip::new(&law.rr, &law.t_over_c)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small uniform load over a tapered chord, mirroring a converged
    /// design's station data (hub → tip).
    fn web_default_like() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            vec![
                0.0060, 0.0129, 0.0198, 0.0267, 0.0336, 0.0404, 0.0473, 0.0542, 0.0611, 0.0680,
            ],
            vec![
                0.0062, 0.0073, 0.0081, 0.0083, 0.0077, 0.0065, 0.0050, 0.0037, 0.0031, 0.0030,
            ],
            vec![
                0.0058, 0.0344, 0.0784, 0.1261, 0.1670, 0.1909, 0.1901, 0.1683, 0.1443, 0.1142,
            ],
        )
    }

    #[test]
    fn closed_form_closes_the_deflection() {
        // Uniform chord and uniform load, no floor binding: the sizes land
        // on a constant curvature, so the trapezoid double integral
        // reproduces the allowed tip deflection (within the tip-panel
        // discretization: the free end has zero moment/curvature, so the
        // last panel loses ~h²/L² of the bend).
        let n = 20;
        let rr: Vec<f64> = (0..n).map(|i| 0.06 + 0.094 * i as f64 / (n - 1) as f64).collect();
        let chords = vec![0.02; n];
        let thrust = vec![0.2; n]; // ~2N total annular
        let law = size_mechanical_thickness(&rr, &chords, &thrust, 3, 3.0e9, 0.05, 0.0).unwrap();
        let rel = (law.tip_deflection - law.deflection_limit).abs() / law.deflection_limit;
        assert!(rel < 2.0e-3, "tip deflection {} vs limit {}", law.tip_deflection, law.deflection_limit);
    }

    #[test]
    fn floor_pins_the_tip() {
        let (rr, c, d) = web_default_like();
        let law = size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.05, 0.06).unwrap();
        // Every station respects the floor, and the tip (zero moment)
        // sits exactly on it.
        for &tc in law.t_over_c.iter() {
            assert!(tc >= 0.06 - 1.0e-12, "t/c {} below floor", tc);
            assert!(tc <= 1.0 + 1.0e-12, "t/c {} above 1", tc);
        }
        assert!((law.t_over_c[law.t_over_c.len() - 1] - 0.06).abs() < 1.0e-12);
        // The predicted deflection stays inside the allowed limit.
        assert!(law.tip_deflection <= law.deflection_limit * (1.0 + 1.0e-9));
        // The root — where the moment peaks — is the thickest section.
        let t_abs: Vec<f64> = law
            .t_over_c
            .iter()
            .zip(c.iter())
            .map(|(tc, ci)| tc * ci)
            .collect();
        assert!(t_abs[0] >= t_abs[t_abs.len() - 1]);
    }

    #[test]
    fn stiffness_scaling_is_cubic() {
        let (rr, c, d) = web_default_like();
        let soft = size_mechanical_thickness(&rr, &c, &d, 3, 2.0e9, 0.05, 0.0).unwrap();
        let stiff = size_mechanical_thickness(&rr, &c, &d, 3, 4.0e9, 0.05, 0.0).unwrap();
        // t ∝ E^(-1/3): doubling E shrinks every loaded section by 2^(-1/3).
        // (The tip station has no moment — the closed form sizes it to
        // zero, a knife edge absorbed by the real law's floor — so compare
        // only inboard of it.)
        let n = soft.t_over_c.len() - 1;
        for (a, b) in soft.t_over_c[..n].iter().zip(stiff.t_over_c[..n].iter()) {
            assert!(
                (a / b - 2.0f64.powf(1.0 / 3.0)).abs() < 1.0e-9,
                "{} vs {}",
                a,
                b
            );
        }
        // t ∝ δ^(-1/3) at fixed E: a looser deflection limit thins the blade.
        let loose = size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.10, 0.0).unwrap();
        let base = size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        for (a, b) in base.t_over_c[..n].iter().zip(loose.t_over_c[..n].iter()) {
            assert!(
                (a / b - 2.0f64.powf(1.0 / 3.0)).abs() < 1.0e-9,
                "{} vs {}",
                a,
                b
            );
        }
    }

    #[test]
    fn load_is_per_blade() {
        // Halving the blade count doubles each blade's load: the blade is
        // 2^(1/3) thicker.
        let (rr, c, d) = web_default_like();
        let three = size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let six = size_mechanical_thickness(&rr, &c, &d, 6, 3.0e9, 0.05, 0.0).unwrap();
        let n = three.t_over_c.len() - 1; // skip the zero-moment tip
        for (a, b) in three.t_over_c[..n].iter().zip(six.t_over_c[..n].iter()) {
            assert!((a / b - 2.0f64.powf(1.0 / 3.0)).abs() < 1.0e-9, "{} vs {}", a, b);
        }
    }

    #[test]
    fn chord_widening_thins_the_section() {
        // I ∝ c t³: doubling the chord at every station shrinks the sized
        // absolute thickness by 2^(-1/3), and (normalised by the doubled
        // chord) the t/c fraction by 2^(-4/3).
        let (rr, _, d) = web_default_like();
        let c1: Vec<f64> = web_default_like().1;
        let c2: Vec<f64> = c1.iter().map(|c| 2.0 * c).collect();
        let narrow = size_mechanical_thickness(&rr, &c1, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let wide = size_mechanical_thickness(&rr, &c2, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n = narrow.t_over_c.len() - 1; // skip the zero-moment tip
        for (a, b) in narrow.t_over_c[..n].iter().zip(wide.t_over_c[..n].iter()) {
            assert!(
                (a / b - 2.0f64.powf(4.0 / 3.0)).abs() < 1.0e-9,
                "{} vs {}",
                a,
                b
            );
        }
        // The absolute thickness still thins with the wider chord.
        for (i, (a, b)) in narrow.t_over_c[..n]
            .iter()
            .zip(wide.t_over_c[..n].iter())
            .enumerate()
        {
            assert!(
                (a * c1[i] / (b * c2[i]) - 2.0f64.powf(1.0 / 3.0)).abs() < 1.0e-9,
                "{} vs {}",
                a * c1[i],
                b * c2[i]
            );
        }
    }

    #[test]
    fn nothing_to_size_returns_none() {
        let (rr, c, _) = web_default_like();
        let d = vec![0.0; rr.len()]; // zero load
        assert!(size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.05, 0.06).is_none());
        assert!(size_mechanical_thickness(&rr[..1], &c[..1], &[1.0], 3, 3.0e9, 0.05, 0.06).is_none());
        assert!(size_mechanical_thickness(&rr, &c, &web_default_like().2, 0, 3.0e9, 0.05, 0.06).is_none());
        assert!(size_mechanical_thickness(&rr, &c, &web_default_like().2, 3, -3.0e9, 0.05, 0.06).is_none());
        assert!(size_mechanical_thickness(&rr, &c, &web_default_like().2, 3, 3.0e9, 0.0, 0.06).is_none());
        assert!(size_mechanical_thickness(&rr, &c, &web_default_like().2, 3, 3.0e9, 0.05, 1.5).is_none());
        // Ragged input.
        assert!(size_mechanical_thickness(&rr, &c[..5], &web_default_like().2, 3, 3.0e9, 0.05, 0.06).is_none());
    }

    #[test]
    fn interp_reproduces_the_law() {
        let (rr, c, d) = web_default_like();
        let law = size_mechanical_thickness(&rr, &c, &d, 3, 3.0e9, 0.05, 0.06).unwrap();
        let interp = law_interpolant(&law);
        for (r, tc) in law.rr.iter().zip(law.t_over_c.iter()) {
            assert!((interp.eval(*r) - tc).abs() < 1.0e-12, "at {}", r);
        }
    }
}
