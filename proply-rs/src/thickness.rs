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
//! The blade bends as an Euler–Bernoulli beam with the section's *real*
//! bending inertia about the rotor-disk chord axis, not the enclosing
//! rectangle.  The section is twisted by the local twist `θ(r)` about the
//! spanwise axis, so the z-bending inertia of the *rotated* section is
//! the sum of its two principal moments, each weighted by the twist (the
//! z-projections of the section's dimensions):
//!
//! ```text
//! I(r) = i_flat(r)·c(r)·t³(r)/12 · cos²θ  +  i_edge(r)·c³(r)·t(r)/12 · sin²θ ,
//! ```
//!
//! where `i_flat` and `i_edge` are the section's actual second moments
//! about the chord line and the thickness axis relative to the rectangle
//! ([`crate::foil::SectionShape`], computed from the foil's shape
//! points).  The chord-line moment contains the **camber**: a curved
//! (cambered) section is stiffer than a flat one because its mean line
//! carries area away from the bending axis (`i_flat ≈ 0.47` for a
//! symmetric NACA section, ≈ 0.62 at 2% camber, ≈ 1.3 at 6% camber);
//! even symmetric, a real section is roughly half as stiff as its
//! rectangle.  The flatwise term is the z-projection of the *thickness*
//! (`t·cos θ`); the chordwise term is the z-projection of the *chord*
//! (`c·sin θ`) — since the chord is much longer than the thickness, at
//! high twist the chord's projection dominates and the needed thickness
//! tends to the chord (the section's total z-extent
//! `c·sin θ + t·cos θ`):
//!
//! ```text
//! w''(r) = M(r) / (E I(r)),    w(r_hub) = 0,  w'(r_hub) = 0 .
//! ```
//!
//! The thickness is laid out for a constant-curvature bend: with `w'' = κ`
//! and the tip deflection closed at the allowed value `δ` over the beam
//! length `L = R - r_hub` (`κ = 2δ/L²`), the section thickness satisfies
//! the cubic
//!
//! ```text
//! i_flat·cos²θ·t³ + i_edge·c²·sin²θ·t = 6 M(r) L² / (E c(r) δ)   (=: K) ,
//! ```
//!
//! which reduces to the closed form `t = (K/(i_flat))^(1/3)` for an
//! untwisted section.  A uniform-curvature bend has no curvature
//! concentration — the blade bends as a smooth arc — and the cubic keeps
//! the sizing stable against the load (a blade with double the thrust
//! needs only 2^(1/3) more thickness).  The result is floored at a
//! minimum fraction of the local chord (the tip, where `M → 0`, would
//! otherwise taper to a knife edge; the default floor is the thinnest
//! ARA-D table, 6%).
//!
//! The hub thickness (`hub_depth`) is deliberately not part of this: it
//! describes the hub mounting, and the blade's airfoil thickness is a
//! separate structural question answered by the load.  The law is returned
//! as a radius → *thickness/chord* curve — the foils carry `t/c`, and the
//! chord is the design's own output, so the sized absolute thickness
//! scales with the final chord exactly.

use crate::foil::SectionShape;
use crate::pchip::Pchip;

/// The sized thickness law: station radii (hub → tip, m) and the
/// thickness-to-chord fraction at each, plus the beam-model bookkeeping
/// (the predicted and allowed tip deflections it was closed on, m).
pub struct ThicknessLaw {
    pub rr: Vec<f64>,
    pub t_over_c: Vec<f64>,
    /// Predicted tip deflection of the sized blade (m), from the
    /// trapezoid double-integration of `M/(E I)` (twist-aware inertia).
    pub tip_deflection: f64,
    /// Allowed tip deflection (m): `deflection_fraction * R`.
    pub deflection_limit: f64,
}

/// Solve the shape-and-twist-aware sizing cubic for the positive
/// thickness: `a·t³ + b·t − K = 0` with `a = i_flat·cos²θ` and
/// `b = i_edge·c²·sin²θ`.  The left side is strictly increasing in
/// `t ≥ 0`, so the root is unique; Newton from a bracketed start (the
/// b = 0 solution is an upper bound) converges monotonically.
fn solve_sized_cubic(a: f64, b: f64, k: f64) -> f64 {
    if a <= 1.0e-12 {
        // Nearly edge-on: the linear (chord-projection) term carries the
        // load — t ≈ K / (i_edge·c²·sin²θ).
        return k / b.max(1.0e-30);
    }
    if b <= 1.0e-30 {
        // Untwisted section: t = (K/(i_flat·cos²θ))^(1/3).
        return (k / a).cbrt();
    }
    let upper = (k / a.max(1.0e-30)).cbrt(); // root < upper (b·t > 0)
    let mut t = 0.5 * upper;
    for _ in 0..40 {
        let f = a * t * t * t + b * t - k;
        let fp = 3.0 * a * t * t + b;
        let dt = f / fp;
        t -= dt;
        if dt.abs() <= 1.0e-14 * t.max(upper) {
            break;
        }
    }
    t.max(0.0)
}

/// Size the blade thickness from the thrust load.
///
/// * `rr` — station radii, hub → tip, strictly increasing (m).
/// * `chords` — the converged design's chord at each station (m).
/// * `twist` — the section twist at each station (rad): the rotation of
///   the foil about the spanwise axis that enters the z-bending inertia
///   (the z-projections of the thickness and of the chord).
/// * `shape` — each station's bending-inertia shape factors relative to
///   the enclosing rectangle ([`crate::foil::SectionShape`]): the real
///   section is roughly half as stiff as its rectangle, and camber raises
///   the flatwise factor — a curved (cambered) section is stiffer than a
///   flat one and sizes thinner.
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
/// ragged input, a non-positive modulus/deflection, a non-positive shape
/// factor, or zero/negative total thrust — a design that produced no load
/// cannot be sized).
#[allow(clippy::too_many_arguments)]
pub fn size_mechanical_thickness(
    rr: &[f64],
    chords: &[f64],
    twist: &[f64],
    shape: &[SectionShape],
    thrust: &[f64],
    n_blades: usize,
    modulus: f64,
    deflection_fraction: f64,
    thickness_floor: f64,
) -> Option<ThicknessLaw> {
    let n = rr.len();
    if n < 2
        || rr.len() != chords.len()
        || rr.len() != twist.len()
        || rr.len() != shape.len()
        || rr.len() != thrust.len()
        || n_blades == 0
        || !modulus.is_finite()
        || modulus <= 0.0
        || !deflection_fraction.is_finite()
        || deflection_fraction <= 0.0
        || !thickness_floor.is_finite()
        || !(0.0..1.0).contains(&thickness_floor)
        || twist.iter().any(|x| !x.is_finite())
        || shape
            .iter()
            .any(|s| !s.i_flat.is_finite() || !s.i_edge.is_finite() || s.i_flat <= 0.0 || s.i_edge <= 0.0)
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

    // Constant-curvature sizing with the shape-and-twist-aware inertia,
    // floored at a minimum fraction of the local chord (the tip moment
    // vanishes, so without the floor the tip would taper to a knife
    // edge).  The cubic `i_flat·cos²θ·t³ + i_edge·c²·sin²θ·t = K`
    // accounts for the real section (a cambered section's higher
    // `i_flat` sizes a thinner foil) and for the twist presenting the
    // chord to the z load.
    let l = rr[n - 1] - rr[0];
    let delta = deflection_fraction * rr[n - 1];
    let mut t = vec![0.0; n];
    for i in 0..n {
        let c = chords[i].abs().max(1.0e-9);
        let floor = thickness_floor * c;
        let sized = if m[i] > 0.0 {
            let k = 6.0 * m[i] * l * l / (modulus * c * delta);
            let (cos2, sin2) = (twist[i].cos().powi(2), twist[i].sin().powi(2));
            solve_sized_cubic(
                shape[i].i_flat * cos2,
                shape[i].i_edge * c * c * sin2,
                k,
            )
        } else {
            0.0
        };
        t[i] = sized.max(floor);
    }

    // The achieved tip deflection (trapezoid double integration of the
    // curvature), with the floor applied — with the floor the blade is
    // stiffer than the sizing everywhere it binds, so the achieved
    // deflection is at most the allowed one.
    let mut wpp = vec![0.0; n];
    for i in 0..n {
        let c = chords[i].abs().max(1.0e-9);
        let (cos2, sin2) = (twist[i].cos().powi(2), twist[i].sin().powi(2));
        let i_sec = shape[i].i_edge * c.powi(3) * t[i] / 12.0 * sin2
            + shape[i].i_flat * c * t[i].powi(3) / 12.0 * cos2;
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

    /// An untwisted blade: the sizing collapses to the closed form.
    fn zero_twist(n: usize) -> Vec<f64> {
        vec![0.0; n]
    }

    /// Rectangle-section shape factors: the model's historical section.
    fn rect(n: usize) -> Vec<SectionShape> {
        vec![SectionShape::default(); n]
    }

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
        // Uniform chord, uniform load, no twist and no floor binding: the
        // sizes land on a constant curvature, so the trapezoid double
        // integral reproduces the allowed tip deflection (within the
        // tip-panel discretization: the free end has zero moment/curvature,
        // so the last panel loses ~h²/L² of the bend).
        let n = 20;
        let rr: Vec<f64> = (0..n).map(|i| 0.06 + 0.094 * i as f64 / (n - 1) as f64).collect();
        let chords = vec![0.02; n];
        let thrust = vec![0.2; n]; // ~2N total annular
        let law = size_mechanical_thickness(
            &rr,
            &chords,
            &zero_twist(n),
            &rect(n),
            &thrust,
            3,
            3.0e9,
            0.05,
            0.0,
        )
        .unwrap();
        let rel = (law.tip_deflection - law.deflection_limit).abs() / law.deflection_limit;
        assert!(
            rel < 2.0e-3,
            "tip deflection {} vs limit {}",
            law.tip_deflection,
            law.deflection_limit
        );
    }

    #[test]
    fn floor_pins_the_tip() {
        let (rr, c, d) = web_default_like();
        let law =
            size_mechanical_thickness(&rr, &c, &zero_twist(rr.len()), &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.06)
                .unwrap();
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
        let tw = zero_twist(rr.len());
        let soft = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 2.0e9, 0.05, 0.0).unwrap();
        let stiff = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 4.0e9, 0.05, 0.0).unwrap();
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
        let loose = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.10, 0.0).unwrap();
        let base = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
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
        let tw = zero_twist(rr.len());
        let three = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let six = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 6, 3.0e9, 0.05, 0.0).unwrap();
        let n = three.t_over_c.len() - 1; // skip the zero-moment tip
        for (a, b) in three.t_over_c[..n].iter().zip(six.t_over_c[..n].iter()) {
            assert!((a / b - 2.0f64.powf(1.0 / 3.0)).abs() < 1.0e-9, "{} vs {}", a, b);
        }
    }

    #[test]
    fn chord_widening_thins_the_section() {
        // I ∝ c t³ (untwisted): doubling the chord at every station
        // shrinks the sized absolute thickness by 2^(-1/3), and
        // (normalised by the doubled chord) the t/c fraction by 2^(-4/3).
        let (rr, _, d) = web_default_like();
        let tw = zero_twist(rr.len());
        let c1: Vec<f64> = web_default_like().1;
        let c2: Vec<f64> = c1.iter().map(|c| 2.0 * c).collect();
        let narrow = size_mechanical_thickness(&rr, &c1, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let wide = size_mechanical_thickness(&rr, &c2, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
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
    fn twist_reduces_the_sized_thickness() {
        // The twisted section presents the z-projection of its chord to
        // the load (the TODO's "twist, as the chord can make the beam
        // stiffer"): at every loaded station the sizing with the local
        // twist is at most the untwisted sizing, and strictly less where
        // the blade is pitched.
        let (rr, c, d) = web_default_like();
        // A real pitch ramp: 31.5° at the hub -> 4.6° at the tip.
        let tw: Vec<f64> = (0..rr.len())
            .map(|i| 0.55 - (0.55 - 0.08) * i as f64 / (rr.len() - 1) as f64)
            .collect();
        let flat =
            size_mechanical_thickness(&rr, &c, &zero_twist(rr.len()), &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0)
                .unwrap();
        let twisted = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n = flat.t_over_c.len() - 1; // skip the zero-moment tip
        for i in 0..n {
            assert!(
                twisted.t_over_c[i] <= flat.t_over_c[i] + 1.0e-12,
                "station {}: {:.6} vs {:.6}",
                i,
                twisted.t_over_c[i],
                flat.t_over_c[i]
            );
        }
        // Inboard, where the twist is substantial, the reduction is real:
        // the wide root presents most of its chord to the z load.
        assert!(
            twisted.t_over_c[0] < 0.5 * flat.t_over_c[0],
            "root {:.4} vs {:.4}",
            twisted.t_over_c[0],
            flat.t_over_c[0]
        );
        // The law never sizes a sub-floor fraction: the interpolant stays
        // within the foil families' supported range.
        for &tc in twisted.t_over_c.iter() {
            assert!((0.0..=1.0).contains(&tc), "t/c {} out of range", tc);
        }
    }

    #[test]
    fn higher_twist_thins_monotonically() {
        let (rr, c, d) = web_default_like();
        let mut prev = vec![f64::INFINITY; rr.len()];
        for deg in [0.0_f64, 5.0, 15.0, 30.0, 60.0, 85.0] {
            let tw = vec![deg.to_radians(); rr.len()];
            let law = size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
            for (i, tc) in law.t_over_c.iter().enumerate().take(law.t_over_c.len() - 1) {
                assert!(
                    *tc <= prev[i] + 1.0e-12,
                    "{}° station {}: {} > {}",
                    deg,
                    i,
                    tc,
                    prev[i]
                );
                prev[i] = *tc;
            }
        }
    }

    #[test]
    fn high_twist_is_governed_by_the_chord_projection() {
        // At ~90° twist the section stands edge-on in the z load: the
        // sizing collapses to the linear term t ≈ K/(c²·sin²θ), so the
        // thickness follows the chord cubed — doubling the chord divides
        // the absolute thickness by eight (16 in t/c).
        let (rr, _, d) = web_default_like();
        let c1: Vec<f64> = web_default_like().1;
        let c2: Vec<f64> = c1.iter().map(|c| 2.0 * c).collect();
        let tw = vec![85.0_f64.to_radians(); rr.len()];
        let n1 = size_mechanical_thickness(&rr, &c1, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n2 = size_mechanical_thickness(&rr, &c2, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n = n1.t_over_c.len() - 1; // skip the zero-moment tip
        for i in 0..n {
            // The cos²θ term never vanishes exactly, so the ratio is close
            // to (not exactly) the pure-linear scalings: 16 (t/c) and 8
            // (absolute).
            let tc_ratio = n1.t_over_c[i] / n2.t_over_c[i];
            assert!((tc_ratio - 16.0).abs() < 5.0e-3, "station {}: {}", i, tc_ratio);
            let t_ratio = n1.t_over_c[i] * c1[i] / (n2.t_over_c[i] * c2[i]);
            assert!((t_ratio - 8.0).abs() < 5.0e-3, "station {} abs: {}", i, t_ratio);
        }
    }

    #[test]
    fn twist_cubic_reduces_to_the_closed_form_when_untwisted() {
        // θ = 0 (and a 360°-equivalent twist) must reproduce the untwisted
        // closed form exactly.
        let (rr, c, d) = web_default_like();
        let tw0 = zero_twist(rr.len());
        let base = size_mechanical_thickness(&rr, &c, &tw0, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let tw_full = vec![2.0 * std::f64::consts::PI; rr.len()];
        let full = size_mechanical_thickness(&rr, &c, &tw_full, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        for (a, b) in base.t_over_c.iter().zip(full.t_over_c.iter()) {
            assert!((a - b).abs() < 1.0e-12, "{} vs {}", a, b);
        }
    }

    #[test]
    fn real_section_factors_thicken_the_untwisted_blade() {
        // A real symmetric section is ≈ 0.47x as stiff as its enclosing
        // rectangle, so the untwisted sizing is t = (K/i_flat)^(1/3):
        // 29% thicker than the rectangle model at i_flat = 0.473.
        let (rr, c, d) = web_default_like();
        let tw = zero_twist(rr.len());
        let rect_shape = rect(rr.len());
        let real_shape: Vec<SectionShape> = (0..rr.len())
            .map(|_| SectionShape {
                i_flat: 0.473,
                i_edge: 0.454,
            })
            .collect();
        let rect_law =
            size_mechanical_thickness(&rr, &c, &tw, &rect_shape, &d, 3, 3.0e9, 0.05, 0.0)
                .unwrap();
        let real_law =
            size_mechanical_thickness(&rr, &c, &tw, &real_shape, &d, 3, 3.0e9, 0.05, 0.0)
                .unwrap();
        let n = rect_law.t_over_c.len() - 1; // skip the zero-moment tip
        for i in 0..n {
            assert!(
                (real_law.t_over_c[i] / rect_law.t_over_c[i] - 0.473f64.powf(-1.0 / 3.0)).abs()
                    < 1.0e-9,
                "station {}: {} vs {}",
                i,
                real_law.t_over_c[i],
                rect_law.t_over_c[i]
            );
        }
    }

    #[test]
    fn cambered_sections_size_thinner_than_flat() {
        // A cambered (curved) section is stiffer than a flat one: its
        // higher i_flat (0.622 vs 0.473 for a 2%-cambered vs symmetric
        // NACA at 12% t/c) sizes it t ∝ i_flat^(-1/3) thinner.
        let (rr, c, d) = web_default_like();
        let tw = zero_twist(rr.len());
        let flat: Vec<SectionShape> = (0..rr.len())
            .map(|_| SectionShape {
                i_flat: 0.473,
                i_edge: 0.454,
            })
            .collect();
        let curved: Vec<SectionShape> = (0..rr.len())
            .map(|_| SectionShape {
                i_flat: 0.622,
                i_edge: 0.454,
            })
            .collect();
        let flat_law =
            size_mechanical_thickness(&rr, &c, &tw, &flat, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let curved_law =
            size_mechanical_thickness(&rr, &c, &tw, &curved, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n = flat_law.t_over_c.len() - 1;
        for i in 0..n {
            let want = (0.473_f64 / 0.622).powf(1.0 / 3.0);
            assert!(
                (curved_law.t_over_c[i] / flat_law.t_over_c[i] - want).abs() < 1.0e-9,
                "station {}: {} vs {}",
                i,
                curved_law.t_over_c[i],
                flat_law.t_over_c[i]
            );
        }
    }

    #[test]
    fn edge_factor_scales_the_chord_projection() {
        // The edgewise shape factor only matters where the twist exposes
        // the chord to the load: at 85° the sizing is nearly the linear
        // term t ≈ K/(i_edge·c²·sin²θ), so i_edge 0.454 vs 1.0 almost
        // exactly doubles the sized thickness.
        let (rr, c, d) = web_default_like();
        let tw = vec![85.0_f64.to_radians(); rr.len()];
        let s1: Vec<SectionShape> = (0..rr.len())
            .map(|_| SectionShape {
                i_flat: 0.473,
                i_edge: 1.0,
            })
            .collect();
        let s2: Vec<SectionShape> = (0..rr.len())
            .map(|_| SectionShape {
                i_flat: 0.473,
                i_edge: 0.454,
            })
            .collect();
        let stiff =
            size_mechanical_thickness(&rr, &c, &tw, &s1, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let soft = size_mechanical_thickness(&rr, &c, &tw, &s2, &d, 3, 3.0e9, 0.05, 0.0).unwrap();
        let n = stiff.t_over_c.len() - 1;
        for i in 0..n {
            assert!(
                (soft.t_over_c[i] / stiff.t_over_c[i] - 1.0 / 0.454).abs() < 1.0e-3,
                "station {}: {} vs {}",
                i,
                soft.t_over_c[i],
                stiff.t_over_c[i]
            );
        }
    }

    #[test]
    fn nothing_to_size_returns_none() {
        let (rr, c, _) = web_default_like();
        let tw = zero_twist(rr.len());
        let d = vec![0.0; rr.len()]; // zero load
        assert!(size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.06).is_none());
        assert!(
            size_mechanical_thickness(&rr[..1], &c[..1], &tw[..1], &rect(1), &[1.0], 3, 3.0e9, 0.05, 0.06)
                .is_none()
        );
        assert!(
            size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &web_default_like().2, 0, 3.0e9, 0.05, 0.06)
                .is_none()
        );
        assert!(
            size_mechanical_thickness(
                &rr,
                &c,
                &tw,
                &rect(rr.len()),
                &web_default_like().2,
                3,
                -3.0e9,
                0.05,
                0.06
            )
            .is_none()
        );
        assert!(
            size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &web_default_like().2, 3, 3.0e9, 0.0, 0.06)
                .is_none()
        );
        assert!(
            size_mechanical_thickness(&rr, &c, &tw, &rect(rr.len()), &web_default_like().2, 3, 3.0e9, 0.05, 1.5)
                .is_none()
        );
        // Ragged input: a twist array of the wrong length.
        assert!(
            size_mechanical_thickness(&rr, &c[..5], &tw, &rect(rr.len()), &web_default_like().2, 3, 3.0e9, 0.05, 0.06)
                .is_none()
        );
        // Non-finite twist.
        let bad_tw = vec![f64::NAN; rr.len()];
        assert!(
            size_mechanical_thickness(&rr, &c, &bad_tw, &rect(rr.len()), &web_default_like().2, 3, 3.0e9, 0.05, 0.06)
                .is_none()
        );
    }

    #[test]
    fn interp_reproduces_the_law() {
        let (rr, c, d) = web_default_like();
        let law =
            size_mechanical_thickness(&rr, &c, &zero_twist(rr.len()), &rect(rr.len()), &d, 3, 3.0e9, 0.05, 0.06)
                .unwrap();
        let interp = law_interpolant(&law);
        for (r, tc) in law.rr.iter().zip(law.t_over_c.iter()) {
            assert!((interp.eval(*r) - tc).abs() < 1.0e-12, "at {}", r);
        }
    }
}
