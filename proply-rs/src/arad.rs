// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! ARA-D propeller airfoil family, ported from `proply/foil_ARA.py`
//! (`ARADFoil`, the class the legacy design loop used via `ARADProp`).
//!
//! The family is table-driven: four ARA-D sections (thickness 6, 10, 13
//! and 20% of chord) as Selig coordinate files under `src/arad/`, embedded
//! with [`include_str!`].  Each surface of each table is smoothed by a
//! degree-12 least-squares polynomial fit (as in the Python), and a section
//! of any thickness t/c is evaluated at 60 uniform chordwise stations by
//! interpolating across the thickness nodes:
//!
//! * t/c in [0.00, 0.04] (7 nodes): the 6% table scaled by t/0.06 — linear
//!   thinning down to a flat plate,
//! * t/c = 0.06, 0.10, 0.13, 0.20: the tables themselves,
//! * t/c in [0.25, 1.00] (9 nodes): the 20% table scaled by t/0.20 —
//!   linear thickening for the root sections.
//!
//! Deviations from the legacy Python (small, deliberate): the blend across
//! thickness is a PCHIP interpolation over those nodes instead of scipy's
//! `CloughTocher2DInterpolator` — identical at the nodes and throughout the
//! linear ramp ranges, C1-smooth between them; the smoothing polynomials
//! are evaluated directly at the output stations instead of through a
//! 90-point cosine resample; and a t/c outside [0, 1] is clamped (the
//! Python produced NaN outside the interpolator's hull).  The family
//! carries its own camber in the tables, so — as in the legacy `ARADProp` —
//! the design loop's camber setting does not apply.

use std::sync::OnceLock;

use crate::foil::{max, min, rotate, Foil, FoilLike};
use crate::pchip::Pchip;
use crate::polyfit::{polyfit, polyval};

/// Thickness (fraction of chord) of the four base tables.
const BASE_T: [f64; 4] = [0.06, 0.10, 0.13, 0.20];

/// Chordwise stations where a section is evaluated: 60 uniform stations
/// over the unit chord (the Python `np.linspace(x0, x1, 60)`).
const N_STATIONS: usize = 60;

/// Parse a Selig-format table: one header line, then `x y` rows running
/// TE -> LE along the upper surface and LE -> TE along the lower (the
/// Python `Foil.load_selig`).  Returns (xl, yl, xu, yu), both LE -> TE.
fn load_selig(data: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for line in data.lines().skip(1) {
        let mut it = line.split_whitespace().filter_map(|t| t.parse::<f64>().ok());
        let (Some(x), Some(y)) = (it.next(), it.next()) else {
            continue;
        };
        xs.push(x);
        ys.push(y);
    }
    // The upper block ends where x stops decreasing; its last row is the
    // leading edge, shared with the lower block's start.
    let split = xs
        .windows(2)
        .position(|w| w[1] >= w[0])
        .map_or(xs.len(), |i| i + 1);
    let xu: Vec<f64> = xs[..split].iter().rev().copied().collect();
    let yu: Vec<f64> = ys[..split].iter().rev().copied().collect();
    (xs[split..].to_vec(), ys[split..].to_vec(), xu, yu)
}

/// The pre-blended family: per-side ordinates of the four base tables at
/// the output stations, over the 20 thickness nodes.  Built once and
/// shared by every [`Arad`] (the Python cached its interpolators in module
/// globals for the same reason).
struct Tables {
    t_nodes: Vec<f64>,
    x: Vec<f64>,
    /// `l_rows[i][j]`: lower-surface y at thickness node i, station j.
    l_rows: Vec<Vec<f64>>,
    u_rows: Vec<Vec<f64>>,
}

fn tables() -> &'static Tables {
    static TABLES: OnceLock<Tables> = OnceLock::new();
    TABLES.get_or_init(|| {
        const FILES: [&str; 4] = [
            include_str!("arad/ara_d_6.dat"),
            include_str!("arad/ara_d_10.dat"),
            include_str!("arad/ara_d_13.dat"),
            include_str!("arad/ara_d_20.dat"),
        ];
        let x: Vec<f64> = (0..N_STATIONS)
            .map(|i| i as f64 / (N_STATIONS - 1) as f64)
            .collect();
        // Degree-12 least-squares smoothing of each surface, evaluated at
        // the output stations (ARADFoil.polyfit + poly1d).
        let mut base_l: Vec<Vec<f64>> = Vec::with_capacity(4);
        let mut base_u: Vec<Vec<f64>> = Vec::with_capacity(4);
        for data in FILES {
            let (xl, yl, xu, yu) = load_selig(data);
            let cl = polyfit(&xl, &yl, 12);
            let cu = polyfit(&xu, &yu, 12);
            base_l.push(x.iter().map(|&xi| polyval(&cl, xi)).collect());
            base_u.push(x.iter().map(|&xi| polyval(&cu, xi)).collect());
        }
        let mut t_nodes: Vec<f64> = Vec::with_capacity(20);
        let mut l_rows: Vec<Vec<f64>> = Vec::with_capacity(20);
        let mut u_rows: Vec<Vec<f64>> = Vec::with_capacity(20);
        let mut add = |t: f64, base: usize| {
            let s = t / BASE_T[base];
            t_nodes.push(t);
            l_rows.push(base_l[base].iter().map(|v| v * s).collect());
            u_rows.push(base_u[base].iter().map(|v| v * s).collect());
        };
        for k in 0..7 {
            add(0.04 * k as f64 / 6.0, 0); // np.linspace(0.0, 0.04, 7)
        }
        for (i, &t) in BASE_T.iter().enumerate() {
            add(t, i); // the tables themselves (scale 1)
        }
        for k in 0..9 {
            add(0.25 + 0.75 * k as f64 / 8.0, 3); // np.linspace(0.25, 1.0, 9)
        }
        Tables {
            t_nodes,
            x,
            l_rows,
            u_rows,
        }
    })
}

/// An ARA-D foil: the family shape at `thickness` (fraction of chord),
/// scaled by the chord.  Construction evaluates the thickness blend once
/// and stores the 60-station section — like the Python stored
/// `self.xl/yl/xu/yu`.
#[derive(Debug, Clone)]
pub struct Arad {
    pub base: Foil,
    x: Vec<f64>,
    yl: Vec<f64>,
    yu: Vec<f64>,
    /// TE gap of the interpolated section (unit chord): the trailing-edge
    /// setting is opened/closed relative to this.
    init_te: f64,
}

impl Arad {
    pub fn new(chord: f64, thickness: f64) -> Self {
        let tb = tables();
        // t/c outside [0, 1] has no table support; clamp to the end shapes.
        let t = thickness.clamp(0.0, 1.0);
        let mut yl = Vec::with_capacity(tb.x.len());
        let mut yu = Vec::with_capacity(tb.x.len());
        for j in 0..tb.x.len() {
            let lc: Vec<f64> = tb.l_rows.iter().map(|r| r[j]).collect();
            let uc: Vec<f64> = tb.u_rows.iter().map(|r| r[j]).collect();
            yl.push(Pchip::new(&tb.t_nodes, &lc).eval(t));
            yu.push(Pchip::new(&tb.t_nodes, &uc).eval(t));
        }
        let init_te = yu[yu.len() - 1] - yl[yl.len() - 1];
        Self {
            base: Foil::new(chord, thickness),
            x: tb.x.clone(),
            yl,
            yu,
            init_te,
        }
    }

    /// Max camber (fraction of chord) of this section: the peak of the mean
    /// line.  The ARA-D tables are inherently cambered; the design loop's
    /// camber setting does not move it.
    pub fn camber(&self) -> f64 {
        self.yl
            .iter()
            .zip(self.yu.iter())
            .map(|(l, u)| 0.5 * (l + u))
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Hash of the (unit-chord) shape: family tag, thickness and normalized
    /// trailing edge, mirroring the Python `"ARAD_I %5.2f,%5.2f"`.
    /// Chord-independent, so a chord scale change reuses the same cached
    /// polar for the same Reynolds bucket; thickness is quantised at 0.01
    /// exactly as in the Python.
    pub fn hash(&self) -> String {
        format!(
            "ARAD {:5.2},{:5.2}",
            self.base.thickness, self.base.trailing_edge
        )
    }

    /// Shape points scaled by the chord: (xl, yl, xu, yu), lower and upper,
    /// each running LE -> TE at cosine-spaced stations (the same convention
    /// as `Naca4::get_shape_points`).  The trailing-edge setting is blended
    /// in as a linear ramp from the LE, so the gap reaches `trailing_edge`
    /// at the TE (`ARADFoil.get_shape_points`).
    pub fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let n5 = n * 5;
        let li = Pchip::new(&self.x, &self.yl);
        let ui = Pchip::new(&self.x, &self.yu);
        let c = self.base.chord;
        let mut xl = Vec::with_capacity(n);
        let mut yl = Vec::with_capacity(n);
        let mut xu = Vec::with_capacity(n);
        let mut yu = Vec::with_capacity(n);
        for i in 0..n5 {
            if i % 5 != 0 {
                continue;
            }
            let beta = std::f64::consts::PI * i as f64 / (n5 - 1) as f64;
            let x = (1.0 - beta.cos()) / 2.0;
            let offset =
                (self.base.trailing_edge - self.init_te) / 2.0 * i as f64 / (n5 - 1) as f64;
            let y_lower = li.eval(x) - offset;
            let y_upper = if i == 0 {
                y_lower // close the leading edge
            } else {
                ui.eval(x) + offset
            };
            xl.push(x * c);
            yl.push(y_lower * c);
            xu.push(x * c);
            yu.push(y_upper * c);
        }
        (xl, yl, xu, yu)
    }
}

impl std::fmt::Display for Arad {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ch={}, te={:4.3}, ARA-D t={:4.2}",
            self.base.chord, self.base.trailing_edge, self.base.thickness
        )
    }
}

impl FoilLike for Arad {
    fn chord(&self) -> f64 {
        self.base.chord
    }
    fn thickness(&self) -> f64 {
        self.base.thickness
    }
    fn modify_chord(&mut self, c: f64) {
        self.base.modify_chord(c);
    }
    fn set_trailing_edge(&mut self, te: f64) {
        self.base.set_trailing_edge(te);
    }
    fn hash(&self) -> String {
        Arad::hash(self)
    }
    fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        Arad::get_shape_points(self, n)
    }
    fn get_points(&self, n: usize, rotation_angle: f64) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
        let (xl, yl, xu, yu) = Arad::get_shape_points(self, n);
        let x0 = 0.67 * (max(&xu) - min(&xu));
        let y0 = 0.0;
        let l = rotate(&xl, &yl, x0, y0, rotation_angle);
        let u = rotate(&xu, &yu, x0, y0, rotation_angle);
        (l, u)
    }
    fn get_bounding_box(&self, theta: f64) -> (f64, f64, f64, f64) {
        let (l, u) = self.get_points(50, theta);
        let mut xs = Vec::with_capacity(100);
        let mut ys = Vec::with_capacity(100);
        for p in l.iter().chain(u.iter()) {
            xs.push(p[0]);
            ys.push(p[1]);
        }
        (min(&xs), max(&xs), min(&ys), max(&ys))
    }
    fn get_max_chord(&self, x_limit: f64, y_limit: f64, theta: f64) -> f64 {
        let (x0, x1, y0, y1) = self.get_bounding_box(theta);
        let dy = y1 - y0;
        let dx = x1 - x0;
        let y_scale = y_limit / dy;
        let x_scale = x_limit / dx;
        self.base.chord * x_scale.min(y_scale)
    }
    fn reynolds(&self, v: f64) -> f64 {
        self.base.reynolds(v)
    }
    fn mach(&self, v: f64) -> f64 {
        self.base.mach(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILES: [&str; 4] = [
        include_str!("arad/ara_d_6.dat"),
        include_str!("arad/ara_d_10.dat"),
        include_str!("arad/ara_d_13.dat"),
        include_str!("arad/ara_d_20.dat"),
    ];

    #[test]
    fn selig_tables_parse() {
        for data in FILES {
            let (xl, _yl, xu, yu) = load_selig(data);
            // 51 rows: 26 upper (including the shared LE) + 25 lower.
            assert_eq!(xu.len(), 26);
            assert_eq!(xl.len(), 25);
            // Both surfaces run LE -> TE over the unit chord.
            assert!(xu[0].abs() < 1e-12 && yu[0].abs() < 1e-12, "LE on upper");
            assert!((xl[0] - 0.003).abs() < 1e-12);
            assert!((xu[25] - 1.0).abs() < 1e-12 && (xl[24] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn base_tables_reproduced_at_nodes() {
        // At the four base thicknesses the blend returns the smoothed
        // table values: max upper-surface y and its position match the
        // .dat files (0.079022 @ 0.40, 0.088961 @ 0.35, 0.100753 @ 0.30,
        // 0.137500 @ 0.25), within the degree-12 fit smoothing.
        let anchors = [
            (0.06, 0.079022, 0.40),
            (0.10, 0.088961, 0.35),
            (0.13, 0.100753, 0.30),
            (0.20, 0.137500, 0.25),
        ];
        for (t, ymax, xat) in anchors {
            let f = Arad::new(1.0, t);
            let mut best = (f64::NEG_INFINITY, 0.0);
            for i in 0..f.x.len() {
                if f.yu[i] > best.0 {
                    best = (f.yu[i], f.x[i]);
                }
            }
            assert!(
                (best.0 - ymax).abs() < 2.0e-3,
                "t={} max y {} want {}",
                t,
                best.0,
                ymax
            );
            assert!(
                (best.1 - xat).abs() < 0.08,
                "t={} x at max {} want {}",
                t,
                best.1,
                xat
            );
        }
    }

    #[test]
    fn thickness_law_reaches_the_section() {
        for t in BASE_T {
            let f = Arad::new(1.0, t);
            let tmax = f
                .yu
                .iter()
                .zip(f.yl.iter())
                .map(|(u, l)| u - l)
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                ((tmax - t) / t).abs() < 0.05,
                "t={} measured {}",
                t,
                tmax
            );
        }
    }

    #[test]
    fn zero_thickness_is_a_flat_plate() {
        let f = Arad::new(0.01, 0.0);
        for (l, u) in f.yl.iter().zip(f.yu.iter()) {
            assert!(l.abs() < 1e-12 && u.abs() < 1e-12);
        }
        assert!(f.camber().abs() < 1e-12);
    }

    #[test]
    fn ramps_scale_linearly() {
        // PCHIP reproduces linear data exactly, so the ramp ranges are
        // exact rescalings of their base table.
        let b06 = Arad::new(1.0, 0.06);
        let low = Arad::new(1.0, 0.02); // node of the low ramp: 0.06/3
        let b20 = Arad::new(1.0, 0.20);
        let high = Arad::new(1.0, 0.50); // inside the high ramp: 0.20 * 2.5
        for i in 0..b06.x.len() {
            assert!((low.yu[i] - b06.yu[i] / 3.0).abs() < 1e-12);
            assert!((low.yl[i] - b06.yl[i] / 3.0).abs() < 1e-12);
            assert!((high.yu[i] - 2.5 * b20.yu[i]).abs() < 1e-12);
            assert!((high.yl[i] - 2.5 * b20.yl[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn thickness_is_clamped_to_the_tables() {
        let over = Arad::new(1.0, 1.5);
        let top = Arad::new(1.0, 1.0);
        let under = Arad::new(1.0, -0.3);
        for i in 0..over.x.len() {
            assert!((over.yu[i] - top.yu[i]).abs() < 1e-12);
            assert!((under.yu[i] - 0.0).abs() < 1e-12);
        }
    }

    #[test]
    fn shape_points_contract() {
        let mut f = Arad::new(1.0, 0.10);
        f.set_trailing_edge(0.001);
        let (xl, yl, xu, yu) = f.get_shape_points(42);
        assert_eq!(xl.len(), 42);
        // Cosine-spaced stations, shared by both surfaces, matching the
        // NACA4 decimation: LE at 0, last sample at x = 0.9990965.
        assert!(xu[0].abs() < 1e-12);
        assert!((xu[41] - 0.9990965).abs() < 1e-6, "TE sample x = {}", xu[41]);
        for i in 0..42 {
            assert!((xu[i] - xl[i]).abs() < 1e-12);
        }
        // Leading edge closed, surfaces on the right sides of zero.
        assert!((yu[0] - yl[0]).abs() < 1e-15);
        assert!(max(&yu) > 0.0 && min(&yl) < 0.0);
    }

    #[test]
    fn trailing_edge_gap_matches_setting() {
        let mut f = Arad::new(0.1, 0.10);
        f.set_trailing_edge(0.001);
        let (_xl, yl, _xu, yu) = f.get_shape_points(42);
        // The decimated last sample sits just before the TE (i = 205 of
        // 209), where the ramp has opened most of the way to the setting;
        // the remaining gap is the section's own init_te tail.
        let gap = yu[41] - yl[41];
        assert!(
            (gap - 0.001 * 205.0 / 209.0).abs() < 2.0e-4,
            "gap {}",
            gap
        );
    }

    #[test]
    fn inherent_camber() {
        // The tables are cambered sections: the 6% table's mean line peaks
        // near 5% of chord.
        let f = Arad::new(1.0, 0.06);
        assert!(
            (f.camber() - 0.05).abs() < 0.01,
            "camber {}",
            f.camber()
        );
    }

    #[test]
    fn hash_is_chord_independent() {
        let a = Arad::new(0.010, 0.12);
        let b = Arad::new(0.020, 0.12);
        assert_eq!(a.hash(), b.hash(), "hash must not depend on chord");
        assert_ne!(Arad::new(0.01, 0.15).hash(), a.hash());
        let mut c = Arad::new(0.01, 0.12);
        c.set_trailing_edge(0.0005);
        assert_ne!(c.hash(), a.hash());
    }
}
