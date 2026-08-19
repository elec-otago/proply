// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Foil geometry, ported from `proply/foil.py`.
//!
//! `Foil` is the base class (a flat plate); `Naca4` generates the NACA
//! 4-series shape with cosine spacing and an optional trailing-edge gap.

/// A foil: chord (m), thickness (fraction of chord) and trailing-edge gap
/// (fraction of chord).
#[derive(Debug, Clone)]
pub struct Foil {
    pub chord: f64,
    pub thickness: f64,
    pub trailing_edge: f64,
}

impl Foil {
    pub fn new(chord: f64, thickness: f64) -> Self {
        Self {
            chord,
            thickness,
            trailing_edge: 0.0,
        }
    }

    pub fn modify_chord(&mut self, c: f64) {
        self.chord = c;
    }

    /// Set the trailing-edge gap from a length in metres (stored normalized).
    pub fn set_trailing_edge(&mut self, te: f64) {
        self.trailing_edge = te / self.chord;
    }

    /// Reynolds number at velocity `v`.
    pub fn reynolds(&self, v: f64) -> f64 {
        1.225 * v * self.chord / 15.11e-6
    }

    /// Mach number at velocity `v` (speed of sound 330 m/s, as in the Python).
    pub fn mach(&self, v: f64) -> f64 {
        v / 330.0
    }

    /// Dynamic-pressure term: 0.5 rho v^2 chord.
    pub fn polar_aux(&self, v: f64) -> f64 {
        0.5 * 1.225 * v * v * self.chord
    }

    pub fn lift_per_unit_length(&self, v: f64, cl: f64) -> f64 {
        self.polar_aux(v) * cl
    }

    pub fn drag_per_unit_length(&self, v: f64, cd: f64) -> f64 {
        self.polar_aux(v) * cd
    }

    /// A unique hash for this foil (used as the polar cache key).
    ///
    /// The polar (from rust-foil, which normalizes to unit chord) depends only
    /// on the *shape* and the Reynolds number — not on the chord magnitude.
    /// `get_shape_points` scales coordinates by `chord`, so the hash is
    /// normalised by `chord` to identify the shape alone.  This lets a chord
    /// scale change reuse the same cached polar for the same Reynolds bucket
    /// (and makes the cache deterministic regardless of the order a design
    /// loop visits chord scales).
    pub fn hash(&self) -> String {
        let (_xl, yl, _xu, yu) = self.get_shape_points(10);
        let s: f64 = (yu.iter().skip(1).sum::<f64>() + yl.iter().skip(1).sum::<f64>())
            / self.chord.max(1.0e-12);
        format!("{}", s)
    }

    /// Shape points scaled by the chord: returns (xl, yl, xu, yu), each of
    /// length `n`.  Lower and upper surfaces both run LE -> TE.
    pub fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        // Flat plate: y = +/- thickness * chord.
        let x: Vec<f64> = (0..n).map(|i| self.chord * i as f64 / (n - 1) as f64).collect();
        let y = self.thickness * self.chord;
        let xl = x.clone();
        let yl = vec![-y; n];
        let xu = x;
        let yu = vec![y; n];
        (xl, yl, xu, yu)
    }

    /// The two profiles (lower, upper) as 2-D point lists, rotated by
    /// `rotation_angle` around the point x0 = 0.67 * chord-span.
    pub fn get_points(
        &self,
        n: usize,
        rotation_angle: f64,
    ) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
        let (xl, yl, xu, yu) = self.get_shape_points(n);
        let x0 = 0.67 * (max(&xu) - min(&xu));
        let y0 = 0.0;
        let l = rotate(&xl, &yl, x0, y0, rotation_angle);
        let u = rotate(&xu, &yu, x0, y0, rotation_angle);
        (l, u)
    }

    /// Lowest/highest extents of the rotated foil at angle `theta`.
    pub fn get_bounding_box(&self, theta: f64) -> (f64, f64, f64, f64) {
        let (l, u) = self.get_points(50, theta);
        let mut xs = Vec::with_capacity(100);
        let mut ys = Vec::with_capacity(100);
        for p in l.iter().chain(u.iter()) {
            xs.push(p[0]);
            ys.push(p[1]);
        }
        (min(&xs), max(&xs), min(&ys), max(&ys))
    }

    /// Largest chord that still fits the rotated foil in the box
    /// (x_limit, y_limit).
    pub fn get_max_chord(&self, x_limit: f64, y_limit: f64, theta: f64) -> f64 {
        let (x0, x1, y0, y1) = self.get_bounding_box(theta);
        let dy = y1 - y0;
        let dx = x1 - x0;
        let y_scale = y_limit / dy;
        let x_scale = x_limit / dx;
        self.chord * x_scale.min(y_scale)
    }

    /// Rotate a profile by `theta` around (x0, y0).
    pub fn rotate_points(points: &[[f64; 2]], x0: f64, y0: f64, theta: f64) -> Vec<[f64; 2]> {
        points
            .iter()
            .map(|p| {
                let x = p[0] - x0;
                let y = p[1] - y0;
                [
                    x * theta.cos() + y * theta.sin(),
                    -x * theta.sin() + y * theta.cos(),
                ]
            })
            .collect()
    }
}

impl std::fmt::Display for Foil {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ch={}, a={}", self.chord, self.thickness)
    }
}

/// Common interface over the foil classes, mirroring the Python duck typing
/// (the simulator and blade element only ever touch these methods).
pub trait FoilLike {
    fn chord(&self) -> f64;
    fn thickness(&self) -> f64;
    fn modify_chord(&mut self, c: f64);
    fn set_trailing_edge(&mut self, te: f64);
    fn hash(&self) -> String;
    /// Shape points scaled by the chord: (xl, yl, xu, yu), lower and upper,
    /// each running LE -> TE.
    fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);
    fn get_points(&self, n: usize, rotation_angle: f64) -> (Vec<[f64; 2]>, Vec<[f64; 2]>);
    fn get_bounding_box(&self, theta: f64) -> (f64, f64, f64, f64);
    fn get_max_chord(&self, x_limit: f64, y_limit: f64, theta: f64) -> f64;
    fn reynolds(&self, v: f64) -> f64;
    fn mach(&self, v: f64) -> f64;
}

impl FoilLike for Foil {
    fn chord(&self) -> f64 {
        self.chord
    }
    fn thickness(&self) -> f64 {
        self.thickness
    }
    fn modify_chord(&mut self, c: f64) {
        Foil::modify_chord(self, c);
    }
    fn set_trailing_edge(&mut self, te: f64) {
        Foil::set_trailing_edge(self, te);
    }
    fn hash(&self) -> String {
        Foil::hash(self)
    }
    fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        Foil::get_shape_points(self, n)
    }
    fn get_points(&self, n: usize, rotation_angle: f64) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
        Foil::get_points(self, n, rotation_angle)
    }
    fn get_bounding_box(&self, theta: f64) -> (f64, f64, f64, f64) {
        Foil::get_bounding_box(self, theta)
    }
    fn get_max_chord(&self, x_limit: f64, y_limit: f64, theta: f64) -> f64 {
        Foil::get_max_chord(self, x_limit, y_limit, theta)
    }
    fn reynolds(&self, v: f64) -> f64 {
        Foil::reynolds(self, v)
    }
    fn mach(&self, v: f64) -> f64 {
        Foil::mach(self, v)
    }
}

impl FoilLike for Naca4 {
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
        Naca4::hash(self)
    }
    fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        Naca4::get_shape_points(self, n)
    }
    fn get_points(&self, n: usize, rotation_angle: f64) -> (Vec<[f64; 2]>, Vec<[f64; 2]>) {
        // NACA4 uses the base Foil's point transforms on its own shape.
        let (xl, yl, xu, yu) = Naca4::get_shape_points(self, n);
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

/// NACA 4-series foil: `m` max camber (fraction of chord), `p` position of
/// max camber (fraction of chord from the LE).
#[derive(Debug, Clone)]
pub struct Naca4 {
    pub base: Foil,
    pub m: f64,
    pub p: f64,
}

impl Naca4 {
    pub fn new(chord: f64, thickness: f64, m: f64, p: f64) -> Self {
        Self {
            base: Foil::new(chord, thickness),
            m,
            p,
        }
    }

    pub fn hash(&self) -> String {
        format!(
            "{:5.2},{:5.2},{:5.2}, {:5.2}",
            self.m, self.p, self.base.thickness, self.base.trailing_edge
        )
    }

    /// Shape points scaled by chord; NACA report 460 with cosine spacing.
    /// Returns (xl, yl, xu, yu), lower and upper, both LE -> TE.
    pub fn get_shape_points(&self, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let n5 = n * 5;
        let t = self.base.thickness;
        let p = self.p;
        let m = self.m;

        let beta: Vec<f64> = (0..n5).map(|i| std::f64::consts::PI * i as f64 / (n5 - 1) as f64).collect();
        let x: Vec<f64> = beta.iter().map(|b| (1.0 - b.cos()) / 2.0).collect();
        let y_offset: Vec<f64> = (0..n5)
            .map(|i| self.base.trailing_edge / 2.0 * i as f64 / (n5 - 1) as f64)
            .collect();

        let yt: Vec<f64> = x
            .iter()
            .zip(y_offset.iter())
            .map(|(xi, yo)| {
                5.0 * t
                    * (0.2969 * xi.sqrt()
                        - 0.1260 * xi
                        - 0.3516 * xi * xi
                        + 0.2843 * xi.powi(3)
                        - 0.1036 * xi.powi(4))
                    + yo
            })
            .collect();

        let yc: Vec<f64> = x
            .iter()
            .map(|xi| {
                if *xi > p {
                    (m / ((1.0 - p) * (1.0 - p))) * (1.0 - 2.0 * p + 2.0 * p * xi - xi * xi)
                } else {
                    (m / (p * p)) * (2.0 * p * xi - xi * xi)
                }
            })
            .collect();

        let dyc: Vec<f64> = x
            .iter()
            .map(|xi| {
                if *xi > p {
                    (2.0 * m * (p - xi)) / ((p - 1.0) * (p - 1.0))
                } else {
                    m * (2.0 * p - 2.0 * xi) / (p * p)
                }
            })
            .collect();

        let theta: Vec<f64> = dyc.iter().map(|d| d.atan()).collect();

        let xu: Vec<f64> = x
            .iter()
            .zip(yt.iter().zip(theta.iter()))
            .map(|(xi, (yti, th))| xi - yti * th.sin())
            .collect();
        let yu: Vec<f64> = yc
            .iter()
            .zip(yt.iter().zip(theta.iter()))
            .map(|(yci, (yti, th))| yci + yti * th.cos())
            .collect();
        let xl: Vec<f64> = x
            .iter()
            .zip(yt.iter().zip(theta.iter()))
            .map(|(xi, (yti, th))| xi + yti * th.sin())
            .collect();
        let yl: Vec<f64> = yc
            .iter()
            .zip(yt.iter().zip(theta.iter()))
            .map(|(yci, (yti, th))| yci - yti * th.cos())
            .collect();

        let c = self.base.chord;
        (
            xl.iter().step_by(5).map(|v| v * c).collect(),
            yl.iter().step_by(5).map(|v| v * c).collect(),
            xu.iter().step_by(5).map(|v| v * c).collect(),
            yu.iter().step_by(5).map(|v| v * c).collect(),
        )
    }
}

impl std::fmt::Display for Naca4 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ch={}, te={:4.3}, NACA{:02}{:02}{:2}",
            self.base.chord,
            self.base.trailing_edge,
            (self.m * 100.0) as i32,
            (self.p * 10.0) as i32,
            (self.base.thickness * 100.0) as i32
        )
    }
}

fn min(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
}

fn max(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Rotate (x, y) by `theta` around (x0, y0) — the Python `Foil.rotate`.
fn rotate(x: &[f64], y: &[f64], x0: f64, y0: f64, theta: f64) -> Vec<[f64; 2]> {
    x.iter()
        .zip(y.iter())
        .map(|(xi, yi)| {
            let xr = xi - x0;
            let yr = yi - y0;
            [
                xr * theta.cos() + yr * theta.sin(),
                -xr * theta.sin() + yr * theta.cos(),
            ]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn naca4_symmetry() {
        // A symmetric NACA 0012: upper and lower are mirror images.
        let f = Naca4::new(1.0, 0.12, 0.0, 0.4);
        let (xl, yl, xu, yu) = f.get_shape_points(42);
        assert_eq!(xl.len(), 42);
        for i in 0..42 {
            assert!((xu[i] - xl[i]).abs() < 1e-12, "x mismatch at {}", i);
            assert!((yu[i] + yl[i]).abs() < 1e-12, "y mismatch at {}", i);
        }
        // LE at x=0.  The decimated profile (n*5 points, every 5th kept) ends
        // at the sample before the TE: beta = pi*205/209 -> x = 0.9990965.
        assert!(xu[0].abs() < 1e-12);
        assert!((xu[41] - 0.9990965).abs() < 1e-6, "TE sample x = {}", xu[41]);
    }

    #[test]
    fn hash_is_chord_independent() {
        // The polar cache key must identify the (unit-chord-normalised) shape
        // only, so that scaling the chord reuses the same polar for the same
        // Reynolds bucket.
        let a = Naca4::new(0.010, 0.12, 0.0, 0.4);
        let b = Naca4::new(0.020, 0.12, 0.0, 0.4);
        assert_eq!(a.hash(), b.hash(), "hash must not depend on chord");
        let mut c = Naca4::new(1.0, 0.12, 0.0, 0.4);
        c.modify_chord(0.005);
        assert_eq!(Naca4::new(0.005, 0.12, 0.0, 0.4).hash(), c.hash());
        // A different thickness (shape) must change the hash.
        assert_ne!(Naca4::new(0.01, 0.15, 0.0, 0.4).hash(), a.hash());
    }

    #[test]
    fn naca4_max_thickness() {
        // NACA 0012 has max HALF-thickness 0.06 (full 0.12) near x=0.3.
        let f = Naca4::new(1.0, 0.12, 0.0, 0.4);
        let (_xl, _yl, xu, yu) = f.get_shape_points(42);
        let max_y = max(&yu);
        assert!((max_y - 0.06).abs() < 0.002, "max half thickness {}", max_y);
        let idx = argmax(&yu);
        let x_at_max = xu[idx];
        assert!((x_at_max - 0.3).abs() < 0.05, "x at max {}", x_at_max);
    }

    #[test]
    fn cambered_foil_lifts() {
        // A cambered foil (m=0.06, p=0.4) has max camber near x=p.
        let f = Naca4::new(1.0, 0.15, 0.06, 0.4);
        let (_xl, yl, _xu, yu) = f.get_shape_points(42);
        let mean: Vec<f64> = (0..42).map(|i| (yu[i] + yl[i]) / 2.0).collect();
        let idx = argmax(&mean);
        let camber = mean[idx];
        assert!((camber - 0.06).abs() < 0.01, "max camber {}", camber);
    }

    fn argmax(v: &[f64]) -> usize {
        let mut bi = 0;
        let mut bv = f64::NEG_INFINITY;
        for (i, x) in v.iter().enumerate() {
            if *x > bv {
                bv = *x;
                bi = i;
            }
        }
        bi
    }

    #[test]
    fn trailing_edge_gap() {
        // With a trailing-edge setting, upper/lower TE are separated by
        // 2*yt*cos(theta) at the last decimated sample (beta = pi*205/209):
        // yt = 0.75*P(0.9991) + te/2*(205/209) = 0.049205, cos(theta)=0.98062.
        let mut f = Naca4::new(0.1, 0.15, 0.06, 0.4);
        f.base.set_trailing_edge(0.01);
        let (_xl, _yl, _xu, yu) = f.get_shape_points(42);
        let (_xl2, yl2, _xu2, _yu2) = f.get_shape_points(42);
        let te_gap = yu[41] - yl2[41];
        assert!((te_gap - 0.0096509).abs() < 1e-5, "te gap {}", te_gap);
    }

    #[test]
    fn hash_format() {
        let f = Naca4::new(0.1, 0.15, 0.06, 0.4);
        f.hash();
        // hash is stable across instances
        let g = Naca4::new(0.2, 0.15, 0.06, 0.4);
        assert_eq!(f.hash(), g.hash());
    }
}
