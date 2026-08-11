//! Airfoil polars via rust-foil, with the flat-plate fallback — a port of
//! `proply/foil_simulator.py` (the `XfoilSimulatedFoil` class).
//!
//! A polar is a 41-point alpha sweep from -20° to 20° at a Reynolds number
//! rounded onto a log-spaced grid; the (cl, cd) data is fitted with a
//! degree-9 polynomial in alpha (radians).  Results are cached per foil per
//! Reynolds number, both in memory and on disk.
//!
//! Note: the Python code passes a Mach number into the cache key but never
//! applies it to the XFOIL run (only `Re` and `max_iter` are set), so the
//! port keys the cache on (foil hash, Reynolds) only — the computed polar
//! is identical.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use crate::cache::{cache_key, PolarStore, StoredPolar};
use crate::foil::FoilLike;
use crate::polyfit::{polyfit, polyval};

const DEG2RAD: f64 = std::f64::consts::PI / 180.0;

/// A simulated foil: returns CL and CD for a velocity and angle of attack,
/// simulating polars with rust-foil on demand.
///
/// The foil is shared with the owning blade element (`Rc<RefCell>`), exactly
/// like the Python `self.fs = FoilSimulator(self.foil)` — so chord changes
/// made through `set_chord` are seen by the simulator (they feed the
/// Reynolds number).
pub struct FoilSimulator<F: FoilLike> {
    foil: Rc<RefCell<F>>,
    store: Arc<Mutex<PolarStore>>,
    /// (hash, reynolds) -> fitted (cl, cd) polynomial coefficients.
    poly_cache: RefCell<HashMap<(String, u64), (Vec<f64>, Vec<f64>)>>,
    /// When set, CL/CD come from the analytic flat-plate model (2 pi alpha,
    /// 1.28 sin alpha) instead of simulated polars — used for testing the
    /// design chain against the Python `PlateSimulatedFoil`.
    plate_mode: bool,
}

impl<F: FoilLike> FoilSimulator<F> {
    pub fn new(foil: Rc<RefCell<F>>, store: Arc<Mutex<PolarStore>>) -> Self {
        Self {
            foil,
            store,
            poly_cache: RefCell::new(HashMap::new()),
            plate_mode: false,
        }
    }

    /// The shared foil's current chord (used by the optimizer as the
    /// initial chord guess, mirroring `foil_simulator.foil.chord`).
    pub fn chord(&self) -> f64 {
        self.foil.borrow().chord()
    }

    /// Switch to the analytic flat-plate polar model.
    pub fn set_plate_mode(&mut self, on: bool) {
        self.plate_mode = on;
    }

    /// Reynolds number rounded onto np.geomspace(30000, 2e6, 20).
    fn get_reynolds(&self, v: f64) -> f64 {
        let re = self.foil.borrow().reynolds(v);
        let mut best = 30000.0;
        let mut best_d = (best - re).abs();
        let geom_ratio: f64 = 2.0e6 / 30000.0;
        for i in 1..20 {
            let cand = 30000.0 * geom_ratio.powf(i as f64 / 19.0);
            let d = (cand - re).abs();
            if d < best_d {
                best_d = d;
                best = cand;
            }
        }
        if best < 30000.0 {
            30000.0
        } else {
            best
        }
    }

    /// Mach number rounded to the nearest 0.05 (np.round(Ma*2, 1) / 2).
    fn get_mach(&self, v: f64) -> f64 {
        let ma = self.foil.borrow().mach(v);
        ((ma * 2.0 * 10.0).round() / 10.0) / 2.0
    }

    pub fn get_cl(&self, v: f64, alpha: f64) -> f64 {
        if self.plate_mode {
            return 2.0 * std::f64::consts::PI * alpha;
        }
        let ma = self.foil.borrow().mach(v);
        if ma > 0.97 || alpha.abs() > 30.0 * DEG2RAD || self.foil.borrow().reynolds(v) < 30000.0 {
            return 2.0 * std::f64::consts::PI * alpha;
        }
        let (cl_poly, _) = self.get_polars(v);
        polyval(&cl_poly, alpha)
    }

    pub fn get_cd(&self, v: f64, alpha: f64) -> f64 {
        if self.plate_mode {
            return 1.28 * alpha.sin();
        }
        let ma = self.foil.borrow().mach(v);
        if ma > 0.97 || alpha.abs() > 30.0 * DEG2RAD || self.foil.borrow().reynolds(v) < 30000.0 {
            return 1.28 * alpha.sin();
        }
        let (_, cd_poly) = self.get_polars(v);
        polyval(&cd_poly, alpha)
    }

    /// The fitted (cl, cd) polynomials for the velocity `v`.
    pub fn get_polars(&self, v: f64) -> (Vec<f64>, Vec<f64>) {
        let reynolds = self.get_reynolds(v);
        let mach = self.get_mach(v);
        let key = (self.foil.borrow().hash(), reynolds.to_bits());

        if let Some(p) = self.poly_cache.borrow().get(&key) {
            return p.clone();
        }

        let ck = cache_key(&key.0, reynolds, mach);
        let cached = {
            let store = self.store.lock().unwrap();
            store.get(&ck).cloned()
        };
        let (cl_poly, cd_poly) = match cached {
            Some(p) if p.alpha.len() > 20 => fit_polar(&p),
            _ => {
                self.xfoil_simulate_polars(reynolds, mach);
                let store = self.store.lock().unwrap();
                match store.get(&ck) {
                    Some(p) if p.alpha.len() > 20 => fit_polar(p),
                    // Degenerate case: too few converged points.  The Python
                    // code recurses forever here (its fallback never writes
                    // to the DB); we return the flat-plate model instead.
                    _ => flat_plate_polys(),
                }
            }
        };

        self.poly_cache
            .borrow_mut()
            .insert(key, (cl_poly.clone(), cd_poly.clone()));
        (cl_poly, cd_poly)
    }

    /// Simulate the polar with rust-foil and store it in the shared cache.
    /// Returns the stored polar (or None for the degenerate case).
    fn xfoil_simulate_polars(&self, reynolds: f64, mach: f64) {
        let (kept_x, kept_y) = {
            let f = self.foil.borrow();
            prepare_airfoil_coords(&*f, 42)
        };

        let mut xf = rust_foil::XFoil::new();
        xf.set_show_output(false);
        xf.set_airfoil(&kept_x, &kept_y);
        xf.set_reynolds(reynolds);
        xf.set_max_iter(80);
        xf.set_n_crit(9.0);

        let polar = xf.aseq(-20.0, 20.0, 80);

        // Prune non-converged points (the Python prunes NaN alpha values).
        let mut alfa = Vec::new();
        let mut cl = Vec::new();
        let mut cd = Vec::new();
        for (a_deg, clv, cdv, _cm, _cp, conv) in polar {
            if conv {
                alfa.push(a_deg * DEG2RAD);
                cl.push(clv);
                cd.push(cdv);
            }
        }

        if alfa.len() < 5 {
            // Foil didn't simulate; nothing is stored (see get_polars).
            return;
        }

        let key = cache_key(&self.foil.borrow().hash(), reynolds, mach);
        let mut store = self.store.lock().unwrap();
        store.insert(
            key,
            StoredPolar {
                alpha: alfa,
                cl,
                cd,
            },
        );
    }
}

/// Degree-9 least-squares fit of cl and cd over alpha (radians).
fn fit_polar(p: &StoredPolar) -> (Vec<f64>, Vec<f64>) {
    (
        polyfit(&p.alpha, &p.cl, 9),
        polyfit(&p.alpha, &p.cd, 9),
    )
}

/// The flat-plate fallback model, fitted like the Python `Foil didn't
/// simulate` branch (degree 4 over a degree grid).
fn flat_plate_polys() -> (Vec<f64>, Vec<f64>) {
    let alpha: Vec<f64> = (-20..=20).map(|a| a as f64).collect();
    let cl: Vec<f64> = alpha.iter().map(|a| 2.0 * std::f64::consts::PI * a).collect();
    let cd: Vec<f64> = alpha.iter().map(|a| 1.28 * (a * DEG2RAD).sin()).collect();
    (polyfit(&alpha, &cl, 4), polyfit(&alpha, &cd, 4))
}

/// Build the airfoil coordinate pair fed to rust-foil for a polar sweep,
/// exactly as the Python's `xfoil_simulate_polars` does — except that the
/// trailing edge is closed.
///
/// rust-foil normalizes the input to unit chord itself (its `set_airfoil`
/// default, matching canonical XFOIL), so the chord-scaled coordinates can
/// be passed through as-is.  The TE close is a deliberate deviation: the
/// proply TE gap (~0.25 mm) has negligible aerodynamic effect but
/// destabilizes the low-Reynolds boundary-layer solve (rust-foil, like the
/// XFOIL family, struggles to converge an open-TE boundary layer below
/// Re ~ 2e5; the converged values there are non-physical).  The design
/// geometry keeps the gap.
fn prepare_airfoil_coords(foil: &impl FoilLike, n_points: usize) -> (Vec<f64>, Vec<f64>) {
    let (xl, yl, xu, yu) = foil.get_shape_points(n_points);

    // xcoords = concat(xl[::-1], xu), ycoords likewise: from the TE,
    // around the LE, back to the TE.
    let mut xc: Vec<f64> = xl.iter().rev().copied().collect();
    let mut yc: Vec<f64> = yl.iter().rev().copied().collect();
    for i in 0..n_points {
        xc.push(xu[i]);
        yc.push(yu[i]);
    }

    // Chop off overhang: keep points with x <= xcoords[0].
    let x0 = xc[0];
    let mut kept_x = Vec::new();
    let mut kept_y = Vec::new();
    for (x, y) in xc.iter().zip(yc.iter()) {
        if *x <= x0 {
            kept_x.push(*x);
            kept_y.push(*y);
        }
    }

    // Close the trailing edge by replacing the two TE endpoints with their
    // mean (see the function docs above).
    let n = kept_x.len();
    let mean_x = 0.5 * (kept_x[0] + kept_x[n - 1]);
    let mean_y = 0.5 * (kept_y[0] + kept_y[n - 1]);
    kept_x[0] = mean_x;
    kept_y[0] = mean_y;
    kept_x[n - 1] = mean_x;
    kept_y[n - 1] = mean_y;

    (kept_x, kept_y)
}

/// Convenience: build a shared polar store, loaded from `path`.
pub fn shared_store(path: &str) -> Arc<Mutex<PolarStore>> {
    Arc::new(Mutex::new(PolarStore::load(path)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::foil::Naca4;

    fn test_store() -> Arc<Mutex<PolarStore>> {
        let dir = std::env::temp_dir().join("proply_rs_sim_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cache.json");
        let _ = std::fs::remove_file(&path);
        Arc::new(Mutex::new(PolarStore::load(path.to_str().unwrap())))
    }

    #[test]
    fn reynolds_grid_rounding() {
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f)), test_store());
        // v = 10 m/s, chord 0.05: Re = 1.225*10*0.05/15.11e-6 = 40536
        let re = fs.get_reynolds(10.0);
        // Nearest geomspace(30000, 2e6, 20) point: index 1 = 30000*(2e6/3e4)^(1/19)
        let ratio: f64 = 2.0e6 / 30000.0;
        let expected = 30000.0 * ratio.powf(1.0 / 19.0);
        assert!((re - expected).abs() < 1e-6, "re = {}", re);
    }

    #[test]
    fn mach_rounding() {
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f)), test_store());
        // v=80 -> Ma = 0.2424 -> rounded to 0.25
        assert!((fs.get_mach(80.0) - 0.25).abs() < 1e-9);
        assert!((fs.get_mach(100.0) - 0.30).abs() < 1e-9);
    }

    #[test]
    fn flat_plate_fallback_outside_envelope() {
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f)), test_store());
        // |alpha| > 30 deg -> flat plate
        let alpha = 35.0 * DEG2RAD;
        let cl = fs.get_cl(10.0, alpha);
        let cd = fs.get_cd(10.0, alpha);
        assert!((cl - 2.0 * std::f64::consts::PI * alpha).abs() < 1e-12);
        assert!((cd - 1.28 * alpha.sin()).abs() < 1e-12);
    }

    #[test]
    fn polar_fit_matches_stored_data() {
        // A smooth synthetic polar is inserted under the exact cache key the
        // simulator computes for a chosen velocity; get_polars must find it
        // and fit degree-9 polynomials that reproduce the data.
        let v = 22.0; // Re = 1.225*22*0.05/15.11e-6 ~ 89179
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());

        // Compute the same key the simulator does.
        let re_raw = f.reynolds(v);
        let geom_ratio: f64 = 2.0e6 / 30000.0;
        let mut re_grid = 30000.0;
        let mut best_d = (re_grid - re_raw).abs();
        for i in 1..20 {
            let cand = 30000.0 * geom_ratio.powf(i as f64 / 19.0);
            let d = (cand - re_raw).abs();
            if d < best_d {
                best_d = d;
                re_grid = cand;
            }
        }
        let key = cache_key(&f.hash(), re_grid, fs.get_mach(v));

        // Smooth degree-5 synthetic polar (in radians).
        let mut alpha = Vec::new();
        let mut cl = Vec::new();
        let mut cd = Vec::new();
        for i in 0..81 {
            let a = (-20.0 + 0.5 * i as f64) * DEG2RAD;
            alpha.push(a);
            cl.push(2.0 * std::f64::consts::PI * a + 3.0 * a * a - 20.0 * a.powi(3));
            cd.push(0.005 + 0.5 * a * a);
        }
        store.lock().unwrap().insert(
            key,
            StoredPolar {
                alpha,
                cl,
                cd,
            },
        );

        let (cl_p, cd_p) = fs.get_polars(v);
        // A degree-9 fit of the degree-5 data reproduces it to ~1e-10.
        let a = 5.0 * DEG2RAD;
        let cl_ref = 2.0 * std::f64::consts::PI * a + 3.0 * a * a - 20.0 * a.powi(3);
        let cl_est = polyval(&cl_p, a);
        assert!((cl_est - cl_ref).abs() < 1e-6, "cl_est {} cl_ref {}", cl_est, cl_ref);
        let cd_ref = 0.005 + 0.5 * a * a;
        let cd_est = polyval(&cd_p, a);
        assert!((cd_est - cd_ref).abs() < 1e-6, "cd_est {} cd_ref {}", cd_est, cd_ref);
    }

    #[test]
    fn airfoil_coords_close_the_trailing_edge() {
        // Regression test for the polar-preparation fix: rust-foil needs a
        // closed-TE airfoil or the low-Reynolds viscous solution breaks.
        // (Normalization is rust-foil's own default now; the coordinates
        // stay chord-scaled here.)
        let mut f = Naca4::new(0.02, 0.10, 0.0, 0.4);
        f.base.set_trailing_edge(0.25 / 1000.0);
        let (x, y) = prepare_airfoil_coords(&f, 42);
        let n = x.len();
        assert_eq!(n, y.len());
        // Closed TE: the first and last points coincide.
        assert!((x[0] - x[n - 1]).abs() < 1e-12, "TE not closed");
        assert!((y[0] - y[n - 1]).abs() < 1e-12, "TE not closed");
        // Coordinates are still chord-scaled (max x ~ chord, not 1).
        let mut mx = f64::NEG_INFINITY;
        for v in &x {
            mx = mx.max(*v);
        }
        assert!((mx - 0.02).abs() < 1e-3, "max x {} (chord-scaled)", mx);
        // 84 points (42 per side), and the LE is a single point.
        assert_eq!(n, 84);
    }
}
