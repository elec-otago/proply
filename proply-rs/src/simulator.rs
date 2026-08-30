// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Airfoil polars via rust-foil, with the flat-plate fallback — a port of
//! `proply/foil_simulator.py` (the `XfoilSimulatedFoil` class).
//!
//! A polar is a 41-point alpha sweep from -20° to 20° at a Reynolds number
//! on a log-spaced grid; the (cl, cd) data is fitted with a degree-9
//! polynomial in alpha (radians).  Results are cached per foil per Reynolds
//! number, both in memory and on disk.
//!
//! Polars are *simulated* at the grid points only, but *evaluated* by
//! log-Re interpolation between the two bracketing buckets' fits (blending
//! the fitted coefficients, which is exact for linear interpolation of the
//! polar functions).  Below the grid floor the same blend runs from the
//! first bucket down to the analytic flat-plate model, so cl/cd — and
//! everything derived from them (best-L/D angles, circulation) — are
//! continuous in Reynolds number instead of stepping at every bucket
//! boundary.
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

/// Reynolds grid the polars are simulated at: np.geomspace(30000, 2e6, 20).
const RE_MIN: f64 = 30000.0;
const RE_MAX: f64 = 2.0e6;
const N_RE: usize = 20;
/// Below the grid floor the polars blend from the first bucket down to the
/// analytic flat-plate model, reaching it at this Reynolds number (XFOIL is
/// not trustworthy there, but a continuous blend beats a hard switch —
/// low-Re hub stations would otherwise jump between two very different
/// cl/cd models across one radius).
const RE_FLAT_PLATE: f64 = 10000.0;
/// Physical floor on the profile drag: at these Reynolds numbers a smooth
/// airfoil's cd never drops below ~0.008 over the working alpha range, but
/// rust-foil occasionally "converges" with cd ~ 0 (a failed viscous solve
/// flagged as converged).  A degree-9 fit of such a sweep produces spurious
/// cl/cd peaks of 200+ that `best_alpha` happily chases, so those buckets
/// are treated as degenerate (flat-plate fallback), exactly like sweeps
/// with too few converged points.
const CD_FLOOR: f64 = 0.004;

/// Flat-plate fallback blend window for high attack angles (radians).
///
/// The simulated polars are fitted over the ±20° alpha sweep; beyond it the
/// degree-9 fit extrapolates uncontrollably and XFOIL's viscous solution is
/// unreliable post-stall, so the legacy code switched to the analytic
/// flat-plate model at |alpha| > 30°.  That hard switch made cl/cd
/// discontinuous (a fitted cl of ~1 just below 30° jumping to 2 pi alpha ~
/// 3.3 just above), and the design optimizers exploited the cliff as a free
/// lift source — the root stations of the browser_demo design settled on
/// the flat-plate branch while the outboard stations used the polars,
/// producing mixed, inconsistent designs.  Instead, blend smoothly from the
/// fitted polar to the flat-plate model over `[ALPHA_BLEND_LO, ALPHA_BLEND_HI]`
/// (cubic smoothstep: C1-continuous at both ends), pure polar below the low
/// end and pure flat plate above the high end.
const ALPHA_BLEND_LO: f64 = 18.0_f64.to_radians();
const ALPHA_BLEND_HI: f64 = 33.0_f64.to_radians();

/// Smoothstep (Hermite) blend weight from the fitted polar to the analytic
/// flat-plate model: 0 for |alpha| <= [`ALPHA_BLEND_LO`], 1 for |alpha| >=
/// [`ALPHA_BLEND_HI`], smooth (C1) in between.
fn stall_blend(alpha: f64) -> f64 {
    let a = alpha.abs();
    if a <= ALPHA_BLEND_LO {
        return 0.0;
    }
    if a >= ALPHA_BLEND_HI {
        return 1.0;
    }
    let t = (a - ALPHA_BLEND_LO) / (ALPHA_BLEND_HI - ALPHA_BLEND_LO);
    t * t * (3.0 - 2.0 * t)
}

/// Reject sweeps whose drag is unphysical (see [`CD_FLOOR`]): the minimum
/// cd over the |alpha| <= 15 deg working window must clear the floor.
fn polar_is_physical(p: &StoredPolar) -> bool {
    let win = 15.0 * DEG2RAD;
    let mut min_cd = f64::INFINITY;
    let mut n_win = 0;
    for (&a, &cd) in p.alpha.iter().zip(p.cd.iter()) {
        if a.abs() <= win {
            if cd.is_nan() {
                return false;
            }
            min_cd = min_cd.min(cd);
            n_win += 1;
        }
    }
    n_win > 0 && min_cd > CD_FLOOR
}

/// Per-key in-flight simulations (process-wide).  A worker pool and the
/// parallel design passes can miss the same cache key concurrently; the
/// gate makes the expensive sweep run once per key, with waiters blocking
/// on the condvar instead of duplicating it.  Lives outside the polar
/// store so a waiter never blocks while holding the store's data lock.
fn sim_keys() -> &'static std::sync::Mutex<std::collections::HashSet<String>> {
    static KEYS: std::sync::OnceLock<std::sync::Mutex<std::collections::HashSet<String>>> =
        std::sync::OnceLock::new();
    KEYS.get_or_init(|| std::sync::Mutex::new(std::collections::HashSet::new()))
}

fn sim_cv() -> &'static std::sync::Condvar {
    static CV: std::sync::OnceLock<std::sync::Condvar> = std::sync::OnceLock::new();
    CV.get_or_init(std::sync::Condvar::new)
}

/// Grid point `i` of the Reynolds grid (log-spaced, `RE_MIN..=RE_MAX`).
fn re_grid(i: usize) -> f64 {
    RE_MIN * (RE_MAX / RE_MIN).powf(i as f64 / (N_RE - 1) as f64)
}

/// Bracket `re >= RE_MIN` on the Reynolds grid: `(lo, hi, w)` with `lo`/`hi`
/// adjacent grid indices and `w` the log-Re blend weight toward `hi`
/// (0 at `re_grid(lo)`, 1 at `re_grid(hi)`).  Clamped at the top of the
/// grid — no extrapolation past `RE_MAX`.
fn re_bracket(re: f64) -> (usize, usize, f64) {
    let mut lo = 0;
    while lo + 2 < N_RE && re_grid(lo + 1) <= re {
        lo += 1;
    }
    let (g_lo, g_hi) = (re_grid(lo), re_grid(lo + 1));
    let w = ((re / g_lo).ln() / (g_hi / g_lo).ln()).clamp(0.0, 1.0);
    (lo, lo + 1, w)
}

/// Linear blend of two polynomial coefficient vectors (highest power first,
/// as `polyfit` returns), padding the shorter with zero high-order
/// coefficients.
fn blend_poly(lo: &[f64], hi: &[f64], w: f64) -> Vec<f64> {
    let n = lo.len().max(hi.len());
    let lift = |p: &[f64]| -> Vec<f64> {
        let mut v = vec![0.0; n - p.len()];
        v.extend_from_slice(p);
        v
    };
    let (a, b) = (lift(lo), lift(hi));
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (1.0 - w) * x + w * y)
        .collect()
}

/// A simulated foil: returns CL and CD for a velocity and angle of attack,
/// simulating polars with rust-foil on demand.
///
/// The foil is shared with the owning blade element (`Rc<RefCell>`), exactly
/// like the Python `self.fs = FoilSimulator(self.foil)` — so chord changes
/// made through `set_chord` are seen by the simulator (they feed the
/// Reynolds number).
/// (hash, reynolds) -> fitted (cl, cd) polynomial coefficients.
type PolarEq = (String, u64);
/// Keyed cache of polar fits, `PolarEq -> (cl, cd) coefficient vectors`.
type PolyCache = RefCell<HashMap<PolarEq, (Vec<f64>, Vec<f64>)>>;
pub struct FoilSimulator<F: FoilLike> {
    foil: Rc<RefCell<F>>,
    store: Arc<Mutex<PolarStore>>,
    poly_cache: PolyCache,
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
        if ma > 0.97 || self.foil.borrow().reynolds(v) < RE_FLAT_PLATE {
            return 2.0 * std::f64::consts::PI * alpha;
        }
        let flat = 2.0 * std::f64::consts::PI * alpha;
        let w = stall_blend(alpha);
        if w >= 1.0 {
            return flat;
        }
        let (cl_poly, _) = self.get_polars(v);
        (1.0 - w) * polyval(&cl_poly, alpha) + w * flat
    }

    pub fn get_cd(&self, v: f64, alpha: f64) -> f64 {
        if self.plate_mode {
            return 1.28 * alpha.sin();
        }
        let ma = self.foil.borrow().mach(v);
        if ma > 0.97 || self.foil.borrow().reynolds(v) < RE_FLAT_PLATE {
            return 1.28 * alpha.sin();
        }
        let flat = 1.28 * alpha.sin();
        let w = stall_blend(alpha);
        if w >= 1.0 {
            return flat;
        }
        let (_, cd_poly) = self.get_polars(v);
        (1.0 - w) * polyval(&cd_poly, alpha) + w * flat
    }

    /// The fitted (cl, cd) polynomials for the velocity `v`, interpolated
    /// across Reynolds buckets (see the module docs): the bracketing grid
    /// buckets' fits are blended in log-Re, and below the grid floor the
    /// blend runs from the first bucket down to the flat-plate model.
    pub fn get_polars(&self, v: f64) -> (Vec<f64>, Vec<f64>) {
        let re = self.foil.borrow().reynolds(v);
        let mach = self.get_mach(v);

        if re < RE_FLAT_PLATE {
            return flat_plate_polys();
        }
        if re < RE_MIN {
            let (fp_cl, fp_cd) = flat_plate_polys();
            let w = (re / RE_FLAT_PLATE).ln() / (RE_MIN / RE_FLAT_PLATE).ln();
            let (cl_hi, cd_hi) = self.bucket_polars(RE_MIN, mach);
            return (blend_poly(&fp_cl, &cl_hi, w), blend_poly(&fp_cd, &cd_hi, w));
        }

        let (lo, hi, w) = re_bracket(re);
        let (cl_lo, cd_lo) = self.bucket_polars(re_grid(lo), mach);
        let (cl_hi, cd_hi) = self.bucket_polars(re_grid(hi), mach);
        (blend_poly(&cl_lo, &cl_hi, w), blend_poly(&cd_lo, &cd_hi, w))
    }

    /// The fitted (cl, cd) polynomials of the single grid bucket at
    /// `reynolds` (a [`re_grid`] point): the memoized degree-9 fit of the
    /// cached polar, simulated with rust-foil on first use.  This is the
    /// bucket-level fetch the interpolation blends.
    fn bucket_polars(&self, reynolds: f64, mach: f64) -> (Vec<f64>, Vec<f64>) {
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
            Some(p) if p.alpha.len() > 20 && polar_is_physical(&p) => fit_polar(&p),
            _ => {
                // Claim the key (blocking while another thread simulates
                // it), then re-check the store once the claim is ours — the
                // previous holder may have stored it in the meantime.
                {
                    let mut keys = sim_keys().lock().unwrap();
                    while keys.contains(&ck) {
                        keys = sim_cv().wait(keys).unwrap();
                    }
                    keys.insert(ck.clone());
                }
                let raced = {
                    let store = self.store.lock().unwrap();
                    store.get(&ck).and_then(|p| {
                        if p.alpha.len() > 20 && polar_is_physical(p) {
                            Some(fit_polar(p))
                        } else {
                            None
                        }
                    })
                };
                let polys = match raced {
                    Some(polys) => polys,
                    None => {
                        self.xfoil_simulate_polars(reynolds, mach);
                        let store = self.store.lock().unwrap();
                        match store.get(&ck) {
                            Some(p) if p.alpha.len() > 20 && polar_is_physical(p) => fit_polar(p),
                            // Degenerate case: too few converged points, or an
                            // unphysical (near-zero drag) sweep.  The Python
                            // code recurses forever here (its fallback never
                            // writes to the DB); we return the flat-plate
                            // model instead.
                            _ => flat_plate_polys(),
                        }
                    }
                };
                let mut keys = sim_keys().lock().unwrap();
                keys.remove(&ck);
                sim_cv().notify_all();
                polys
            }
        };

        self.poly_cache
            .borrow_mut()
            .insert(key, (cl_poly.clone(), cd_poly.clone()));
        (cl_poly, cd_poly)
    }

    /// The `(Reynolds-grid bucket, Mach)` warm targets an evaluation at
    /// speed `v` will consult: the two grid buckets bracketing the foil's
    /// raw Reynolds number there (the low-Re blend zone needs only the
    /// first bucket; below the blend floor nothing is simulated).
    /// Split out so a worker pool can deduplicate its queue on these —
    /// adjacent stations share bracketing buckets, and per-station tasks
    /// would otherwise pile every worker onto the same first bucket.
    pub fn warm_plan(&self, v: f64) -> Vec<(f64, f64)> {
        if self.plate_mode {
            return Vec::new();
        }
        let re = self.foil.borrow().reynolds(v);
        let mach = self.get_mach(v);
        if re < RE_FLAT_PLATE {
            return Vec::new(); // analytic branch: nothing is simulated
        }
        if re < RE_MIN {
            return vec![(RE_MIN, mach)];
        }
        let (lo, hi, _) = re_bracket(re);
        vec![(re_grid(lo), mach), (re_grid(hi), mach)]
    }

    /// Simulate and cache one [`warm_plan`] target (a Reynolds-grid bucket
    /// at a Mach number).  Cheap when the bucket is already cached; this
    /// is the expensive first-touch path of [`get_polars`], split out so a
    /// worker pool can populate the shared store one distinct bucket per
    /// task.  No-op in plate mode (no polars are consulted).
    pub fn warm_bucket(&self, reynolds: f64, mach: f64) {
        if self.plate_mode {
            return;
        }
        self.bucket_polars(reynolds, mach);
    }

    /// Simulate and cache the polar buckets an evaluation at speed `v`
    /// would need ([`warm_plan`] + [`warm_bucket`]).
    pub fn warm_polars(&self, v: f64) {
        for (re, mach) in self.warm_plan(v) {
            self.warm_bucket(re, mach);
        }
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
    (polyfit(&p.alpha, &p.cl, 9), polyfit(&p.alpha, &p.cd, 9))
}

/// The flat-plate fallback model, fitted like the Python `Foil didn't
/// simulate` branch (degree 4 over a degree grid).
fn flat_plate_polys() -> (Vec<f64>, Vec<f64>) {
    let alpha: Vec<f64> = (-20..=20).map(|a| a as f64).collect();
    let cl: Vec<f64> = alpha
        .iter()
        .map(|a| 2.0 * std::f64::consts::PI * a)
        .collect();
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
    fn reynolds_bracket_on_log_grid() {
        // v = 10 m/s, chord 0.05: Re = 1.225*10*0.05/15.11e-6 = 40536,
        // which lies between grid points 1 and 2.
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let _fs = FoilSimulator::new(Rc::new(RefCell::new(f)), test_store());
        let re = 1.225 * 10.0 * 0.05 / 15.11e-6;
        let (lo, hi, w) = re_bracket(re);
        assert_eq!((lo, hi), (1, 2), "re = {}", re);
        assert!(
            w > 0.0 && w < 1.0,
            "w = {} for re = {} between {} and {}",
            w,
            re,
            re_grid(lo),
            re_grid(hi)
        );
        // Exactly on a grid point: weight zero, bracket starts there.
        let (lo, _hi, w) = re_bracket(re_grid(5));
        assert_eq!(lo, 5);
        assert!(w.abs() < 1.0e-12);
        // Above the grid: clamped to the top pair, no extrapolation.
        let (lo, hi, w) = re_bracket(RE_MAX * 10.0);
        assert_eq!((lo, hi), (N_RE - 2, N_RE - 1));
        assert!((w - 1.0).abs() < 1.0e-12);
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
    fn warm_plan_lists_distinct_bucket_targets() {
        // The warm targets are the bracketing Reynolds-grid buckets at the
        // rounded Mach — what a worker pool deduplicates its queue on.
        let chord = 0.05;
        let v = 22.0; // Re ~ 89179, between grid points 1 and 2
        let f = Naca4::new(chord, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());
        let (lo, hi, _) = re_bracket(f.reynolds(v));
        let plan = fs.warm_plan(v);
        assert_eq!(plan.len(), 2, "{:?}", plan);
        assert_eq!(plan[0], (re_grid(lo), fs.get_mach(v)));
        assert_eq!(plan[1], (re_grid(hi), fs.get_mach(v)));

        // The low-Re blend zone needs only the first bucket; below the
        // blend floor nothing is simulated at all.
        let v_for = |re: f64| re * 15.11e-6 / (1.225 * chord);
        assert_eq!(
            fs.warm_plan(v_for(RE_MIN * 0.9)),
            vec![(RE_MIN, fs.get_mach(v_for(RE_MIN * 0.9)))]
        );
        assert!(fs.warm_plan(v_for(RE_FLAT_PLATE * 0.9)).is_empty());
    }

    #[test]
    fn flat_plate_fallback_outside_envelope() {
        let f = Naca4::new(0.05, 0.12, 0.0, 0.4);
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f)), test_store());
        // |alpha| beyond the stall-blend window (>= 33 deg) is pure flat plate.
        let alpha = 35.0 * DEG2RAD;
        let cl = fs.get_cl(10.0, alpha);
        let cd = fs.get_cd(10.0, alpha);
        assert!((cl - 2.0 * std::f64::consts::PI * alpha).abs() < 1e-12);
        assert!((cd - 1.28 * alpha.sin()).abs() < 1e-12);
    }

    #[test]
    fn stall_blend_is_continuous_across_thirty_degrees() {
        // Regression (browser_demo plate-vs-full-polar discrepancy): the
        // legacy hard switch at |alpha| > 30° made cl/cd jump (a fitted cl of
        // ~1 just below 30° vs 2 pi alpha ~ 3.3 just above), and the design
        // optimizers exploited the cliff as a free lift source.  The blend
        // must be continuous — and smooth — across the old switch point, so
        // an optimizer cannot gain a free force jump by crossing it.
        let chord = 0.05;
        let v_for = |re: f64| re * 15.11e-6 / (1.225 * chord);
        let v = 22.0; // Re ~ 89179, between grid points 1 and 2
        let f = Naca4::new(chord, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());

        // Distinctive seeded polars (cl = 0.5 * 2 pi alpha, cd quadratic) on
        // the bracketing buckets: the blend mixes the polar with the
        // flat-plate model, so the result must move continuously through the
        // 30° boundary with no step between the two models.
        let re_raw = f.reynolds(v);
        let (lo, hi, _) = re_bracket(re_raw);
        seed_bucket(&f, &store, &fs, re_grid(lo), v_for(re_grid(lo)), 0.5);
        seed_bucket(&f, &store, &fs, re_grid(hi), v_for(re_grid(hi)), 0.5);

        let eps = 1.0e-9;
        let a30 = 30.0 * DEG2RAD;
        let cl_lo = fs.get_cl(v, a30 - eps);
        let cl_hi = fs.get_cl(v, a30 + eps);
        let cd_lo = fs.get_cd(v, a30 - eps);
        let cd_hi = fs.get_cd(v, a30 + eps);
        // No jump: the 30°-crossing step must be gone entirely (the old code
        // stepped by most of the difference between the polar and the plate).
        let scale = (2.0 * std::f64::consts::PI * a30).abs().max(1.0);
        assert!(
            (cl_lo - cl_hi).abs() < 1.0e-6 * scale,
            "cl jump at 30°: {} vs {}",
            cl_lo,
            cl_hi
        );
        assert!(
            (cd_lo - cd_hi).abs() < 1.0e-6,
            "cd jump at 30°: {} vs {}",
            cd_lo,
            cd_hi
        );
        // The blend lies between the weak polar and the flat plate inside the
        // window, and the flat plate alone governs beyond it.
        let a40 = 40.0 * DEG2RAD;
        assert!(
            (fs.get_cl(v, a40) - 2.0 * std::f64::consts::PI * a40).abs() < 1e-12,
            "pure flat plate above the window"
        );
        let a10 = 10.0 * DEG2RAD;
        assert!(
            (fs.get_cl(v, a10) - 0.5 * 2.0 * std::f64::consts::PI * a10).abs() < 1e-6,
            "pure polar below the window"
        );
    }

    /// Seed `k`-scaled synthetic polar data (cl = k * 2 pi alpha, exact under
    /// the degree-9 fit) under the store keys for the given grid Reynolds
    /// numbers, at the snapped Mach the simulator will query for `v_for(re)`.
    fn seed_bucket(
        f: &Naca4,
        store: &Arc<Mutex<PolarStore>>,
        fs: &FoilSimulator<Naca4>,
        g: f64,
        v_for_g: f64,
        k: f64,
    ) {
        let mut alpha = Vec::new();
        let mut cl = Vec::new();
        let mut cd = Vec::new();
        for i in 0..81 {
            let a = (-20.0 + 0.5 * i as f64) * DEG2RAD;
            alpha.push(a);
            cl.push(k * 2.0 * std::f64::consts::PI * a);
            cd.push(0.005 + 0.3 * a * a);
        }
        let key = cache_key(&f.hash(), g, fs.get_mach(v_for_g));
        store
            .lock()
            .unwrap()
            .insert(key, StoredPolar { alpha, cl, cd });
    }

    #[test]
    fn polar_fit_matches_stored_data() {
        // get_polars must find the stored bucket data and reproduce it with
        // the degree-9 fit.  The interpolation fetches both bracketing
        // buckets, so seed them with the same polar.
        let chord = 0.05;
        let v_for = |re: f64| re * 15.11e-6 / (1.225 * chord);
        let v = 22.0; // Re = 1.225*22*0.05/15.11e-6 ~ 89179
        let f = Naca4::new(chord, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());

        let re_raw = f.reynolds(v);
        let (lo, hi, _w) = re_bracket(re_raw);
        seed_bucket(&f, &store, &fs, re_grid(lo), v_for(re_grid(lo)), 1.0);
        seed_bucket(&f, &store, &fs, re_grid(hi), v_for(re_grid(hi)), 1.0);

        let (cl_p, cd_p) = fs.get_polars(v);
        // A degree-9 fit of linear/quadratic data reproduces it to ~1e-10.
        let a = 5.0 * DEG2RAD;
        let cl_ref = 2.0 * std::f64::consts::PI * a;
        let cl_est = polyval(&cl_p, a);
        assert!(
            (cl_est - cl_ref).abs() < 1e-6,
            "cl_est {} cl_ref {}",
            cl_est,
            cl_ref
        );
        let cd_ref = 0.005 + 0.3 * a * a;
        let cd_est = polyval(&cd_p, a);
        assert!(
            (cd_est - cd_ref).abs() < 1e-6,
            "cd_est {} cd_ref {}",
            cd_est,
            cd_ref
        );
    }

    #[test]
    fn polars_interpolate_across_reynolds() {
        // Adjacent buckets carry different polars; cl must move continuously
        // from one to the other in log-Re, with no step at the bucket
        // boundary (the snapped behaviour jumped by the full
        // bucket-to-bucket difference).
        let chord = 0.05;
        let v_for = |re: f64| re * 15.11e-6 / (1.225 * chord);
        let f = Naca4::new(chord, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());

        // cl_k(alpha) = k * 2 pi alpha; buckets 1..3 get k = 1, 3, 5.
        seed_bucket(&f, &store, &fs, re_grid(1), v_for(re_grid(1)), 1.0);
        seed_bucket(&f, &store, &fs, re_grid(2), v_for(re_grid(2)), 3.0);
        seed_bucket(&f, &store, &fs, re_grid(3), v_for(re_grid(3)), 5.0);

        let a = 5.0 * DEG2RAD;
        let tau = 2.0 * std::f64::consts::PI;
        let cl_at = |re: f64| fs.get_cl(v_for(re), a);

        // On grid points the bucket's own polar is returned.
        assert!((cl_at(re_grid(1)) - tau * a).abs() < 1e-9);
        assert!((cl_at(re_grid(2)) - 3.0 * tau * a).abs() < 1e-9);
        // Log-midpoint: exactly the mean of the bracketing polars.
        let re_mid = (re_grid(1) * re_grid(2)).sqrt();
        assert!((cl_at(re_mid) - 2.0 * tau * a).abs() < 1e-9);
        // Continuity across a bucket boundary (this is the regression: the
        // snapped code stepped by the full (3-1) * tau * a here).
        let eps = 1.0e-6;
        let below = cl_at(re_grid(2) * (1.0 - eps));
        let above = cl_at(re_grid(2) * (1.0 + eps));
        assert!(
            (below - above).abs() < 1.0e-4 * tau * a,
            "bucket-boundary step: {} vs {}",
            below,
            above
        );
    }

    #[test]
    fn polar_physicality_guard() {
        // rust-foil sometimes "converges" with cd ~ 0 (a failed viscous
        // solve); the degree-9 fit of such a sweep gives spurious cl/cd
        // peaks of 200+ that best_alpha chases.  The guard must reject those
        // sweeps (and NaNs) while accepting healthy drag levels.  (A
        // rejected cached sweep is re-simulated; the flat-plate fallback only
        // applies when that fails too, which cannot be forced here.)
        let mk = |cd_val: f64| {
            let mut alpha = Vec::new();
            let mut cl = Vec::new();
            let mut cd = Vec::new();
            for i in 0..41 {
                let a = (-20.0 + i as f64) * DEG2RAD;
                alpha.push(a);
                cl.push(2.0 * std::f64::consts::PI * a);
                cd.push(cd_val.max(0.01 + 0.3 * a.abs()));
            }
            StoredPolar { alpha, cl, cd }
        };
        // Healthy drag passes.
        assert!(polar_is_physical(&mk(0.02)));
        // Near-zero drag inside the working window is rejected...
        let mut garbage = mk(0.02);
        for (i, cdv) in garbage.cd.iter_mut().enumerate() {
            if garbage.alpha[i].abs() <= 10.0 * DEG2RAD {
                *cdv = 1.0e-9;
            }
        }
        assert!(!polar_is_physical(&garbage), "zero-drag sweep accepted");
        // ...as is NaN drag.
        let mut nan_p = mk(0.02);
        nan_p.cd[20] = f64::NAN;
        assert!(!polar_is_physical(&nan_p), "NaN-drag sweep accepted");
    }

    #[test]
    fn low_re_blends_flat_plate_with_first_bucket() {
        // Below the grid floor the polar blends from the first bucket down
        // to the flat-plate model over [RE_FLAT_PLATE, RE_MIN); below
        // RE_FLAT_PLATE the analytic branch is exact.
        let chord = 0.05;
        let v_for = |re: f64| re * 15.11e-6 / (1.225 * chord);
        let f = Naca4::new(chord, 0.12, 0.0, 0.4);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());

        // A distinctive first bucket: cl = 0.5 * 2 pi alpha.
        seed_bucket(&f, &store, &fs, RE_MIN, v_for(RE_MIN), 0.5);

        let a = 5.0 * DEG2RAD;
        let tau = 2.0 * std::f64::consts::PI;
        // Just below the grid floor: nearly the first bucket's polar.
        let cl = fs.get_cl(v_for(RE_MIN * (1.0 - 1.0e-9)), a);
        assert!((cl - 0.5 * tau * a).abs() < 1e-6, "cl {}", cl);
        // Near the bottom of the blend: nearly the flat plate.
        let cl = fs.get_cl(v_for(RE_FLAT_PLATE * (1.0 + 1.0e-9)), a);
        assert!((cl - tau * a).abs() < 1e-6, "cl {}", cl);
        // Below the blend floor the analytic flat plate is exact.
        let cl = fs.get_cl(v_for(RE_FLAT_PLATE * 0.999), a);
        assert!((cl - tau * a).abs() < 1e-12, "cl {}", cl);
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
