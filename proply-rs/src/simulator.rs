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
use std::sync::{Arc, Mutex, OnceLock};

use crate::cache::{PolarStore, StoredPolar};
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

/// Transition-turbulence level for the e^n model.  Above [`N_CRIT_HI_RE`]
/// the standard XFOIL free-flight value (Ncrit 9) is used; at and below it
/// the environment is noisier (low-Re propellers/rotors run in their own
/// turbulent wake), and Ncrit 9 makes the boundary-layer solve fail
/// wholesale — measured on thin cambered sections (NACA 2406-class):
/// 0 of 81 sweep points converge at Re 15k/30k with Ncrit 9, while
/// Ncrit 5 converges 76-80 of 81 at Re 10k-30k with physical drag.
const N_CRIT_HI_RE: f64 = 100_000.0;
const N_CRIT_LOW_RE: f64 = 5.0;
const N_CRIT_HIGH_RE: f64 = 9.0;

/// The cache key of one simulated polar bucket.  Below [`N_CRIT_HI_RE`] the
/// key is versioned: those buckets were previously simulated with Ncrit 9
/// (or cached as degenerate markers when that failed wholesale), and the
/// Ncrit-5 policy below only takes effect if the old entries are never
/// found again.
fn bucket_cache_key(hash: &str, reynolds: f64, mach: f64) -> String {
    let k = crate::cache::cache_key(hash, reynolds, mach);
    if reynolds <= N_CRIT_HI_RE {
        format!("{k}|n5")
    } else {
        k
    }
}

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

/// Cache keys whose sweep came back degenerate this session: too few
/// converged points, or an unphysical near-zero-drag polar (see
/// [`polar_is_physical`] and [`CD_FLOOR`]).  rust-foil's viscous solve
/// fails *deterministically* on those (foil, Reynolds, Mach) targets, so
/// re-running the sweep reproduces the same garbage — yet every fresh
/// design pass (new simulators, empty fit caches) used to re-run it, once
/// per pass, for each stubborn bucket.  A key lands here the first time a
/// sweep proves it degenerate, and [`bucket_fits`] then falls straight
/// back to the flat-plate model — numerically identical to what the
/// repeated sweeps returned, without the sweeps.
fn bad_keys() -> &'static std::sync::Mutex<std::collections::HashSet<String>> {
    static BAD: std::sync::OnceLock<std::sync::Mutex<std::collections::HashSet<String>>> =
        std::sync::OnceLock::new();
    BAD.get_or_init(|| std::sync::Mutex::new(std::collections::HashSet::new()))
}

/// The sidecar file that carries the degenerate-key markers across runs:
/// the polar-cache path with `.bad.json` appended (`foil_cache.bad.json`
/// beside `foil_cache.json`).  The CLI loads it at startup and saves it at
/// exit, so keys proven degenerate in one run (or one prop of a `make`
/// sweep) are not swept again by the next.
pub fn bad_keys_sidecar(cache_path: &str) -> String {
    format!("{cache_path}.bad.json")
}

/// Load the persisted degenerate-key markers into the session set (the
/// native CLI calls this at startup, next to loading the polar cache).
/// Missing or unreadable files are fine — the markers only suppress
/// re-sweeping keys already proven degenerate.  Deleting the sidecar
/// clears them (e.g. after a code update that might simulate them).
pub fn load_persisted_bad_keys(cache_path: &str) {
    let Ok(text) = std::fs::read_to_string(bad_keys_sidecar(cache_path)) else {
        return;
    };
    if let Ok(keys) = serde_json::from_str::<Vec<String>>(&text) {
        let mut set = bad_keys().lock().unwrap();
        for key in keys {
            set.insert(key);
        }
    }
}

/// Persist the session's degenerate-key markers (the native CLI calls this
/// at exit, next to saving the polar cache).  An empty set removes any
/// stale sidecar.  Only meaningful on hosts with a filesystem — the
/// WebAssembly build never calls it.
pub fn save_persisted_bad_keys(cache_path: &str) {
    let mut keys: Vec<String> = bad_keys().lock().unwrap().iter().cloned().collect();
    keys.sort();
    let path = bad_keys_sidecar(cache_path);
    if keys.is_empty() {
        let _ = std::fs::remove_file(&path);
        return;
    }
    if let Ok(json) = serde_json::to_string(&keys) {
        let _ = std::fs::write(&path, json);
    }
}

/// The Reynolds grid the polars are simulated at, computed once
/// (`np.geomspace(30000, 2e6, 20)`); [`re_grid`] indexes it.
fn re_grid_values() -> &'static [f64; N_RE] {
    static GRID: OnceLock<[f64; N_RE]> = OnceLock::new();
    GRID.get_or_init(|| {
        std::array::from_fn(|i| RE_MIN * (RE_MAX / RE_MIN).powf(i as f64 / (N_RE - 1) as f64))
    })
}

/// Grid point `i` of the Reynolds grid (log-spaced, `RE_MIN..=RE_MAX`).
fn re_grid(i: usize) -> f64 {
    re_grid_values()[i]
}

/// Bracket `re >= RE_MIN` on the Reynolds grid: `(lo, hi, w)` with `lo`/`hi`
/// adjacent grid indices and `w` the log-Re blend weight toward `hi`
/// (0 at `re_grid(lo)`, 1 at `re_grid(hi)`).  Clamped at the top of the
/// grid — no extrapolation past `RE_MAX`.
///
/// The grid is geometric (`re_grid(i) = RE_MIN * q^i`), so the bracket
/// index is a single logarithm; the correction loops only absorb
/// floating-point rounding at exact grid points.  (The old linear scan
/// recomputed grid points with `powf` on every step — up to a microsecond
/// of the hot lookup path at high Reynolds numbers.)
fn re_bracket(re: f64) -> (usize, usize, f64) {
    let grid = re_grid_values();
    let spread = (RE_MAX / RE_MIN).ln();
    let mut lo = (((re / RE_MIN).ln() / spread) * (N_RE - 1) as f64) as usize;
    lo = lo.min(N_RE - 2);
    while lo < N_RE - 2 && grid[lo + 1] <= re {
        lo += 1;
    }
    while lo > 0 && grid[lo] > re {
        lo -= 1;
    }
    let (g_lo, g_hi) = (grid[lo], grid[lo + 1]);
    let w = ((re / g_lo).ln() / (g_hi / g_lo).ln()).clamp(0.0, 1.0);
    (lo, lo + 1, w)
}

/// Zero-pad a polynomial coefficient vector (highest power first, up to
/// [`N_FIT`] entries) into a [`Fit`].
fn into_fit(coeffs: &[f64]) -> Fit {
    debug_assert!(coeffs.len() <= N_FIT);
    let mut fit = [0.0; N_FIT];
    fit[N_FIT - coeffs.len()..].copy_from_slice(coeffs);
    fit
}

/// Elementwise log-Re blend of two fixed fits into `out`
/// (`out = (1 - w) * lo + w * hi`).
fn blend_fit(out: &mut Fit, lo: &Fit, hi: &Fit, w: f64) {
    for i in 0..N_FIT {
        out[i] = (1.0 - w) * lo[i] + w * hi[i];
    }
}

/// A simulated foil: returns CL and CD for a velocity and angle of attack,
/// simulating polars with rust-foil on demand.
///
/// The foil is shared with the owning blade element (`Rc<RefCell>`), exactly
/// like the Python `self.fs = FoilSimulator(self.foil)` — so chord changes
/// made through `set_chord` are seen by the simulator (they feed the
/// Reynolds number).
/// One fitted (cl or cd) polynomial: degree 9, coefficients highest power
/// first, zero-padded to ten entries.  Fixed size so the hot lookup path
/// (bucket blend + evaluation, run millions of times per design) never
/// touches the heap.
type Fit = [f64; N_FIT];
const N_FIT: usize = 10;
/// Keyed cache of polar fits.  The cache belongs to one [`FoilSimulator`]
/// whose foil shape is fixed for its lifetime (chord changes do not alter
/// the foil hash), so the Reynolds bucket alone identifies the fit — the
/// old `(hash String, Reynolds)` key cost a string hash on every lookup.
type PolyCache = RefCell<HashMap<u64, (Fit, Fit)>>;
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
        let (cl_fit, _) = self.blended_fits(v);
        (1.0 - w) * polyval(&cl_fit, alpha) + w * flat
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
        let (_, cd_fit) = self.blended_fits(v);
        (1.0 - w) * polyval(&cd_fit, alpha) + w * flat
    }

    /// The fitted (cl, cd) polynomials for the velocity `v`, interpolated
    /// across Reynolds buckets (see the module docs): the bracketing grid
    /// buckets' fits are blended in log-Re, and below the grid floor the
    /// blend runs from the first bucket down to the flat-plate model.
    pub fn get_polars(&self, v: f64) -> (Vec<f64>, Vec<f64>) {
        let (cl_fit, cd_fit) = self.blended_fits(v);
        (cl_fit.to_vec(), cd_fit.to_vec())
    }

    /// The blended (cl, cd) fits for the velocity `v` — exactly the
    /// [`FoilSimulator::get_polars`] coefficient blend, computed into
    /// fixed arrays with no per-call heap allocation (the design loops run
    /// this millions of times per pass).
    fn blended_fits(&self, v: f64) -> (Fit, Fit) {
        let re = self.foil.borrow().reynolds(v);
        let mach = self.get_mach(v);

        if re < RE_FLAT_PLATE {
            return *flat_plate_fit();
        }
        if re < RE_MIN {
            let fp = flat_plate_fit();
            let (cl_hi, cd_hi) = self.bucket_fits(RE_MIN, mach);
            let w = (re / RE_FLAT_PLATE).ln() / (RE_MIN / RE_FLAT_PLATE).ln();
            let (mut cl, mut cd) = ([0.0; N_FIT], [0.0; N_FIT]);
            blend_fit(&mut cl, &fp.0, &cl_hi, w);
            blend_fit(&mut cd, &fp.1, &cd_hi, w);
            return (cl, cd);
        }

        let (lo, hi, w) = re_bracket(re);
        let (cl_lo, cd_lo) = self.bucket_fits(re_grid(lo), mach);
        let (cl_hi, cd_hi) = self.bucket_fits(re_grid(hi), mach);
        let (mut cl, mut cd) = ([0.0; N_FIT], [0.0; N_FIT]);
        blend_fit(&mut cl, &cl_lo, &cl_hi, w);
        blend_fit(&mut cd, &cd_lo, &cd_hi, w);
        (cl, cd)
    }

    /// The fitted (cl, cd) polynomials of the single grid bucket at
    /// `reynolds` (a [`re_grid`] point): the memoized degree-9 fit of the
    /// cached polar, simulated with rust-foil on first use.  This is the
    /// bucket-level fetch the interpolation blends.
    fn bucket_fits(&self, reynolds: f64, mach: f64) -> (Fit, Fit) {
        // The fit cache is per simulator and its foil shape is fixed for
        // the simulator's lifetime (chord changes do not alter the foil
        // hash), so the Reynolds bucket alone keys the fit.
        let key = reynolds.to_bits();
        if let Some(fits) = self.poly_cache.borrow().get(&key) {
            return *fits;
        }

        let hash = self.foil.borrow().hash();
        let ck = bucket_cache_key(&hash, reynolds, mach);
        // A key this session (or a previous run — see
        // load_persisted_bad_keys) already proved degenerate: skip the
        // store and the sweep entirely — the flat-plate fallback is what
        // the sweep would have returned anyway.
        if bad_keys().lock().unwrap().contains(&ck) {
            let fp = *flat_plate_fit();
            self.poly_cache.borrow_mut().insert(key, fp);
            return fp;
        }
        let cached = {
            let store = self.store.lock().unwrap();
            store.get(&ck).cloned()
        };
        let fits = match &cached {
            Some(p) if p.alpha.len() > 20 && polar_is_physical(p) => fit_stored(p),
            _ => {
                // Claim the key (blocking while another thread simulates
                // it), then inspect the store once the claim is ours — the
                // previous holder may have stored a good polar in the
                // meantime.
                {
                    let mut keys = sim_keys().lock().unwrap();
                    while keys.contains(&ck) {
                        keys = sim_cv().wait(keys).unwrap();
                    }
                    keys.insert(ck.clone());
                }
                let fits = {
                    let store = self.store.lock().unwrap();
                    let entry = store.get(&ck).cloned();
                    match entry {
                        // A good cached polar: fit it.
                        Some(p) if p.alpha.len() > 20 && polar_is_physical(&p) => fit_stored(&p),
                        // A *stored but degenerate* sweep.  rust-foil fails
                        // deterministically on these targets, so the entry
                        // itself is the proof: re-sweeping reproduces the
                        // same garbage.  Blacklist the key for the session
                        // and take the flat-plate fallback immediately —
                        // the value the doomed sweep would have produced.
                        Some(_) => {
                            drop(store);
                            bad_keys().lock().unwrap().insert(ck.clone());
                            *flat_plate_fit()
                        }
                        // Nothing stored: this is the cold first touch (or
                        // the previous holder's sweep stored nothing, which
                        // the blacklist check below catches).
                        None => {
                            drop(store);
                            if bad_keys().lock().unwrap().contains(&ck) {
                                *flat_plate_fit()
                            } else {
                                self.xfoil_simulate_polars(reynolds, mach);
                                let store = self.store.lock().unwrap();
                                match store.get(&ck) {
                                    Some(p) if p.alpha.len() > 20 && polar_is_physical(p) => {
                                        fit_stored(p)
                                    }
                                    // Degenerate case: too few converged
                                    // points, or an unphysical (near-zero
                                    // drag) sweep.  The Python code recurses
                                    // forever here (its fallback never writes
                                    // to the DB); we return the flat-plate
                                    // model instead — and remember the key so
                                    // the design never pays for this doomed
                                    // sweep again.
                                    _ => {
                                        drop(store);
                                        bad_keys().lock().unwrap().insert(ck.clone());
                                        *flat_plate_fit()
                                    }
                                }
                            }
                        }
                    }
                };
                let mut keys = sim_keys().lock().unwrap();
                keys.remove(&ck);
                sim_cv().notify_all();
                fits
            }
        };

        self.poly_cache.borrow_mut().insert(key, fits);
        fits
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
        self.bucket_fits(reynolds, mach);
    }

    /// Simulate and cache the polar buckets an evaluation at speed `v`
    /// would need ([`warm_plan`] + [`warm_bucket`]).
    pub fn warm_polars(&self, v: f64) {
        for (re, mach) in self.warm_plan(v) {
            self.warm_bucket(re, mach);
        }
    }

    /// Simulate the polar with rust-foil and store it in the shared cache:
    /// the polar itself when the sweep converged, a degenerate marker when
    /// it failed (so the failure is cached too and never re-attempted).
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
        // Low-Re buckets trip transition early (see [`N_CRIT_HI_RE`]): the
        // e^n critical amplification is lowered there so the viscous solve
        // converges at all.
        if reynolds <= N_CRIT_HI_RE {
            xf.set_n_crit(N_CRIT_LOW_RE);
        } else {
            xf.set_n_crit(N_CRIT_HIGH_RE);
        }

        // The first operating point can diverge wholesale (a sweep that
        // "converged" 0 of 81 points in a fraction of a second: the BL
        // state at the sweep's starting alpha fails and the whole sequence
        // aborts, while a fresh solve of the same foil converges 79/81).
        // Retry once with a brand-new solve before declaring the bucket
        // degenerate — a doomed key still ends in a marker after both
        // attempts, so the cost is one extra sweep per truly failing bucket.
        let polar = xf.aseq(-20.0, 20.0, 80);
        let polar = if polar.iter().filter(|p| p.5).count() < 5 {
            let mut xf2 = rust_foil::XFoil::new();
            xf2.set_show_output(false);
            xf2.set_airfoil(&kept_x, &kept_y);
            xf2.set_reynolds(reynolds);
            xf2.set_max_iter(80);
            if reynolds <= N_CRIT_HI_RE {
                xf2.set_n_crit(N_CRIT_LOW_RE);
            } else {
                xf2.set_n_crit(N_CRIT_HIGH_RE);
            }
            xf2.aseq(-20.0, 20.0, 80)
        } else {
            polar
        };

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

        let key = bucket_cache_key(&self.foil.borrow().hash(), reynolds, mach);
        let mut store = self.store.lock().unwrap();
        if alfa.len() < 5 {
            // The sweep failed — nothing usable to fit (see get_polars).
            // Cache the *failure* as a marker polar: a marker never passes
            // the polar checks, so any later fetch — this run or a future
            // one, since markers persist like any other polar — treats the
            // key as degenerate and never runs the doomed sweep again.
            // Without this, every fresh run re-attempted each failing
            // bucket once (the seeding warm-up on thick mechanical-law
            // root sections failed ~half its sweeps).
            crate::dprintln!(
                "polar sweep FAILED ({} pts): {} -> cached as degenerate marker",
                alfa.len(),
                key
            );
            store.insert(
                key,
                StoredPolar {
                    alpha: Vec::new(),
                    cl: Vec::new(),
                    cd: Vec::new(),
                },
            );
            return;
        }

        crate::dprintln!("polar sweep: {} ({} pts)", key, alfa.len());
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

/// [`fit_polar`] into fixed [`Fit`] arrays.
fn fit_stored(p: &StoredPolar) -> (Fit, Fit) {
    let (cl, cd) = fit_polar(p);
    (into_fit(&cl), into_fit(&cd))
}

/// The flat-plate fallback fits, computed once: degree 4 over the same
/// degree-abscissa grid as [`flat_plate_polys`], zero-padded to [`Fit`]s.
fn flat_plate_fit() -> &'static (Fit, Fit) {
    static FP: OnceLock<(Fit, Fit)> = OnceLock::new();
    FP.get_or_init(|| {
        let (cl, cd) = flat_plate_polys();
        (into_fit(&cl), into_fit(&cd))
    })
}

/// The flat-plate fallback model, fitted like the Python `Foil didn't
/// simulate` branch (degree 4 over the ±20° alpha sweep).
///
/// The fit's abscissa is the angle in *radians*, exactly as the fitted
/// polars' (a degree-grid abscissa was a units bug: evaluating the degree
/// fit at radian arguments shrank the drag ~57x, so the "flat-plate"
/// fallback advertised L/D ~ 281 instead of ~4.9 and best-L/D scans rode
/// its corner).
fn flat_plate_polys() -> (Vec<f64>, Vec<f64>) {
    let alpha: Vec<f64> = (-20..=20).map(|a| (a as f64) * DEG2RAD).collect();
    let cl: Vec<f64> = alpha.iter().map(|a| 2.0 * std::f64::consts::PI * a).collect();
    let cd: Vec<f64> = alpha.iter().map(|a| 1.28 * a.sin()).collect();
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
        let key = bucket_cache_key(&f.hash(), g, fs.get_mach(v_for_g));
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
    fn stored_degenerate_polar_is_not_resimulated() {
        // A stored sweep with unphysical (near-zero) drag must fall back to
        // the flat-plate model WITHOUT re-running the sweep: rust-foil
        // fails deterministically on those targets, and re-sweeping them on
        // every design pass used to dominate warm runs (the honda cache
        // held 29 such polars of 348; a re-run spent ~90% of its time
        // re-simulating them).  If the code regressed and re-swept, the
        // *real* polar (cd ~ 0.008 at 5°) would come back instead of the
        // fitted flat-plate fallback, and this assertion fails.
        let chord = 0.06;
        let v = 22.0; // Re ~ 107k, between buckets 4 and 5
        // A foil configuration no other test uses, so the session-wide
        // degenerate-key blacklist below cannot leak into their keys.
        let f = Naca4::new(chord, 0.14, 0.01, 0.35);
        let store = test_store();
        let fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), store.clone());
        let re_raw = f.reynolds(v);
        let (lo, hi, _) = re_bracket(re_raw);
        let tau = 2.0 * std::f64::consts::PI;
        for g in [re_grid(lo), re_grid(hi)] {
            let mut alpha = Vec::new();
            let mut cl = Vec::new();
            let mut cd = Vec::new();
            for i in 0..41 {
                let a = (-20.0 + i as f64) * DEG2RAD;
                alpha.push(a);
                cl.push(tau * a);
                cd.push(1.0e-9); // the degenerate signature: ~zero drag
            }
            let vg = g * 15.11e-6 / (1.225 * chord);
            let ma = f.mach(vg);
            let mach = ((ma * 20.0).round()) / 20.0;
            let key = bucket_cache_key(&f.hash(), g, mach);
            store
                .lock()
                .unwrap()
                .insert(key, StoredPolar { alpha, cl, cd });
        }
        let a = 5.0 * DEG2RAD;
        let cd = fs.get_cd(v, a);
        // The fallback is the fitted flat-plate model (its own quirk: the
        // fit's abscissa is degrees, so at a radian alpha it returns the
        // near-zero-drag value — the point is the value must come from the
        // fallback, never from a fresh sweep of the degenerate bucket).
        let (_, cd_fp) = flat_plate_polys();
        let cd_expected = polyval(&cd_fp, a);
        assert!(
            (cd - cd_expected).abs() < 1.0e-12,
            "cd {} should be the flat-plate fallback {} (a stored degenerate polar must not be re-simulated)",
            cd,
            cd_expected
        );
    }

    #[test]
    fn bad_key_sidecar_round_trips() {
        // The CLI persists the session's degenerate-key markers so future
        // runs skip the doomed sweeps.  The probe key is unique to this
        // test, so the shared session set is only ever *extended* by it
        // (never cleared) and no other test can collide with it.
        let dir = std::env::temp_dir().join(format!("proply-bad-keys-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cache = dir.join("cache.json");
        let cs = cache.to_str().unwrap().to_string();
        let side = bad_keys_sidecar(&cs);
        let probe = "proply-probe-bad-key|12345.6|0.05";

        // A hand-written sidecar loads into the session set...
        std::fs::write(&side, format!(r#"["{probe}"]"#)).unwrap();
        load_persisted_bad_keys(&cs);
        assert!(bad_keys().lock().unwrap().contains(probe));
        // ...and saving persists it (with whatever else the session holds).
        save_persisted_bad_keys(&cs);
        let text = std::fs::read_to_string(&side).expect("sidecar written");
        let parsed: Vec<String> = serde_json::from_str(&text).expect("sidecar parses");
        assert!(parsed.contains(&probe.to_string()), "probe key persisted");

        bad_keys().lock().unwrap().remove(probe);
        let _ = std::fs::remove_dir_all(&dir);
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