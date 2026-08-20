// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! The propeller design loop, ported from `proply/prop.py`.
//!
//! `full_optimize` sweeps radial stations from tip to hub, sizing each
//! blade element with `optimize_all`, then smooths the twist and chord
//! distributions and re-evaluates the forces with `get_forces`.

use std::cell::RefCell;
use std::collections::HashSet;
use std::rc::Rc;
use std::sync::{Arc, Mutex};
use std::thread;

use crate::blade_element::BladeElement;
use crate::cache::PolarStore;
use crate::design_parameters::DesignParameters;
use crate::foil::{Cst, FoilFamily, FoilLike, Naca4};
use crate::lift_line::{self, Station};
use crate::optimize;
use crate::pchip::Pchip;
use crate::polyfit::{polyfit, polyval};
use crate::simulator::FoilSimulator;
use crate::smooth::smooth;
use indicatif::{ProgressBar, ProgressStyle};

/// Camber values scanned by the lifting-line design when no explicit camber
/// is given (the JSON `camber` key / `--camber`): each candidate gets a full
/// design pass, plus one composed per-station distribution (the best of
/// these candidates at each radius, smoothed — see [`compose_camber`]);
/// the best torque at matched thrust wins.  A small fixed set, not a
/// continuous variable — the polar cache hashes camber at 0.01 granularity,
/// so a continuous camber would turn the objective into a staircase of
/// discrete polar families.
pub const CAMBER_CANDIDATES: [f64; 3] = [0.0, 0.02, 0.04];

/// Minimum chord (m) anywhere on the blade.
const CH_FLOOR: f64 = 0.002;
/// Bounds on the common attack-angle offset `da` over the best-L/D angles
/// (rad).  `DA_MAX` sits ~7 deg above best-L/D, below deep stall; `DA_MIN`
/// may also trim *through* zero lift: thrust is monotone in `da` across the
/// whole range, so the inner bisection can match a low thrust target at any
/// chord — without it, a full-chord blade whose thrust already exceeds the
/// target could never converge, as `da > 0` only adds thrust.
const DA_MIN: f64 = -0.45;
const DA_MAX: f64 = 0.12;

/// Allowed chord as a function of radius (m): `tip_chord*R^2/r^2`, capped by
/// the blade-count spacing (and the local twist).  The parameter form shared
/// by the design passes (worker threads cannot borrow `Prop`).
fn max_chord(param: &DesignParameters, n_blades: usize, r: f64, twist: f64) -> f64 {
    let k = param.tip_chord * param.radius.powi(2);
    let c = k / r.powi(2);
    let upper_limit = (2.0 * std::f64::consts::PI * r / (n_blades as f64 + 2.0)) / twist.cos();
    c.min(upper_limit)
}

/// Free-chord law `tip_chord*R^2/r^2`, capped by the geometric limits
/// ([`max_chord`]).  `scale < 1` thins the blade for a higher aspect ratio.
/// The operation order matches the original method exactly — the optimizer
/// trajectories are sensitive to last-ulp changes.
fn chord_law(param: &DesignParameters, n_blades: usize, r: f64, twist: f64, scale: f64) -> f64 {
    let k = param.tip_chord * param.radius.powi(2);
    let cap = max_chord(param, n_blades, r, twist);
    (scale * k / r.powi(2)).min(cap).max(1.0e-6)
}

/// Allowed foil thickness (m) as a function of radius: a p = 0.3 power law
/// between the hub depth and a tenth of it at the tip.
fn foil_thickness(param: &DesignParameters, r: f64) -> f64 {
    let thickness_root = param.hub_depth * 1.0;
    let thickness_end = param.hub_depth * 0.1;
    let p = 0.3;
    let k = (thickness_root - thickness_end)
        / (param.hub_radius.powf(p) - param.radius.powf(p));
    let s = thickness_end - k * param.radius.powf(p);
    s + k * r.powf(p)
}

/// The foil of one station: the thickness law at radius `r`, chord `c`,
/// camber fraction `camber` and the parameter trailing edge.  The family is
/// [`DesignParameters::cst`]: the default NACA 4-series, or a CST (Kulfan)
/// section — the default 18-parameter shape re-thicknessed and cambered to
/// the same laws.  Plain data (clonable, sendable) so worker threads can
/// build their own.
fn station_foil(param: &DesignParameters, r: f64, c: f64, camber: f64) -> FoilFamily {
    let thickness = foil_thickness(param, r);
    if param.cst {
        let mut f = Cst::default(c);
        f.set_thickness(thickness / c);
        f.set_camber(camber);
        f.set_trailing_edge(param.trailing_edge / 1000.0);
        FoilFamily::Cst(f)
    } else {
        let mut f = Naca4::new(c, thickness / c, camber, 0.4);
        f.base.set_trailing_edge(param.trailing_edge / 1000.0);
        FoilFamily::Naca4(f)
    }
}

/// A propeller: a collection of blade elements plus the design parameters.
pub struct Prop {
    pub param: DesignParameters,
    pub radial_resolution: f64,
    pub radial_steps: usize,
    pub n_blades: usize,
    pub blade_elements: Vec<BladeElement<FoilFamily>>,
    store: Arc<Mutex<PolarStore>>,
    scimitar_interpolator: Option<Pchip>,
    max_depth_interpolator: Option<(Vec<f64>, Vec<f64>)>,
    plate_mode: bool,
}

impl Prop {
    pub fn new(param: DesignParameters, resolution: f64, store: Arc<Mutex<PolarStore>>) -> Self {
        let radial_steps = (param.radius / resolution) as usize;
        Self {
            param,
            radial_resolution: resolution,
            radial_steps,
            n_blades: 2,
            blade_elements: Vec::new(),
            store,
            scimitar_interpolator: None,
            max_depth_interpolator: None,
            plate_mode: false,
        }
    }

    /// Use the analytic flat-plate polar model for every blade element
    /// (testing only — mirrors the Python `PlateSimulatedFoil`).
    pub fn set_plate_mode(&mut self, on: bool) {
        self.plate_mode = on;
    }

    /// Create a blade element with a foil sized for the station (the
    /// configured family, default NACA4).
    fn new_blade_element(&mut self, r: f64, rpm: f64, twist: f64) -> BladeElement<FoilFamily> {
        let y_limit = self.get_max_depth(r);
        let x_limit = self.get_max_chord(r, twist);

        let foil = Rc::new(RefCell::new(station_foil(
            &self.param,
            r,
            x_limit,
            self.param.camber.unwrap_or(0.0),
        )));

        let c_max = {
            let f = foil.borrow();
            f.get_max_chord(x_limit, y_limit, twist)
        };
        println!("Max Chord {}", c_max);
        foil.borrow_mut().modify_chord(c_max);

        let mut be = BladeElement::new(
            r,
            self.radial_resolution,
            foil,
            twist,
            rpm,
            self.param.forward_airspeed,
            self.store.clone(),
        );
        if self.plate_mode {
            be.set_plate_mode(true);
        }
        be
    }

    /// Create a blade element with an explicit chord `c` (the lifting-line
    /// design sizes the chord directly instead of the geometric fit) and
    /// camber fraction `camber` (NACA 4-series `m`).
    #[cfg(test)]
    fn new_blade_element_with_chord(
        &mut self,
        r: f64,
        rpm: f64,
        twist: f64,
        c: f64,
        camber: f64,
    ) -> BladeElement<FoilFamily> {
        let foil = Rc::new(RefCell::new(station_foil(&self.param, r, c, camber)));

        let mut be = BladeElement::new(
            r,
            self.radial_resolution,
            foil,
            twist,
            rpm,
            self.param.forward_airspeed,
            self.store.clone(),
        );
        if self.plate_mode {
            be.set_plate_mode(true);
        }
        be
    }

    /// Allowed chord as a function of radius (m): k / r^2, capped by the
    /// blade-count spacing.  See [`max_chord`] (the parameter form the
    /// design passes use).
    pub fn get_max_chord(&self, r: f64, twist: f64) -> f64 {
        max_chord(&self.param, self.n_blades, r, twist)
    }

    /// Scimitar offset (m) as a function of radius.
    pub fn get_scimitar_offset(&mut self, r: f64) -> f64 {
        if self.scimitar_interpolator.is_none() {
            let hub_r = self.param.hub_radius;
            let max_r = self.param.radius * 0.8;
            let max_c = self.param.radius * (self.param.scimitar_percent / 100.0);
            let x = vec![0.0, hub_r, max_r, self.param.radius];
            let y = vec![0.0, 1.1 * 0.0, max_c, 0.0];
            self.scimitar_interpolator = Some(Pchip::new(&x, &y));
        }
        self.scimitar_interpolator.as_ref().unwrap().eval(r)
    }

    /// Allowed foil thickness (m) as a function of radius: a p = 0.3 power
    /// law between the hub depth and a tenth of it at the tip.  See
    /// [`foil_thickness`] (the parameter form the design passes use).
    pub fn get_foil_thickness(&self, r: f64) -> f64 {
        foil_thickness(&self.param, r)
    }

    /// Allowed depth of the prop as a function of radius (m), from the
    /// fixed environment points in the Python (linear interpolation).
    pub fn get_max_depth(&mut self, r: f64) -> f64 {
        if self.max_depth_interpolator.is_none() {
            let hub_r = self.param.hub_radius;
            let hub_depth = self.param.hub_depth;
            let max_depth = self.param.hub_depth * 3.0;
            let max_r = self.param.radius / 3.0;
            let end_depth = self.param.hub_depth * 2.0;
            let x = vec![
                0.0,
                hub_r / 2.0,
                hub_r,
                1.1 * hub_r,
                1.5 * hub_r,
                max_r,
                0.9 * self.param.radius,
                self.param.radius,
            ];
            let y = vec![
                hub_depth,
                hub_depth,
                hub_depth,
                hub_depth,
                1.1 * hub_depth,
                max_depth,
                1.2 * end_depth,
                end_depth,
            ];
            self.max_depth_interpolator = Some((x, y));
        }
        let (x, y) = self.max_depth_interpolator.as_ref().unwrap();
        lin_interp(x, y, r)
    }

    /// Prandtl-style tip and hub loss factor.
    pub fn tip_loss(&self, r: f64, phi: f64) -> f64 {
        let f = (self.n_blades as f64 * (self.param.radius - r * 0.96))
            / (2.0 * r * phi.sin());
        let tip_loss = 2.0 * (-f).exp().acos() / std::f64::consts::PI;
        let f = (self.n_blades as f64 * (r - self.param.hub_radius * 0.95))
            / (2.0 * r * phi.sin());
        let hub_loss = 2.0 * (-f).exp().acos() / std::f64::consts::PI;
        tip_loss * hub_loss
    }

    /// Total torque and thrust at a given RPM, re-running the BEM on each
    /// element (mirrors `get_forces`).
    pub fn get_forces(&mut self, rpm: f64) -> (f64, f64) {
        let mut torque = 0.0;
        let mut thrust = 0.0;
        for be in self.blade_elements.iter_mut() {
            be.rpm = rpm;
            be.omega = optimize::rpm2omega(rpm);
            let dv_goal = be.dv;
            let (dv, a_prime, err) = be.bem(self.n_blades);

            if err < 0.01 {
                let dt = be.d_t();
                let dm = be.d_m();
                thrust += dt;
                torque += dm;
                println!(
                    "r={}, theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={}",
                    be.r,
                    be.get_twist().to_degrees(),
                    dv,
                    a_prime,
                    dt,
                    dm,
                    dt / dm
                );
            } else {
                eprintln!(
                    "r={}: BEM did not converge {} {} {}",
                    be.r, be.dv, dv_goal, a_prime
                );
                if err > 0.5 {
                    be.dv = 1.0;
                    be.a_prime = 0.0;
                }
            }
        }
        (torque, thrust)
    }

    /// Design every blade element for the target thrust at the optimum RPM,
    /// smooth the twist/chord distributions, and return (torque, thrust).
    pub fn full_optimize(&mut self, optimum_rpm: f64, thrust: f64) -> (f64, f64) {
        self.blade_elements.clear();
        let u_0 = self.param.forward_airspeed;

        let dv_goal = optimize::dv_from_thrust(thrust, self.param.radius, u_0);
        let radial_points: Vec<f64> = (0..self.radial_steps)
            .map(|i| {
                self.param.radius
                    - (self.param.radius - self.param.hub_radius) * i as f64
                        / (self.radial_steps - 1) as f64
            })
            .collect();
        // radial_points runs tip -> hub, as in the Python.

        let mut total_thrust = 0.0;
        let mut total_torque = 0.0;
        let omega = optimize::rpm2omega(optimum_rpm);
        let dr = (radial_points[0] - radial_points[1]).abs();

        let mut twist_angles: Vec<f64> = Vec::new();
        let mut chords: Vec<f64> = Vec::new();
        let mut prev_twist = 0.0;

        // Station bar: each iteration sizes one blade element, triggering the
        // lazy first-touch polar sweeps on a cold cache.
        let pb = ProgressBar::new(radial_points.len() as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] {pos}/{len} BEM stations (eta {eta})",
            )
            .unwrap(),
        );

        for r in &radial_points {
            let u = u_0 + dv_goal;
            let v = omega * r;
            let phi = (u / v).atan();

            let dv_modified = dv_goal * self.tip_loss(*r, phi);
            let mut be = self.new_blade_element(*r, optimum_rpm, prev_twist);

            let y_limit = self.get_max_depth(*r);
            let x_limit = self.get_max_chord(*r, prev_twist);
            let maxchord = be.foil.borrow().get_max_chord(x_limit, y_limit, prev_twist);

            let x = optimize::optimize_all(
                &be.fs,
                dv_modified,
                optimum_rpm,
                *r,
                dr,
                u_0,
                self.n_blades as f64,
                maxchord,
            )
            .0;
            let (theta, dv, a_prime, chord) = (x[0], x[1], x[2], x[3]);
            be.set_chord(chord);
            be.set_twist(theta);
            be.set_bem(dv, a_prime);
            twist_angles.push(theta);
            chords.push(be.foil.borrow().chord());

            let d_t = be.d_t();
            let d_m = be.d_m();
            total_thrust += d_t;
            total_torque += d_m;

            pb.println(format!(
                "r={} theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={} ",
                r,
                theta.to_degrees(),
                dv,
                a_prime,
                d_t,
                d_m,
                d_t / d_m
            ));
            pb.println(be.to_string());

            self.blade_elements.push(be);
            prev_twist = theta;
            pb.inc(1);
        }
        pb.finish();

        self.blade_elements.reverse();
        twist_angles.reverse();
        chords.reverse();

        // Smooth the twist (degree-4 polyfit) and the chords (PCHIP over the
        // smoothed distribution), as in `full_optimize`.
        let radial_hub_to_tip: Vec<f64> = radial_points.iter().rev().copied().collect();
        let twist_poly = polyfit(&radial_hub_to_tip, &twist_angles, 4);

        let mut c_points: Vec<f64> =
            vec![0.0, self.param.hub_radius / 2.0, 0.9 * self.param.hub_radius];
        c_points.extend(radial_hub_to_tip.iter());
        let mut extra_chords: Vec<f64> = vec![
            0.9 * self.param.hub_depth,
            0.9 * self.param.hub_depth,
            0.9 * self.param.hub_depth,
        ];
        extra_chords.extend(chords.iter());
        let smoothed = smooth(&extra_chords, 11, "hanning");
        let chord_poly = Pchip::new(&c_points, &smoothed);

        println!("Smoothed Blade Form");
        for be in self.blade_elements.iter_mut() {
            let c = chord_poly.eval(be.r);
            let t = polyval(&twist_poly, be.r);
            be.set_chord(c);
            be.set_twist(t);
            println!("{}", be);
        }

        let (torque, thrust_final) = self.get_forces(optimum_rpm);
        println!("Total Thrust: {:5.2}, Torque: {:5.3}", thrust_final, torque);
        let _ = (total_thrust, total_torque);
        (torque, thrust_final)
    }

    /// The maximum CL/CD at speed `v` and the angle of attack (rad) that
    /// gives it, found by scanning the polar and refining around the winner.
    /// Uses the cached polar, so repeated calls are cheap after the first.
    fn best_ld<F: FoilLike>(fs: &FoilSimulator<F>, v: f64) -> (f64, f64) {
        let ld = |alpha: f64| -> f64 {
            let cl = fs.get_cl(v, alpha);
            let cd = fs.get_cd(v, alpha).max(1.0e-6);
            cl / cd
        };
        let mut best = 0.0;
        let mut best_ld = f64::NEG_INFINITY;
        for a in -8..=16 {
            let alpha = (a as f64).to_radians();
            let l = ld(alpha);
            if l > best_ld {
                best_ld = l;
                best = alpha;
            }
        }
        // Degenerate/analytic models (flat-plate: cl/cd constant) tie at every
        // alpha; fall back to a working moderate attack angle instead of the
        // argmax tie (which would pick the scan start).
        if best_ld.is_finite() && best_ld <= 4.9 + 1.0e-6 {
            return (0.10, best_ld);
        }
        // The scan quantises the angle to whole degrees; refine to a
        // continuous maximum with a golden-section search in the +/- 1 deg
        // neighbourhood of the winner (cl/cd is a smooth polynomial ratio of
        // alpha there).
        const INV_PHI: f64 = 0.618_033_988_749_895;
        let step = 1.0_f64.to_radians();
        let (mut lo, mut hi) = (best - step, best + step);
        let (mut c, mut d) = (hi - INV_PHI * (hi - lo), lo + INV_PHI * (hi - lo));
        let (mut fc, mut fd) = (ld(c), ld(d));
        for _ in 0..40 {
            if fc > fd {
                hi = d;
                d = c;
                fd = fc;
                c = hi - INV_PHI * (hi - lo);
                fc = ld(c);
            } else {
                lo = c;
                c = d;
                fc = fd;
                d = lo + INV_PHI * (hi - lo);
                fd = ld(d);
            }
        }
        (0.5 * (lo + hi), fc.max(fd))
    }

    /// Populate the polar cache for `(foil, v)` work items on all available
    /// cores: the first-touch rust-foil sweeps dominate a cold design run,
    /// and the pool overlaps them.  The caller must deduplicate items
    /// (workers do not coordinate).  No-op in plate mode.  Progress is
    /// reported on a bar labelled `what` (hidden when stderr is not a TTY).
    fn warm_polar_pool(&self, work: &[(FoilFamily, f64)], what: &str) {
        if work.is_empty() || self.plate_mode {
            return;
        }
        let pb = ProgressBar::new(work.len() as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} foil polars ({msg}, eta {eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.set_message(what.to_string());
        let queue = Mutex::new(work.to_vec());
        let n_workers = thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
            .min(work.len());
        let store = self.store.clone();
        thread::scope(|s| {
            for _ in 0..n_workers {
                let queue = &queue;
                let store = store.clone();
                let pb = pb.clone();
                s.spawn(move || {
                    while let Some((foil, v)) = queue.lock().unwrap().pop() {
                        let fs = FoilSimulator::new(Rc::new(RefCell::new(foil)), store.clone());
                        fs.warm_polars(v);
                        pb.inc(1);
                    }
                });
            }
        });
        pb.finish();
    }

    /// One camber candidate's design pass: build the station elements from
    /// `camber_dist`, prescribe the attack angles `alpha_base` (plus a common
    /// offset `da` matched to the thrust target by an inner monotone
    /// bisection) and optimise the chord-spline controls with three
    /// warm-chained Nelder-Mead starts.  Self-contained (own elements and
    /// warm state, parameter-only geometry laws) so the camber candidates can
    /// run in parallel on scoped threads; returns the measured geometry, or
    /// `None` when no thrust match was found.
    #[allow(clippy::too_many_arguments)]
    fn run_design_pass(
        param: &DesignParameters,
        n_blades: usize,
        plate_mode: bool,
        store: &Arc<Mutex<PolarStore>>,
        radial_res: f64,
        rr: &[f64],
        ref_r: &[f64],
        cap_ctl: &[f64],
        infl: &lift_line::Influence,
        rpm: f64,
        omega: f64,
        u_0: f64,
        thrust: f64,
        n_ctrl: usize,
        s_cap: f64,
        label: &str,
        camber_dist: &[f64],
        alpha_base: &[f64],
    ) -> Option<PassOutcome> {
        let m = rr.len();

        // Station elements: the seed chord law at each radius, twisted to a
        // first guess (the circulation solve sets the final twist).
        let mut elements: Vec<BladeElement<FoilFamily>> = Vec::with_capacity(m);
        for (i, &r) in rr.iter().enumerate() {
            let phi0 = u_0.atan2(omega * r).max(1.0e-3);
            let c = chord_law(param, n_blades, r, 0.0, s_cap);
            let foil = Rc::new(RefCell::new(station_foil(param, r, c, camber_dist[i])));
            let mut be = BladeElement::new(r, radial_res, foil, phi0 + 0.12, rpm, u_0, store.clone());
            if plate_mode {
                be.set_plate_mode(true);
            }
            elements.push(be);
        }

        // Design variables `x ∈ [0,1]^n_ctrl` log-map to the control chords;
        // the common attack offset `da ∈ [DA_MIN, DA_MAX]` over the best-L/D
        // angles is matched to the thrust target by an inner monotone
        // bisection (`meet_thrust` below), not optimized.  `alpha_i =
        // best_L/D_i + da` is prescribed (`twist = phi + alpha`), so there is
        // no twist<->induction feedback; the circulation solve is converged
        // by damped Newton (warm-started from the previous evaluation).  The
        // chord at any radius is the smooth PCHIP spline through the
        // controls at `ref_r` (clipped to the local geometric cap for
        // safety), so there are no taper/cap kinks.
        let state = RefCell::new(elements);
        let pg = RefCell::new(vec![0.0; m]); // warm-start gamma
        let da_hint = RefCell::new(None::<f64>);

        // eval(controls, da): chord = smooth shape-preserving cubic (PCHIP)
        // spline through the N control values at `ref_r` (kink-free, clipped to
        // the geometrically-allowed chord for safety), attack angle prescribed.
        let eval = |controls: &[f64],
                    da: f64,
                    elems: &mut Vec<BladeElement<FoilFamily>>,
                    seed: &[f64]|
         -> (f64, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
            let alphas: Vec<f64> = alpha_base
                .iter()
                .map(|&ab| (ab + da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX))
                .collect();
            let spl = Pchip::new(ref_r, controls);
            for i in 0..m {
                let cap_i = max_chord(param, n_blades, rr[i], alphas[i]).max(CH_FLOOR);
                let c = spl.eval(rr[i]).clamp(CH_FLOOR, cap_i);
                elems[i].set_chord(c);
            }
            let stations: Vec<Station> = (0..m)
                .map(|i| Station {
                    r: rr[i],
                    c: elems[i].foil.borrow().chord(),
                    alpha: alphas[i],
                })
                .collect();
            let fs_refs: Vec<&FoilSimulator<FoilFamily>> = elems.iter().map(|be| &be.fs).collect();
            let res = lift_line::solve_with_influence(&stations, n_blades, omega, u_0, &fs_refs, seed, infl);
            for i in 0..m {
                elems[i].set_twist(res.phi[i] + alphas[i]);
            }
            let chords: Vec<f64> = elems.iter().map(|be| be.foil.borrow().chord()).collect();
            (res.thrust, res.torque, alphas, chords, res.gamma.clone(), res.phi.clone())
        };

        // Outer optimizer over the N control values: for each candidate shape,
        // an inner **monotone `da` bisection** (below) reliably matches the
        // thrust target, so the outer NelderMead only needs to minimise the
        // torque (efficiency) among shapes that meet the target.  Seeded at the
        // full chord — the thrust-capable geometry — so it cannot be trapped in
        // a thin/low-thrust region.
        // (da, err, thrust, torque, alpha, chord, phi) — one design candidate.
        type MeetOutcome = (f64, f64, f64, f64, Vec<f64>, Vec<f64>, Vec<f64>);
        let meet_thrust = |controls: &[f64],
                           elems: &mut Vec<BladeElement<FoilFamily>>,
                           pg_r: &RefCell<Vec<f64>>,
                           hint_da: Option<f64>|
         -> (MeetOutcome, bool) {
            // returns (da, err, thrust, torque, alpha, chord, phi, reachable)
            let mut cur: Vec<f64> = pg_r.borrow().clone();
            let mut bst: Option<MeetOutcome> = None; // da, err, t, q, alpha, chord, phi
            let mut prev: Option<(f64, f64)> = None; // (da, thrust)
            let mut bracket: Option<(f64, f64)> = None;
            let better = |err: f64,
                          q: f64,
                          b: &Option<MeetOutcome>|
             -> bool {
                match b {
                    None => q > 0.0,
                    Some((_, be, _, bq, _, _, _)) => {
                        (q > 0.0)
                            && (err < *be - 1.0e-9 || ((err - *be).abs() < 1.0e-9 && q < *bq))
                    }
                }
            };

            // Warm start: evaluate the previous `da` first; if it already meets
            // the target, return immediately (a single circulation solve).
            if let Some(hd) = hint_da {
                if (DA_MIN..=DA_MAX).contains(&hd) {
                    let (t, q, alphas, chords, g, phis) = eval(controls, hd, elems, &cur);
                    cur = g;
                    let err = (t - thrust).abs() / thrust.max(1.0e-9);
                    if err <= 0.03 {
                        *pg_r.borrow_mut() = cur.clone();
                        println!(
                            "lift-line [{}] ctrl=[{}] da={:.3} (warm): T={:.4} Q={:.4}",
                            label,
                            controls.iter().map(|c| format!("{:.4}", c)).collect::<Vec<_>>().join(" "),
                            hd, t, q
                        );
                        return ((hd, err, t, q, alphas, chords, phis), true);
                    }
                    bst = Some((hd, err, t, q, alphas, chords, phis));
                    prev = Some((hd, t));
                }
            }

            // Full bounded scan (including the hint if not already scanned),
            // then bisection on the bracketing interval.
            let mut cands: Vec<f64> = vec![DA_MIN, 0.5 * (DA_MIN + DA_MAX), DA_MAX];
            if let Some(hd) = hint_da {
                if (DA_MIN..=DA_MAX).contains(&hd) && !cands.contains(&hd) {
                    cands.push(hd);
                }
            }
            cands.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for &da in &cands {
                let (t, q, alphas, chords, g, phis) = eval(controls, da, elems, &cur);
                cur = g;
                let err = (t - thrust).abs() / thrust.max(1.0e-9);
                if better(err, q, &bst) {
                    bst = Some((da, err, t, q, alphas, chords, phis));
                }
                if let Some((pd, pt)) = prev {
                    if (pt < thrust && t >= thrust) || (pt >= thrust && t < thrust) {
                        bracket = Some((pd, da));
                    }
                }
                prev = Some((da, t));
            }
            if let Some((mut lo, mut hi)) = bracket {
                for _ in 0..6 {
                    let mid = 0.5 * (lo + hi);
                    let (t, q, alphas, chords, g, phis) = eval(controls, mid, elems, &cur);
                    cur = g;
                    let err = (t - thrust).abs() / thrust.max(1.0e-9);
                    if better(err, q, &bst) {
                        bst = Some((mid, err, t, q, alphas, chords, phis));
                    }
                    if t < thrust {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
            }
            *pg_r.borrow_mut() = cur.clone();
            match bst {
                Some((da, err, t, q, alphas, chords, phis)) => {
                    ((da, err, t, q, alphas, chords, phis), err <= 0.05)
                }
                None => (
                    (0.0, 1.0, 0.0, 0.0, Vec::new(), Vec::new(), Vec::new()),
                    false,
                ),
            }
        };

        let obj = |x: &[f64]| -> f64 {
            let controls: Vec<f64> = (0..n_ctrl)
                .map(|k| CH_FLOOR * (cap_ctl[k] / CH_FLOOR).powf(x[k].clamp(0.0, 1.0)))
                .collect();
            let hint = *da_hint.borrow();
            let ((da, err, t, q, _al, _ch, _ph), _reach) = {
                let mut e = state.borrow_mut();
                meet_thrust(&controls, &mut e, &pg, hint)
            };
            *da_hint.borrow_mut() = Some(da);
            if hint.is_none() {
                println!(
                    "lift-line [{}] ctrl=[{}] da={:.3}: T={:.4} Q={:.4}",
                    label,
                    controls.iter().map(|c| format!("{:.4}", c)).collect::<Vec<_>>().join(" "),
                    da, t, q
                );
            }
            // Meet the target (50*err dominates) and then minimise torque.
            q + 50.0 * err
        };

        let nm = optimize::NelderMead {
            maxiter: 120,
            ..Default::default()
        };
        let bounds: Vec<Option<(f64, f64)>> = vec![Some((0.0, 1.0)); n_ctrl];
        // First start at the full chord (the thrust-capable geometry); each
        // later start warm-starts from the running best point so the simplex
        // keeps working the promising region (and the shared gamma/da warm
        // state — pg / da_hint — is already near the solution).
        let mut best_x: Vec<f64> = vec![1.0; n_ctrl];
        let mut best_f = f64::INFINITY;
        for i in 0..3 {
            let x0: Vec<f64> = if i == 0 {
                vec![1.0; n_ctrl]
            } else {
                let scale = if i == 1 { 0.85 } else { 0.7 };
                best_x.iter().map(|&v| (scale * v).clamp(0.0, 1.0)).collect()
            };
            let (x, f) = nm.minimize(obj, &x0, &bounds);
            if f < best_f {
                best_f = f;
                best_x = x;
            }
            if best_f < 0.03 {
                break;
            }
        }
        let _ = best_f;

        // Final measured candidate, seeded with the optimizer's converged
        // `da` so the reported geometry is exactly the design the outer loop
        // evaluated as its best.  A cold full scan (hint=None) re-solves the
        // candidates from a different Newton path and can land on another
        // branch of the nonlinear system, reporting a *worse* thrust match
        // than the optimizer actually achieved; the warm check still falls
        // back to the full scan if it does not meet the target.
        let mut elements = std::mem::take(&mut *state.borrow_mut());
        let controls: Vec<f64> = (0..n_ctrl)
            .map(|k| CH_FLOOR * (cap_ctl[k] / CH_FLOOR).powf(best_x[k].clamp(0.0, 1.0)))
            .collect();
        let ((da, err, r_t, r_q, _al, chords, phis), _reach) = {
            let hint = *da_hint.borrow();
            meet_thrust(&controls, &mut elements, &pg, hint)
        };
        if phis.len() != m {
            println!(
                "lift-line camber {}: no thrust match (T={:.4}, err {:.3})",
                label, r_t, err
            );
            return None;
        }
        for i in 0..m {
            elements[i].set_twist(
                phis[i] + (alpha_base[i] + da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX),
            );
        }
        let f = r_q + 50.0 * err;
        println!(
            "lift-line camber {}: T={:.4} Q={:.4} (obj {:.4})",
            label, r_t, r_q, f
        );
        Some(PassOutcome {
            f,
            q: r_q,
            t: r_t,
            label: label.to_string(),
            da,
            alpha_base: alpha_base.to_vec(),
            phis,
            chords,
            camber_dist: camber_dist.to_vec(),
        })
    }

    /// Design a blade with the coupled lifting line, targeting `thrust`.
    ///
    /// `ar` optionally forces a minimum blade aspect ratio
    /// `(R - hub) / mean_chord`, which thins the blade (lowers induced loss
    /// and raises efficiency).  Returns `(torque, thrust)` of the converged
    /// design.
    pub fn lift_line_design(&mut self, rpm: f64, thrust: f64, ar: Option<f64>) -> (f64, f64) {
        let u_0 = self.param.forward_airspeed;
        let omega = optimize::rpm2omega(rpm);
        let r_hub = self.param.hub_radius;
        let r_tip = self.param.radius;
        let n_blades = self.n_blades;
        let m = self.radial_steps.max(2);
        let rr: Vec<f64> = (0..m)
            .map(|i| r_hub + (r_tip - r_hub) * i as f64 / (m - 1) as f64)
            .collect();

        // --- Chord is a smooth (shape-preserving cubic / PCHIP) spline ---
        // through `chord_spline_n` control values at radii spread hub->tip.
        // Those control values are the design variables, so the optimum
        // *shape* (not just a global scale) is found.  `--ar` only caps the
        // upper bound each control may take.
        let shape: f64 = rr
            .iter()
            .map(|&r| chord_law(&self.param, n_blades, r, 0.0, 1.0))
            .sum::<f64>()
            / m as f64;
        let ar_shape = (r_tip - r_hub) / shape;
        const S_FLOOR: f64 = 0.2;
        let s_cap = ar.map(|a| (ar_shape / a).clamp(S_FLOOR, 1.0)).unwrap_or(1.0);
        let n_ctrl = self.param.chord_spline_n.max(2);
        let ref_r: Vec<f64> = (0..n_ctrl)
            .map(|k| r_hub + (r_tip - r_hub) * k as f64 / (n_ctrl.max(2) - 1) as f64)
            .collect();
        // The trailed-wake influence depends only on the station radii and the
        // blade count, so build it once and share it across the passes.
        let infl = lift_line::influence(&rr, n_blades);
        // Per-control upper bound: the geometrically-allowed chord (taper
        // capped by blade spacing) at each reference radius, thinned by `s_cap`
        // when a minimum aspect ratio was requested.
        let cap_ctl: Vec<f64> = ref_r
            .iter()
            .map(|&r| (s_cap * self.get_max_chord(r, 0.0)).max(CH_FLOOR))
            .collect();

        // --- Camber: one full design pass per candidate, best wins -------
        // An explicit `camber` parameter pins a single value; otherwise the
        // [`CAMBER_CANDIDATES`] set is scanned (see its docs for why a fixed
        // set rather than a continuous variable) plus a composed per-station
        // distribution: the best section L/D at each radius, smoothed.
        let cambers: Vec<f64> = match self.param.camber {
            Some(m) => vec![m],
            None => CAMBER_CANDIDATES.to_vec(),
        };

        // Candidate seeds as plain per-station data (the passes build their
        // own blade elements: `Rc`-based elements cannot cross threads).  The
        // station foils are collected for every candidate first so a single
        // worker pool can warm all the polar buckets the seeding queries
        // (the first-touch rust-foil sweeps dominate a cold run).
        let mut rows: Vec<Vec<(FoilFamily, f64)>> = Vec::with_capacity(cambers.len());
        let mut work: Vec<(FoilFamily, f64)> = Vec::new();
        let mut seen: HashSet<(String, u64)> = HashSet::new();
        for &m_c in &cambers {
            let row: Vec<(FoilFamily, f64)> = rr
                .iter()
                .map(|&r| {
                    let vv = ((u_0 * u_0) + (omega * r).powi(2)).sqrt();
                    let c = chord_law(&self.param, n_blades, r, 0.0, s_cap);
                    (station_foil(&self.param, r, c, m_c), vv)
                })
                .collect();
            for (f, vv) in &row {
                let key = (f.hash(), f.reynolds(*vv).to_bits());
                if seen.insert(key) {
                    work.push((f.clone(), *vv));
                }
            }
            rows.push(row);
        }
        self.warm_polar_pool(&work, "seeding warm-up");

        // Seed the best-L/D attack angles from the (now cached) polars.
        // Each seed is (label, per-station camber, smoothed attack angles).
        let mut seeds: Vec<(String, Vec<f64>, Vec<f64>)> = Vec::with_capacity(cambers.len() + 1);
        let mut raw_alphas: Vec<Vec<f64>> = Vec::with_capacity(cambers.len());
        let mut lds: Vec<Vec<f64>> = Vec::with_capacity(cambers.len());
        for (row, &m_c) in rows.iter().zip(cambers.iter()) {
            let mut alpha_raw: Vec<f64> = Vec::with_capacity(m);
            let mut ld_row: Vec<f64> = Vec::with_capacity(m);
            for (f, vv) in row {
                let mut fs = FoilSimulator::new(Rc::new(RefCell::new(f.clone())), self.store.clone());
                if self.plate_mode {
                    fs.set_plate_mode(true);
                }
                let (a, l) = Self::best_ld::<FoilFamily>(&fs, *vv);
                alpha_raw.push(a);
                ld_row.push(l);
            }
            // The raw best-L/D angles are jagged station-to-station: quantised to
            // whole degrees by the scan, and stepped wherever the Reynolds number
            // crosses a polar bucket or the low-Re flat-plate cutoff.  Re and Mach
            // vary smoothly with r, so the underlying optimum does too: use a
            // least-squares fit (as the BEM path does for the twist) everywhere
            // downstream (`eval`, `meet_thrust`, the final twist assembly).
            let alpha_base = smooth_alpha_curve(&rr, &alpha_raw);
            seeds.push((format!("m={:.2}", m_c), vec![m_c; m], alpha_base));
            raw_alphas.push(alpha_raw);
            lds.push(ld_row);
        }

        // The composed per-station distribution: best candidate at each
        // radius (thick root sections prefer little or no camber, thin
        // outboard sections benefit from it), smoothed and quantised to the
        // 0.01 polar-hash grid, competing against the uniform candidates on
        // the same objective.
        if self.param.camber.is_none() {
            let (m_dist, winners) = compose_camber(&rr, &lds, &cambers);
            let mut work: Vec<(FoilFamily, f64)> = Vec::new();
            for (i, &r) in rr.iter().enumerate() {
                let vv = ((u_0 * u_0) + (omega * r).powi(2)).sqrt();
                let c = chord_law(&self.param, n_blades, r, 0.0, s_cap);
                let f = station_foil(&self.param, r, c, m_dist[i]);
                let key = (f.hash(), f.reynolds(vv).to_bits());
                if seen.insert(key) {
                    work.push((f, vv));
                }
            }
            self.warm_polar_pool(&work, "composed-camber warm-up");
            // Each station keeps its winning candidate's attack angle.
            let alpha_raw: Vec<f64> = (0..m).map(|i| raw_alphas[winners[i]][i]).collect();
            let alpha_base = smooth_alpha_curve(&rr, &alpha_raw);
            println!(
                "lift-line composed camber m(r): [{}]",
                m_dist
                    .iter()
                    .map(|&m| format!("{:.2}", m))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            seeds.push(("m(r) per-station".into(), m_dist, alpha_base));
        }

        // --- One full design pass per candidate, in parallel --------------
        // The passes are independent (separate elements, attack angles and
        // warm state), so they run on scoped threads; the best objective
        // (`q + 50*err` — the same one the optimizer minimises) wins, so
        // camber competes with chord shape on equal terms.
        let param = &self.param;
        let store = self.store.clone();
        let radial_res = self.radial_resolution;
        let plate_mode = self.plate_mode;
        let mut outcomes: Vec<PassOutcome> = Vec::new();
        thread::scope(|s| {
            let handles: Vec<_> = seeds
                .into_iter()
                .map(|(label, camber_dist, alpha_base)| {
                    let (store, rr, ref_r, cap_ctl, infl) =
                        (&store, &rr, &ref_r, &cap_ctl, &infl);
                    s.spawn(move || {
                        Self::run_design_pass(
                            param, n_blades, plate_mode, store, radial_res, rr, ref_r, cap_ctl,
                            infl, rpm, omega, u_0, thrust, n_ctrl, s_cap, &label, &camber_dist,
                            &alpha_base,
                        )
                    })
                })
                .collect();
            for h in handles {
                if let Some(o) = h.join().expect("design pass panicked") {
                    outcomes.push(o);
                }
            }
        });

        let win = match outcomes
            .into_iter()
            .min_by(|a, b| a.f.partial_cmp(&b.f).unwrap())
        {
            Some(w) => w,
            // No candidate produced a usable design (every pass failed to
            // meet the thrust target); nothing to report or export.
            None => return (0.0, 0.0),
        };
        println!(
            "Lifting-line design: camber {} da={:.3} T={:.4} Q={:.4}",
            win.label, win.da, win.t, win.q
        );

        // Rebuild the winning elements on this thread from the pass's plain
        // geometry outcome (chord and twist per station).
        println!("Lifting-line blade stations");
        let mut elements: Vec<BladeElement<FoilFamily>> = Vec::with_capacity(m);
        for (i, &ri) in rr.iter().enumerate() {
            let alpha =
                (win.alpha_base[i] + win.da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX);
            let phi0 = u_0.atan2(omega * ri).max(1.0e-3);
            let c = chord_law(&self.param, n_blades, ri, 0.0, s_cap);
            let foil =
                Rc::new(RefCell::new(station_foil(&self.param, ri, c, win.camber_dist[i])));
            let mut be = BladeElement::new(
                ri,
                self.radial_resolution,
                foil,
                phi0 + 0.12,
                rpm,
                u_0,
                self.store.clone(),
            );
            if self.plate_mode {
                be.set_plate_mode(true);
            }
            be.set_chord(win.chords[i]);
            be.set_twist(win.phis[i] + alpha);
            println!(
                "r={} camber={} alpha_base={} alpha={} phi={} twist={} chord={} ",
                ri,
                win.camber_dist[i],
                win.alpha_base[i].to_degrees(),
                alpha.to_degrees(),
                win.phis[i].to_degrees(),
                (win.phis[i] + alpha).to_degrees(),
                win.chords[i]
            );
            elements.push(be);
        }
        self.blade_elements = elements;
        (win.q, win.t)
    }

}

/// Piecewise-linear interpolation (scipy `interp1d(x, y, "linear")`).
fn lin_interp(x: &[f64], y: &[f64], t: f64) -> f64 {
    if t <= x[0] {
        return y[0];
    }
    let n = x.len();
    if t >= x[n - 1] {
        return y[n - 1];
    }
    let mut i = 0;
    while x[i + 1] < t {
        i += 1;
    }
    let f = (t - x[i]) / (x[i + 1] - x[i]);
    y[i] + f * (y[i + 1] - y[i])
}

/// Least-squares smoothing (degree-4 polynomial) of a per-station angle
/// distribution: absorbs quantisation and polar-bucket steps while still
/// tracking the data (the BEM path applies the same idea to the twist).
fn smooth_alpha_curve(rr: &[f64], alpha: &[f64]) -> Vec<f64> {
    let deg = 4.min(rr.len().saturating_sub(1));
    let poly = polyfit(rr, alpha, deg);
    rr.iter().map(|&r| polyval(&poly, r)).collect()
}

/// Compose a per-station camber distribution from `lds[c][i]`, the section
/// L/D of candidate camber `c` at station `i`: each station takes its best
/// candidate's camber, the resulting step distribution is smoothed with a
/// low-order polynomial fit (like the attack angles), clamped to the
/// candidate range, and quantised to the 0.01 polar-hash grid — the built
/// foil then hashes to exactly the polars it was selected on.  Returns the
/// camber per station and each station's winning candidate index (for its
/// attack-angle seed).
fn compose_camber(rr: &[f64], lds: &[Vec<f64>], candidates: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let n = rr.len();
    let mut winners = vec![0usize; n];
    let mut raw = vec![candidates[0]; n];
    for i in 0..n {
        let mut best = f64::NEG_INFINITY;
        for (c, l) in lds.iter().enumerate() {
            if l[i] > best {
                best = l[i];
                winners[i] = c;
                raw[i] = candidates[c];
            }
        }
    }
    let deg = 3.min(n.saturating_sub(1));
    let poly = polyfit(rr, &raw, deg);
    let hi = candidates.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let m_dist = rr
        .iter()
        .map(|&r| {
            let v = polyval(&poly, r).clamp(0.0, hi);
            (100.0 * v).round() / 100.0
        })
        .collect();
    (m_dist, winners)
}

/// One camber candidate's converged design pass: the objective the
/// candidates compete on (`f = q + 50*err`), the measured operating point
/// and the per-station geometry.  Plain data so the passes can run on
/// worker threads; the caller rebuilds the winning blade elements from it.
struct PassOutcome {
    f: f64,
    q: f64,
    t: f64,
    label: String,
    da: f64,
    alpha_base: Vec<f64>,
    phis: Vec<f64>,
    chords: Vec<f64>,
    camber_dist: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_prop() -> Prop {
        let param = DesignParameters::default();
        let store = Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json")));
        Prop::new(param, 0.002, store)
    }

    #[test]
    fn max_chord_decreases_with_radius() {
        let p = test_prop();
        let c_root = p.get_max_chord(0.01, 0.0);
        let c_tip = p.get_max_chord(0.05, 0.0);
        assert!(c_root > c_tip, "{} vs {}", c_root, c_tip);
    }

    #[test]
    fn thickness_law_bounds() {
        let p = test_prop();
        let t_root = p.get_foil_thickness(p.param.hub_radius);
        let t_tip = p.get_foil_thickness(p.param.radius);
        // thickness_root = hub_depth, thickness_end = 0.1*hub_depth
        assert!((t_root - p.param.hub_depth).abs() < 1e-6, "root {}", t_root);
        assert!((t_tip - 0.1 * p.param.hub_depth).abs() < 1e-6, "tip {}", t_tip);
    }

    #[test]
    fn scimitar_offset_zero_at_ends() {
        let mut p = test_prop();
        let at_hub = p.get_scimitar_offset(p.param.hub_radius);
        let at_tip = p.get_scimitar_offset(p.param.radius);
        assert!(at_hub.abs() < 1e-9, "hub {}", at_hub);
        assert!(at_tip.abs() < 1e-9, "tip {}", at_tip);
    }

    #[test]
    fn tip_loss_is_one_at_mid_radius() {
        let p = test_prop();
        let l = p.tip_loss(0.03, 0.2);
        assert!(l > 0.9 && l <= 1.0, "tip loss {}", l);
    }

    #[test]
    fn smooth_alpha_curve_removes_staircase() {
        // A jagged best-L/D distribution: a low-Re flat-plate region at the
        // hub, then whole-degree quantised angles with neighbour flips.
        let rr: Vec<f64> = (0..40).map(|i| 0.006 + 0.069 * i as f64 / 39.0).collect();
        let raw: Vec<f64> = rr
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if i < 6 {
                    0.10 // low-Re flat-plate fallback region
                } else {
                    (5.0 + ((i * 7) % 3) as f64 * 0.5).to_radians()
                }
            })
            .collect();
        let s = smooth_alpha_curve(&rr, &raw);

        // Smooth: small second differences station-to-station.
        let mut max_d2 = 0.0_f64;
        for i in 1..s.len() - 1 {
            let d2 = (s[i + 1] - 2.0 * s[i] + s[i - 1]).abs();
            max_d2 = max_d2.max(d2);
        }
        assert!(max_d2 < 1.0e-4, "second differences too large: {}", max_d2);

        // Still tracks the data it was fit to.
        let dev: f64 = s
            .iter()
            .zip(raw.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / s.len() as f64;
        assert!(dev < 1.0_f64.to_radians(), "fit drifted from data: {} rad", dev);
    }

    #[test]
    fn compose_camber_smooths_and_quantises() {
        // A span where the thick root prefers no camber and the outboard
        // sections prefer 0.04: the composed distribution must rise
        // smoothly from ~0 at the root to ~0.04 outboard, quantised to the
        // 0.01 polar-hash grid and clamped to the candidate range.
        let n = 30;
        let rr: Vec<f64> = (0..n).map(|i| 0.006 + 0.062 * i as f64 / (n - 1) as f64).collect();
        let candidates = [0.0_f64, 0.02, 0.04];
        // lds[c][i]: the flat candidate wins the inner ~20% of the span,
        // the cambered one the rest (the middle candidate never wins).
        let lds: Vec<Vec<f64>> = [0, 1, 2]
            .iter()
            .map(|&c| {
                rr.iter()
                    .map(|&r| {
                        let x = (r - 0.006) / 0.062;
                        match c {
                            0 => 55.0 - 25.0 * x,
                            1 => 30.0,
                            _ => 45.0 + 25.0 * x,
                        }
                    })
                    .collect()
            })
            .collect();

        let (m_dist, winners) = compose_camber(&rr, &lds, &candidates);
        // Root region takes the flat candidate, outboard the cambered one.
        assert_eq!(winners[0], 0);
        assert_eq!(winners[n - 1], 2);
        // Quantised to the 0.01 grid, clamped to [0, 0.04].
        for &m in &m_dist {
            let on_grid = (m * 100.0).fract().abs() < 1.0e-9;
            assert!(on_grid && (0.0..=0.04).contains(&m), "m = {}", m);
        }
        // Monotone rise (small tolerance for fit wiggle at the ends).
        assert!(m_dist[n - 1] > m_dist[0] + 0.02, "no rise: {} -> {}", m_dist[0], m_dist[n - 1]);
    }

    #[test]
    fn blade_element_camber_is_applied() {
        // The camber candidate must reach the foil (it keys the polar
        // cache, so a wrong m silently designs a different blade).
        let mut p = test_prop();
        let be = p.new_blade_element_with_chord(0.02, 1000.0, 0.1, 0.008, 0.04);
        let f = be.foil.borrow();
        let m = match &*f {
            FoilFamily::Naca4(n) => n.m,
            FoilFamily::Cst(_) => panic!("expected a NACA4 foil (default family)"),
        };
        assert!((m - 0.04).abs() < 1.0e-12, "camber not applied");
    }

    #[test]
    fn blade_element_cst_family_is_applied() {
        // With the CST family selected, stations carry Cst foils re-sized to
        // the station laws: the camber candidate mapped onto the LEM weight,
        // and the thickness law (thickness/chord at this station) reaching
        // the foil exactly.
        let mut p = test_prop();
        p.param.cst = true;
        p.param.hub_depth = 0.006; // nonzero root depth: thickness law > 0
        let be = p.new_blade_element_with_chord(0.02, 1000.0, 0.1, 0.008, 0.04);
        let f = be.foil.borrow();
        match &*f {
            FoilFamily::Cst(c) => {
                assert!(
                    (c.params.leading_edge_weight).abs() > 0.1,
                    "LEM weight {}",
                    c.params.leading_edge_weight
                );
                let t_law = p.get_foil_thickness(0.02) / 0.008;
                assert!(
                    (c.thickness() - t_law).abs() < 1e-6,
                    "thickness {} vs law {}",
                    c.thickness(),
                    t_law
                );
            }
            FoilFamily::Naca4(_) => panic!("expected a CST foil with param.cst set"),
        }
    }

    #[test]
    fn lift_line_design_finite_and_thrust_matching() {
        // Plate polar (no rust-foil) so this is fast; checks the coupled
        // lifting-line produces finite thrust that converges to the target.
        let mut prop = test_prop();
        // coarser grid for speed
        prop.radial_steps = 12;
        prop.set_plate_mode(true);
        let (q, t) = prop.lift_line_design(12000.0, 2.0, Some(4.0));
        assert!(t.is_finite() && t > 0.0, "thrust {} not finite/positive", t);
        assert!(q.is_finite() && q > 0.0, "torque {} not finite/positive", q);
        // The target must actually be met (a u/v-swap in the force projection
        // once made the model thrust ~v/u times too small, so the design
        // could never converge onto the target).
        assert!(
            (t - 2.0).abs() / 2.0 < 0.05,
            "thrust {} does not match the 2.0 N target",
            t
        );
        assert!(q < 0.5, "torque {} not bounded", q);
    }
}
