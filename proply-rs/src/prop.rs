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
use crate::foil::{Arad, Cst, FoilFamily, FoilLike, Naca4};
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
/// may also trim *through* zero lift: torque is monotone in `da` across the
/// whole range, so the inner bisection can match a low torque target at any
/// chord — without it, a full-chord blade whose torque already exceeds the
/// target could never converge, as `da > 0` only adds torque.
const DA_MIN: f64 = -0.45;
const DA_MAX: f64 = 0.12;

/// Relative torque tolerance of the direct lifting-line operating-point
/// match ([`Prop::lift_line_design`]): above this the geometry could not
/// absorb the demanded torque and the result carries a warning.
const TOL_LIFT_LINE: f64 = 1.0e-2;

/// Worker threads available to the design's scoped-thread pools.
/// WebAssembly has no threads (`thread::scope`'s spawn panics there), and a
/// single-core host gains nothing from spawning: both fall back to running
/// the work lists inline.
fn worker_count() -> usize {
    if cfg!(target_arch = "wasm32") {
        1
    } else {
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    }
}

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
    let k = (thickness_root - thickness_end) / (param.hub_radius.powf(p) - param.radius.powf(p));
    let s = thickness_end - k * param.radius.powf(p);
    s + k * r.powf(p)
}

/// The foil of one station: the thickness law at radius `r` (the geometric
/// power law by default, or the mechanical [`crate::thickness`] law when
/// `law` is given — it carries thickness/chord directly, so the chord
/// scaling never distorts the sized section), chord `c`, camber fraction
/// `camber` and the parameter trailing edge.  The family is
/// [`DesignParameters::arad`] / [`DesignParameters::cst`]: the default NACA
/// 4-series, a CST (Kulfan) section — the default 18-parameter shape
/// re-thicknessed and cambered to the same laws — or the table-driven
/// ARA-D family, which carries its own camber (the `camber` argument does
/// not apply, as in the Python `ARADProp`).  Plain data (clonable,
/// sendable) so worker threads can build their own.
fn station_foil(param: &DesignParameters, law: Option<&Pchip>, r: f64, c: f64, camber: f64) -> FoilFamily {
    let t_frac = match law {
        Some(s) => s.eval(r).clamp(0.0, 1.0),
        None => foil_thickness(param, r) / c,
    };
    if param.arad {
        let mut f = Arad::new(c, t_frac);
        f.base.set_trailing_edge(param.trailing_edge / 1000.0);
        FoilFamily::Arad(f)
    } else if param.cst {
        let mut f = Cst::default(c);
        f.set_thickness(t_frac);
        f.set_camber(camber);
        f.set_trailing_edge(param.trailing_edge / 1000.0);
        FoilFamily::Cst(f)
    } else {
        let mut f = Naca4::new(c, t_frac, camber, 0.4);
        f.base.set_trailing_edge(param.trailing_edge / 1000.0);
        FoilFamily::Naca4(f)
    }
}

/// The outcome of a force evaluation over the blade: the totals summed over
/// the stations whose BEM solve converged, plus the station convergence
/// coverage — the design loop and the output summary report how much of the
/// blade actually reached a momentum equilibrium rather than silently
/// treating failed stations as zeros.
#[derive(Debug, Clone, Copy, Default)]
pub struct Forces {
    pub torque: f64,
    pub thrust: f64,
    pub converged: usize,
    pub total: usize,
}

/// The converged design: the absorbed torque and produced thrust at the
/// design RPM, the station convergence coverage, and — when the geometry
/// cannot absorb the design torque — an explicit warning describing the
/// closest design that was reached.
#[derive(Debug, Clone)]
pub struct DesignResult {
    pub torque: f64,
    pub thrust: f64,
    /// Set when the operating point could not be matched: explains the gap
    /// and the closest design returned.
    pub warning: Option<String>,
    pub converged_stations: usize,
    pub total_stations: usize,
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
    /// The mechanical thickness law (radius → thickness/chord) installed by
    /// [`Prop::size_mechanical_thickness`]; `None` keeps the geometric
    /// power law ([`foil_thickness`]).
    pub mech_thickness_law: Option<Pchip>,
    /// Predicted tip deflection of the sized mechanical law (m),
    /// reported in the design summary.
    pub mech_tip_deflection: Option<f64>,
    /// The winning design of the most recent lifting-line pass: the next
    /// torque-match iteration seeds its incumbent candidate from this, so
    /// consecutive matches stay on the same geometry branch (fewer, less
    /// oscillatory torque-match iterations).
    prev_win: Option<PassOutcome>,
}

impl Prop {
    pub fn new(param: DesignParameters, element_width: f64, store: Arc<Mutex<PolarStore>>) -> Self {
        let radial_steps = (param.radius / element_width) as usize;
        Self {
            param,
            radial_resolution: element_width,
            radial_steps,
            n_blades: 2,
            blade_elements: Vec::new(),
            store,
            scimitar_interpolator: None,
            max_depth_interpolator: None,
            plate_mode: false,
            mech_thickness_law: None,
            mech_tip_deflection: None,
            prev_win: None,
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
            self.mech_thickness_law.as_ref(),
            r,
            x_limit,
            self.param.camber.unwrap_or(0.0),
        )));

        let c_max = {
            let f = foil.borrow();
            f.get_max_chord(x_limit, y_limit, twist)
        };
        crate::dprintln!("Max Chord {}", c_max);
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
        let foil = Rc::new(RefCell::new(station_foil(
            &self.param,
            self.mech_thickness_law.as_ref(),
            r,
            c,
            camber,
        )));

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

    /// The *geometric* thickness law's foil thickness (m) as a function of
    /// radius: a p = 0.3 power law between the hub depth and a tenth of it
    /// at the tip.  See [`foil_thickness`] (the parameter form the design
    /// passes use).  This is the fallback law; when the mechanical law is
    /// active the design loop instead evaluates
    /// [`Prop::mech_thickness_law`] (a radius → thickness/chord curve).
    pub fn get_foil_thickness(&self, r: f64) -> f64 {
        foil_thickness(&self.param, r)
    }

    /// Size the blade thickness from the converged design's station loads
    /// ([`crate::thickness::size_mechanical_thickness`]) and install the
    /// resulting radius → thickness/chord curve as the design's thickness
    /// law.  The sizing uses each station's twist (the section's z-bending
    /// inertia mixes the z-projections of the thickness and of the chord,
    /// so a twisted section is stiffer against the thrust).  The hub
    /// thickness is deliberately not involved — the airfoil thickness
    /// follows the beam-deflection sizing, not the hub geometry.
    /// Returns `false` (and leaves the geometric law active) when the
    /// design has nothing to size: no stations, or no load on the blade.
    pub fn size_mechanical_thickness(&mut self) -> bool {
        if self.blade_elements.len() < 2 {
            return false;
        }
        let mut rows: Vec<(f64, f64, f64, crate::foil::SectionShape, f64)> = self
            .blade_elements
            .iter()
            .map(|be| {
                let f = be.foil.borrow();
                let c = f.chord();
                let t = be.thrust_n.unwrap_or_else(|| be.d_t());
                // The real section's bending inertia relative to its
                // enclosing rectangle (camber raises the flatwise factor —
                // a curved section is stiffer than a flat one).
                let shape = f.section_shape_factors(100);
                (be.r, c, be.get_twist(), shape, t)
            })
            .collect();
        // The elements are hub → tip, but sort defensively anyway.
        rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let rr: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let chords: Vec<f64> = rows.iter().map(|r| r.1).collect();
        let twist: Vec<f64> = rows.iter().map(|r| r.2).collect();
        let shape: Vec<crate::foil::SectionShape> = rows.iter().map(|r| r.3).collect();
        let thrust: Vec<f64> = rows.iter().map(|r| r.4).collect();
        let Some(law) = crate::thickness::size_mechanical_thickness(
            &rr,
            &chords,
            &twist,
            &shape,
            &thrust,
            self.n_blades,
            self.param.modulus,
            self.param.deflection_fraction,
            self.param.thickness_floor,
        ) else {
            crate::deprintln!(
                "proply: mechanical thickness could not be sized (no station loads) — keeping the geometric law"
            );
            return false;
        };
        crate::dprintln!(
            "Mechanical thickness: sized from the station loads — predicted tip deflection {:.3} mm (allowed {:.3} mm, E = {} Pa)",
            law.tip_deflection * 1000.0,
            law.deflection_limit * 1000.0,
            self.param.modulus
        );
        self.mech_thickness_law = Some(crate::thickness::law_interpolant(&law));
        self.mech_tip_deflection = Some(law.tip_deflection);
        true
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
        let f = (self.n_blades as f64 * (self.param.radius - r * 0.96)) / (2.0 * r * phi.sin());
        let tip_loss = 2.0 * (-f).exp().acos() / std::f64::consts::PI;
        let f = (self.n_blades as f64 * (r - self.param.hub_radius * 0.95)) / (2.0 * r * phi.sin());
        let hub_loss = 2.0 * (-f).exp().acos() / std::f64::consts::PI;
        tip_loss * hub_loss
    }

    /// Total torque and thrust at a given RPM, re-running the BEM on each
    /// element (mirrors `get_forces`).  Stations whose fixed-point solve
    /// fails are flagged (`BladeElement::converged` = false) and excluded
    /// from the totals, but keep their last induction state — previously
    /// they were silently reset to (dv=1, a_prime=0), which turned a
    /// struggling design into a zero-torque one without any explanation.
    pub fn get_forces(&mut self, rpm: f64) -> Forces {
        let mut torque = 0.0;
        let mut thrust = 0.0;
        let mut converged = 0usize;
        for be in self.blade_elements.iter_mut() {
            be.rpm = rpm;
            be.omega = optimize::rpm2omega(rpm);
            let dv_goal = be.dv;
            let (dv, a_prime, err) = be.bem(self.n_blades);
            be.converged = err < 0.01;

            if be.converged {
                converged += 1;
                let dt = be.d_t();
                let dm = be.d_m();
                thrust += dt;
                torque += dm;
                crate::dprintln!(
                    "r={}, theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={}, bem_err={}",
                    be.r,
                    be.get_twist().to_degrees(),
                    dv,
                    a_prime,
                    dt,
                    dm,
                    dt / dm,
                    err
                );
            } else {
                crate::deprintln!(
                    "r={}: BEM did not converge (bem_err {:.6}, dv_goal {}) dv={} a_prime={}",
                    be.r, err, dv_goal, dv, a_prime
                );
            }
        }
        Forces {
            torque,
            thrust,
            converged,
            total: self.blade_elements.len(),
        }
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
            // A degenerate station (zero chord extent makes both bounding-box
            // scales infinite/undefined) must bound the chord to zero, not
            // NaN: a NaN bound panics the optimizer's clamp.  This happens
            // with a zero hub_depth (zero-thickness blade law).
            let maxchord = if maxchord.is_finite() { maxchord } else { 0.0 };

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

        let mut c_points: Vec<f64> = vec![
            0.0,
            self.param.hub_radius / 2.0,
            0.9 * self.param.hub_radius,
        ];
        c_points.extend(radial_hub_to_tip.iter());
        let mut extra_chords: Vec<f64> = vec![
            0.9 * self.param.hub_depth,
            0.9 * self.param.hub_depth,
            0.9 * self.param.hub_depth,
        ];
        extra_chords.extend(chords.iter());
        let smoothed = smooth(&extra_chords, 11, "hanning");
        let chord_poly = Pchip::new(&c_points, &smoothed);

        crate::dprintln!("Smoothed Blade Form");
        for be in self.blade_elements.iter_mut() {
            let c = chord_poly.eval(be.r);
            let t = polyval(&twist_poly, be.r);
            be.set_chord(c);
            be.set_twist(t);
            crate::dprintln!("{}", be);
        }

        let f = self.get_forces(optimum_rpm);
        crate::dprintln!(
            "Total Thrust: {:5.2}, Torque: {:5.3} ({} of {} stations converged)",
            f.thrust, f.torque, f.converged, f.total
        );
        let _ = (total_thrust, total_torque);
        (f.torque, f.thrust)
    }

    /// Converge the design onto the operating point: at `rpm` the blade
    /// absorbs exactly `q_target` of torque (the motor's torque there), and
    /// the geometry maximises efficiency — equivalently thrust, since the
    /// shaft power Q·ω is fixed once (torque, RPM) is.
    ///
    /// The lifting-line path reaches the operating point *directly*: every
    /// design evaluation matches a common attack offset `da` so the
    /// absorbed torque equals `q_target` (torque is monotone in `da` over
    /// the attached range) and the chord/camber search then maximises the
    /// resulting thrust in a single pass.  No thrust target is iterated —
    /// the achieved thrust is an output, not a constraint.
    ///
    /// The BEM path sizes each station from momentum theory for a target
    /// induced velocity and has no equivalent single global offset, so it
    /// keeps the damped multiplicative iteration over the thrust target
    /// (torque rises monotonically with it); the iteration starts from a
    /// momentum-theory estimate of the thrust the demanded torque can
    /// produce (see [`Prop::bem_thrust_seed`]).
    ///
    /// If no geometry can absorb the demanded torque (chord/depth caps
    /// limit the loading), the closest design is returned with a warning.
    pub fn design_for_torque(&mut self, rpm: f64, q_target: f64, ar: Option<f64>) -> DesignResult {
        if self.param.lifting_line {
            let (q, t, err) = self.lift_line_design(rpm, q_target, ar);
            let (converged, total) = self.station_coverage();
            let warning = (err > TOL_LIFT_LINE).then(|| {
                let msg = format!(
                    "design torque {:.4} Nm at {:.0} rpm not achievable: the closest design absorbs {:.4} Nm ({:.1}% of the target) at {:.2} N thrust; {}/{} BEM stations converged",
                    q_target, rpm, q, 100.0 * err, t, converged, total
                );
                crate::dprintln!("proply: WARNING {}", msg);
                msg
            });
            return DesignResult {
                torque: q,
                thrust: t,
                warning,
                converged_stations: converged,
                total_stations: total,
            };
        }
        let seed = self.bem_thrust_seed(rpm, q_target);
        self.bem_design_for_torque(rpm, q_target, seed)
    }

    /// A starting thrust target for the BEM fixed-point iteration,
    /// estimated from the demanded torque by momentum theory: the
    /// iteration's absorbed torque is monotone in the target, so it
    /// converges from either side and the seed only needs the right
    /// magnitude.  In forward flight the thrust is bounded by P/u_0
    /// (with the ideal Froude efficiency at the disk ~P/(u_0+dv/2));
    /// hovering, the ideal disk gives T = (2 ρ A P²)^{1/3} and a real
    /// rotor achieves roughly half of that (figure of merit ~0.3-0.5), so
    /// the estimates sit within a factor ~2-3 of any sane answer.
    fn bem_thrust_seed(&self, rpm: f64, q_target: f64) -> f64 {
        let p = (q_target * optimize::rpm2omega(rpm)).abs();
        let u = self.param.forward_airspeed;
        let area = std::f64::consts::PI
            * (self.param.radius.powi(2) - self.param.hub_radius.powi(2)).max(1.0e-6);
        // Ideal-disk thrust for the shaft power P: T = (2 ρ A P²)^{1/3}.
        let hover = (2.0 * crate::optimize::RHO * area * p * p).cbrt();
        let seed = if u > 0.5 {
            // Forward flight: T ≤ P/u_0 (with η ≤ 1); blend the cruise
            // bound with the hover estimate so slow flight stays sane.
            let cruise = p / u;
            0.5 * (cruise * cruise * hover).cbrt()
        } else {
            0.5 * hover
        };
        seed.max(1.0e-3)
    }

    /// The BEM path of [`design_for_torque`]: iterate the thrust target
    /// until a full BEM design absorbs `q_target` at `rpm` (each iteration
    /// re-optimises every station for the induced velocity of that target).
    fn bem_design_for_torque(
        &mut self,
        rpm: f64,
        q_target: f64,
        thrust_seed: f64,
    ) -> DesignResult {
        const MAX_ITERS: usize = 30;
        /// Relative torque tolerance for a converged match.
        const TOL: f64 = 1.0e-2;
        /// Baseline damping exponent for the target update: ideally
        /// T ∝ Q^(2/3) in hover (T ∝ P^(2/3)) and ~Q^1 in cruise, so 0.8 is
        /// a sound starting point.  The exponent is then adapted to the
        /// locally measured response slope (see below).
        const DAMPING: f64 = 0.8;
        /// The measured thrust → torque response is steep with the
        /// mechanical thickness law (elasticity up to ~5, from the design
        /// logs) and slightly noisy, so an unbounded update overshoots and
        /// oscillates; cap every step at ±25% of the current target.
        const STEP_CAP: f64 = 0.25;

        let mut thrust = thrust_seed.max(1.0e-3);
        // (relative error, thrust target) of the closest design so far —
        // the fallback if the iteration cannot close on the target.
        let mut best: Option<(f64, f64)> = None;
        // (thrust, torque) of the previous sample: the local slope adapts
        // the damping exponent to the measured elasticity.  A fixed 0.8
        // exponent under-damped the mechanical-law phase (the log shows a
        // response elasticity of 2-5 and an 11-match oscillation; damping
        // above 1/elasticity is linearly unstable there).
        let mut prev: Option<(f64, f64)> = None;
        for iter in 0..MAX_ITERS {
            let (q, t) = self.full_optimize(rpm, thrust);
            let err = (q - q_target).abs() / q_target;
            crate::dprintln!(
                "Operating point match {:2}: thrust target {:6.3} N -> T={:6.3} N, Q={:6.4} Nm (design Q={:6.4}, err {:4.2}%)",
                iter + 1,
                thrust,
                t,
                q,
                q_target,
                100.0 * err
            );
            if err < TOL {
                let (converged, total) = self.station_coverage();
                return DesignResult {
                    torque: q,
                    thrust: t,
                    warning: None,
                    converged_stations: converged,
                    total_stations: total,
                };
            }
            if best.is_none_or(|(e, _)| err < e) {
                best = Some((err, thrust));
            }
            if !q.is_finite() || q <= 1.0e-9 {
                crate::deprintln!(
                    "proply: torque match stalled (absorbed Q={:.4}, design Q={:.4})",
                    q, q_target
                );
                break;
            }
            // Next target: the damped power law, with the exponent reduced
            // towards 1 / (measured elasticity) so the update stays inside
            // the linearly-stable range (elasticity from the last two
            // samples, when they are usable).
            let mut damp = DAMPING;
            if let Some((t0, q0)) = prev {
                let rel_dt = (thrust - t0) / t0;
                let rel_dq = (q - q0) / q0;
                if rel_dt.abs() > 1.0e-3 && rel_dq.abs() > 5.0e-3 {
                    let elasticity = rel_dq / rel_dt;
                    if (0.2..=6.0).contains(&elasticity) && elasticity.is_finite() {
                        damp = (0.9 / elasticity).clamp(0.25, DAMPING);
                    }
                }
            }
            let mut next = thrust * (q_target / q).powf(damp);
            let rel = (next - thrust) / thrust;
            if rel.abs() > STEP_CAP {
                next = thrust * (1.0 + STEP_CAP * rel.signum());
            }
            if (next - thrust).abs() / thrust < 1.0e-3 {
                crate::deprintln!(
                    "proply: torque match stalled at Q={:.4} Nm (design Q={:.4})",
                    q, q_target
                );
                break;
            }
            prev = Some((thrust, q));
            thrust = next;
        }
        // Not converged: re-run the closest design so the blade in place is
        // the one being reported, and hand back an explicit warning (a
        // silently "successful" but unmatched design is how the browser_demo
        // discrepancy surfaced: the tool reported a broken 10%-thrust design
        // as if nothing was wrong).
        if let Some((_, t_best)) = best {
            if (t_best - thrust).abs() > 1.0e-9 {
                crate::deprintln!(
                    "proply: falling back to the closest design (thrust target {:.3} N)",
                    t_best
                );
                return self.finish_unconverged(rpm, t_best, q_target);
            }
        }
        self.finish_unconverged(rpm, thrust, q_target)
    }

    /// Re-run the BEM design at `thrust` (the best match found by
    /// [`bem_design_for_torque`]) and report it with a warning describing
    /// the gap to the demanded operating point.
    fn finish_unconverged(&mut self, rpm: f64, thrust: f64, q_target: f64) -> DesignResult {
        let (q, t) = self.full_optimize(rpm, thrust);
        let (converged, total) = self.station_coverage();
        let err = (q - q_target).abs() / q_target;
        let warning = Some(format!(
            "design torque {:.4} Nm at {:.0} rpm not achievable: the closest design absorbs {:.4} Nm ({:.1}% of the target) at {:.2} N thrust; {}/{} BEM stations converged",
            q_target, rpm, q, 100.0 * err, t, converged, total
        ));
        crate::dprintln!("proply: WARNING {}", warning.as_ref().unwrap());
        DesignResult {
            torque: q,
            thrust: t,
            warning,
            converged_stations: converged,
            total_stations: total,
        }
    }

    /// How many of the current blade elements reached a converged BEM state:
    /// the YAML summary and the warning report the coverage.  Lifting-line
    /// elements never run `get_forces`, so their flag stays at the default
    /// (converged) — the count only carries meaning for the BEM loop.
    fn station_coverage(&self) -> (usize, usize) {
        let total = self.blade_elements.len();
        let converged = self.blade_elements.iter().filter(|be| be.converged).count();
        (converged, total)
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
    /// and the pool overlaps them.  The work list is flattened to one task
    /// per *distinct* `(foil, Reynolds bucket, Mach)` warm target —
    /// adjacent stations share bracketing buckets, so per-station tasks
    /// would pile every worker onto the same first bucket and serialize
    /// behind the per-key claim gate; distinct tasks keep every worker on a
    /// distinct simulation.  No-op in plate mode.  Progress is reported on
    /// a bar labelled `what` (hidden when stderr is not a TTY).
    fn warm_polar_pool(&self, work: &[(FoilFamily, f64)], what: &str) {
        if work.is_empty() || self.plate_mode {
            return;
        }
        let mut seen = HashSet::new();
        let mut tasks: Vec<(FoilFamily, f64, f64)> = Vec::new();
        for (foil, v) in work {
            let fs = FoilSimulator::new(Rc::new(RefCell::new(foil.clone())), self.store.clone());
            for (re, mach) in fs.warm_plan(*v) {
                if seen.insert((foil.hash(), re.to_bits(), mach.to_bits())) {
                    tasks.push((foil.clone(), re, mach));
                }
            }
        }
        crate::dprintln!("warm-up {}: {} polar tasks", what, tasks.len());
        let pb = ProgressBar::new(tasks.len() as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} foil polars ({msg}, eta {eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.set_message(what.to_string());
        let n_tasks = tasks.len();
        let queue = Mutex::new(tasks);
        let n_workers = worker_count().min(n_tasks);
        let store = self.store.clone();

        // One worker's share of the queue.  Pop through a `match` (not
        // `while let`): the guard of `while let ... = queue.lock().unwrap().pop()`
        // lives to the end of the loop body, so the first worker would hold
        // the queue lock for its whole simulation and serialize the entire
        // pool.
        fn worker(
            queue: &Mutex<Vec<(FoilFamily, f64, f64)>>,
            store: &Arc<Mutex<PolarStore>>,
            pb: &ProgressBar,
        ) {
            loop {
                let (foil, re, mach) = match queue.lock().unwrap().pop() {
                    Some(task) => task,
                    None => break,
                };
                let fs = FoilSimulator::new(Rc::new(RefCell::new(foil)), store.clone());
                fs.warm_bucket(re, mach);
                pb.inc(1);
            }
        }

        if n_workers > 1 {
            thread::scope(|s| {
                for _ in 0..n_workers {
                    let (queue, store, pb) = (&queue, store.clone(), pb.clone());
                    s.spawn(move || worker(queue, &store, &pb));
                }
            });
        } else {
            // No thread support (wasm) or a single core: drain inline.
            worker(&queue, &store, &pb);
        }
        pb.finish();
        crate::dprintln!("warm-up {}: done", what);
    }

    /// One camber candidate's design pass: build the station elements from
    /// `camber_dist`, prescribe the attack angles `alpha_base` (plus a common
    /// offset `da` matched so the blade absorbs exactly `q_target` of torque,
    /// by an inner monotone bisection) and optimise the chord-spline controls
    /// with warm-chained Nelder-Mead starts on the resulting thrust.  At a
    /// fixed (torque, RPM) operating point the shaft power Q·ω is fixed, so
    /// maximum thrust *is* maximum efficiency.  Self-contained (own elements
    /// and warm state, parameter-only geometry laws) so the camber candidates
    /// can run in parallel on scoped threads; returns the measured geometry,
    /// or `None` when no usable design was found.  `law` is the active
    /// thickness
    /// law (the mechanical radius → t/c curve, or `None` for the geometric
    /// power law); it is shared read-only, so the passes stay independent.
    #[allow(clippy::too_many_arguments)]
    fn run_design_pass(
        param: &DesignParameters,
        law: Option<&Pchip>,
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
        q_target: f64,
        n_ctrl: usize,
        s_cap: f64,
        label: &str,
        camber_dist: &[f64],
        alpha_base: &[f64],
        // Warm start for the first Nelder-Mead run: the previous torque
        // match's winning controls (`None` = start from the full chord).
        // The incumbent's pass re-optimizes from its previous optimum, so
        // consecutive matches' designs (and torque samples) stay close.
        x0_seed: Option<&[f64]>,
        pb: &ProgressBar,
    ) -> Option<PassOutcome> {
        let m = rr.len();

        // Station elements: the seed chord law at each radius, twisted to a
        // first guess (the circulation solve sets the final twist).
        let mut elements: Vec<BladeElement<FoilFamily>> = Vec::with_capacity(m);
        for (i, &r) in rr.iter().enumerate() {
            let phi0 = u_0.atan2(omega * r).max(1.0e-3);
            let c = chord_law(param, n_blades, r, 0.0, s_cap);
            let foil =
                Rc::new(RefCell::new(station_foil(param, law, r, c, camber_dist[i])));
            let mut be =
                BladeElement::new(r, radial_res, foil, phi0 + 0.12, rpm, u_0, store.clone());
            if plate_mode {
                be.set_plate_mode(true);
            }
            elements.push(be);
        }

        // Design variables `x ∈ [0,1]^n_ctrl` log-map to the control chords;
        // the common attack offset `da ∈ [DA_MIN, DA_MAX]` over the best-L/D
        // angles is matched so the absorbed torque equals `q_target` by an
        // inner monotone bisection (`meet_torque` below), not optimized.
        // `alpha_i =
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
        #[allow(clippy::type_complexity)]
        let eval = |controls: &[f64],
                    da: f64,
                    elems: &mut Vec<BladeElement<FoilFamily>>,
                    seed: &[f64]|
         -> (f64, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
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
            let res = lift_line::solve_with_influence(
                &stations, n_blades, omega, u_0, &fs_refs, seed, infl,
            );
            for i in 0..m {
                elems[i].set_twist(res.phi[i] + alphas[i]);
            }
            let chords: Vec<f64> = elems.iter().map(|be| be.foil.borrow().chord()).collect();
            (
                res.thrust,
                res.torque,
                alphas,
                chords,
                res.gamma.clone(),
                res.phi.clone(),
                res.u_i.clone(),
            )
        };

        // Outer optimizer over the N control values: for each candidate shape,
        // an inner **monotone `da` bisection** (below) matches the absorbed
        // torque to `q_target`, so the outer NelderMead only needs to
        // maximise the thrust among shapes that absorb exactly the target
        // torque.  Seeded at the full chord — the torque-capable geometry —
        // so it cannot be trapped in a thin/low-torque region.
        // (da, err, thrust, torque, alpha, chord, phi) — one design candidate.
        type MeetOutcome = (f64, f64, f64, f64, Vec<f64>, Vec<f64>, Vec<f64>);
        // A converged solve is usable only while every station carries
        // non-negative circulation: the section models (and their polar
        // fits) can support spurious flow branches with an outboard
        // *brake* — negative gamma contributing negative torque and thrust.
        // Such states must never seed the running best: with the old
        // thrust-matched objective a brake dumps torque, so the garbage
        // branch *won* the competition (low Q at the matched T) — the
        // torque cliff between consecutive matches.  At a torque target a
        // brake wastes torque budget, but only if it is never adopted here.
        let physical = |q: f64, g: &[f64], ui: &[f64]| -> bool {
            q.is_finite()
                && q > 0.0
                && g.iter().all(|&x| x.is_finite() && x >= -1.0e-9)
                && circulation_smooth(g)
                && flow_ui_sane(ui, rr, u_0, omega)
        };
        let meet_torque = |controls: &[f64],
                           elems: &mut Vec<BladeElement<FoilFamily>>,
                           pg_r: &RefCell<Vec<f64>>,
                           hint_da: Option<f64>|
         -> (MeetOutcome, bool) {
            // returns (da, err, thrust, torque, alpha, chord, phi, reachable)
            let mut cur: Vec<f64> = pg_r.borrow().clone();
            let mut bst: Option<MeetOutcome> = None; // da, err, t, q, alpha, chord, phi
            let mut prev: Option<(f64, f64)> = None; // (da, torque)
            let mut bracket: Option<(f64, f64)> = None;
            // Prefer the sample closest to the target torque; at equal error
            // prefer the larger thrust (the efficiency objective).
            let better = |err: f64, t: f64, b: &Option<MeetOutcome>| -> bool {
                match b {
                    None => true,
                    Some((_, be, bt, _, _, _, _)) => {
                        err < *be - 1.0e-9 || ((err - *be).abs() < 1.0e-9 && t > *bt)
                    }
                }
            };

            // Warm start: evaluate the previous `da` first; if it already
            // meets the target, return immediately (a single circulation
            // solve).  The acceptance requires a physical result: a
            // warm-started Newton solve can land on a brake-tip branch, and
            // the objective must not reward that state.
            if let Some(hd) = hint_da {
                if (DA_MIN..=DA_MAX).contains(&hd) {
                    let (t, q, alphas, chords, g, phis, u_i) = eval(controls, hd, elems, &cur);
                    let err = (q - q_target).abs() / q_target.max(1.0e-9);
                    cur = g;
                    if err <= 0.03 && physical(q, &cur, &u_i) {
                        *pg_r.borrow_mut() = cur.clone();
                        crate::dprintln!(
                            "lift-line [{}] ctrl=[{}] da={:.3} (warm): T={:.4} Q={:.4}",
                            label,
                            controls
                                .iter()
                                .map(|c| format!("{:.4}", c))
                                .collect::<Vec<_>>()
                                .join(" "),
                            hd,
                            t,
                            q
                        );
                        return ((hd, err, t, q, alphas, chords, phis), true);
                    }
                    // The sample may still bracket the target (its torque is
                    // meaningful), but a non-physical state must not seed the
                    // running best: it would block the honest samples.
                    if physical(q, &cur, &u_i) {
                        bst = Some((hd, err, t, q, alphas, chords, phis));
                    }
                    prev = Some((hd, q));
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
                let (t, q, alphas, chords, g, phis, u_i) = eval(controls, da, elems, &cur);
                cur = g;
                let err = (q - q_target).abs() / q_target.max(1.0e-9);
                if physical(q, &cur, &u_i) && better(err, t, &bst) {
                    bst = Some((da, err, t, q, alphas, chords, phis));
                }
                if let Some((pd, pq)) = prev {
                    if (pq < q_target && q >= q_target) || (pq >= q_target && q < q_target) {
                        bracket = Some((pd, da));
                    }
                }
                prev = Some((da, q));
            }
            if let Some((mut lo, mut hi)) = bracket {
                for _ in 0..6 {
                    let mid = 0.5 * (lo + hi);
                    let (t, q, alphas, chords, g, phis, u_i) = eval(controls, mid, elems, &cur);
                    cur = g;
                    let err = (q - q_target).abs() / q_target.max(1.0e-9);
                    if physical(q, &cur, &u_i) && better(err, t, &bst) {
                        bst = Some((mid, err, t, q, alphas, chords, phis));
                    }
                    if q < q_target {
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
                meet_torque(&controls, &mut e, &pg, hint)
            };
            *da_hint.borrow_mut() = Some(da);
            pb.inc(1);
            pb.set_message(format!("{}: T={:.3} Q={:.4}", label, t, q));
            if hint.is_none() {
                crate::dprintln!(
                    "lift-line [{}] ctrl=[{}] da={:.3}: T={:.4} Q={:.4}",
                    label,
                    controls
                        .iter()
                        .map(|c| format!("{:.4}", c))
                        .collect::<Vec<_>>()
                        .join(" "),
                    da,
                    t,
                    q
                );
            }
            // Absorb the target torque exactly (1000*err dominates: a design
            // that cannot absorb the demanded torque is not the operating
            // point) and then maximise the resulting thrust — at the fixed
            // shaft power Q·ω, thrust ⇔ efficiency.
            1000.0 * err - t
        };

        let nm = optimize::NelderMead {
            maxiter: 120,
            ..Default::default()
        };
        let bounds: Vec<Option<(f64, f64)>> = vec![Some((0.0, 1.0)); n_ctrl];
        // First run: warm start from the previous operating point's controls
        // when this is the incumbent design (the torque target is the same
        // and only the thickness law moved between the pipeline's two runs),
        // else the full chord — the torque-capable geometry.  When a seed
        // was used, a second run from the full chord anchors reachability,
        // so the search cannot be trapped in a thin region.
        let full_chord: Vec<f64> = vec![1.0; n_ctrl];
        let seeded: Vec<f64> = x0_seed
            .filter(|s| s.len() == n_ctrl)
            .map(|s| s.to_vec())
            .unwrap_or_else(|| full_chord.clone());
        let mut runs: Vec<Vec<f64>> = vec![seeded];
        if runs[0]
            .iter()
            .zip(&full_chord)
            .any(|(a, b)| (a - b).abs() > 1.0e-12)
        {
            runs.push(full_chord);
        }
        let mut best_x: Vec<f64> = vec![1.0; n_ctrl];
        let mut best_f = f64::INFINITY;
        for x0 in runs {
            let (x, f) = nm.minimize(obj, &x0, &bounds);
            if f < best_f {
                best_f = f;
                best_x = x;
            }
        }
        // Warm-chained scaled restarts around the running best (they keep
        // working the promising region, and the shared gamma/da warm state
        // — pg / da_hint — is already near the solution).  A restart that
        // cannot beat the running best stops the chain: the simplex is
        // already in the solution's basin, and a scaled copy would only
        // re-grind it (the log showed ~half of all evaluations repeating
        // the previous one).
        for scale in [0.85, 0.7] {
            let x0: Vec<f64> = best_x
                .iter()
                .map(|&v| (scale * v).clamp(0.0, 1.0))
                .collect();
            let before = best_f;
            let (x, f) = nm.minimize(obj, &x0, &bounds);
            if f < best_f {
                best_f = f;
                best_x = x;
            }
            if f >= before {
                break;
            }
        }
        let _ = best_f;

        // Final measured candidate, seeded with the optimizer's converged
        // `da` so the reported geometry is exactly the design the operating
        // point evaluated as its best.  A cold full scan (hint=None) re-solves
        // the candidates from a different Newton path and can land on another
        // branch of the nonlinear system, reporting a *worse* torque match
        // than the optimizer actually achieved; the warm check still falls
        // back to the full scan if it does not meet the target.
        let mut elements = std::mem::take(&mut *state.borrow_mut());
        let controls: Vec<f64> = (0..n_ctrl)
            .map(|k| CH_FLOOR * (cap_ctl[k] / CH_FLOOR).powf(best_x[k].clamp(0.0, 1.0)))
            .collect();
        let ((mut da, mut err, mut r_t, mut r_q, _al, mut chords, mut phis), _reach) = {
            let hint = *da_hint.borrow();
            meet_torque(&controls, &mut elements, &pg, hint)
        };
        if phis.len() != m {
            crate::dprintln!(
                "lift-line camber {}: no usable design (T={:.4}, err {:.3})",
                label, r_t, err
            );
            return None;
        }
        // The warm hint's acceptance (err <= 3%) leaves the measured
        // operating point up to a few percent off the target torque.
        // Candidates must compete on the thrust at the *exact* target
        // torque, so refine `da` with a bounded bisection around the warm
        // value: the absorbed torque is monotone in `da`, and every solve
        // warm-starts from the previous circulation, so the Newton path
        // stays on the branch the optimizer converged on.
        if err > 1.0e-3 && r_q > 0.0 {
            let mut cur: Vec<f64> = pg.borrow().clone();
            let (mut lo, mut hi) = (da, da);
            let (mut q_lo, mut q_hi) = (r_q, r_q);
            // Expand a bracket around the warm da until the absorbed torque
            // straddles the target.
            for _ in 0..7 {
                if (q_lo - q_target) * (q_hi - q_target) <= 0.0 && lo < hi {
                    break;
                }
                if q_hi < q_target && hi < DA_MAX {
                    lo = hi;
                    q_lo = q_hi;
                    hi = (hi + 0.1).min(DA_MAX);
                    let (_t2, q2, _a, _c, g2, _p, _u2) = eval(&controls, hi, &mut elements, &cur);
                    cur = g2;
                    q_hi = q2;
                } else if q_lo > q_target && lo > DA_MIN {
                    hi = lo;
                    q_hi = q_lo;
                    lo = (lo - 0.1).max(DA_MIN);
                    let (_t2, q2, _a, _c, g2, _p, _u2) = eval(&controls, lo, &mut elements, &cur);
                    cur = g2;
                    q_lo = q2;
                } else {
                    break;
                }
            }
            if (q_lo - q_target) * (q_hi - q_target) <= 0.0 {
                let (mut l, mut h) = (lo, hi);
                for _ in 0..14 {
                    let mid = 0.5 * (l + h);
                    let (t2, q2, _a, c2, g2, p2, u2) = eval(&controls, mid, &mut elements, &cur);
                    cur = g2;
                    let e2 = (q2 - q_target).abs() / q_target.max(1.0e-9);
                    // Adopt the tighter sample only when it is physical
                    // (positive torque, every station lifting, attached
                    // flow): a warm-started Newton solve at some bisection
                    // point can land on a brake-tip or unattached branch
                    // whose total torque still matches, and the objective
                    // must not reward that state.
                    if e2 < err && physical(q2, &cur, &u2) {
                        err = e2;
                        da = mid;
                        r_t = t2;
                        r_q = q2;
                        chords = c2;
                        phis = p2;
                    }
                    if q2 <= q_target {
                        l = mid;
                    } else {
                        h = mid;
                    }
                }
            }
        }
        for i in 0..m {
            elements[i].set_twist(
                phis[i] + (alpha_base[i] + da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX),
            );
        }
        let f = 1000.0 * err - r_t;
        crate::dprintln!(
            "lift-line camber {}: T={:.4} Q={:.4} (err {:.2}%, obj {:.4})",
            label, r_t, r_q, 100.0 * err, f
        );
        // The adopted state's own circulation: re-solve it on the warm
        // branch (the last bisection sample moved `pg` past it) so the
        // outcome carries the design branch's seed for the exported
        // blade's verification.
        let gamma = {
            let cur = pg.borrow().clone();
            let (_t3, _q3, _a3, _c3, g3, _p3, _u3) = eval(&controls, da, &mut elements, &cur);
            g3
        };
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
            gamma,
            design_x: best_x,
        })
    }

    /// Design a blade with the coupled lifting line that absorbs exactly
    /// `q_target` of torque at `rpm` and maximises the resulting thrust
    /// (i.e. the efficiency: the shaft power Q·ω is fixed by the operating
    /// point).  One design pass — no iteration over a thrust target.
    ///
    /// `ar` optionally forces a minimum blade aspect ratio
    /// `(R - hub) / mean_chord`, which thins the blade (lowers induced loss
    /// and raises efficiency).  Returns `(torque, thrust, relative torque
    /// error)` of the converged design: the error is small when the
    /// geometry can absorb the demanded torque, and larger (with the
    /// closest design returned) when chord/depth caps limit the loading.
    pub fn lift_line_design(&mut self, rpm: f64, q_target: f64, ar: Option<f64>) -> (f64, f64, f64) {
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
        let s_cap = ar
            .map(|a| (ar_shape / a).clamp(S_FLOOR, 1.0))
            .unwrap_or(1.0);
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
                    (station_foil(&self.param, self.mech_thickness_law.as_ref(), r, c, m_c), vv)
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
        // The per-station best-L/D scans mostly hit the warmed polars, but
        // on a cold cache they simulate their own (serially) and the phase
        // runs minutes with no other output: show a bar (hidden when
        // stderr is not a TTY).
        let mut seeds: Vec<PassSeed> =
            Vec::with_capacity(cambers.len() + 2);
        let mut raw_alphas: Vec<Vec<f64>> = Vec::with_capacity(cambers.len());
        let mut lds: Vec<Vec<f64>> = Vec::with_capacity(cambers.len());
        let seed_pb = ProgressBar::new((rows.len() * m) as u64);
        seed_pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] {pos}/{len} station seeds ({msg}, eta {eta})",
            )
            .unwrap(),
        );
        for (row, &m_c) in rows.iter().zip(cambers.iter()) {
            seed_pb.set_message(format!("m={:.2}", m_c));
            let mut alpha_raw: Vec<f64> = Vec::with_capacity(m);
            let mut ld_row: Vec<f64> = Vec::with_capacity(m);
            for (f, vv) in row {
                let mut fs =
                    FoilSimulator::new(Rc::new(RefCell::new(f.clone())), self.store.clone());
                if self.plate_mode {
                    fs.set_plate_mode(true);
                }
                let (a, l) = Self::best_ld::<FoilFamily>(&fs, *vv);
                alpha_raw.push(a);
                ld_row.push(l);
                seed_pb.inc(1);
            }
            // The raw best-L/D angles are jagged station-to-station: quantised to
            // whole degrees by the scan, and stepped wherever the Reynolds number
            // crosses a polar bucket or the low-Re flat-plate cutoff.  Re and Mach
            // vary smoothly with r, so the underlying optimum does too: use a
            // least-squares fit (as the BEM path does for the twist) everywhere
            // downstream (`eval`, `meet_torque`, the final twist assembly).  The
            // smoothed angles are floored at zero: a *negative* best-L/D (or a
            // polynomial edge overshoot below it) would prescribe a negative-lift
            // section — a tip brake — and the polar fits do produce spurious
            // negative-alpha optima on thin low-Re outboard sections.  A real
            // section's best L/D sits at positive lift; the floor keeps the
            // design out of the negative-lift flow branch entirely.
            let alpha_base: Vec<f64> =
                smooth_alpha_curve(&rr, &alpha_raw).iter().map(|&a| a.max(0.0)).collect();
            seeds.push((format!("m={:.2}", m_c), vec![m_c; m], alpha_base, None));
            raw_alphas.push(alpha_raw);
            lds.push(ld_row);
        }
        seed_pb.finish();

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
                let f = station_foil(&self.param, self.mech_thickness_law.as_ref(), r, c, m_dist[i]);
                let key = (f.hash(), f.reynolds(vv).to_bits());
                if seen.insert(key) {
                    work.push((f, vv));
                }
            }
            self.warm_polar_pool(&work, "composed-camber warm-up");
            // Each station keeps its winning candidate's attack angle.
            let alpha_raw: Vec<f64> = (0..m).map(|i| raw_alphas[winners[i]][i]).collect();
            let alpha_base: Vec<f64> =
                smooth_alpha_curve(&rr, &alpha_raw).iter().map(|&a| a.max(0.0)).collect();
            crate::dprintln!(
                "lift-line composed camber m(r): [{}]",
                m_dist
                    .iter()
                    .map(|&m| format!("{:.2}", m))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            seeds.push(("m(r) per-station".into(), m_dist, alpha_base, None));
        }

        // The previous operating point's winning design re-competes as the
        // incumbent, warm-started from its own controls: the pipeline's
        // second run (the mechanical-thickness law sized on the first
        // run's loads) then starts from the first run's geometry instead of
        // the full chord — the two designs stay on one branch.
        if let Some(prev) = self.prev_win.as_ref() {
            if prev.camber_dist.len() == m
                && prev.alpha_base.len() == m
                && prev.design_x.len() == n_ctrl
            {
                let alpha_base: Vec<f64> = prev.alpha_base.iter().map(|&a| a.max(0.0)).collect();
                seeds.insert(
                    0,
                    (
                        "prev".to_string(),
                        prev.camber_dist.clone(),
                        alpha_base,
                        Some(prev.design_x.clone()),
                    ),
                );
            }
        }

        // --- One full design pass per candidate, in parallel --------------
        // The passes are independent (separate elements, attack angles and
        // warm state), so they run on scoped threads; the best objective
        // (`1000*err - thrust` — the same one the optimizer minimises) wins,
        // so camber competes with chord shape on equal terms.
        let param = &self.param;
        let law = self.mech_thickness_law.as_ref();
        let store = self.store.clone();
        let radial_res = self.radial_resolution;
        let plate_mode = self.plate_mode;
        let mut outcomes: Vec<PassOutcome> = Vec::new();
        // Each pass's Nelder-Mead evaluations carry full circulation solves
        // and, on a cold cache, fresh polar simulations (the moved chords
        // hit new Reynolds buckets), so minutes can pass between printed
        // improvements: one shared evaluation counter across the parallel
        // passes keeps the phase visibly moving (hidden when stderr is not
        // a TTY).
        let eval_pb = ProgressBar::new_spinner();
        eval_pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] {pos} design evaluations ({msg})",
            )
            .unwrap(),
        );
        eval_pb.set_message("camber candidate passes");
        if worker_count() > 1 {
            thread::scope(|s| {
                let handles: Vec<_> = seeds
                    .into_iter()
                    .map(|(label, camber_dist, alpha_base, x0_seed)| {
                        let (store, rr, ref_r, cap_ctl, infl) =
                            (&store, &rr, &ref_r, &cap_ctl, &infl);
                        let eval_pb = &eval_pb;
                        s.spawn(move || {
                            Self::run_design_pass(
                                param,
                                law,
                                n_blades,
                                plate_mode,
                                store,
                                radial_res,
                                rr,
                                ref_r,
                                cap_ctl,
                                infl,
                                rpm,
                                omega,
                                u_0,
                                q_target,
                                n_ctrl,
                                s_cap,
                                &label,
                                &camber_dist,
                                &alpha_base,
                                x0_seed.as_deref(),
                                eval_pb,
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
        } else {
            // No thread support (wasm) or a single core: the candidate
            // passes run inline, in seed order.  Identical results — the
            // passes are independent by construction.
            for (label, camber_dist, alpha_base, x0_seed) in seeds {
                if let Some(o) = Self::run_design_pass(
                    param,
                    law,
                    n_blades,
                    plate_mode,
                    &store,
                    radial_res,
                    &rr,
                    &ref_r,
                    &cap_ctl,
                    &infl,
                    rpm,
                    omega,
                    u_0,
                    q_target,
                    n_ctrl,
                    s_cap,
                    &label,
                    &camber_dist,
                    &alpha_base,
                    x0_seed.as_deref(),
                    &eval_pb,
                ) {
                    outcomes.push(o);
                }
            }
        }
        eval_pb.finish();

        let mut win = match outcomes
            .into_iter()
            .min_by(|a, b| a.f.partial_cmp(&b.f).unwrap())
        {
            Some(w) => w,
            // No candidate produced a usable design (every pass failed to
            // match the torque target); nothing to report or export.
            None => return (0.0, 0.0, 1.0),
        };
        crate::dprintln!(
            "Lifting-line design: camber {} da={:.3} T={:.4} Q={:.4}",
            win.label, win.da, win.t, win.q
        );
        // The exported camber is a solved smooth spline too: the composed
        // m(r) candidate is quantised to the 0.01 polar-hash grid, leaving
        // small steps between adjacent stations.  Fit a low-order
        // polynomial over the span (least squares — the solved spline
        // parameters) and clamp it to the data's own range so the fit
        // cannot overshoot at the ends.
        let camber_smooth: Vec<f64> = {
            let lo = win.camber_dist.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = win
                .camber_dist
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max);
            let deg = 3.min(m.saturating_sub(1));
            let fit = polyfit(&rr, &win.camber_dist, deg);
            (0..m)
                .map(|i| polyval(&fit, rr[i]).clamp(lo, hi))
                .collect()
        };

        // Rebuild the winning elements on this thread from the pass's plain
        // geometry outcome (chord and twist per station).
        crate::dprintln!("Lifting-line blade stations");
        let mut elements: Vec<BladeElement<FoilFamily>> = Vec::with_capacity(m);
        for (i, &ri) in rr.iter().enumerate() {
            let phi0 = u_0.atan2(omega * ri).max(1.0e-3);
            let c = chord_law(&self.param, n_blades, ri, 0.0, s_cap);
            let foil = Rc::new(RefCell::new(station_foil(
                &self.param,
                self.mech_thickness_law.as_ref(),
                ri,
                c,
                camber_smooth[i],
            )));
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
            elements.push(be);
        }
        // The exported blade is a smooth realisation of the design: the
        // chord's geometric cap can kink where the optimizer's spline
        // meets it, so the exported chord re-fits the winning chord over
        // the span and the cold re-match below settles the operating point
        // on this smoothed blade (the analysis stations then use these
        // chords).
        let chord_shape: Vec<f64>;
        {
            // The exported chord is a solved smooth spline: least-squares
            // polynomial parameters over the span (degree 5), clamped to
            // the data's range so the fit cannot overshoot at the ends.
            // This rounds the geometric-cap knee that the design spline
            // otherwise kinks on, and smooths the unloaded hub station
            // (the optimizer pinches it).  Its scale is left free: the
            // operating-point solve below picks it so the smoothed blade
            // absorbs exactly the target torque.
            let raw: Vec<f64> = elements
                .iter()
                .map(|be| be.foil.borrow().chord())
                .collect();
            let lo = raw.iter().cloned().fold(f64::INFINITY, f64::min);
            let hi = raw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let deg = 5.min(m.saturating_sub(1));
            let fit = polyfit(&rr, &raw, deg);
            chord_shape = (0..m).map(|i| polyval(&fit, rr[i]).clamp(lo, hi)).collect();
        }
        // The exported blade's own operating point.  An independent cold
        // solve (starting from zero circulation, not the optimizer's warm
        // state) can settle on a neighbouring branch of the nonlinear flow
        // and absorb a few percent more or less torque than the warm state
        // the pass measured — and the exported blade is what the summary
        // reports.  Re-match `da` with cold starts so the cold-absorbed
        // torque hits the target, and export/report *that* state: the
        // design point is then exactly what an independent solve of the
        // exported geometry measures.  The same cold solve provides the
        // per-station flow diagnostics and the station loads the mechanical
        // law sizes on.
        let eval_at =
            |da: f64,
             elems: &mut Vec<BladeElement<FoilFamily>>,
             seed: &[f64]|
             -> (lift_line::LiftLineResult, Vec<f64>) {
                let alphas: Vec<f64> = win
                    .alpha_base
                    .iter()
                    .map(|&ab| (ab + da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX))
                    .collect();
                let stations: Vec<lift_line::Station> = (0..m)
                    .map(|i| lift_line::Station {
                        r: rr[i],
                        c: elems[i].foil.borrow().chord(),
                        alpha: alphas[i],
                    })
                    .collect();
                let fs_refs: Vec<&FoilSimulator<FoilFamily>> =
                    elems.iter().map(|be| &be.fs).collect();
                let res = lift_line::solve_with_influence(
                    &stations,
                    self.n_blades,
                    omega,
                    u_0,
                    &fs_refs,
                    seed,
                    &infl,
                );
                (res, alphas)
            };
        // The exported flow must also be *attached*: outside the inner hub
        // zone the inflow angle and the induced axial velocity cannot blow
        // up (see [`flow_phi_sane`], [`flow_ui_sane`]) — the wake model's
        // small-perturbation assumption is gone past those bounds and the
        // solved state (huge phi, wild exported twists) is not a physical
        // loading.
        let physical = |r: &lift_line::LiftLineResult| -> bool {
            r.torque.is_finite()
                && r.torque > 0.0
                && r.gamma.iter().all(|&x| x.is_finite() && x >= -1.0e-9)
                && circulation_smooth(&r.gamma)
                && flow_ui_sane(&r.u_i, &rr, u_0, omega)
        };
        let target = q_target.max(1.0e-9);
        // Solve the exported spline's scale so the smoothed blade absorbs
        // exactly the target torque: the scale is the spline's one free
        // parameter and Q scales almost linearly with chord, so a short
        // multiplicative iteration on the measured Q (each step a full
        // cold da re-match) converges in a few rounds.  The reported
        // operating point is then the smoothed, scaled blade's own state.
        let set_chords = |elems: &mut Vec<BladeElement<FoilFamily>>, scale: f64| {
            for (i, be) in elems.iter_mut().enumerate() {
                be.set_chord(chord_shape[i] * scale);
            }
        };
        let rematch =
            |elems: &mut Vec<BladeElement<FoilFamily>>, seed: &[f64]| -> (lift_line::LiftLineResult, Vec<f64>, f64, f64, bool) {
            let mut da = win.da;
            let (mut res, mut alphas) = eval_at(da, elems, seed);
            let mut err = (res.torque - q_target).abs() / target;
            let mut exported = physical(&res);
            // A physical solve whose circulation seeds continuation retries
            // for samples the primary seed lands off-branch on (the solve
            // can flip to a spurious branch over a small da step).
            let mut cont: Option<Vec<f64>> = if exported {
                Some(res.gamma.clone())
            } else {
                None
            };
            let mut sample =
                |sda: f64, elems: &mut Vec<BladeElement<FoilFamily>>| -> (lift_line::LiftLineResult, Vec<f64>, bool) {
                    let (r2, a2) = eval_at(sda, elems, seed);
                    if physical(&r2) {
                        cont = Some(r2.gamma.clone());
                        return (r2, a2, true);
                    }
                    if let Some(g) = &cont {
                        let (r3, a3) = eval_at(sda, elems, g);
                        if physical(&r3) {
                            cont = Some(r3.gamma.clone());
                            return (r3, a3, true);
                        }
                    }
                    (r2, a2, false)
                };
            if exported && err > 2.0e-3 {
                // Expand a bracket around the warm da until the cold-absorbed
                // torque straddles the target, then bisect.  Samples that fall
                // onto a non-physical branch are not adopted but still bracket.
                let mut best: Option<(f64, f64)> = None; // (err, da)
                let (mut lo, mut hi) = (da, da);
                let (mut q_lo, mut q_hi) = (res.torque, res.torque);
                // Expand in small steps so the continuation seed tracks the
                // physical branch: a 0.1 rad jump can leave it behind and
                // every further sample lands off-branch.  A failed sample
                // marks its end non-physical (INF); the expansion keeps
                // walking that direction, since the branch may reappear or
                // the bracketing torque may still lie past the dead spot.
                for _ in 0..10 {
                    let bracketed = q_lo.is_finite()
                        && q_hi.is_finite()
                        && (q_lo - q_target) * (q_hi - q_target) <= 0.0
                        && lo < hi;
                    if bracketed {
                        break;
                    }
                    if ((q_hi.is_finite() && q_hi < q_target) || !q_hi.is_finite())
                        && hi < DA_MAX
                    {
                        lo = hi;
                        q_lo = q_hi;
                        hi = (hi + 0.02).min(DA_MAX);
                        let (r2, _a2, ok) = sample(hi, elems);
                        q_hi = if ok { r2.torque } else { f64::INFINITY };
                    } else if ((q_lo.is_finite() && q_lo > q_target) || !q_lo.is_finite())
                        && lo > DA_MIN
                    {
                        hi = lo;
                        q_hi = q_lo;
                        lo = (lo - 0.02).max(DA_MIN);
                        let (r2, _a2, ok) = sample(lo, elems);
                        q_lo = if ok { r2.torque } else { f64::NEG_INFINITY };
                    } else {
                        break;
                    }
                }
                // Bisect only a true bracket: both ends physical and the
                // target between their torques.
                if q_lo.is_finite()
                    && q_hi.is_finite()
                    && (q_lo - q_target) * (q_hi - q_target) <= 0.0
                {
                    let (mut l, mut h) = (lo, hi);
                    for _ in 0..12 {
                        let mid = 0.5 * (l + h);
                        let (r2, _a2, ok) = sample(mid, elems);
                        if ok {
                            let e2 = (r2.torque - q_target).abs() / target;
                            if best.is_none_or(|(be, _)| e2 < be) {
                                best = Some((e2, mid));
                            }
                            if r2.torque <= q_target {
                                l = mid;
                            } else {
                                h = mid;
                            }
                        }
                    }
                }
                if let Some((e2, da2)) = best {
                    if e2 < err {
                        let (r3, a3, ok3) = sample(da2, elems);
                        if ok3 {
                            da = da2;
                            err = (r3.torque - q_target).abs() / target;
                            (res, alphas) = (r3, a3);
                            exported = true;
                        }
                    }
                }
            }
            (res, alphas, da, err, exported)
        };
        let (mut res, mut alphas, mut da, mut err, mut exported);
        let mut verified: &str = "cold";
        {
            let mut scale = 1.0;
            // The smoothed chord shape is the exported blade's chord: apply
            // it before the first solve so the cold re-match below measures
            // (and reports) the blade that is actually exported — even when
            // the raw design already sits on the target torque.
            set_chords(&mut elements, scale);
            let zeros = vec![0.0; m];
            (res, alphas, da, err, exported) = rematch(&mut elements, &zeros);
            for _ in 0..6 {
                if !exported || err <= 2.0e-3 {
                    break;
                }
                let next = (scale * q_target / res.torque.max(1.0e-9)).clamp(0.5, 2.0);
                if (next - scale).abs() / scale < 2.0e-3 {
                    break;
                }
                let prev_scale = scale;
                scale = next;
                set_chords(&mut elements, scale);
                // Continue the branch the last physical state lives on:
                // a fresh zero-seed solve at the scaled chord can settle on
                // a different root and defeat the monotone Q(scale) the
                // multiplicative update relies on.
                let seed = res.gamma.clone();
                let (r2, a2, d2, e2, x2) = rematch(&mut elements, &seed);
                if x2 && e2 < err {
                    (res, alphas, da, err, exported) = (r2, a2, d2, e2, true);
                } else {
                    // The scaled state was no better than the one the report
                    // keeps: restore its chords so the exported geometry and
                    // the reported operating point still describe one blade.
                    set_chords(&mut elements, prev_scale);
                    break;
                }
            }
        }
        // The zero-seed solve can die before the target torque: the cold
        // branch it settles on may end just short of it (or the scaled
        // continuation above may fail), while the design branch — the one
        // the optimizer converged on — does absorb it.  Verify the
        // exported blade on that branch too, seeded with the pass's own
        // converged circulation, and report whichever verification reached
        // the target; the cold state, when it does, stays the reported one.
        if win.gamma.len() == m && (!exported || err > 2.0e-3) {
            let (r4, a4, d4, e4, x4) = rematch(&mut elements, &win.gamma);
            if x4 && e4 < err {
                (res, alphas, da, err, exported) = (r4, a4, d4, e4, true);
                verified = "branch";
            }
        }
        // The exported blade: twist follows the solved inflow at the
        // matched `da` (cold- or design-branch-verified above), so the
        // geometry, the trace and the reported operating point all
        // describe the same state.  When no physical state of the exported
        // blade reaches the target torque (positive torque, every station
        // lifting, no circulation spike) — the pass itself only adopted
        // physical samples, so this should be rare — export the warm state
        // the optimizer measured instead: its twists (win.phis) at its own
        // matched `da`, and its matched operating point.  A re-solve can
        // settle on a spurious circulation-spike branch that the
        // optimizer's warm-started path never visits; that branch must not
        // become the exported blade.
        let warm_twist = |da: f64| -> Vec<f64> {
            win.alpha_base
                .iter()
                .map(|&ab| (ab + da).clamp(-lift_line::ALPHA_MAX, lift_line::ALPHA_MAX))
                .collect()
        };
        let (q, t) = if !exported {
            let warm_alphas = warm_twist(win.da);
            for (i, be) in elements.iter_mut().enumerate() {
                be.set_twist(win.phis[i] + warm_alphas[i]);
            }
            da = win.da;
            let err_warm = (win.q - q_target).abs() / target;
            crate::deprintln!(
                "proply: no physical re-solve of the winning blade (cold or design branch); exporting the warm state (T={:.4} Q={:.4}, err {:.2}%)",
                win.t,
                win.q,
                100.0 * err_warm
            );
            crate::dprintln!(
                "Lifting-line blade stations (warm state: T={:.4} Q={:.4}, err {:.2}%)",
                win.t,
                win.q,
                100.0 * err_warm
            );
            err = err_warm;
            // The station rows below describe the cold re-solve (the last
            // one attempted); label them as the reference state, not the
            // exported geometry's flow.
            (win.q, win.t)
        } else {
            for (i, be) in elements.iter_mut().enumerate() {
                be.set_twist(res.phi[i] + alphas[i]);
            }
            if self.param.mech_thickness {
                for (i, be) in elements.iter_mut().enumerate() {
                    be.thrust_n = Some(res.d_thrust[i]);
                }
            }
            crate::dprintln!(
                "Lifting-line blade stations ({}-verified: T={:.4} Q={:.4}, err {:.2}%)",
                verified,
                res.thrust,
                res.torque,
                100.0 * err
            );
            (res.torque, res.thrust)
        };
        // Smooth the exported twist too: the per-station flow solve adds
        // station-to-station noise to the inflow angles (twist = phi +
        // alpha), so re-fit it over the span with a solved least-squares
        // polynomial (the same idea the BEM path applies to its twist).
        // The flow state (and the reported operating point) does not depend
        // on the twist — the analysis prescribes the attack angles directly
        // — so this only straightens the manufactured geometry.
        {
            // The exported twist is a solved smooth spline too: among the
            // degree-4..6 least-squares polynomials fit to the per-station
            // phi + alpha, keep the one with the least spanwise curvature
            // (see [`smooth_twist_spline`]).  The fit may overshoot the
            // data's range slightly — clamping would put a corner exactly
            // where the curve turns.
            let raw: Vec<f64> = elements.iter().map(|be| be.get_twist()).collect();
            let smooth = smooth_twist_spline(&rr, &raw);
            for (i, be) in elements.iter_mut().enumerate() {
                be.set_twist(smooth[i]);
            }
        }
        for (i, &ri) in rr.iter().enumerate() {
            crate::dprintln!(
                "r={} camber={} alpha_base={} alpha={} phi={} twist={} chord={} ",
                ri,
                camber_smooth[i],
                win.alpha_base[i].to_degrees(),
                alphas[i].to_degrees(),
                res.phi[i].to_degrees(),
                elements[i].get_twist().to_degrees(),
                elements[i].foil.borrow().chord()
            );
            crate::dprintln!(
                "st{:2}: gamma={:8.3} ui={:7.2} vi={:7.2} dT={:8.4} dQ={:8.5}",
                i,
                res.gamma[i],
                res.u_i[i],
                res.v_i[i],
                res.d_thrust[i],
                res.d_torque[i]
            );
        }
        self.blade_elements = elements;
        // Remember the winning design so the pipeline's second run (the
        // mechanical-thickness law) can warm-start its incumbent pass from it.
        win.da = da;
        win.q = q;
        win.t = t;
        self.prev_win = Some(win);
        (q, t, err)
    }
}

/// True when the circulation profile is smooth: no station carries a
/// spurious spike.  A converged discrete solve can sit on a root where one
/// station holds several times its neighbours' circulation (seen up to ~9x
/// on low-aspect-ratio blades at 30+ stations): the largest bound
/// circulation may not exceed [`CIRC_SPIKE_RATIO`] times the second
/// largest, which is far above any smooth loading (those peak at ~1.5x).
/// Such states distort the induced inflow and the exported twist, so they
/// are rejected wherever the brake branches are.
fn circulation_smooth(gamma: &[f64]) -> bool {
    if gamma.len() < 3 {
        return true;
    }
    let mut s: Vec<f64> = gamma.iter().copied().filter(|g| g.is_finite()).collect();
    if s.len() < 3 {
        return true;
    }
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let largest = *s.last().unwrap();
    let second = s[s.len() - 2];
    largest <= CIRC_SPIKE_RATIO * second.max(1.0e-9) + 1.0e-6
}

/// [`circulation_smooth`]'s largest/second-largest ratio bound.
const CIRC_SPIKE_RATIO: f64 = 4.0;

/// Largest attached-flow ratio of the induced axial velocity to the local
/// speed `u_0 + ωr` outside the hub zone: beyond ~this the wake model's
/// small-perturbation assumption fails and the solved state is not a
/// physical loading.  The spurious states it rejects are far past any
/// sane value — flywoo's unattached hover state reached ~1.2-1.9, and a
/// diverged solve (turnigy's mechanical re-design) reached ~8600; healthy
/// designs sit well below 1.0.
const UI_RATIO_MAX: f64 = 1.0;

/// Hub-zone width as a fraction of the blade span: inside it the root /
/// hub-loss region legitimately sees reversed or extreme inflow, so the
/// attached-flow checks start only past this radius.
const HUB_ZONE_FRAC: f64 = 0.15;

/// Attached-flow check on the solved induced axial velocity (m/s): outside
/// the hub zone, `|u_i|` may not approach the local rotational speed.
fn flow_ui_sane(u_i: &[f64], rr: &[f64], u_0: f64, omega: f64) -> bool {
    if u_i.len() != rr.len() {
        return true;
    }
    let span = rr[rr.len() - 1] - rr[0];
    let hub_zone = rr[0] + HUB_ZONE_FRAC * span;
    u_i.iter().zip(rr.iter()).all(|(&ui, &ri)| {
        ri < hub_zone || ui.abs() <= UI_RATIO_MAX * (u_0 + omega * ri) + 1.0e-6
    })
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

/// Solved smooth spline for the exported twist: the per-station flow solve
/// adds station-to-station noise to the inflow angles, so fit degree-4..6
/// least-squares polynomials and keep the one with the least spanwise
/// curvature (largest second difference of the fit).  The winning
/// parameters come straight from the fit — no post-smoothing — and the
/// fit may overshoot the data slightly rather than being clamped into a
/// corner at its extremes.
fn smooth_twist_spline(rr: &[f64], twist: &[f64]) -> Vec<f64> {
    let mut best: Option<(f64, Vec<f64>)> = None;
    for deg in 4..=6 {
        if deg >= twist.len() {
            continue;
        }
        let fit = polyfit(rr, twist, deg);
        let vals: Vec<f64> = rr.iter().map(|&r| polyval(&fit, r)).collect();
        let d2 = vals
            .windows(3)
            .map(|w| ((w[2] - w[1]) - (w[1] - w[0])).abs())
            .fold(0.0, f64::max);
        if best.as_ref().is_none_or(|(b, _)| d2 < *b) {
            best = Some((d2, vals));
        }
    }
    best.expect("a degree-4..6 fit always exists for 5+ stations")
        .1
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

/// One lifting-line camber candidate seed: (label, per-station camber,
/// per-station smoothed best-L/D attack angles floored at zero lift,
/// warm-start controls for the first Nelder-Mead run — `Some` only for the
/// previous operating point's incumbent design, see [`Prop::prev_win`]).
type PassSeed = (String, Vec<f64>, Vec<f64>, Option<Vec<f64>>);

/// One camber candidate's converged design pass: the objective the
/// candidates compete on (`f = 1000*err - thrust`, i.e. absorb the target
/// torque exactly — err is the relative torque error — then maximise the
/// thrust at that torque), the measured operating point and the per-station
/// geometry.  Plain data so the passes can run on worker threads; the
/// caller rebuilds the winning blade elements from it.
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
    /// The converged circulation of the pass's own final state: the design
    /// branch's seed for the exported blade's verification, when the cold
    /// (zero-seed) solve settles off that branch and cannot reach the
    /// target torque.
    gamma: Vec<f64>,
    /// The winning chord-spline controls (`x ∈ [0, 1]`, the Nelder-Mead
    /// design variables) — kept so the next torque match can warm-start
    /// from this design instead of the full chord.
    design_x: Vec<f64>,
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
        assert!(
            (t_tip - 0.1 * p.param.hub_depth).abs() < 1e-6,
            "tip {}",
            t_tip
        );
    }

    #[test]
    fn mechanical_thickness_law_is_sized_from_the_loads() {
        // A synthetic converged blade (fixed induction, as a BEM design
        // leaves behind): the sizing installs a radius -> t/c law, and the
        // hub depth does not influence it — the load, chord and stiffness
        // decide, not the hub geometry (the TODO's requirement).
        let make = |hub_depth: f64| -> Prop {
            let param = DesignParameters {
                hub_depth,
                ..Default::default()
            };
            let store = Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json")));
            let mut p = Prop::new(param, 0.002, store);
            p.n_blades = 2;
            p.set_plate_mode(true);
            for r in [0.006, 0.02, 0.04, 0.06, 0.0625] {
                let foil = Rc::new(RefCell::new(FoilFamily::Naca4(crate::foil::Naca4::new(
                    0.012, 0.12, 0.06, 0.4,
                ))));
                let mut be = BladeElement::new(r, 0.002, foil, 0.4, 10000.0, 1.0, p.store.clone());
                // A believable induction profile: dv falling toward the tip.
                be.set_bem(6.0 * (1.0 - r / 0.0625) + 1.0, 0.05);
                p.blade_elements.push(be);
            }
            p
        };
        let mut p1 = make(0.006);
        assert!(p1.size_mechanical_thickness());
        assert!(p1.mech_thickness_law.is_some());
        let tip_def = p1.mech_tip_deflection.unwrap();
        // The sizing closes on the allowed deflection (5% of R by default),
        // so the predicted deflection is at most the limit.
        assert!(tip_def <= 0.05 * p1.param.radius * (1.0 + 1.0e-9), "tip {tip_def}");

        // The mechanical law is installed for the design loop: the stations
        // whose elements are rebuilt with it carry exactly the sized t/c.
        let law = p1.mech_thickness_law.as_ref().unwrap();
        for r in [0.006, 0.02, 0.04, 0.06, 0.0625] {
            let tc = law.eval(r);
            assert!((0.06 - 1.0e-12..=1.0).contains(&tc), "t/c {tc} outside law range");
        }

        // hub_depth is not involved: a different hub depth sizes the same
        // (radius, chord, load) blade identically.
        let mut p2 = make(0.012);
        assert!(p2.size_mechanical_thickness());
        let law2 = p2.mech_thickness_law.unwrap();
        for r in [0.006, 0.02, 0.04, 0.06, 0.0625] {
            assert!(
                (law.eval(r) - law2.eval(r)).abs() < 1.0e-12,
                "hub_depth leaked: {} vs {} at {}",
                law.eval(r),
                law2.eval(r),
                r
            );
        }
    }

    #[test]
    fn mechanical_thickness_nothing_to_size_keeps_geometric_law() {
        let mut p = test_prop();
        p.set_plate_mode(true);
        // No blade elements at all: sizing is a no-op, the geometric law
        // stays active.
        assert!(!p.size_mechanical_thickness());
        assert!(p.mech_thickness_law.is_none());

        // Elements with zero induction (lifting-line geometry before its
        // loads are stored): zero thrust, nothing to size.
        for r in [0.02, 0.04, 0.06] {
            let foil = Rc::new(RefCell::new(FoilFamily::Naca4(crate::foil::Naca4::new(
                0.012, 0.12, 0.06, 0.4,
            ))));
            let be = BladeElement::new(r, 0.002, foil, 0.4, 10000.0, 1.0, p.store.clone());
            p.blade_elements.push(be);
        }
        assert!(!p.size_mechanical_thickness());
        assert!(p.mech_thickness_law.is_none());
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
        assert!(
            dev < 1.0_f64.to_radians(),
            "fit drifted from data: {} rad",
            dev
        );
    }

    #[test]
    fn compose_camber_smooths_and_quantises() {
        // A span where the thick root prefers no camber and the outboard
        // sections prefer 0.04: the composed distribution must rise
        // smoothly from ~0 at the root to ~0.04 outboard, quantised to the
        // 0.01 polar-hash grid and clamped to the candidate range.
        let n = 30;
        let rr: Vec<f64> = (0..n)
            .map(|i| 0.006 + 0.062 * i as f64 / (n - 1) as f64)
            .collect();
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
        assert!(
            m_dist[n - 1] > m_dist[0] + 0.02,
            "no rise: {} -> {}",
            m_dist[0],
            m_dist[n - 1]
        );
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
            FoilFamily::Cst(_) | FoilFamily::Arad(_) => {
                panic!("expected a NACA4 foil (default family)")
            }
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
            FoilFamily::Arad(_) => panic!("expected a CST foil with param.cst set"),
        }
    }

    #[test]
    fn blade_element_arad_family_is_applied() {
        // With the ARA-D family selected, stations carry Arad foils sized to
        // the station thickness law (thickness/chord reaching the foil
        // exactly).  The camber candidate does not apply: the ARA-D tables
        // carry their own camber, so every candidate yields the same foil.
        let mut p = test_prop();
        p.param.arad = true;
        p.param.hub_depth = 0.006; // nonzero root depth: thickness law > 0
        let be = p.new_blade_element_with_chord(0.02, 1000.0, 0.1, 0.008, 0.04);
        let f = be.foil.borrow();
        match &*f {
            FoilFamily::Arad(a) => {
                let t_law = p.get_foil_thickness(0.02) / 0.008;
                assert!(
                    (a.thickness() - t_law).abs() < 1e-12,
                    "thickness {} vs law {}",
                    a.thickness(),
                    t_law
                );
                assert!(a.camber() > 0.01, "inherent camber {}", a.camber());
            }
            _ => panic!("expected an ARA-D foil with param.arad set"),
        }
    }

    #[test]
    fn lift_line_design_finite_and_torque_matching() {
        // Plate polar (no rust-foil) so this is fast; checks the coupled
        // lifting-line produces a finite thrust while absorbing the target
        // torque at the design RPM.
        let mut prop = test_prop();
        // coarser grid for speed
        prop.radial_steps = 12;
        prop.set_plate_mode(true);
        // The plate test propeller's full-chord torque ceiling is ~0.06 Nm
        // at 12000 rpm (see the design's saturated ceiling below), so pick
        // a target inside the absorbable range.
        let q_target = 0.03;
        let (q, t, err) = prop.lift_line_design(12000.0, q_target, Some(4.0));
        assert!(t.is_finite() && t > 0.0, "thrust {} not finite/positive", t);
        assert!(q.is_finite() && q > 0.0, "torque {} not finite/positive", q);
        // The target must actually be met (a u/v-swap in the force projection
        // once made the model thrust ~v/u times too small, so the design
        // could never converge onto the target).
        assert!(
            err < 0.05,
            "torque error {:.4} vs the {:.3} Nm target",
            err,
            q_target
        );
        assert!(q < 0.5, "torque {} not bounded", q);
    }

    #[test]
    fn bem_design_matches_target_torque() {
        // Plate polars (no rust-foil) so this is fast: the operating-point
        // match must land the absorbed torque on the target at the design
        // RPM, converging from the internal momentum-theory seed.
        let mut p = test_prop();
        p.param.hub_depth = 0.003; // nonzero thickness law, as real props have
        p.radial_steps = 10;
        p.set_plate_mode(true);
        let (q0, _t0) = p.full_optimize(12000.0, 2.0);
        assert!(q0 > 0.0, "reference torque {} not positive", q0);
        let q_target = 0.8 * q0;
        let res = p.design_for_torque(12000.0, q_target, None);
        assert!(
            (res.torque - q_target).abs() / q_target < 0.01,
            "Q {} vs target {}",
            res.torque,
            q_target
        );
        assert!(
            res.thrust.is_finite() && res.thrust > 0.0,
            "thrust {}",
            res.thrust
        );
        // A feasible design converges onto the operating point without a
        // warning (the unachievable-operating-point warning must stay off).
        assert!(res.warning.is_none(), "unexpected warning {:?}", res.warning);
        assert_eq!(res.converged_stations, res.total_stations);
        // The converged geometry stays in place for the STEP/YAML writers.
        assert!(!p.blade_elements.is_empty());
    }

    #[test]
    fn lifting_line_design_matches_target_torque() {
        // The direct operating-point design through the coupled lifting
        // line: `design_for_torque` must absorb the target torque at the
        // design RPM in a single pass (no thrust-target iteration).
        let mut p = test_prop();
        p.param.lifting_line = true;
        p.param.hub_depth = 0.003; // nonzero thickness law, as real props have
        p.radial_steps = 10;
        p.set_plate_mode(true);
        // Reference torque at a working scale: the legacy BEM design for
        // the same propeller at 2 N (the lifting line's absorbed-torque
        // range covers this operating point).
        let (q0, _t0) = p.full_optimize(12000.0, 2.0);
        assert!(q0 > 0.0, "reference torque {} not positive", q0);
        let q_target = 0.7 * q0;
        let res = p.design_for_torque(12000.0, q_target, None);
        assert!(
            (res.torque - q_target).abs() / q_target < 0.01,
            "Q {} vs target {}",
            res.torque,
            q_target
        );
        assert!(
            res.thrust.is_finite() && res.thrust > 0.0,
            "thrust {}",
            res.thrust
        );
        assert!(res.warning.is_none(), "unexpected warning {:?}", res.warning);
    }

    #[test]
    fn unreachable_torque_target_reports_warning() {
        // Regression (browser_demo plate-vs-full-polar discrepancy): when
        // the geometry cannot absorb the design torque, `design_for_torque`
        // must say so — a warning with the closest design and the station
        // coverage — instead of silently reporting a broken (near-zero
        // torque) design as if it had converged.  A zero `hub_depth` makes
        // every station's max chord zero, so no station can carry any load:
        // the loop stalls at Q = 0 on the first pass.
        let mut p = test_prop();
        p.param.hub_depth = 0.0; // degenerate depth law: no chord allowed
        p.radial_steps = 10;
        p.set_plate_mode(true);
        let res = p.design_for_torque(12000.0, 0.05, None);
        assert!(res.torque == 0.0, "torque {}", res.torque);
        assert!(res.thrust == 0.0, "thrust {}", res.thrust);
        let w = res
            .warning
            .expect("unreachable operating point must carry a warning");
        assert!(w.contains("not achievable"), "warning: {}", w);
        assert!(
            res.converged_stations <= res.total_stations && res.total_stations > 0,
            "coverage {} / {}",
            res.converged_stations,
            res.total_stations
        );
    }
}
