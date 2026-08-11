//! The propeller design loop, ported from `proply/prop.py`.
//!
//! `full_optimize` sweeps radial stations from tip to hub, sizing each
//! blade element with `optimize_all`, then smooths the twist and chord
//! distributions and re-evaluates the forces with `get_forces`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use crate::blade_element::BladeElement;
use crate::cache::PolarStore;
use crate::design_parameters::DesignParameters;
use crate::foil::{FoilLike, Naca4};
use crate::optimize;
use crate::pchip::Pchip;
use crate::polyfit::{polyfit, polyval};
use crate::smooth::smooth;

/// A propeller: a collection of blade elements plus the design parameters.
pub struct Prop {
    pub param: DesignParameters,
    pub radial_resolution: f64,
    pub radial_steps: usize,
    pub n_blades: usize,
    pub blade_elements: Vec<BladeElement<Naca4>>,
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

    /// Create a blade element with a NACA4 foil sized for the station.
    fn new_blade_element(&mut self, r: f64, rpm: f64, twist: f64) -> BladeElement<Naca4> {
        let y_limit = self.get_max_depth(r);
        let x_limit = self.get_max_chord(r, twist);
        let thickness = self.get_foil_thickness(r);

        let foil = Rc::new(RefCell::new(Naca4::new(x_limit, thickness / x_limit, 0.0, 0.4)));
        foil.borrow_mut().set_trailing_edge(self.param.trailing_edge / 1000.0);

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

    /// Allowed chord as a function of radius (m): k / r^2, capped by the
    /// blade-count spacing.
    pub fn get_max_chord(&self, r: f64, twist: f64) -> f64 {
        let k = self.param.tip_chord * self.param.radius.powi(2);
        let c = k / r.powi(2);
        let upper_limit = (2.0 * std::f64::consts::PI * r / (self.n_blades as f64 + 2.0)) / twist.cos();
        c.min(upper_limit)
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
    /// law between the hub depth and a tenth of it at the tip.
    pub fn get_foil_thickness(&self, r: f64) -> f64 {
        let thickness_root = self.param.hub_depth * 1.0;
        let thickness_end = self.param.hub_depth * 0.1;
        let p = 0.3;
        let k = (thickness_root - thickness_end)
            / (self.param.hub_radius.powf(p) - self.param.radius.powf(p));
        let s = thickness_end - k * self.param.radius.powf(p);
        s + k * r.powf(p)
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

            println!(
                "r={} theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={} ",
                r,
                theta.to_degrees(),
                dv,
                a_prime,
                d_t,
                d_m,
                d_t / d_m
            );
            println!("{}", be);

            self.blade_elements.push(be);
            prev_twist = theta;
        }

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
}
