// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Blade element: a single radial slice of the propeller blade, ported from
//! `proply/blade_element.py`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use crate::cache::PolarStore;
use crate::foil::FoilLike;
use crate::optimize;
use crate::simulator::FoilSimulator;

/// One blade element at radius `r` with span `dr`, holding a foil and twist.
pub struct BladeElement<F: FoilLike> {
    pub r: f64,
    pub dr: f64,
    pub foil: Rc<RefCell<F>>,
    pub fs: FoilSimulator<F>,
    pub zero_lift_angle: Option<f64>,
    pub dv: f64,
    pub a_prime: f64,
    pub velocity: f64,
    pub rpm: f64,
    pub omega: f64,
    pub u_0: f64,
    /// Whether the last BEM solve (`bem` / `get_forces`) converged on a
    /// momentum-equilibrium state.  A station whose solve did not converge
    /// keeps its last induction state and is reported (and counted
    /// separately) instead of being silently dropped.
    pub converged: bool,
    /// This station's annular thrust (N, all blades) from the design loop.
    /// BEM elements compute it from the converged induction (`d_t()`); the
    /// lifting-line loop leaves no induction on its elements, so
    /// [`Prop::lift_line_design`] stores the circulation solve's element
    /// thrust here for the mechanical-thickness sizing.
    pub thrust_n: Option<f64>,
    twist: f64,
}

impl<F: FoilLike> BladeElement<F> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        r: f64,
        dr: f64,
        foil: Rc<RefCell<F>>,
        twist: f64,
        rpm: f64,
        u_0: f64,
        store: Arc<Mutex<PolarStore>>,
    ) -> Self {
        let fs = FoilSimulator::new(foil.clone(), store);
        Self {
            r,
            dr,
            foil,
            fs,
            zero_lift_angle: None,
            dv: 0.0,
            a_prime: 0.0,
            velocity: 0.0,
            rpm,
            omega: 2.0 * std::f64::consts::PI * rpm / 60.0,
            u_0,
            converged: true,
            thrust_n: None,
            twist,
        }
    }

    pub fn set_chord(&mut self, c: f64) {
        self.foil.borrow_mut().modify_chord(c);
    }

    /// Switch the element's polar model to the analytic flat plate.
    pub fn set_plate_mode(&mut self, on: bool) {
        self.fs.set_plate_mode(on);
    }

    pub fn set_twist(&mut self, twist: f64) {
        self.twist = twist;
    }

    pub fn get_twist(&self) -> f64 {
        self.twist
    }

    /// Element thrust (newtons).
    pub fn d_t(&self) -> f64 {
        optimize::d_t(self.dv, self.r, self.dr, self.u_0)
    }

    /// Element torque (N m).
    pub fn d_m(&self) -> f64 {
        optimize::d_m(self.dv, self.a_prime, self.r, self.dr, self.omega, self.u_0)
    }

    /// Run one BEM iteration against the current induced velocity.
    pub fn bem(&mut self, n_blades: usize) -> (f64, f64, f64) {
        let (dv, a_prime, err) = optimize::bem_iterate(
            &self.fs,
            self.dv,
            self.twist,
            self.rpm,
            self.r,
            self.dr,
            self.u_0,
            n_blades as f64,
        );
        self.set_bem(dv, a_prime);
        (dv, a_prime, err)
    }

    /// Store the converged (dv, a_prime) and the resulting flow velocity.
    pub fn set_bem(&mut self, dv: f64, a_prime: f64) {
        self.dv = dv;
        self.a_prime = a_prime;
        let u = self.u_0 + dv;
        let v = self.omega * self.r * (1.0 - a_prime);
        self.velocity = (u * u + v * v).sqrt();
    }

    /// The 3-D profile points at this station: returns (lower_line,
    /// upper_line), each an `n`-element array of [x, y, z] points.  The
    /// profile is wrapped onto the cylinder of radius `r` and offset by the
    /// scimitar angle (mirrors `get_foil_points`), then re-centred
    /// vertically: the twist rotation pivots on a point of the chord line
    /// (see `get_points`), which swings a twisted section to one side of
    /// z = 0 — without the re-centring the root sections ride high and
    /// protrude above the hub (the hub is centred on z = 0).
    pub fn get_foil_points(
        &self,
        n: usize,
        scimitar_offset: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let (pl, pu) = {
            let f = self.foil.borrow();
            f.get_points(n, self.twist)
        };
        let r = self.r;
        let scimitar_angle = (scimitar_offset / r).atan();
        let circumference = 2.0 * std::f64::consts::PI * r;

        // Centre the station on its own vertical extent, so the stacking
        // axis (z = 0) passes through the rotated section's midpoint.
        let mut zmin = f64::INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        for p in pl.iter().chain(pu.iter()) {
            zmin = zmin.min(p[1]);
            zmax = zmax.max(p[1]);
        }
        let zc = 0.5 * (zmin + zmax);

        let wrap = |y: f64| {
            let theta = 2.0 * std::f64::consts::PI * y / circumference + scimitar_angle;
            [r * theta.cos(), r * theta.sin()]
        };

        let mut lower = Vec::with_capacity(n);
        let mut upper = Vec::with_capacity(n);
        for i in 0..n {
            // pl = (yl = foil x/chord dir, zl = foil y/thickness dir)
            let [lx, ly] = wrap(pl[i][0]);
            lower.push([lx, ly, pl[i][1] - zc]);
            let [ux, uy] = wrap(pu[i][0]);
            upper.push([ux, uy, pu[i][1] - zc]);
        }
        (lower, upper)
    }
}

impl<F: FoilLike> std::fmt::Display for BladeElement<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dt = self.d_t();
        let dm = self.d_m();
        let eff = if dm != 0.0 { dt / dm } else { f64::INFINITY };
        write!(
            f,
            "BladeElement(r={:5.3}, twist={:5.2}, dv={:4.1}, eff={:4.1})",
            self.r,
            self.get_twist().to_degrees(),
            self.dv,
            eff
        )
    }
}
