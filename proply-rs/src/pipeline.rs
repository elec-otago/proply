// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! The portable design pipeline shared by the CLI and the WebAssembly
//! build.
//!
//! Everything between "parsed parameters" and "output text" lives here so
//! the browser entry point runs exactly the code the CLI does: no file
//! I/O, just the design compute and the STEP/YAML text.

use std::sync::{Arc, Mutex};

use crate::cache::PolarStore;
use crate::design_parameters::DesignParameters;
use crate::prop::Prop;
use crate::step_out;
use crate::yaml_out;

/// A finished design: everything the CLI writes to files, as text.
pub struct DesignOutcome {
    /// The propeller as a STEP (AP242) document.
    pub step: String,
    /// The design summary (`yaml_out::summary`).
    pub yaml: String,
    /// Absorbed torque at the design point (N m).
    pub torque: f64,
    /// Produced thrust at the design point (N).
    pub thrust: f64,
    /// Design rotation rate (rpm).
    pub rpm: f64,
    /// Power at the motor operating point (W).
    pub power: f64,
    /// Set when the design could not absorb the demanded torque: an explicit
    /// note describing the closest design that was reached (CLI, YAML and
    /// browser demo all surface it).  `None` for a converged operating point.
    pub warning: Option<String>,
}

/// Run one full propeller design from parsed parameters against `store`,
/// returning the STEP and YAML output text.  This is exactly the pipeline
/// the CLI runs after parsing: build the [`Prop`], converge the design onto
/// the motor operating point, then serialize.
pub fn run_design(
    param: &DesignParameters,
    store: Arc<Mutex<PolarStore>>,
) -> Result<DesignOutcome, String> {
    let element_width = (param.radius - param.hub_radius) / param.element_count as f64;
    let mut p = Prop::new(param.clone(), element_width, store);
    p.n_blades = param.blades;
    p.set_plate_mode(param.plate);

    // An explicitly specified operating point (motor_torque + motor_RPM in
    // the design file, e.g. an engine) overrides the electric motor model's
    // maximum-efficiency point.
    let (optimum_torque, optimum_rpm, power) = param.motor_operating_point();

    // The design converges onto the motor operating point: it absorbs the
    // design torque at the design RPM and maximises efficiency (thrust, at
    // the fixed shaft power) there — directly for the lifting-line loop,
    // by iterating the thrust target for the legacy BEM loop.  When the
    // geometry cannot absorb the demanded torque, the result carries a
    // warning describing the closest design.
    let mut res = p.design_for_torque(optimum_rpm, optimum_torque, param.ar);

    // The mechanical thickness law is sized on the converged design's
    // station loads (the blade as a cantilever beam under its own thrust),
    // and the design then re-runs on that law so the reported operating
    // point and geometry match the mechanically sized blade.
    if param.mech_thickness && p.size_mechanical_thickness() {
        let res_mech = p.design_for_torque(optimum_rpm, optimum_torque, param.ar);
        if res_mech.torque > 0.0 || res_mech.warning.is_none() {
            // The mechanical re-design reached the operating point (or at
            // least produced a real closest design): report it.
            res = res_mech;
        } else {
            // The mechanical re-design found no usable state at all (the
            // sized blade cannot be solved at this operating point, e.g.
            // turnigy_CA_120: only unattached/diverged flow states exist).
            // Keep the first design and say so, instead of reporting a
            // blade that was never built.
            let base = res.warning.clone().unwrap_or_default();
            let note = if base.is_empty() {
                "the mechanical-thickness re-design at this operating point was not achievable; this design uses the geometric thickness law".to_string()
            } else {
                format!("{base} (the mechanical-thickness re-design at this operating point was not achievable; this design uses the geometric thickness law)")
            };
            res.warning = Some(note);
        }
    }
    let (q, t) = (res.torque, res.thrust);

    let step = step_out::write_prop(&mut p, param.n)?;
    let motor_info = yaml_out::MotorInfo {
        kv_rpm_per_volt: param.motor_Kv,
        voltage_v: param.motor_volts,
        winding_resistance_ohm: param.motor_winding_resistance,
        no_load_current_a: param.motor_no_load_current,
        optimum_rpm,
        optimum_torque_nm: optimum_torque,
        max_power_w: power,
    };
    let yaml = yaml_out::summary(&p, optimum_rpm, t, q, &motor_info, res.warning.as_deref());

    Ok(DesignOutcome {
        step,
        yaml,
        torque: q,
        thrust: t,
        rpm: optimum_rpm,
        power,
        warning: res.warning,
    })
}
