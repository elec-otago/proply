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
    let resolution_m = (param.radius - param.hub_radius) / param.resolution as f64;
    let mut p = Prop::new(param.clone(), resolution_m, store);
    p.n_blades = param.blades;
    p.set_plate_mode(param.plate);

    // An explicitly specified operating point (motor_torque + motor_RPM in
    // the design file, e.g. an engine) overrides the electric motor model's
    // maximum-efficiency point.
    let (optimum_torque, optimum_rpm, power) = param.motor_operating_point();

    // The design converges onto the motor operating point: the thrust
    // target is iterated until the blade absorbs the design torque at the
    // design RPM.  The inner loops maximise efficiency at a matched thrust,
    // so the converged design is the maximum-efficiency design at that
    // operating point.  When the geometry cannot absorb the demanded torque,
    // the result carries a warning describing the closest design.
    let res = p.design_for_torque(optimum_rpm, optimum_torque, param.thrust, param.ar);
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
