// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! End-to-end test of the portable design pipeline — the exact code path
//! the WebAssembly build exposes (`proply_rs::wasm::PropSession::design`):
//! parsed parameters in, STEP + YAML text out.

use std::sync::{Arc, Mutex};

use proply_rs::cache::PolarStore;
use proply_rs::design_parameters::DesignParameters;
use proply_rs::pipeline;

/// A small, fast design: the known-converging default geometry and motor,
/// with plate polars (no rust-foil sweeps), coarse stations and a coarse
/// STEP loft.
fn fast_params() -> DesignParameters {
    let json = r#"{
        "name": "pipeline_test",
        "altitude": 0.0,
        "forward_airspeed": 1.0,
        "motor_Kv": 1200,
        "motor_volts": 11.0,
        "motor_no_load_current": 0.5,
        "motor_winding_resistance": 0.206,
        "thrust": 2.0,
        "radius": 0.0625,
        "tip_chord": 0.007,
        "hub_radius": 0.005,
        "hub_depth": 0.003,
        "blades": 2,
        "bem": true,
        "plate": true,
        "resolution": 10,
        "n": 12
    }"#;
    DesignParameters::from_json(json).expect("valid parameters")
}

#[test]
fn pipeline_produces_step_and_yaml() {
    let param = fast_params();
    let store: Arc<Mutex<PolarStore>> =
        Arc::new(Mutex::new(PolarStore::in_memory()));
    let outcome = pipeline::run_design(&param, store).expect("design converges");

    // The STEP text is a valid Part 21 document with a blade and a hub.
    assert!(outcome.step.starts_with("ISO-10303-21;"), "wrong STEP header");
    assert!(outcome.step.contains("'blade'"), "no blade part");
    assert!(outcome.step.contains("'hub'"), "no hub part");

    // The YAML summary describes the finished design.
    assert!(outcome.yaml.contains("pipeline_test"), "no prop name");
    assert!(outcome.yaml.contains("rpm:"), "no rpm entry");

    // The design matched a real operating point.
    assert!(outcome.thrust.is_finite() && outcome.thrust > 0.0, "thrust");
    assert!(outcome.torque.is_finite() && outcome.torque > 0.0, "torque");
    assert!(outcome.rpm > 0.0, "rpm");
    // A feasible design converges onto the operating point without warning.
    assert!(outcome.warning.is_none(), "warning: {:?}", outcome.warning);
}
