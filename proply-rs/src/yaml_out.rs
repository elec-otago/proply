// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! YAML summary of a finished prop design.
//!
//! Every design run writes `<propname>.yml` beside the STEP output (or
//! beside an explicit `--step-file`): a machine-readable record of how the
//! prop performs (RPM, thrust, torque, shaft power, efficiencies) at the
//! motor operating point it was designed for, plus the per-station section
//! list (radius, chord, twist, camber, thickness) of the blade that was
//! built.

use serde::Serialize;

use crate::foil::{FoilFamily, FoilLike};
use crate::lift_line::RHO;
use crate::optimize::rpm2omega;
use crate::prop::Prop;

/// The motor operating point the design was run at (from `Motor::get_qmax`
/// and `Motor::get_pmax` at the design voltage).
#[derive(Debug, Serialize)]
pub struct MotorInfo {
    pub kv_rpm_per_volt: f64,
    pub voltage_v: f64,
    pub winding_resistance_ohm: f64,
    pub no_load_current_a: f64,
    pub optimum_rpm: f64,
    pub optimum_torque_nm: f64,
    pub max_power_w: f64,
}

/// One blade station, hub -> tip.  The induction-dependent fields (induced
/// velocity and the element loads) are emitted only for BEM designs: that
/// design loop leaves the converged induction on the elements, while the
/// lifting-line loop solves the coupled system internally and reports only
/// the totals.
#[derive(Debug, Serialize)]
#[allow(non_snake_case)] // r/R is the standard propeller notation
struct Section {
    r_m: f64,
    r_over_R: f64,
    chord_m: f64,
    twist_deg: f64,
    camber: f64,
    thickness_m: f64,
    thickness_fraction: f64,
    dr_m: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    dv_mps: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    a_prime: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    thrust_n: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    torque_nm: Option<f64>,
    /// Whether this station's BEM solve reached a momentum equilibrium (BEM
    /// designs only): the design loop and the summary report the coverage
    /// instead of silently treating failed stations as zeros.
    #[serde(skip_serializing_if = "Option::is_none")]
    converged: Option<bool>,
}

#[derive(Debug, Serialize)]
struct DesignInfo {
    mode: &'static str,
    blades: usize,
    radius_m: f64,
    hub_radius_m: f64,
    forward_airspeed_mps: f64,
    altitude_m: f64,
    stations: usize,
    element_count: usize,
    loft_points: usize,
    foil_family: &'static str,
    /// "geometric" (the p = 0.3 power law on the hub depth) or "mechanical"
    /// (the beam-deflection sizing — see [`crate::thickness`]).
    thickness_law: &'static str,
    /// Predicted tip deflection of the mechanical thickness law (mm),
    /// present only for mechanical-law designs.
    #[serde(skip_serializing_if = "Option::is_none")]
    tip_deflection_mm: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    camber: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    min_aspect_ratio: Option<f64>,
    chord_spline_n: usize,
    auto: bool,
}

#[derive(Debug, Serialize)]
struct Performance {
    rpm: f64,
    thrust_n: f64,
    torque_nm: f64,
    shaft_power_w: f64,
    /// `T u_0 / P_shaft` — the useful-work fraction in forward flight
    /// (zero in hover).
    propulsive_efficiency: f64,
    /// Hover figure of merit `T^{3/2} / (sqrt(2 rho A) P_shaft)` against
    /// the ideal actuator disk.
    figure_of_merit: f64,
    tip_speed_mps: f64,
}

#[derive(Debug, Serialize)]
struct DesignSummary {
    name: String,
    generated_by: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    warning: Option<String>,
    design: DesignInfo,
    motor: MotorInfo,
    performance: Performance,
    sections: Vec<Section>,
}

/// Max camber of a station foil (fraction of chord).
fn family_camber(f: &FoilFamily) -> f64 {
    match f {
        FoilFamily::Naca4(n) => n.m,
        FoilFamily::Cst(c) => c.camber(),
        FoilFamily::Arad(a) => a.camber(),
    }
}

/// Round to 6 decimal places: metre and degree values in a summary do not
/// need full float precision, and fixed rounding keeps the file diffable.
fn r6(x: f64) -> f64 {
    (x * 1.0e6).round() / 1.0e6
}

/// Build the YAML summary of a finished design: `rpm`/`thrust`/`torque` are
/// the design operating point and totals (as printed at the end of the run).
/// `warning`, when set, records that the design could not reach its
/// operating point (an explicit note instead of a silently unmatched one).
pub fn summary(
    p: &Prop,
    rpm: f64,
    thrust: f64,
    torque: f64,
    motor: &MotorInfo,
    warning: Option<&str>,
) -> String {
    let param = &p.param;
    let mode = if param.lifting_line {
        "lifting-line"
    } else {
        "bem"
    };
    // Element induction state is only meaningful for the BEM design loop.
    let bem_aero = mode == "bem";

    let sections: Vec<Section> = p
        .blade_elements
        .iter()
        .map(|be| {
            let (chord, t_frac, camber) = {
                let f = be.foil.borrow();
                (f.chord(), f.thickness(), family_camber(&f))
            };
            Section {
                r_m: r6(be.r),
                r_over_R: r6(be.r / param.radius),
                chord_m: r6(chord),
                twist_deg: r6(be.get_twist().to_degrees()),
                camber: r6(camber),
                thickness_m: r6(chord * t_frac),
                thickness_fraction: r6(t_frac),
                dr_m: r6(be.dr),
                dv_mps: bem_aero.then(|| r6(be.dv)),
                a_prime: bem_aero.then(|| r6(be.a_prime)),
                thrust_n: bem_aero.then(|| r6(be.d_t())),
                torque_nm: bem_aero.then(|| r6(be.d_m())),
                converged: bem_aero.then_some(be.converged),
            }
        })
        .collect();

    let omega = rpm2omega(rpm);
    let shaft_power = torque * omega;
    let eta_prop = if shaft_power > 0.0 {
        thrust * param.forward_airspeed / shaft_power
    } else {
        0.0
    };
    let fom = if shaft_power > 0.0 {
        thrust.powf(1.5) / ((2.0 * RHO * param.area()).sqrt() * shaft_power)
    } else {
        0.0
    };

    let s = DesignSummary {
        name: param.name.clone(),
        generated_by: format!("proply-rs {}", env!("CARGO_PKG_VERSION")),
        warning: warning.map(|w| w.to_string()),
        design: DesignInfo {
            mode,
            blades: p.n_blades,
            radius_m: r6(param.radius),
            hub_radius_m: r6(param.hub_radius),
            forward_airspeed_mps: r6(param.forward_airspeed),
            altitude_m: r6(param.altitude),
            stations: sections.len(),
            element_count: param.element_count,
            loft_points: param.n,
            foil_family: if param.arad {
                "arad"
            } else if param.cst {
                "cst"
            } else {
                "naca4"
            },
            thickness_law: if p.mech_thickness_law.is_some() {
                "mechanical"
            } else {
                "geometric"
            },
            tip_deflection_mm: p.mech_tip_deflection.map(|d| r6(d * 1000.0)),
            camber: param.camber.map(r6),
            min_aspect_ratio: param.ar,
            chord_spline_n: param.chord_spline_n,
            auto: param.auto,
        },
        motor: MotorInfo {
            kv_rpm_per_volt: r6(motor.kv_rpm_per_volt),
            voltage_v: r6(motor.voltage_v),
            winding_resistance_ohm: r6(motor.winding_resistance_ohm),
            no_load_current_a: r6(motor.no_load_current_a),
            optimum_rpm: r6(motor.optimum_rpm),
            optimum_torque_nm: r6(motor.optimum_torque_nm),
            max_power_w: r6(motor.max_power_w),
        },
        performance: Performance {
            rpm: r6(rpm),
            thrust_n: r6(thrust),
            torque_nm: r6(torque),
            shaft_power_w: r6(shaft_power),
            propulsive_efficiency: r6(eta_prop),
            figure_of_merit: r6(fom),
            tip_speed_mps: r6(omega * param.radius),
        },
        sections,
    };
    let body = serde_yaml::to_string(&s).expect("summary serializes");
    format!("# proply-rs design summary: {}\n{}", param.name, body)
}

/// Write a summary file, creating parent directories as needed (mirrors
/// `step_out::write_step_file`).
pub fn write_yaml_file(path: &str, text: &str) -> std::io::Result<()> {
    if let Some(parent) = std::path::Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    std::fs::write(path, text)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::cell::RefCell;
    use std::rc::Rc;
    use std::sync::{Arc, Mutex};

    use crate::blade_element::BladeElement;
    use crate::cache::PolarStore;
    use crate::design_parameters::DesignParameters;

    fn synthetic_prop() -> Prop {
        // Same shape as tests/step_output.rs: default parameters, four
        // stations with plain NACA foils, no design loop.
        let param = DesignParameters::default();
        let store: Arc<Mutex<PolarStore>> =
            Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json")));
        let mut p = Prop::new(param, 0.002, store.clone());
        p.n_blades = 2;
        for r in [0.006, 0.02, 0.04, 0.06] {
            let foil = Rc::new(RefCell::new(FoilFamily::Naca4(crate::foil::Naca4::new(
                0.012, 0.12, 0.06, 0.4,
            ))));
            let mut be = BladeElement::new(r, 0.002, foil, 0.4, 10000.0, 1.0, store.clone());
            // Converged-looking induction so the BEM aero fields exercise
            // their formatting path.
            be.set_bem(2.0, 0.05);
            p.blade_elements.push(be);
        }
        p
    }

    fn motor_info() -> MotorInfo {
        MotorInfo {
            kv_rpm_per_volt: 1200.0,
            voltage_v: 11.0,
            winding_resistance_ohm: 0.206,
            no_load_current_a: 0.5,
            optimum_rpm: 5000.0,
            optimum_torque_nm: 0.15,
            max_power_w: 150.0,
        }
    }

    #[test]
    fn bem_summary_parses_with_expected_content() {
        let p = synthetic_prop();
        let text = summary(&p, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");

        assert_eq!(y["name"], "hello world");
        assert!(text.starts_with("# proply-rs design summary:"));
        assert_eq!(y["design"]["mode"], "bem");
        assert_eq!(y["design"]["foil_family"], "naca4");
        assert_eq!(y["design"]["stations"], 4);
        assert_eq!(y["design"]["blades"], 2);
        assert_eq!(y["design"]["camber"], serde_yaml::Value::Null);
        // No warning on a converged design: the key is absent, not empty.
        assert!(y.get("warning").is_none(), "unexpected warning field");

        assert_eq!(y["motor"]["optimum_rpm"], 5000.0);
        assert_eq!(y["motor"]["voltage_v"], 11.0);

        assert_eq!(y["performance"]["rpm"], 5000.0);
        assert_eq!(y["performance"]["thrust_n"], 2.25);
        assert_eq!(y["performance"]["torque_nm"], 0.05);
        let omega = rpm2omega(5000.0);
        assert_eq!(y["performance"]["shaft_power_w"], r6(0.05 * omega));
        assert_eq!(y["performance"]["tip_speed_mps"], r6(omega * 0.0625));

        let sections = y["sections"].as_sequence().unwrap();
        assert_eq!(sections.len(), 4);
        let s0 = &sections[0];
        assert_eq!(s0["r_m"], 0.006);
        assert_eq!(s0["r_over_R"], r6(0.006 / 0.0625));
        assert_eq!(s0["chord_m"], 0.012);
        assert_eq!(s0["twist_deg"], r6(0.4f64.to_degrees()));
        assert_eq!(s0["camber"], 0.06);
        assert_eq!(s0["thickness_fraction"], 0.12);
        assert_eq!(s0["thickness_m"], r6(0.012 * 0.12));
        assert_eq!(s0["dr_m"], 0.002);
        // BEM designs carry the element induction and loads.
        assert_eq!(s0["dv_mps"], 2.0);
        assert_eq!(s0["a_prime"], 0.05);
        assert!(s0["thrust_n"].as_f64().unwrap() > 0.0);
        assert!(s0["torque_nm"].as_f64().unwrap() > 0.0);
        // BEM sections carry the station convergence flag (default: converged).
        assert_eq!(s0["converged"], true);

        // Radii run hub -> tip.
        let rs: Vec<f64> = sections
            .iter()
            .map(|s| s["r_m"].as_f64().unwrap())
            .collect();
        assert!(
            rs.windows(2).all(|w| w[0] < w[1]),
            "radii not ordered: {:?}",
            rs
        );
    }

    #[test]
    fn mechanical_thickness_law_is_reported_in_the_summary() {
        // The summary names the active thickness law and, for the mechanical
        // law, the predicted tip deflection it was sized to close on.
        let mut p = synthetic_prop();
        let text = summary(&p, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");
        assert_eq!(y["design"]["thickness_law"], "geometric");
        assert!(y["design"].get("tip_deflection_mm").is_none());

        p.mech_thickness_law = Some(crate::pchip::Pchip::new(
            &[0.006, 0.06],
            &[0.28, 0.06],
        ));
        p.mech_tip_deflection = Some(0.0034);
        let text = summary(&p, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");
        assert_eq!(y["design"]["thickness_law"], "mechanical");
        assert_eq!(y["design"]["tip_deflection_mm"], r6(3.4));
    }

    #[test]
    fn lifting_line_sections_omit_induction_fields() {        let mut p = synthetic_prop();
        p.param.lifting_line = true;
        let text = summary(&p, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");

        assert_eq!(y["design"]["mode"], "lifting-line");
        let sections = y["sections"].as_sequence().unwrap();
        for s in sections {
            assert!(s.get("dv_mps").is_none() || s["dv_mps"].is_null());
            assert!(s.get("a_prime").is_none() || s["a_prime"].is_null());
            assert!(s.get("thrust_n").is_none() || s["thrust_n"].is_null());
            assert!(s.get("torque_nm").is_none() || s["torque_nm"].is_null());
            // Geometry is always present.
            assert!(s.get("chord_m").is_some());
            assert!(s.get("twist_deg").is_some());
            // Station convergence is BEM-only.
            assert!(s.get("converged").is_none(), "converged leaked into lifting-line");
        }
    }

    #[test]
    fn warning_is_emitted_when_the_operating_point_was_not_reached() {
        let p = synthetic_prop();
        let w = "design torque 0.0500 Nm not achievable: the closest design absorbs 0.0000 Nm (100.0% of the target) at 0.00 N thrust; 0/4 BEM stations converged";
        let text = summary(&p, 5000.0, 0.0, 0.0, &motor_info(), Some(w));
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");
        assert_eq!(y["warning"], w);
        // The performance block still reports the closest design's numbers.
        assert_eq!(y["performance"]["thrust_n"], 0.0);
        assert_eq!(y["performance"]["torque_nm"], 0.0);
    }

    #[test]
    fn efficiencies_match_their_definitions() {
        let mut p = synthetic_prop();
        p.param.forward_airspeed = 0.0;
        // Hover (u_0 = 0): propulsive efficiency is zero, FoM positive.
        let text = summary(&p, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");
        assert_eq!(y["performance"]["propulsive_efficiency"], 0.0);
        assert!(y["performance"]["figure_of_merit"].as_f64().unwrap() > 0.0);

        // Forward flight: eta = T u_0 / P.
        let mut p2 = synthetic_prop();
        p2.param.forward_airspeed = 10.0;
        let text = summary(&p2, 5000.0, 2.25, 0.05, &motor_info(), None);
        let y: serde_yaml::Value = serde_yaml::from_str(&text).expect("valid yaml");
        let omega = rpm2omega(5000.0);
        let want = r6(2.25 * 10.0 / (0.05 * omega));
        assert_eq!(y["performance"]["propulsive_efficiency"], want);
    }
}
