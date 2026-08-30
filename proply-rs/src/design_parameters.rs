// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Propeller design parameters, parsed from a JSON file.
//!
//! Mirrors `proply/design_parameters.py`: the same JSON schema, with the
//! Python class attributes as defaults for missing keys.  `center_hole` is
//! optional and defaults to half the `hub_radius` (the Python used a fixed
//! 5 mm).  Physical quantities also accept unit-suffixed strings
//! (`"6 mm"`, `"6.8cm"`, `"500g"`); a bare number keeps its historical
//! unit — metres for lengths, newtons for thrust, millimetres for
//! `trailing_edge` (see [`crate::units`]).

use serde::Deserialize;

fn d(v: f64) -> f64 {
    v
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
#[allow(non_snake_case)] // field names match the JSON schema and the Python attributes
pub struct DesignParameters {
    pub name: String,
    #[serde(deserialize_with = "crate::units::de_length_m")]
    pub radius: f64,
    #[serde(deserialize_with = "crate::units::de_force_n")]
    pub thrust: f64,
    pub blades: usize,
    /// Mounting bore radius (m).  Absent in the JSON: half the hub radius
    /// (see [`DesignParameters::center_hole`]).
    #[serde(default, deserialize_with = "crate::units::de_opt_length_m")]
    pub center_hole: Option<f64>,
    #[serde(deserialize_with = "crate::units::de_length_m")]
    pub tip_chord: f64,
    #[serde(deserialize_with = "crate::units::de_length_m")]
    pub hub_radius: f64,
    #[serde(deserialize_with = "crate::units::de_length_m")]
    pub hub_depth: f64,
    /// Trailing-edge gap, in millimetres when bare (its historical unit —
    /// divided by 1000 at the use site).
    #[serde(deserialize_with = "crate::units::de_length_mm")]
    pub trailing_edge: f64,
    pub forward_airspeed: f64,
    pub altitude: f64,
    pub motor_volts: f64,
    pub motor_Kv: f64,
    pub motor_winding_resistance: f64,
    pub motor_no_load_current: f64,
    pub scimitar_percent: f64,
    /// Number of control points for the smooth (cubic-spline) blade chord
    /// distribution used by the lifting-line design (default 3).
    #[serde(default = "default_chord_spline_n")]
    pub chord_spline_n: usize,
    /// Foil camber (max camber as a fraction of chord, NACA 4-series `m`).
    /// When set, every station uses exactly this camber; when absent, the
    /// lifting-line design scans the [`crate::prop::CAMBER_CANDIDATES`] set
    /// plus a composed per-station distribution and keeps the
    /// best-performing one.
    #[serde(default)]
    pub camber: Option<f64>,
    /// Explicit motor operating point (engine-style): when *both* are given,
    /// the design runs at this torque (N m) and RPM instead of the electric
    /// motor model's maximum-efficiency point derived from `motor_Kv` & co.
    #[serde(default)]
    pub motor_torque: Option<f64>,
    #[serde(default)]
    pub motor_RPM: Option<f64>,
    /// Use the CST (Kulfan) foil family instead of the NACA 4-series: every
    /// station foil is the default 18-parameter section, re-thicknessed and
    /// cambered to the same radial laws.
    #[serde(default)]
    pub cst: bool,
    /// Use the ARA-D foil family (the table-driven propeller sections of
    /// the legacy proply, blended over the radial thickness law) instead of
    /// the NACA 4-series.  The family carries its own camber, so the
    /// `camber` setting/scan does not apply.  Mutually exclusive with `cst`.
    #[serde(default)]
    pub arad: bool,
    /// Use the mechanical blade-thickness law ([`crate::thickness`]):
    /// treat the blade as a cantilever beam anchored at the hub and size
    /// the section thickness from the deflection the thrust load causes,
    /// instead of the geometric power law.  The hub thickness
    /// (`hub_depth`) is deliberately not involved — it describes the hub
    /// mounting, not the blade's airfoil sections.  The law is sized on
    /// the converged design's station loads before the final design run.
    #[serde(default)]
    pub mech_thickness: bool,
    /// Elastic modulus of the blade material (Pa), the `E` of the
    /// mechanical-thickness beam model.  Bare number in pascals, or a
    /// unit-suffixed string (`"3 GPa"`, `"3000 MPa"`).  Defaults to
    /// 3 GPa, typical of a moulded nylon/ABS propeller.
    #[serde(
        default = "default_modulus",
        deserialize_with = "crate::units::de_pressure_pa"
    )]
    pub modulus: f64,
    /// Allowed tip deflection of the mechanical-thickness law as a
    /// fraction of the prop radius `R` (default 0.05 = 5% of R).
    #[serde(default = "default_deflection_fraction")]
    pub deflection_fraction: f64,
    /// Minimum section thickness of the mechanical-thickness law as a
    /// fraction of the local chord (default 0.06, the thinnest ARA-D
    /// table): the sized thickness never goes below this.
    #[serde(default = "default_thickness_floor")]
    pub thickness_floor: f64,
    // ---- run / design options (mirror the CLI flags so a JSON file can
    // carry the whole design) ----
    pub bem: bool,
    pub lifting_line: bool,
    /// Implied since the operating-point match: the design loop always
    /// iterates the thrust target until the blade absorbs the design
    /// torque at the design RPM.  Kept (parsed but unused) so old JSON
    /// files and scripts passing `--auto` keep working.
    #[serde(default)]
    pub auto: bool,
    pub resolution: usize,
    pub n: usize,
    pub ar: Option<f64>,
    pub plate: bool,
    pub dir: String,
    pub step_file: String,
}

fn default_chord_spline_n() -> usize {
    3
}

fn default_modulus() -> f64 {
    3.0e9
}

fn default_deflection_fraction() -> f64 {
    0.05
}

fn default_thickness_floor() -> f64 {
    0.06
}

impl Default for DesignParameters {
    fn default() -> Self {
        Self {
            name: "hello world".into(),
            radius: d(0.0625),
            thrust: d(2.0),
            blades: 2,
            center_hole: None,
            tip_chord: d(7.0 / 1000.0),
            hub_radius: d(5.0 / 1000.0),
            hub_depth: d(0.0),
            // A millimetre field: 0.5 mm (the old 0.5/1000.0 value was
            // divided by 1000 again at the use site, a 0.5 um edge).
            trailing_edge: d(0.5),
            forward_airspeed: d(1.0),
            altitude: d(0.0),
            motor_volts: d(11.0),
            motor_Kv: d(1200.0),
            motor_winding_resistance: d(0.206),
            motor_no_load_current: d(0.5),
            scimitar_percent: d(0.0),
            chord_spline_n: 3,
            camber: None,
            motor_torque: None,
            motor_RPM: None,
            cst: false,
            arad: false,
            mech_thickness: false,
            modulus: 3.0e9,
            deflection_fraction: 0.05,
            thickness_floor: 0.06,
            bem: true,
            lifting_line: false,
            auto: false,
            resolution: 40,
            n: 40,
            ar: None,
            plate: false,
            dir: ".".into(),
            step_file: String::new(),
        }
    }
}

impl DesignParameters {
    /// The mounting bore radius (m): the JSON `center_hole` when given,
    /// else half the hub radius.
    pub fn center_hole(&self) -> f64 {
        self.center_hole.unwrap_or(0.5 * self.hub_radius)
    }

    /// Load from a JSON file, exactly like `DesignParameters(filename)`.
    pub fn from_file(filename: &str) -> Result<Self, String> {
        let data = std::fs::read_to_string(filename)
            .map_err(|e| format!("Cannot read {}: {}", filename, e))?;
        Self::from_json(&data)
    }

    /// Parse from a JSON string, mirroring `from_json`.
    pub fn from_json(data: &str) -> Result<Self, String> {
        serde_json::from_str(data).map_err(|e| format!("Invalid design parameters: {}", e))
    }

    /// Disc area (kept for parity with the Python `area()` method).
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius
    }

    /// The design's motor operating point: `(torque N m, RPM, power W)`.
    /// An explicitly specified `motor_torque` + `motor_RPM` pair (an
    /// engine with a known operating point) overrides the electric motor
    /// model's maximum-efficiency derivation; the power is then the shaft
    /// power at the specified point.  One without the other is ignored.
    pub fn motor_operating_point(&self) -> (f64, f64, f64) {
        match (self.motor_torque, self.motor_RPM) {
            (Some(q), Some(rpm)) => (q, rpm, q * crate::optimize::rpm2omega(rpm)),
            _ => {
                let m = crate::motor::Motor::new(
                    self.motor_Kv,
                    self.motor_no_load_current,
                    self.motor_winding_resistance,
                );
                let (q, rpm) = m.get_qmax(self.motor_volts);
                (q, rpm, m.get_pmax(self.motor_volts))
            }
        }
    }
}

impl std::fmt::Display for DesignParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Design Parameters: r={:5.3}, u_0={}",
            self.radius, self.forward_airspeed
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_example_json() {
        let json = r#"{
            "name": "test",
            "altitude": 0.0,
            "forward_airspeed": 0.0,
            "motor_Kv": 11500,
            "motor_volts": 3.7,
            "motor_no_load_current": 0.075,
            "motor_winding_resistance": 0.035,
            "thrust": 0.5,
            "blades": 3,
            "radius": 0.02,
            "tip_chord": 0.003,
            "center_hole": 1.5e-3,
            "scimitar_percent": -5.0,
            "trailing_edge": 0.25,
            "hub_radius": 0.004,
            "hub_depth": 0.003
        }"#;
        let p = DesignParameters::from_json(json).unwrap();
        assert_eq!(p.name, "test");
        assert_eq!(p.blades, 3);
        assert!((p.radius - 0.02).abs() < 1e-12);
        assert!((p.center_hole() - 1.5e-3).abs() < 1e-12);
        assert!((p.motor_Kv - 11500.0).abs() < 1e-9);
        assert!((p.scimitar_percent + 5.0).abs() < 1e-12);
    }

    #[test]
    fn cli_options_roundtrip_via_json() {
        // Every CLI option must be settable in the JSON design file and read
        // back with defaults when absent.
        let json = r#"{
            "name": "x",
            "radius": 0.05,
            "thrust": 1.0,
            "blades": 2,
            "bem": false,
            "lifting_line": true,
            "auto": true,
            "resolution": 60,
            "n": 24,
            "ar": 4.5,
            "plate": true,
            "dir": "out",
            "step_file": "x.step",
            "chord_spline_n": 5,
            "camber": 0.03,
            "cst": true,
            "arad": true
        }"#;
        let p = DesignParameters::from_json(json).unwrap();
        assert!(!p.bem);
        assert!(p.lifting_line);
        assert!(p.auto);
        assert_eq!(p.resolution, 60);
        assert_eq!(p.n, 24);
        assert_eq!(p.ar, Some(4.5));
        assert!(p.plate);
        assert_eq!(p.dir, "out");
        assert_eq!(p.step_file, "x.step");
        assert_eq!(p.chord_spline_n, 5);
        assert!((p.camber.unwrap() - 0.03).abs() < 1e-12);
        assert!(p.cst);
        assert!(p.arad);

        // Absent keys fall back to the defaults.
        let p2 = DesignParameters::from_json(
            r#"{
            "name": "y", "radius": 0.05, "thrust": 1.0, "blades": 2
        }"#,
        )
        .unwrap();
        assert!(p2.bem);
        assert!(!p2.lifting_line);
        assert!(!p2.auto);
        assert_eq!(p2.resolution, 40);
        assert_eq!(p2.n, 40);
        assert_eq!(p2.ar, None);
        assert!(!p2.plate);
        assert_eq!(p2.chord_spline_n, 3);
        assert!(p2.camber.is_none());
        assert!(!p2.cst);
        assert!(!p2.arad);
    }

    #[test]
    fn specified_motor_operating_point_overrides_electric_model() {
        // Kv=1900, I0=0.5, Rm=0.405, V=11 -> the electric max-efficiency
        // point (matches the motor model test).
        let base = r#"{
            "name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
            "motor_Kv": 1900, "motor_volts": 11.0,
            "motor_no_load_current": 0.5, "motor_winding_resistance": 0.405
        }"#;
        let p = DesignParameters::from_json(base).unwrap();
        assert!(p.motor_torque.is_none() && p.motor_RPM.is_none());
        let (q, rpm, _) = p.motor_operating_point();
        assert!((q - 0.016006).abs() < 1.0e-3, "q = {}", q);
        assert!((rpm - 18064.6).abs() < 5.0, "rpm = {}", rpm);

        // Both given: the specified point wins, power = shaft power there.
        let specified = r#"{
            "name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
            "motor_Kv": 1900, "motor_volts": 11.0,
            "motor_no_load_current": 0.5, "motor_winding_resistance": 0.405,
            "motor_torque": 334.8, "motor_RPM": 1950.0
        }"#;
        let p = DesignParameters::from_json(specified).unwrap();
        let (q, rpm, power) = p.motor_operating_point();
        assert!((q - 334.8).abs() < 1.0e-9);
        assert!((rpm - 1950.0).abs() < 1.0e-9);
        assert!(
            (power - 334.8 * crate::optimize::rpm2omega(1950.0)).abs() < 1.0e-6,
            "power = {}",
            power
        );

        // One without the other: ignored, the electric model is used.
        let half = r#"{
            "name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
            "motor_Kv": 1900, "motor_volts": 11.0,
            "motor_no_load_current": 0.5, "motor_winding_resistance": 0.405,
            "motor_torque": 334.8
        }"#;
        let p = DesignParameters::from_json(half).unwrap();
        let (q, rpm, _) = p.motor_operating_point();
        assert!((q - 0.016006).abs() < 1.0e-3, "q = {}", q);
        assert!((rpm - 18064.6).abs() < 5.0, "rpm = {}", rpm);
    }

    #[test]
    fn mechanical_thickness_parameters_parse_and_default() {
        // Set everything explicitly, with modulus as a unit-suffixed string.
        let json = r#"{
            "name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
            "mech_thickness": true,
            "modulus": "4 GPa",
            "deflection_fraction": 0.02,
            "thickness_floor": 0.10
        }"#;
        let p = DesignParameters::from_json(json).unwrap();
        assert!(p.mech_thickness);
        assert!((p.modulus - 4.0e9).abs() < 1e-6, "modulus {}", p.modulus);
        assert!((p.deflection_fraction - 0.02).abs() < 1e-12);
        assert!((p.thickness_floor - 0.10).abs() < 1e-12);

        // A bare number for modulus means pascals.
        let json2 = r#"{
            "name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
            "modulus": 3.0e9
        }"#;
        let p2 = DesignParameters::from_json(json2).unwrap();
        assert!(!p2.mech_thickness);
        assert!((p2.modulus - 3.0e9).abs() < 1e-6);

        // Absent keys fall back to the defaults.
        let p3 = DesignParameters::default();
        assert!(!p3.mech_thickness);
        assert!((p3.modulus - 3.0e9).abs() < 1e-6);
        assert!((p3.deflection_fraction - 0.05).abs() < 1e-12);
        assert!((p3.thickness_floor - 0.06).abs() < 1e-12);

        // A force unit is not a modulus.
        let err = DesignParameters::from_json(
            r#"{"name": "x", "radius": 0.05, "thrust": 1.0, "blades": 2,
                "modulus": "500g"}"#,
        )
        .unwrap_err();
        assert!(err.contains("Pa, kPa, MPa, GPa"), "{}", err);
    }

    #[test]
    fn missing_center_hole_defaults() {
        let json = r#"{
            "name": "x",
            "thrust": 1.0,
            "blades": 2,
            "radius": 0.05,
            "tip_chord": 0.003,
            "hub_radius": 0.004,
            "hub_depth": 0.003,
            "scimitar_percent": 0.0,
            "trailing_edge": 0.25,
            "forward_airspeed": 0.0,
            "altitude": 0.0,
            "motor_Kv": 1900.0,
            "motor_volts": 11.0,
            "motor_winding_resistance": 0.4,
            "motor_no_load_current": 0.5
        }"#;
        let p = DesignParameters::from_json(json).unwrap();
        assert!(
            (p.center_hole() - 0.5 * 0.004).abs() < 1e-12,
            "bore {} should be half the hub radius",
            p.center_hole()
        );
        // The Default (hub_radius 5 mm) follows the same rule.
        assert!(
            (DesignParameters::default().center_hole() - 0.0025).abs() < 1e-12,
            "default bore {}",
            DesignParameters::default().center_hole()
        );
    }

    #[test]
    fn unit_suffixed_json_matches_bare_numbers() {
        // The same design written with unit suffixes must parse to the
        // same parameters as the bare-number (SI) form.  trailing_edge is
        // the exception: its bare unit is millimetres, so 0.25 == "0.25mm".
        let base = r#"{
            "name": "units", "blades": 2,
            "radius": 0.068, "tip_chord": 0.005, "hub_radius": 0.006,
            "hub_depth": 0.006, "center_hole": 0.0015, "trailing_edge": 0.25,
            "thrust": 4.0
        }"#;
        let suffixed = r#"{
            "name": "units", "blades": 2,
            "radius": "6.8cm", "tip_chord": "5mm", "hub_radius": "6 mm",
            "hub_depth": "6mm", "center_hole": "1.5mm", "trailing_edge": "0.25mm",
            "thrust": "4N"
        }"#;
        let a = DesignParameters::from_json(base).unwrap();
        let b = DesignParameters::from_json(suffixed).unwrap();
        for (what, x, y) in [
            ("radius", a.radius, b.radius),
            ("tip_chord", a.tip_chord, b.tip_chord),
            ("hub_radius", a.hub_radius, b.hub_radius),
            ("hub_depth", a.hub_depth, b.hub_depth),
            ("thrust", a.thrust, b.thrust),
        ] {
            assert!((x - y).abs() < 1e-15, "{}: {} vs {}", what, x, y);
        }
        assert_eq!(a.center_hole, b.center_hole);
        assert!((a.trailing_edge - b.trailing_edge).abs() < 1e-15);
        assert!((a.trailing_edge - 0.25).abs() < 1e-12, "bare mm field");
    }

    #[test]
    fn thrust_in_kilograms_and_grams() {
        let json = |t: &str| {
            format!(r#"{{"name": "x", "radius": 0.05, "blades": 2, "thrust": "{t}"}}"#)
        };
        let kg = DesignParameters::from_json(&json("0.5kg")).unwrap();
        let g = DesignParameters::from_json(&json("500g")).unwrap();
        assert!(
            (kg.thrust - g.thrust).abs() < 1e-9,
            "kg {} vs g {}",
            kg.thrust,
            g.thrust
        );
        assert!((kg.thrust - 0.5 * 9.80665).abs() < 1e-9, "kgf {}", kg.thrust);
    }

    #[test]
    fn wrong_unit_kind_errors() {
        let err = DesignParameters::from_json(
            r#"{"name": "x", "radius": "500g", "blades": 2, "thrust": 1.0}"#,
        )
        .unwrap_err();
        assert!(err.contains("m, cm, mm"), "{}", err);
        let err = DesignParameters::from_json(
            r#"{"name": "x", "radius": 0.05, "blades": 2, "thrust": "5 furlongs"}"#,
        )
        .unwrap_err();
        assert!(err.contains("N, kg, g"), "{}", err);
    }
}
