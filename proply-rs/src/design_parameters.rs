//! Propeller design parameters, parsed from a JSON file.
//!
//! Mirrors `proply/design_parameters.py`: the same JSON schema, with the
//! Python class attributes as defaults for missing keys (`center_hole` was
//! already optional in the Python version).

use serde::Deserialize;

fn d(v: f64) -> f64 {
    v
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
#[allow(non_snake_case)] // field names match the JSON schema and the Python attributes
pub struct DesignParameters {
    pub name: String,
    pub radius: f64,
    pub thrust: f64,
    pub blades: usize,
    pub center_hole: f64,
    pub tip_chord: f64,
    pub hub_radius: f64,
    pub hub_depth: f64,
    pub trailing_edge: f64,
    pub forward_airspeed: f64,
    pub altitude: f64,
    pub motor_volts: f64,
    pub motor_Kv: f64,
    pub motor_winding_resistance: f64,
    pub motor_no_load_current: f64,
    pub scimitar_percent: f64,
}

impl Default for DesignParameters {
    fn default() -> Self {
        Self {
            name: "hello world".into(),
            radius: d(0.0625),
            thrust: d(2.0),
            blades: 2,
            center_hole: d(5.0 / 1000.0),
            tip_chord: d(7.0 / 1000.0),
            hub_radius: d(5.0 / 1000.0),
            hub_depth: d(0.0),
            trailing_edge: d(0.5 / 1000.0),
            forward_airspeed: d(1.0),
            altitude: d(0.0),
            motor_volts: d(11.0),
            motor_Kv: d(1200.0),
            motor_winding_resistance: d(0.206),
            motor_no_load_current: d(0.5),
            scimitar_percent: d(0.0),
        }
    }
}

impl DesignParameters {
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
        3.1415 * self.radius * self.radius
    }
}

impl std::fmt::Display for DesignParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Design Parameters: r={:5.3}, u_0={}", self.radius, self.forward_airspeed)
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
        assert!((p.center_hole - 1.5e-3).abs() < 1e-12);
        assert!((p.motor_Kv - 11500.0).abs() < 1e-9);
        assert!((p.scimitar_percent + 5.0).abs() < 1e-12);
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
        assert!((p.center_hole - 0.005).abs() < 1e-12);
    }
}
