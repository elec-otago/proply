// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Unit-suffixed quantities in the design JSON.
//!
//! A quantity field accepts a bare JSON number — meaning exactly what it
//! has always meant (metres for lengths, newtons for thrust, millimetres
//! for `trailing_edge`, pascals for `modulus`) — or a quoted string
//! carrying a unit suffix: `"5mm"`, `"6 mm"`, `"6.8cm"`, `"500g"`,
//! `"0.5kg"`, `"3.2N"`, `"3 GPa"`.  JSON has no unquoted `5mm` scalar, so
//! suffixed values must be strings.
//!
//! Length units: `m`, `cm`, `mm`.  Force units: `N`, `kg` (kilogram-force,
//! g0 = 9.80665 N), `g` (gram-force).  Pressure units: `Pa`, `kPa`, `MPa`,
//! `GPa`.  Units match case-insensitively and may be separated from the
//! number by whitespace.  A string without a unit uses the field's
//! bare-number unit, and a unit from the wrong kind (`"tip_chord":
//! "500g"`) is an error.

use serde::Deserialize;
use serde::Deserializer;

/// Standard gravity (m/s^2): the kgf/gf -> N conversion.
const G0: f64 = 9.80665;

/// Split a quantity string into (number, unit): the unit is the trailing
/// run of ASCII letters (a leading scan would cut exponent forms like
/// `"2e-3m"` at the 'e').  `"1.5 mm"` -> (1.5, "mm"); `"5"` -> (5.0, None).
fn split_quantity(text: &str) -> Result<(f64, Option<&str>), String> {
    let t = text.trim();
    // The unit run is all-ASCII, so its byte length equals its char count.
    let unit_len = t
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_alphabetic())
        .count();
    let (num, unit) = t.split_at(t.len() - unit_len);
    let value: f64 = num
        .trim()
        .parse()
        .map_err(|_| format!("invalid quantity '{text}' (cannot read a number)"))?;
    let unit = unit.trim();
    if unit.is_empty() {
        Ok((value, None))
    } else {
        Ok((value, Some(unit)))
    }
}

/// Scale factor of a length unit to metres.
fn length_scale(unit: &str) -> Option<f64> {
    match unit.to_ascii_lowercase().as_str() {
        "m" => Some(1.0),
        "cm" => Some(1.0e-2),
        "mm" => Some(1.0e-3),
        _ => None,
    }
}

/// Scale factor of a force unit to newtons.
fn force_scale(unit: &str) -> Option<f64> {
    match unit.to_ascii_lowercase().as_str() {
        "n" => Some(1.0),
        "kg" => Some(G0),
        "g" => Some(G0 * 1.0e-3),
        _ => None,
    }
}

/// Scale factor of a pressure unit to pascals.
fn pressure_scale(unit: &str) -> Option<f64> {
    match unit.to_ascii_lowercase().as_str() {
        "pa" => Some(1.0),
        "kpa" => Some(1.0e3),
        "mpa" => Some(1.0e6),
        "gpa" => Some(1.0e9),
        _ => None,
    }
}

/// Parse a length for a field whose storage unit is `storage_m` metres
/// (1.0 for metre fields, 1e-3 for the millimetre `trailing_edge`).  A
/// bare number or unitless string is already in the storage unit and
/// passes through; a suffixed string converts into the storage unit.
fn length(value: &serde_json::Value, storage_m: f64) -> Result<f64, String> {
    match value {
        serde_json::Value::Number(n) => Ok(n.as_f64().unwrap_or(f64::NAN)),
        serde_json::Value::String(s) => {
            let (v, unit) = split_quantity(s)?;
            match unit {
                None => Ok(v),
                Some(u) => length_scale(u)
                    .map(|scale| v * scale / storage_m)
                    .ok_or_else(|| format!("invalid quantity '{s}' for a length (expected m, cm, mm)")),
            }
        }
        other => Err(format!(
            "invalid quantity '{other}' (expected a number or a unit-suffixed string)"
        )),
    }
}

/// Parse a force whose bare unit is newtons.
fn force(value: &serde_json::Value) -> Result<f64, String> {
    match value {
        serde_json::Value::Number(n) => Ok(n.as_f64().unwrap_or(f64::NAN)),
        serde_json::Value::String(s) => {
            let (v, unit) = split_quantity(s)?;
            match unit {
                None => Ok(v),
                Some(u) => force_scale(u)
                    .map(|scale| v * scale)
                    .ok_or_else(|| format!("invalid quantity '{s}' for a force (expected N, kg, g)")),
            }
        }
        other => Err(format!(
            "invalid quantity '{other}' (expected a number or a unit-suffixed string)"
        )),
    }
}

/// Parse a pressure whose bare unit is pascals (the `modulus` field).
fn pressure(value: &serde_json::Value) -> Result<f64, String> {
    match value {
        serde_json::Value::Number(n) => Ok(n.as_f64().unwrap_or(f64::NAN)),
        serde_json::Value::String(s) => {
            let (v, unit) = split_quantity(s)?;
            match unit {
                None => Ok(v),
                Some(u) => pressure_scale(u)
                    .map(|scale| v * scale)
                    .ok_or_else(|| {
                        format!("invalid quantity '{s}' for a pressure (expected Pa, kPa, MPa, GPa)")
                    }),
            }
        }
        other => Err(format!(
            "invalid quantity '{other}' (expected a number or a unit-suffixed string)"
        )),
    }
}

/// `deserialize_with` helper: a length in metres (bare number = metres).
pub fn de_length_m<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let v = serde_json::Value::deserialize(d)?;
    length(&v, 1.0).map_err(serde::de::Error::custom)
}

/// `deserialize_with` helper: a length in millimetres (bare number = mm,
/// the historical `trailing_edge` unit).
pub fn de_length_mm<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let v = serde_json::Value::deserialize(d)?;
    length(&v, 1.0e-3).map_err(serde::de::Error::custom)
}

/// `deserialize_with` helper: a force in newtons (bare number = N).
pub fn de_force_n<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let v = serde_json::Value::deserialize(d)?;
    force(&v).map_err(serde::de::Error::custom)
}

/// `deserialize_with` helper: a pressure in pascals (bare number = Pa).
pub fn de_pressure_pa<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let v = serde_json::Value::deserialize(d)?;
    pressure(&v).map_err(serde::de::Error::custom)
}

/// Parse a pressure value from a CLI argument (number or unit-suffixed
/// string), for options like `--modulus`.
pub fn parse_pressure(text: &str) -> Result<f64, String> {
    let (v, unit) = split_quantity(text)?;
    match unit {
        None => Ok(v),
        Some(u) => pressure_scale(u)
            .map(|scale| v * scale)
            .ok_or_else(|| format!("invalid quantity '{text}' for a pressure (expected Pa, kPa, MPa, GPa)")),
    }
}

/// `deserialize_with` helper: an optional length in metres (missing or
/// null = None, as for `center_hole`).
pub fn de_opt_length_m<'de, D>(d: D) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    match Option::<serde_json::Value>::deserialize(d)? {
        None | Some(serde_json::Value::Null) => Ok(None),
        Some(v) => length(&v, 1.0).map(Some).map_err(serde::de::Error::custom),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Lengths with a metre bare unit.
    #[test]
    fn lengths_in_metres() {
        assert_eq!(length(&serde_json::json!(0.068), 1.0).unwrap(), 0.068);
        assert_eq!(length(&serde_json::json!("0.068"), 1.0).unwrap(), 0.068);
        assert_eq!(length(&serde_json::json!("6.8cm"), 1.0).unwrap(), 0.068);
        assert_eq!(length(&serde_json::json!("68 mm"), 1.0).unwrap(), 0.068);
        assert_eq!(length(&serde_json::json!("5mm"), 1.0).unwrap(), 0.005);
        assert_eq!(length(&serde_json::json!("2e-3m"), 1.0).unwrap(), 0.002);
        assert_eq!(length(&serde_json::json!("1.5 MM"), 1.0).unwrap(), 0.0015);
    }

    /// The millimetre bare unit (trailing_edge).
    #[test]
    fn lengths_in_millimetres() {
        assert_eq!(length(&serde_json::json!(0.25), 1.0e-3).unwrap(), 0.25);
        assert_eq!(length(&serde_json::json!("0.25"), 1.0e-3).unwrap(), 0.25);
        assert_eq!(length(&serde_json::json!("0.25mm"), 1.0e-3).unwrap(), 0.25);
        assert_eq!(length(&serde_json::json!("0.025cm"), 1.0e-3).unwrap(), 0.25);
        assert_eq!(length(&serde_json::json!("1m"), 1.0e-3).unwrap(), 1000.0);
    }

    #[test]
    fn forces_in_newtons() {
        assert_eq!(force(&serde_json::json!(4.0)).unwrap(), 4.0);
        assert_eq!(force(&serde_json::json!("4N")).unwrap(), 4.0);
        assert_eq!(force(&serde_json::json!("4 n")).unwrap(), 4.0);
        assert!((force(&serde_json::json!("0.5kg")).unwrap() - 0.5 * G0).abs() < 1e-12);
        assert!((force(&serde_json::json!("500g")).unwrap() - 0.5 * G0).abs() < 1e-9);
        assert!((force(&serde_json::json!("2KG")).unwrap() - 2.0 * G0).abs() < 1e-12);
    }

    #[test]
    fn pressures_in_pascals() {
        assert_eq!(pressure(&serde_json::json!(3.0e9)).unwrap(), 3.0e9);
        assert_eq!(pressure(&serde_json::json!("3e9")).unwrap(), 3.0e9);
        assert_eq!(pressure(&serde_json::json!("3GPa")).unwrap(), 3.0e9);
        assert_eq!(pressure(&serde_json::json!("3 GPa")).unwrap(), 3.0e9);
        assert_eq!(pressure(&serde_json::json!("3000 MPa")).unwrap(), 3.0e9);
        assert_eq!(pressure(&serde_json::json!("2.5kpa")).unwrap(), 2500.0);
        assert_eq!(pressure(&serde_json::json!("1Pa")).unwrap(), 1.0);
        // A force unit is not a pressure.
        let err = pressure(&serde_json::json!("500g")).unwrap_err();
        assert!(err.contains("Pa, kPa, MPa, GPa"), "{}", err);
        // The CLI parser follows the same rules.
        assert_eq!(parse_pressure("3 GPa").unwrap(), 3.0e9);
        assert_eq!(parse_pressure("3e9").unwrap(), 3.0e9);
        assert!(parse_pressure("5 furlongs").is_err());
    }

    #[test]
    fn bad_quantities_error() {
        // A force unit on a length field.
        let err = length(&serde_json::json!("500g"), 1.0).unwrap_err();
        assert!(err.contains("500g") && err.contains("m, cm, mm"), "{}", err);
        // A length unit on a force field.
        let err = force(&serde_json::json!("5mm")).unwrap_err();
        assert!(err.contains("N, kg, g"), "{}", err);
        // Unknown units.
        assert!(force(&serde_json::json!("5 furlongs")).is_err());
        // Non-numeric value.
        assert!(length(&serde_json::json!("amm"), 1.0).is_err());
        // Wrong JSON type.
        assert!(length(&serde_json::json!(true), 1.0).is_err());
        assert!(force(&serde_json::json!([1.0])).is_err());
    }
}
