// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Motor model: optimum torque and RPM at maximum efficiency.
//!
//! Ported 1:1 from `proply/motor_model.py`.  `Kv` is in RPM per volt,
//! `Rm` the winding resistance in ohms, `I0` the no-load current in amps.

/// Torque constant derived from `Kv`: `Kq = 30 / (pi * Kv)`.
fn torque_constant(kv: f64) -> f64 {
    30.0 / (std::f64::consts::PI * kv)
}

pub struct Motor {
    pub kv: f64,
    pub i0: f64,
    pub rm: f64,
    pub kq: f64,
}

impl Motor {
    pub fn new(kv: f64, i0: f64, rm: f64) -> Self {
        Self {
            kv,
            i0,
            rm,
            kq: torque_constant(kv),
        }
    }

    /// Torque (N m) at current `I`.
    pub fn get_torque(&self, i: f64) -> f64 {
        self.kq * (i - self.i0)
    }

    /// RPM at torque `q_in`.
    pub fn get_rpm(&self, q_in: f64) -> f64 {
        std::f64::consts::PI * self.kv.powi(2) * q_in / 30.0
    }

    /// Efficiency at voltage `V` and current `I`.
    pub fn get_efficiency(&self, v: f64, i: f64) -> f64 {
        (i - self.i0) * (-i * self.rm + v) / (i * v)
    }

    /// Current at maximum efficiency (amps).
    pub fn get_imax(&self, v: f64) -> f64 {
        (v * self.i0 / self.rm).sqrt()
    }

    /// (Torque, RPM) at maximum efficiency.
    pub fn get_qmax(&self, v: f64) -> (f64, f64) {
        let imax = self.get_imax(v);
        let qmax = self.kq * (imax - self.i0);
        let rpm = self.kv * (v - imax * self.rm);
        (qmax, rpm)
    }

    /// Power (watts) at maximum efficiency.
    pub fn get_pmax(&self, v: f64) -> f64 {
        let (qmax, rpm) = self.get_qmax(v);
        2.0 * std::f64::consts::PI * qmax * (rpm / 60.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_python_motor_model() {
        // Kv=1900, I0=0.5, Rm=0.405, V=11.  Hand-computed; confirmed against
        // the Python Motor class by the golden tests.
        let m = Motor::new(1900.0, 0.5, 0.405);
        assert!((m.get_imax(11.0) - 3.6847).abs() < 1e-3, "Imax = {}", m.get_imax(11.0));
        let (q, rpm) = m.get_qmax(11.0);
        assert!((q - 0.016006).abs() < 1e-4, "Qmax = {}", q);
        assert!((rpm - 18064.6).abs() < 5.0, "RPMmax = {}", rpm);
        let p = m.get_pmax(11.0);
        assert!((p - 30.28).abs() < 0.1, "Pmax = {}", p);
    }
}
