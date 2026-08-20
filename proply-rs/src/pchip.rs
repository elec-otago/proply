// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! PCHIP — piecewise cubic Hermite interpolating polynomial with monotone
//! slopes (Fritsch & Carlson 1980), matching `scipy.interpolate.PchipInterpolator`.
//!
//! Used for the scimitar offset and the chord smoothing in the design loop.

#[derive(Debug, Clone)]
pub struct Pchip {
    x: Vec<f64>,
    y: Vec<f64>,
    d: Vec<f64>,
}

impl Pchip {
    /// Build the interpolant through (x, y); `x` must be strictly increasing.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        assert_eq!(x.len(), y.len(), "x and y must have equal length");
        assert!(x.len() >= 2, "at least 2 points required");
        let n = x.len();

        let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
        let delta: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();

        // Two points collapse to a single linear segment.
        if n == 2 {
            return Self {
                x: x.to_vec(),
                y: y.to_vec(),
                d: vec![delta[0], delta[0]],
            };
        }

        let mut d = vec![0.0; n];
        // Interior points
        for i in 1..n - 1 {
            if delta[i - 1] * delta[i] <= 0.0 {
                d[i] = 0.0;
            } else {
                let w1 = 2.0 * h[i] + h[i - 1];
                let w2 = h[i] + 2.0 * h[i - 1];
                d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
            }
        }
        // Endpoints (one-sided, scipy's non-centered scheme)
        d[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
        if d[0].signum() != delta[0].signum() {
            d[0] = 0.0;
        } else if delta[0].signum() != delta[1].signum() && d[0].abs() > (3.0 * delta[0]).abs() {
            d[0] = 3.0 * delta[0];
        }
        d[n - 1] = ((2.0 * h[n - 2] + h[n - 3]) * delta[n - 2] - h[n - 2] * delta[n - 3])
            / (h[n - 2] + h[n - 3]);
        if d[n - 1].signum() != delta[n - 2].signum() {
            d[n - 1] = 0.0;
        } else if delta[n - 2].signum() != delta[n - 3].signum()
            && d[n - 1].abs() > (3.0 * delta[n - 2]).abs()
        {
            d[n - 1] = 3.0 * delta[n - 2];
        }

        Self {
            x: x.to_vec(),
            y: y.to_vec(),
            d,
        }
    }

    /// Evaluate at `t`.  Out-of-range `t` is clamped to the domain end
    /// (scipy raises; proply never evaluates out of range).
    pub fn eval(&self, t: f64) -> f64 {
        let n = self.x.len();
        if t <= self.x[0] {
            return self.y[0];
        }
        if t >= self.x[n - 1] {
            return self.y[n - 1];
        }
        // Locate the segment (x is strictly increasing).
        let mut i = 0;
        while self.x[i + 1] < t {
            i += 1;
        }
        let h = self.x[i + 1] - self.x[i];
        let s = (t - self.x[i]) / h;
        let h00 = 2.0 * s.powi(3) - 3.0 * s.powi(2) + 1.0;
        let h10 = s.powi(3) - 2.0 * s.powi(2) + s;
        let h01 = -2.0 * s.powi(3) + 3.0 * s.powi(2);
        let h11 = s.powi(3) - s.powi(2);
        h00 * self.y[i] + h10 * h * self.d[i] + h01 * self.y[i + 1] + h11 * h * self.d[i + 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_points_are_a_linear_segment() {
        let x = vec![0.0, 1.0];
        let y = vec![1.0, 3.0];
        let p = Pchip::new(&x, &y);
        for (xi, yi) in x.iter().zip(y.iter()) {
            assert!((p.eval(*xi) - yi).abs() < 1e-12);
        }
        // Linear interpolation at the midpoint (and elsewhere).
        assert!((p.eval(0.5) - 2.0).abs() < 1e-12);
        assert!((p.eval(0.25) - 1.5).abs() < 1e-12);
        // Clamped outside the domain.
        assert!((p.eval(-1.0) - 1.0).abs() < 1e-12);
        assert!((p.eval(2.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn monotone_data_is_preserved() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 1.0, 1.5, 1.6, 2.0];
        let p = Pchip::new(&x, &y);
        // Knots are interpolated exactly.
        for (xi, yi) in x.iter().zip(y.iter()) {
            assert!((p.eval(*xi) - yi).abs() < 1e-12);
        }
        // Monotonicity: sample densely and check no overshoot.
        let mut prev = p.eval(0.0);
        let mut t = 0.01;
        while t < 4.0 {
            let v = p.eval(t);
            assert!(v >= prev - 1e-12, "not monotone at {}", t);
            prev = v;
            t += 0.01;
        }
    }

    #[test]
    fn matches_scipy_classic_case() {
        // Hand-computed PCHIP through (0,0), (1,sin1), (2,sin2) at x=0.1:
        // interior slope d1 = 0.12555, endpoint d0 = 1.2283 -> y = 0.1219.
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0_f64.sin(), 2.0_f64.sin()];
        let p = Pchip::new(&x, &y);
        let v = p.eval(0.1);
        assert!((v - 0.121923).abs() < 1e-6, "pchip(0.1) = {}", v);
    }

    #[test]
    fn non_monotone_flat_spot() {
        // Local extremum forces a zero slope (PCHIP property).
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, 1.0, 0.0];
        let p = Pchip::new(&x, &y);
        assert!((p.eval(1.0) - 1.0).abs() < 1e-12);
        assert!((p.eval(2.0) - 1.0).abs() < 1e-12);
        // between the plateau the interpolant stays within [0,1]
        let mut t = 0.0;
        while t <= 3.0 {
            let v = p.eval(t);
            assert!(
                (-1e-12..=1.0 + 1e-12).contains(&v),
                "overshoot at {}: {}",
                t,
                v
            );
            t += 0.05;
        }
    }
}
