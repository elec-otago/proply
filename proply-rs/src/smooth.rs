// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Windowed smoothing of a 1-D signal, ported from `proply/smooth.py`
//! (which in turn is the classic numpy cookbook recipe).
//!
//! The signal is extended by reflected copies, convolved with a normalized
//! window, and trimmed back to the input length.

/// Smooth `x` with a window of `window_len` samples (must be odd) using the
/// named window (`"hanning"` or `"flat"`).
///
/// Panics if `x.len() < window_len` or the window name is unknown, matching
/// the Python `ValueError` behaviour.
pub fn smooth(x: &[f64], window_len: usize, window: &str) -> Vec<f64> {
    assert!(x.len() >= window_len, "Input vector needs to be bigger than window size");
    if window_len < 3 {
        return x.to_vec();
    }
    let w: Vec<f64> = match window {
        "flat" => vec![1.0; window_len],
        "hanning" => (0..window_len)
            .map(|i| 0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / (window_len - 1) as f64).cos()))
            .collect(),
        _ => panic!("Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"),
    };
    let wsum: f64 = w.iter().sum();
    let w: Vec<f64> = w.iter().map(|v| v / wsum).collect();

    // Reflected extension: x[w-1:0:-1] ++ x ++ x[-1:-w:-1]
    let mut s: Vec<f64> = Vec::with_capacity(x.len() + 2 * (window_len - 1));
    for i in (1..window_len).rev() {
        s.push(x[i]);
    }
    s.extend_from_slice(x);
    for k in 1..window_len {
        s.push(x[x.len() - k]);
    }

    // np.convolve(w, s, 'valid')
    let n = s.len() - w.len() + 1;
    let mut y = vec![0.0; n];
    for i in 0..n {
        for (j, wj) in w.iter().enumerate() {
            y[i] += wj * s[i + j];
        }
    }

    // y[w//2 : -(w//2)]
    y[window_len / 2..y.len() - window_len / 2].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_window_is_moving_average() {
        let x = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        let y = smooth(&x, 3, "flat");
        // constant signal is unchanged
        for v in y {
            assert!((v - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn preserves_length() {
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y = smooth(&x, 11, "hanning");
        assert_eq!(y.len(), x.len());
    }

    #[test]
    fn impulse_response_is_window() {
        // Smoothing a Kronecker delta gives the normalized window (the
        // moving-average kernel).  With window 5 the reflected extension
        // shifts the response: y[i] averages s[i+2..i+7], and the delta
        // lands at s[11], so y[5..=9] = 0.2.
        let mut x = vec![0.0; 15];
        x[7] = 1.0;
        let y = smooth(&x, 5, "flat");
        for i in 5..=9 {
            assert!((y[i] - 0.2).abs() < 1e-12, "y[{}] = {}", i, y[i]);
        }
        assert!(y[4].abs() < 1e-12);
        assert!(y[10].abs() < 1e-12);
    }
}
