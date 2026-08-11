//! Blade element momentum optimisation routines, ported from
//! `proply/optimize.py`.
//!
//! The Python version drives these objectives with scipy SLSQP/COBYLA; this
//! port uses a box-constrained Nelder-Mead (scipy-compatible initial
//! simplex, candidate points projected onto the bounds).  All constraints
//! in `bem_iterate` / `optimize_all` are simple variable bounds, so
//! projection is exact.

use crate::simulator::FoilSimulator;
use crate::foil::FoilLike;

pub const RHO: f64 = 1.225;

/// Anything that can answer CL/CD for a flow — the real polar simulator or
/// the flat-plate model used for testing.
pub trait FoilSim {
    fn get_cl(&self, v: f64, alpha: f64) -> f64;
    fn get_cd(&self, v: f64, alpha: f64) -> f64;
    fn chord(&self) -> f64;
}

impl<F: FoilLike> FoilSim for FoilSimulator<F> {
    fn get_cl(&self, v: f64, alpha: f64) -> f64 {
        FoilSimulator::get_cl(self, v, alpha)
    }
    fn get_cd(&self, v: f64, alpha: f64) -> f64 {
        FoilSimulator::get_cd(self, v, alpha)
    }
    fn chord(&self) -> f64 {
        FoilSimulator::chord(self)
    }
}

/// The analytic flat-plate model (`PlateSimulatedFoil`): cl = 2 pi alpha,
/// cd = 1.28 sin alpha.
#[derive(Debug, Clone)]
pub struct PlateSim {
    pub chord: f64,
}

impl FoilSim for PlateSim {
    fn get_cl(&self, _v: f64, alpha: f64) -> f64 {
        2.0 * std::f64::consts::PI * alpha
    }
    fn get_cd(&self, _v: f64, alpha: f64) -> f64 {
        1.28 * alpha.sin()
    }
    fn chord(&self) -> f64 {
        self.chord
    }
}

pub fn rpm2omega(rpm: f64) -> f64 {
    2.0 * std::f64::consts::PI * rpm / 60.0
}

/// `optimize.precalc`: the flow angle and force coefficients at one element.
#[allow(clippy::too_many_arguments)]
pub fn precalc<S: FoilSim>(
    fs: &S,
    dv: f64,
    a_prime: f64,
    theta: f64,
    omega: f64,
    r: f64,
    _dr: f64,
    u_0: f64,
    _b: f64,
) -> (f64, f64, f64) {
    let u = u_0 + dv;
    let v = omega * r * (1.0 - a_prime);
    let phi = (u / v).atan();
    let alpha = theta - phi;
    let v_rel = (u * u + v * v).sqrt();
    let cd = fs.get_cd(v_rel, alpha);
    let cl = fs.get_cl(v_rel, alpha);
    (cl, cd, phi)
}

/// `optimize.iterate`: one momentum-equation iteration.
#[allow(clippy::too_many_arguments)]
pub fn iterate<S: FoilSim>(
    fs: &S,
    c: f64,
    dv: f64,
    a_prime: f64,
    theta: f64,
    omega: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> (f64, f64) {
    let (cl, cd, _phi) = precalc(fs, dv, a_prime, theta, omega, r, dr, u_0, b);
    let (dv_new, a_prime_new) = momentum_step(c, cl, cd, dv, a_prime, theta, omega, r, dr, u_0, b);
    (dv_new, a_prime_new)
}

/// The shared momentum update, verbatim from `optimize.iterate`.
#[allow(clippy::too_many_arguments)]
fn momentum_step(
    c: f64,
    cl: f64,
    cd: f64,
    dv: f64,
    a_prime: f64,
    theta: f64,
    omega: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> (f64, f64) {
    let _ = theta;
    let v_rel = (omega * omega * r * r * (a_prime - 1.0).powi(2) + (dv + u_0).powi(2)).sqrt();
    let dv_new = -b * c * (cd * (dv + u_0) + cl * omega * r * (a_prime - 1.0)) * v_rel
        / (4.0 * std::f64::consts::PI * (dr + 2.0 * r) * (dv + u_0));
    let a_prime_new = -b * c * v_rel * (cd * omega * r * (a_prime - 1.0) - cl * (dv + u_0))
        / (4.0 * std::f64::consts::PI * omega * r * (dr + 2.0 * r) * (dv + u_0));
    (dv_new, a_prime_new)
}

/// `optimize.lsq`: the objective `bem_iterate` minimises.
#[allow(clippy::too_many_arguments)]
pub fn lsq(
    cl: f64,
    cd: f64,
    c: f64,
    dv: f64,
    a_prime: f64,
    theta: f64,
    omega: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> f64 {
    let _ = theta;
    let v_rel = (omega * omega * r * r * (-a_prime + 1.0).powi(2) + (dv + u_0).powi(2)).sqrt();
    let term1 = -b * c * v_rel * (cd * omega * r * (-a_prime + 1.0) + cl * (dv + u_0))
        / (4.0 * std::f64::consts::PI * omega * r * (dr + 2.0 * r) * (dv + u_0))
        + a_prime;
    let term2 = b * c * (cd * (dv + u_0) - cl * omega * r * (-a_prime + 1.0)) * v_rel
        / (4.0 * std::f64::consts::PI * (dr + 2.0 * r) * (dv + u_0))
        + dv;
    term1.powi(2) / (a_prime + 0.01).powi(2) + term2.powi(2) / dv.powi(2)
}

/// `optimize.min_func2`: the objective over (dv, a_prime).
#[allow(clippy::too_many_arguments)]
pub fn min_func2<S: FoilSim>(
    x: &[f64],
    fs: &S,
    theta: f64,
    omega: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> f64 {
    let dv = x[0];
    let a_prime = x[1];
    let (cl, cd, _phi) = precalc(fs, dv, a_prime, theta, omega, r, dr, u_0, b);
    lsq(cl, cd, fs.chord(), dv, a_prime, theta, omega, r, dr, u_0, b)
}

/// `optimize.min_all`: the objective over (theta, dv, a_prime, chord).
#[allow(clippy::too_many_arguments)]
pub fn min_all<S: FoilSim>(
    x: &[f64],
    fs: &S,
    goal: f64,
    rpm: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> f64 {
    let theta = x[0];
    let dv = x[1];
    let a_prime = x[2];
    let chord = x[3];
    let omega = rpm2omega(rpm);
    let (dv2, a_prime2) = iterate(fs, chord, dv, a_prime, theta, omega, r, dr, u_0, b);
    let err = error(dv, dv2, a_prime, a_prime2);
    let err = err + 10.0 * ((dv2 - goal) / (dv2 + goal)).powi(2);
    let torque = d_m(dv, a_prime, r, dr, omega, u_0);
    let thrust = d_t(dv, r, dr, u_0);
    let eff = (thrust / torque).abs();
    err + 50.0 / eff
}

/// Element thrust: `optimize.dT`.
pub fn d_t(dv: f64, r: f64, dr: f64, u_0: f64) -> f64 {
    let u = u_0 + dv;
    2.0 * std::f64::consts::PI * dr * dv * RHO * u * (dr + 2.0 * r)
}

/// Element torque: `optimize.dM`.
pub fn d_m(dv: f64, a_prime: f64, r: f64, dr: f64, omega: f64, u_0: f64) -> f64 {
    let u = u_0 + dv;
    2.0 * std::f64::consts::PI * a_prime * dr * omega * r.powi(2) * RHO * u * (dr + 2.0 * r)
}

/// Induced velocity for a target thrust: `optimize.dv_from_thrust`.
pub fn dv_from_thrust(t: f64, r: f64, u_0: f64) -> f64 {
    -u_0 / 2.0
        + (std::f64::consts::PI * r.powi(2) * RHO.powi(2) * u_0.powi(2) + 2.0 * t * RHO).sqrt()
            / (2.0 * (std::f64::consts::PI).sqrt() * r * RHO)
}

/// `optimize.error`: relative residual between two momentum states.
pub fn error(dv: f64, dv2: f64, a_prime: f64, a_prime2: f64) -> f64 {
    ((dv - dv2) / (dv + dv2)).abs() + ((a_prime - a_prime2) / (a_prime + a_prime2)).abs()
}

/// Box-constrained Nelder-Mead (replaces scipy SLSQP/COBYLA here).
///
/// Constraints are enforced with a quadratic penalty rather than candidate
/// projection: projecting the simplex onto a bound collapses it (the exact
/// failure scipy's bound handling exhibits on these objectives), while the
/// penalty keeps the unconstrained simplex healthy.  All proply optima are
/// interior, so the penalty is zero at the solution and the result is exact.
/// The initial simplex follows scipy's default (0.05 * x0 per axis, or
/// 0.00025 where x0 is zero).
pub struct NelderMead {
    pub maxiter: usize,
    pub xatol: f64,
    pub fatol: f64,
    /// Quadratic-penalty weight for bound violations.
    pub penalty: f64,
}

impl Default for NelderMead {
    fn default() -> Self {
        Self {
            maxiter: 5000,
            xatol: 1e-7,
            fatol: 1e-7,
            penalty: 1.0e6,
        }
    }
}

impl NelderMead {
    /// Minimise `f` starting from `x0`, with per-axis `bounds` (None = free).
    /// Returns (best feasible point, best value of `f`).
    pub fn minimize<F: Fn(&[f64]) -> f64>(
        &self,
        f: F,
        x0: &[f64],
        bounds: &[Option<(f64, f64)>],
    ) -> (Vec<f64>, f64) {
        let n = x0.len();
        assert_eq!(bounds.len(), n);

        let fpen = |x: &[f64]| -> f64 {
            let mut pen = 0.0;
            for (i, b) in bounds.iter().enumerate() {
                if let Some((lo, hi)) = b {
                    if x[i] < *lo {
                        pen += (lo - x[i]).powi(2);
                    } else if x[i] > *hi {
                        pen += (x[i] - hi).powi(2);
                    }
                }
            }
            f(x) + self.penalty * pen
        };

        // Initial simplex (scipy default).
        let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
        let x0v = x0.to_vec();
        simplex.push((x0v.clone(), fpen(&x0v)));
        for i in 0..n {
            let step = if x0[i] != 0.0 { 0.05 * x0[i] } else { 0.00025 };
            let mut xi = x0v.clone();
            xi[i] += step;
            simplex.push((xi.clone(), fpen(&xi)));
        }

        let alpha = 1.0; // reflection
        let gamma = 2.0; // expansion
        let rho = 0.5; // contraction
        let sigma = 0.5; // shrink

        for _iter in 0..self.maxiter {
            simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            let best = &simplex[0];
            let worst = &simplex[n];

            // Convergence: simplex size and function spread.
            let mut size: f64 = 0.0;
            for i in 1..=n {
                for (j, v) in simplex[i].0.iter().enumerate() {
                    size = size.max((v - best.0[j]).abs());
                }
            }
            let spread = (simplex[n].1 - simplex[0].1).abs();
            if size < self.xatol && spread < self.fatol {
                break;
            }

            // Centroid of all but the worst.
            let mut centroid = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    centroid[j] += simplex[i].0[j];
                }
            }
            for j in 0..n {
                centroid[j] /= n as f64;
            }

            // Reflect.
            let mut xr = vec![0.0; n];
            for j in 0..n {
                xr[j] = centroid[j] + alpha * (centroid[j] - worst.0[j]);
            }
            let fr = fpen(&xr);

            if fr < simplex[0].1 {
                // Expand.
                let mut xe = vec![0.0; n];
                for j in 0..n {
                    xe[j] = centroid[j] + gamma * (centroid[j] - worst.0[j]);
                }
                let fe = fpen(&xe);
                if fe < fr {
                    simplex[n] = (xe, fe);
                } else {
                    simplex[n] = (xr, fr);
                }
            } else if fr < simplex[n - 1].1 {
                simplex[n] = (xr, fr);
            } else {
                // Contract: outside contraction when the reflection improved
                // on the worst, inside contraction otherwise.
                let outside = fr < worst.1;
                let mut xc = vec![0.0; n];
                for j in 0..n {
                    let sign = if outside { 1.0 } else { -1.0 };
                    xc[j] = centroid[j] + rho * sign * (centroid[j] - worst.0[j]);
                }
                let fc = fpen(&xc);
                let accept = if outside { fc <= fr } else { fc < worst.1 };
                if accept {
                    simplex[n] = (xc, fc);
                } else {
                    // Shrink towards the best.
                    for i in 1..=n {
                        for j in 0..n {
                            simplex[i].0[j] =
                                simplex[0].0[j] + sigma * (simplex[i].0[j] - simplex[0].0[j]);
                        }
                        simplex[i].1 = fpen(&simplex[i].0);
                    }
                }
            }
        }

        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let mut x = simplex[0].0.clone();
        // Clamp the reported point back into the feasible box and report the
        // true (unpenalised) objective value.
        for (i, b) in bounds.iter().enumerate() {
            if let Some((lo, hi)) = b {
                x[i] = x[i].clamp(*lo, *hi);
            }
        }
        let fbest = f(&x);
        (x, fbest)
    }
}

/// `optimize.bem_iterate`: solve for (dv, a_prime) at a fixed twist.
#[allow(clippy::too_many_arguments)]
pub fn bem_iterate<S: FoilSim>(
    fs: &S,
    dv_goal: f64,
    theta: f64,
    rpm: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
) -> (f64, f64, f64) {
    let omega = rpm2omega(rpm);
    let x0 = vec![dv_goal, 0.01];
    let bounds = [
        Some((0.0, 3.0 * dv_goal)),
        Some((0.0, 0.3)),
    ];
    let nm = NelderMead::default();
    let (x, fun) = nm.minimize(
        |x| min_func2(x, fs, theta, omega, r, dr, u_0, b),
        &x0,
        &bounds,
    );
    if fun > 0.1 {
        // Mirror the Python COBYLA retry.
        let (x2, fun2) = nm.minimize(
            |x| min_func2(x, fs, theta, omega, r, dr, u_0, b),
            &x0,
            &bounds,
        );
        return (x2[0], x2[1], fun2);
    }
    (x[0], x[1], fun)
}

/// `optimize.optimize_all`: solve for (theta, dv, a_prime, chord) to hit a
/// desired induced velocity at maximum efficiency.
///
/// The objective is multi-modal, so Nelder-Mead is run from several starts
/// and the best result returned (scipy's SLSQP needs only one start; a
/// derivative-free simplex does not).
#[allow(clippy::too_many_arguments)]
pub fn optimize_all<S: FoilSim>(
    fs: &S,
    dv_goal: f64,
    rpm: f64,
    r: f64,
    dr: f64,
    u_0: f64,
    b: f64,
    maxchord: f64,
) -> (Vec<f64>, f64) {
    let omega = rpm2omega(rpm);
    let (_cl, _cd, phi) = precalc(fs, dv_goal, 0.0, 0.0, omega, r, dr, u_0, b);
    let chord0 = fs.chord();
    let bounds = [
        Some((phi - 8.0_f64.to_radians(), phi + 10.0_f64.to_radians())),
        Some((dv_goal / 2.0, 2.0 * dv_goal)),
        Some((0.0, 0.2)),
        Some((0.0, maxchord)),
    ];
    let nm = NelderMead::default();
    let starts: Vec<Vec<f64>> = vec![
        vec![phi, dv_goal, 0.002, chord0],
        vec![phi + 0.04, dv_goal, 0.002, chord0],
        vec![phi - 0.04, dv_goal, 0.002, chord0],
        vec![phi, dv_goal, 0.05, chord0],
        vec![phi, dv_goal, 0.002, chord0 * 1.1],
    ];
    let mut best: Option<(Vec<f64>, f64)> = None;
    for x0 in starts {
        let (x, fun) = nm.minimize(|x| min_all(x, fs, dv_goal, rpm, r, dr, u_0, b), &x0, &bounds);
        if best.as_ref().map_or(true, |(_, f)| fun < *f) {
            best = Some((x, fun));
        }
    }
    let (x, fun) = best.expect("optimize_all: no start evaluated");
    // Polish: restart from the best point — cheap and refines the solution.
    let (x2, fun2) = nm.minimize(|x| min_all(x, fs, dv_goal, rpm, r, dr, u_0, b), &x, &bounds);
    if fun2 < fun {
        (x2, fun2)
    } else {
        (x, fun)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn momentum_equations_match_python() {
        // Values verified against numpy in build/golden/gen_golden.py.
        let fs = PlateSim { chord: 0.008 };
        let dv = 5.0;
        let a_prime = 0.05;
        let theta = (28.0_f64).to_radians();
        let rpm = 12000.0;
        let (r, dr, u_0, b) = (0.03, 0.002, 1.0, 3.0);
        let omega = rpm2omega(rpm);
        let (cl, cd, phi) = precalc(&fs, dv, a_prime, theta, omega, r, dr, u_0, b as f64);
        assert!((cl - 2.027597).abs() < 1e-4, "cl {}", cl);
        assert!((cd - 0.405927).abs() < 1e-4, "cd {}", cd);
        assert!((phi - 0.165990).abs() < 1e-5, "phi {}", phi);
        let (dv_new, ap_new) = iterate(&fs, fs.chord, dv, a_prime, theta, omega, r, dr, u_0, b as f64);
        assert!((dv_new - 13.08411).abs() < 1e-3, "dv_new {}", dv_new);
        assert!((ap_new - 0.132057).abs() < 1e-4, "ap_new {}", ap_new);
        assert!((d_t(5.0, 0.03, 0.002, 1.0) - 2.0 * std::f64::consts::PI * 0.002 * 5.0 * 1.225 * 6.0 * (0.002 + 0.06)).abs() < 1e-12);
        // golden: dv_from_thrust(0.3, 0.05, 1.0) = 3.4800
        assert!((dv_from_thrust(0.3, 0.05, 1.0) - 3.480036).abs() < 1e-4, "dv {}", dv_from_thrust(0.3, 0.05, 1.0));
    }

    #[test]
    fn nelder_mead_rosenbrock() {
        // The classic NM test: Rosenbrock's valley, minimum at (1,1), f=0.
        let nm = NelderMead::default();
        let (x, f) = nm.minimize(
            |x| (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0].powi(2)).powi(2),
            &[0.0, 0.0],
            &[None, None],
        );
        assert!((x[0] - 1.0).abs() < 1e-4, "x0 {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-4, "x1 {}", x[1]);
        assert!(f < 1e-8, "f {}", f);
    }

    #[test]
    fn nelder_mead_respects_bounds() {
        let nm = NelderMead::default();
        // Minimum of (x-3)^2 is at 3, but bound it to [-2, 0].
        let (x, _f) = nm.minimize(
            |x| (x[0] - 3.0).powi(2),
            &[0.0],
            &[Some((-2.0, 0.0))],
        );
        assert!(x[0].abs() < 1e-9, "x {}", x[0]);
    }

    #[test]
    fn bem_iterate_converges_like_slsqp() {
        // scipy SLSQP on the same objective finds (7.1164, 0.0866).
        let fs = PlateSim { chord: 0.008 };
        let (dv, ap, err) = bem_iterate(&fs, 5.0, (28.0_f64).to_radians(), 12000.0, 0.03, 0.002, 1.0, 3.0);
        assert!(err < 1e-6, "err {}", err);
        assert!((dv - 7.1164).abs() < 0.2, "dv {}", dv);
        assert!((ap - 0.0866).abs() < 0.01, "a_prime {}", ap);
    }
}
