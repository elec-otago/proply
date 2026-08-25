# Copyright (c) Tim Molteno tim@elec.ac.nz 2026
#!/usr/bin/env python3
"""Generate golden reference values for the proply-rs tests.

Each section replicates the corresponding Python proply / numpy / scipy
formula exactly and writes a JSON file under proply-rs/tests/golden/.

Run with the project venv:  build/venv/bin/python build/golden/gen_golden.py
"""
import json
import os
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

OUT = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)


def naca4_shape_points(m, p, t, chord, te, n):
    """foil.py NACA4.get_shape_points, verbatim."""
    n = n * 5
    beta = np.linspace(0, np.pi, n)
    x = (1.0 - np.cos(beta)) / 2
    y_offset = np.linspace(0, te / 2, n)
    yt = (5.0 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                     + 0.2843 * x ** 3 - 0.1036 * x ** 4) + y_offset)
    yc = (m / (p ** 2)) * (2.0 * p * x - x ** 2)
    yc2 = (m / ((1.0 - p) ** 2)) * (1.0 - 2.0 * p + 2 * p * x - x ** 2)
    yc[x > p] = yc2[x > p]
    dyc = m * (2.0 * p - 2 * x) / p ** 2
    dyc[x > p] = (2 * m * (p - x) / (p - 1.0) ** 2)[x > p]
    theta = np.arctan(dyc)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    c = chord
    return (xl[::5] * c, yl[::5] * c, xu[::5] * c, yu[::5] * c)


def foil_rotate(x, y, x0, y0, theta):
    x2 = (x - x0) * np.cos(theta) + (y - y0) * np.sin(theta)
    y2 = -(x - x0) * np.sin(theta) + (y - y0) * np.cos(theta)
    return x2, y2


def naca4_bounding_box(m, p, t, chord, te, theta):
    """foil.py Foil.get_bounding_box via get_points(50, theta)."""
    xl, yl, xu, yu = naca4_shape_points(m, p, t, chord, te, 50)
    x0 = 0.67 * (np.max(xu) - np.min(xu))
    y0 = 0.0
    xl, yl = foil_rotate(xl, yl, x0, y0, theta)
    xu, yu = foil_rotate(xu, yu, x0, y0, theta)
    x = np.concatenate([xl, xu])
    y = np.concatenate([yl, yu])
    return float(np.min(x)), float(np.max(x)), float(np.min(y)), float(np.max(y))


def smooth(x, window_len=11, window="hanning"):
    """smooth.py, verbatim."""
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == "flat":
        w = np.ones(window_len, "d")
    else:
        w = eval("np." + window + "(window_len)")
    y = np.convolve(w / w.sum(), s, mode="valid")
    return y[(window_len // 2):-(window_len // 2)]


def bem_precalc(foil_simulator, dv, a_prime, theta, omega, r, dr, u_0, B):
    """optimize.py precalc (with a flat-plate simulator)."""
    u = u_0 + dv
    v = omega * r * (1.0 - a_prime)
    phi = np.arctan(u / v)
    alpha = theta - phi
    v_rel = np.sqrt(u ** 2 + v ** 2)
    C_D = 1.28 * np.sin(alpha)  # PlateSimulatedFoil.get_cd
    C_L = 2.0 * np.pi * alpha   # PlateSimulatedFoil.get_cl
    return float(C_L), float(C_D), float(phi)


def bem_iterate(c, dv, a_prime, theta, omega, r, dr, u_0, B):
    C_L, C_D, phi = bem_precalc(None, dv, a_prime, theta, omega, r, dr, u_0, B)
    dv_new = (-B * c * (C_D * (dv + u_0) + C_L * omega * r * (a_prime - 1))
              * np.sqrt(omega ** 2 * r ** 2 * (a_prime - 1) ** 2 + (dv + u_0) ** 2)
              / (4 * np.pi * (dr + 2 * r) * (dv + u_0)))
    a_prime_new = (-B * c * np.sqrt(omega ** 2 * r ** 2 * (a_prime - 1) ** 2 + (dv + u_0) ** 2)
                   * (C_D * omega * r * (a_prime - 1) - C_L * (dv + u_0))
                   / (4 * np.pi * omega * r * (dr + 2 * r) * (dv + u_0)))
    return float(dv_new), float(a_prime_new)


def lsq(C_L, C_D, c, dv, a_prime, theta, omega, r, dr, u_0, B):
    """optimize.py lsq, verbatim."""
    minfun = (-B * c * np.sqrt(omega ** 2 * r ** 2 * (-a_prime + 1) ** 2 + (dv + u_0) ** 2)
              * (C_D * omega * r * (-a_prime + 1) + C_L * (dv + u_0))
              / (4 * np.pi * omega * r * (dr + 2 * r) * (dv + u_0)) + a_prime) ** 2 / (a_prime + 0.01) ** 2 + (
        B * c * (C_D * (dv + u_0) - C_L * omega * r * (-a_prime + 1))
        * np.sqrt(omega ** 2 * r ** 2 * (-a_prime + 1) ** 2 + (dv + u_0) ** 2)
        / (4 * np.pi * (dr + 2 * r) * (dv + u_0)) + dv) ** 2 / dv ** 2
    return float(minfun)


def min_func2(x, theta, omega, r, dr, u_0, B, chord):
    dv, a_prime = x
    C_L, C_D, phi = bem_precalc(None, dv, a_prime, theta, omega, r, dr, u_0, B)
    return lsq(C_L, C_D, chord, dv, a_prime, theta, omega, r, dr, u_0, B)


def min_all(x, goal, rpm, r, dr, u_0, B, chord0, maxchord):
    theta, dv, a_prime, chord = x
    omega = 2 * np.pi * rpm / 60
    dv2, a_prime2 = bem_iterate(chord, dv, a_prime, theta, omega, r, dr, u_0, B)
    err = abs((dv - dv2) / (dv + dv2)) + abs((a_prime - a_prime2) / (a_prime + a_prime2))
    err += 10 * ((dv2 - goal) / (dv2 + goal)) ** 2
    u = u_0 + dv
    torque = 2 * np.pi * a_prime * dr * omega * r ** 2 * 1.225 * u * (dr + 2 * r)
    thrust = 2 * np.pi * dr * dv * 1.225 * u * (dr + 2 * r)
    eff = abs(thrust / torque)
    err += 50.0 / eff
    return float(err)


# ---------------------------------------------------------------------------
# 1. NACA4 shape points
# ---------------------------------------------------------------------------
def gen_naca4():
    out = {}
    for name, (m, p, t, chord, te) in {
        "symmetric_12": (0.0, 0.4, 0.12, 1.0, 0.0),
        "cambered_15": (0.06, 0.4, 0.15, 0.1, 0.01 / 0.1),  # te set via set_trailing_edge
    }.items():
        xl, yl, xu, yu = naca4_shape_points(m, p, t, chord, te, 42)
        out[name] = {
            "chord": chord,
            "xl": xl.tolist(), "yl": yl.tolist(),
            "xu": xu.tolist(), "yu": yu.tolist(),
        }
    # bounding box at a twist angle
    x0, x1, y0, y1 = naca4_bounding_box(0.06, 0.4, 0.15, 0.02, 0.0, 0.3)
    out["bounding_box"] = {"x0": x0, "x1": x1, "y0": y0, "y1": y1}
    with open(os.path.join(OUT, "naca4.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 2. PCHIP
# ---------------------------------------------------------------------------
def gen_pchip():
    out = {}
    # scimitar-style dataset from prop.py get_scimitar_offset
    radius, hub_r, scimitar_percent = 0.0625, 0.005, -5.0
    x = np.array([0, hub_r, radius * 0.8, radius])
    y = np.array([0.0, 1.1 * 0.0, radius * scimitar_percent / 100.0, 0.0])
    p = PchipInterpolator(x, y)
    q = np.linspace(hub_r, radius, 23)
    out["scimitar"] = {"x": q.tolist(), "y": p(q).tolist()}
    # chord-smoothing style dataset (c_points / extra_chords from prop.py)
    hub_depth = 0.003
    chords = np.linspace(0.003, 0.0012, 40) + 0.0004 * np.sin(np.linspace(0, 8, 40))
    extra = np.concatenate((0.9 * np.array([hub_depth, hub_depth, hub_depth]), chords))
    sm = smooth(extra)
    c_pts = np.concatenate((np.array([0, hub_r / 2, 0.9 * hub_r]),
                            np.linspace(hub_r, radius, 40)))
    pc = PchipInterpolator(c_pts, sm)
    rq = np.linspace(hub_r, radius, 19)
    out["chord_smooth"] = {
        "smooth_input": extra.tolist(),
        "smooth_output": sm.tolist(),
        "r": rq.tolist(),
        "chord": pc(rq).tolist(),
    }
    with open(os.path.join(OUT, "pchip.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 3. polyfit (degree 9 on polar-like data, degree 4 on twist-like data)
# ---------------------------------------------------------------------------
def gen_polyfit():
    out = {}
    alpha = np.radians(np.arange(-20, 20.5, 0.5))
    # smooth, realistic cl/cd curves (cambered foil-ish)
    cl = 0.62 * np.sin(2.0 * alpha) + 0.05 * alpha + 0.4 * alpha ** 2 - 3.0 * alpha ** 3
    cd = 0.008 + 0.2 * alpha ** 2 + 1.5 * alpha ** 4
    cl9 = np.polyfit(alpha, cl, 9)
    cd9 = np.polyfit(alpha, cd, 9)
    out["polar_cl9"] = cl9.tolist()
    out["polar_cd9"] = cd9.tolist()
    # evaluation points (radians)
    ev = np.radians([-15.0, -5.0, 0.0, 5.0, 12.0])
    out["eval_rad"] = ev.tolist()
    out["eval_cl"] = np.poly1d(cl9)(ev).tolist()
    out["eval_cd"] = np.poly1d(cd9)(ev).tolist()
    # degree-4 twist fit (prop.py full_optimize)
    r = np.linspace(0.0625, 0.005, 40)
    twist = np.radians(25.0 * (r / 0.0625) ** 0.6 + 5.0 * np.sin(r * 40))
    c4 = np.polyfit(r[::-1], twist, 4)
    out["twist4"] = c4.tolist()
    out["twist_r"] = r[::-1].tolist()
    out["twist_eval"] = np.poly1d(c4)(r[::-1]).tolist()
    with open(os.path.join(OUT, "polyfit.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 4. motor model
# ---------------------------------------------------------------------------
def gen_motor():
    Kv, I0, Rm, V = 1900.0, 0.5, 0.405, 11.0
    Kq = 30.0 / (np.pi * Kv)
    Imax = np.sqrt(V * I0 / Rm)
    Qmax = Kq * (Imax - I0)
    RPMmax = Kv * (V - Imax * Rm)
    Pmax = 2.0 * np.pi * Qmax * (RPMmax / 60)
    out = {
        "Kq": float(Kq),
        "Imax": float(Imax),
        "Qmax": float(Qmax),
        "RPMmax": float(RPMmax),
        "Pmax": float(Pmax),
        "torque_at_3A": float(Kq * (3.0 - I0)),
        "rpm_at_0_01": float(np.pi * Kv ** 2 * 0.01 / 30),
    }
    with open(os.path.join(OUT, "motor.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 5. BEM equations + optimizer reference (flat-plate simulator)
# ---------------------------------------------------------------------------
def gen_bem():
    out = {}
    # a typical mid-blade station
    dv, a_prime, theta, rpm, r, dr, u_0, B, c = 5.0, 0.05, np.radians(28.0), 12000.0, 0.03, 0.002, 1.0, 3, 0.008
    omega = 2.0 * np.pi * rpm / 60.0
    C_L, C_D, phi = bem_precalc(None, dv, a_prime, theta, omega, r, dr, u_0, B)
    dv_new, a_prime_new = bem_iterate(c, dv, a_prime, theta, omega, r, dr, u_0, B)
    u = u_0 + dv
    dT = 2 * np.pi * dr * dv * 1.225 * u * (dr + 2 * r)
    dM = 2 * np.pi * a_prime * dr * omega * r ** 2 * 1.225 * u * (dr + 2 * r)
    dv_from_thrust = -u_0 / 2 + np.sqrt(np.pi * 0.05 ** 2 * 1.225 ** 2 * u_0 ** 2
                                        + 2 * 0.3 * 1.225) / (2 * np.sqrt(np.pi) * 0.05 * 1.225)
    out["precalc"] = {"CL": C_L, "CD": C_D, "phi": phi, "omega": omega}
    out["iterate"] = {"dv_new": dv_new, "a_prime_new": a_prime_new}
    out["forces"] = {"dT": float(dT), "dM": float(dM),
                     "dv_from_thrust": float(dv_from_thrust)}
    # SLSQP reference for bem_iterate (min_func2), the Python optimizer path
    x0 = np.array([dv, 0.01])
    cons = [{"type": "ineq", "fun": lambda x: x[0]},
            {"type": "ineq", "fun": lambda x: 3 * dv - x[0]},
            {"type": "ineq", "fun": lambda x: x[1]},
            {"type": "ineq", "fun": lambda x: 0.3 - x[1]}]
    res = minimize(min_func2, x0, args=(theta, omega, r, dr, u_0, B, c),
                   method="SLSQP", constraints=cons,
                   options={"disp": False, "maxiter": 1000})
    out["slsqp_bem"] = {"x": res.x.tolist(), "fun": float(res.fun)}
    # SLSQP reference for optimize_all (min_all)
    x0 = np.array([phi, dv, 0.002, c])
    cons = [{"type": "ineq", "fun": lambda x: x[0] - (phi - np.radians(8))},
            {"type": "ineq", "fun": lambda x: (phi + np.radians(15)) - x[0]},
            {"type": "ineq", "fun": lambda x: x[1] - dv / 2},
            {"type": "ineq", "fun": lambda x: 2 * dv - x[1]},
            {"type": "ineq", "fun": lambda x: x[2]},
            {"type": "ineq", "fun": lambda x: 0.2 - x[2]},
            {"type": "ineq", "fun": lambda x: x[3]},
            {"type": "ineq", "fun": lambda x: 0.012 - x[3]}]
    res = minimize(min_all, x0, args=(dv, rpm, r, dr, u_0, B, c, 0.012),
                   method="SLSQP", constraints=cons,
                   options={"disp": False, "maxiter": 1000})
    out["slsqp_all"] = {"x": res.x.tolist(), "fun": float(res.fun)}
    with open(os.path.join(OUT, "bem.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 6. Buhl (2005) turbulent-wake CT(a) relation
# ---------------------------------------------------------------------------
def gen_buhl():
    """Buhl (2005), NREL/TP-500-36834, Eqs. 1 + 18 (decelerating-disk
    convention, a in [0, 1], F the tip/hub loss factor):
        a <= 0.4:  CT = 4 F a (1 - a)
        a >  0.4:  CT = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2
    """
    out = {}
    a_grid = np.linspace(0.0, 0.999, 25)
    for name, F in {"F1": 1.0, "F08": 0.8}.items():
        ct = np.where(
            a_grid <= 0.4,
            4.0 * F * a_grid * (1.0 - a_grid),
            8.0 / 9.0
            + (4.0 * F - 40.0 / 9.0) * a_grid
            + (50.0 / 9.0 - 4.0 * F) * a_grid ** 2,
        )
        out[name] = {"a": a_grid.tolist(), "ct": ct.tolist()}
    # inverse: a(CT) for both loss factors
    ct_pts = np.array([0.2, 0.5, 0.9, 0.96, 1.0, 1.2, 1.5, 1.9])
    inv = {}
    for name, F in {"F1": 1.0, "F08": 0.8}.items():
        a_inv = []
        for ct in ct_pts:
            q = ct / F
            if q <= 0.96:
                a = 0.5 * (1.0 - np.sqrt(1.0 - q))
            else:
                c2 = 50.0 / 9.0 - 4.0 * F
                c1 = 4.0 * F - 40.0 / 9.0
                c0 = 8.0 / 9.0 - ct
                a = (-c1 + np.sqrt(c1 * c1 - 4.0 * c2 * c0)) / (2.0 * c2)
            a_inv.append(float(a))
        inv[name] = a_inv
    out["invert"] = {"ct": ct_pts.tolist(), "F1": inv["F1"], "F08": inv["F08"]}
    with open(os.path.join(OUT, "buhl.json"), "w") as f:
        json.dump(out, f, indent=1)


# ---------------------------------------------------------------------------
# 7. ARA-D family (src/arad.rs, ported from legacy foil_ARA.py ARADFoil)
# ---------------------------------------------------------------------------
def arad_load_selig(path):
    """The Selig parser of src/arad.rs (header line, upper TE->LE, lower
    LE->TE); returns (xl, yl, xu, yu), both LE->TE."""
    xs, ys = [], []
    with open(path) as f:
        for line in f.read().splitlines()[1:]:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                x, y = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            xs.append(x)
            ys.append(y)
    split = len(xs)
    for i in range(len(xs) - 1):
        if xs[i + 1] >= xs[i]:
            split = i + 1
            break
    return xs[split:], ys[split:], xs[:split][::-1], ys[:split][::-1]


def arad_section(t, n_stations=60):
    """The 60-station section at thickness t: degree-12 polyfit smoothing
    per surface, then the PCHIP thickness blend over the 20 nodes."""
    dat = os.path.join(OUT, "..", "..", "src", "arad")
    base_t = [0.06, 0.10, 0.13, 0.20]
    files = ["ara_d_6.dat", "ara_d_10.dat", "ara_d_13.dat", "ara_d_20.dat"]
    x = np.array([j / (n_stations - 1) for j in range(n_stations)])
    base_l, base_u = [], []
    for fname in files:
        xl, yl, xu, yu = arad_load_selig(os.path.join(dat, fname))
        base_l.append(np.poly1d(np.polyfit(xl, yl, 12))(x))
        base_u.append(np.poly1d(np.polyfit(xu, yu, 12))(x))
    low = np.linspace(0.0, 0.04, 7)
    high = np.linspace(0.25, 1.0, 9)
    t_nodes = np.concatenate([low, base_t, high])
    rows_l = [base_l[0] * (t / 0.06) for t in low] + list(base_l) \
        + [base_l[3] * (t / 0.20) for t in high]
    rows_u = [base_u[0] * (t / 0.06) for t in low] + list(base_u) \
        + [base_u[3] * (t / 0.20) for t in high]
    tt = float(np.clip(t, 0.0, 1.0))
    yl = np.array([float(PchipInterpolator(t_nodes, [r[j] for r in rows_l])(tt))
                   for j in range(n_stations)])
    yu = np.array([float(PchipInterpolator(t_nodes, [r[j] for r in rows_u])(tt))
                   for j in range(n_stations)])
    return x, yl, yu


def arad_shape_points(t, chord, te_m, n=42):
    """Arad::get_shape_points (the ported ARADFoil.get_shape_points)."""
    x, yl_st, yu_st = arad_section(t)
    init_te = yu_st[-1] - yl_st[-1]
    n5 = n * 5
    beta = np.linspace(0.0, np.pi, n5)
    xx = (1.0 - np.cos(beta)) / 2.0
    off = np.linspace(0.0, (te_m / chord - init_te) / 2.0, n5)
    lo = PchipInterpolator(x, yl_st)(xx) - off
    up = PchipInterpolator(x, yu_st)(xx) + off
    up[0] = lo[0]
    return xx[::5] * chord, lo[::5] * chord, xx[::5] * chord, up[::5] * chord


def gen_arad():
    out = {}
    for name, t in {"t06_node": 0.06, "t09_blend": 0.09, "t30_ramp": 0.30}.items():
        chord, te_m = 0.1, 0.001
        xl, yl, xu, yu = arad_shape_points(t, chord, te_m)
        out[name] = {
            "thickness": t, "chord": chord, "te": te_m,
            "xl": xl.tolist(), "yl": yl.tolist(),
            "xu": xu.tolist(), "yu": yu.tolist(),
        }
    with open(os.path.join(OUT, "arad.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    gen_naca4()
    gen_pchip()
    gen_polyfit()
    gen_motor()
    gen_bem()
    gen_buhl()
    gen_arad()
    print("golden files written to", os.path.abspath(OUT))
