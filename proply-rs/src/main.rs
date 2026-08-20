// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! proply-rs: propeller design in Rust.
//!
//! Port of the Python `proply` command: design a propeller blade with blade
//! element momentum theory (airfoil polars from rust-foil) and write the
//! propeller as a STEP (AP242) file.

use std::process::exit;
use std::sync::{Arc, Mutex};

use proply_rs::cache::PolarStore;
use proply_rs::design_parameters::DesignParameters;
use proply_rs::motor::Motor;
use proply_rs::optimize;
use proply_rs::prop::Prop;
use proply_rs::step_out;

/// CLI options.  Value/flag options are `Option` so an explicitly-passed flag
/// can be distinguished from "not given" (allowing the JSON parameter file to
/// carry the defaults, with the CLI overriding only what is passed).
struct Args {
    param: String,
    n: Option<usize>,
    bem: Option<bool>,
    lifting_line: Option<bool>,
    auto: Option<bool>,
    resolution: Option<usize>,
    dir: Option<String>,
    step_file: Option<String>,
    plate: Option<bool>,
    ar: Option<f64>,
    chord_spline_n: Option<usize>,
    camber: Option<f64>,
    cst: Option<bool>,
    help: bool,
}

fn print_help() {
    println!(
        "proply-rs: propeller design in Rust.

USAGE:
    proply-rs --bem [OPTIONS]

The BEM (or --lifting-line) design loop must be selected to produce a
propeller.  The propeller is written as a STEP (AP242) file.

DESIGN OPTIONS:
    --bem                  Use the blade-element momentum design loop.
    --lifting-line         Use the coupled vortex lifting-line design loop
                           (spanwise-induced losses; smooth chord).
    --param <FILE>         JSON propeller parameter file
                           (default: prop_design.json).
    --n <N>                STEP loft resolution (default: 40).
    --resolution <MM>      Radial (spanwise) resolution in millimetres;
                           2 x (radius - hub_radius) stations
                           (default: 40).
    --naca                 NACA airfoil family (default).
    --cst                  CST (Kulfan) airfoil family: every station uses
                           the default 18-parameter section, re-thicknessed
                           and cambered to the design's radial laws.
    --auto                 Re-run the design loop, reducing the thrust
                           target until the torque drops below ~1.5 x Qmax.
    --ar <N>               Minimum blade aspect ratio (R - hub)/<mean chord>;
                           caps the chord (thinner blade) in the lifting line.
    --chord-spline-n <N>   Number of control points for the smooth spline
                           chord (lifting line; default: 3).
    --camber <M>           Fixed foil camber (fraction of chord).  Without
                           it, the lifting-line design scans {{0, 0.02, 0.04}}
                           plus a per-station distribution composed from them,
                           keeping the best performer.
    --plate                Use analytic flat-plate polars (testing only).

OUTPUT OPTIONS:
    --dir <DIR>            Directory for the output STEP (created if needed;
                           default: .).
    --step-file <FILE>     Explicit output STEP path (overrides --dir +
                           <param name>.step).

ALL OPTIONS IN JSON:
    Every design/run option above can instead be set in the --param JSON
    file (keys: bem, lifting_line, auto, resolution, n, ar, plate, cst,
    dir, step_file, chord_spline_n, camber).  An explicit CLI flag
    overrides the JSON value, which overrides the built-in default.

OTHER:
    --help, -h             Print this help and exit.
"
    );
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        param: "prop_design.json".into(),
        n: None,
        bem: None,
        auto: None,
        resolution: None,
        dir: None,
        step_file: None,
        plate: None,
        lifting_line: None,
        ar: None,
        chord_spline_n: None,
        camber: None,
        cst: None,
        help: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        // Support both "--flag value" and "--flag=value".
        let (name, inline) = match arg.split_once('=') {
            Some((n, v)) => (n.to_string(), Some(v.to_string())),
            None => (arg.clone(), None),
        };
        let mut value = || -> Result<String, String> {
            if let Some(v) = &inline {
                Ok(v.clone())
            } else {
                it.next().ok_or_else(|| format!("{} needs a value", name))
            }
        };
        match name.as_str() {
            "--param" => a.param = value()?,
            "--n" => a.n = Some(value()?.parse().map_err(|_| "bad --n".to_string())?),
            "--bem" => a.bem = Some(true),
            "--auto" => a.auto = Some(true),
            "--naca" => {} // NACA 4-series family (the default)
            "--cst" => a.cst = Some(true), // CST (Kulfan) foil family
            "--resolution" => a.resolution = Some(value()?.parse().map_err(|_| "bad --resolution".to_string())?),
            "--dir" => a.dir = Some(value()?),
            "--step-file" => a.step_file = Some(value()?),
            "--plate" => a.plate = Some(true), // testing: analytic flat-plate polars
            "--lifting-line" => a.lifting_line = Some(true), // coupled vortex-lattice design
            "--ar" => a.ar = Some(value()?.parse().map_err(|_| "bad --ar".to_string())?),
            "--chord-spline-n" => {
                a.chord_spline_n = Some(value()?.parse().map_err(|_| "bad --chord-spline-n".to_string())?)
            }
            "--camber" => {
                a.camber = Some(value()?.parse().map_err(|_| "bad --camber".to_string())?)
            }
            "--help" | "-h" => a.help = true,
            "--mesh" => return Err("--mesh (GMSH) is not yet ported".into()),
            "--arad" => return Err("--arad (ARA-D foils) is not yet ported".into()),
            other => return Err(format!("unknown argument: {}", other)),
        }
    }
    Ok(a)
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("proply-rs: {}", e);
            exit(1);
        }
    };

    if args.help {
        print_help();
        exit(0);
    }

    let mut param = match DesignParameters::from_file(&args.param) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("proply-rs: {}", e);
            exit(1);
        }
    };
    // Every CLI option can equally be set in the JSON design file; an
    // explicitly-passed CLI flag overrides the JSON value (which overrides the
    // default).
    if let Some(v) = args.bem {
        param.bem = v;
    }
    if let Some(v) = args.lifting_line {
        param.lifting_line = v;
    }
    if let Some(v) = args.auto {
        param.auto = v;
    }
    if let Some(v) = args.resolution {
        param.resolution = v;
    }
    if let Some(v) = args.n {
        param.n = v;
    }
    if let Some(v) = args.ar {
        param.ar = Some(v);
    }
    if let Some(v) = args.chord_spline_n {
        param.chord_spline_n = v;
    }
    if let Some(v) = args.camber {
        param.camber = Some(v);
    }
    if let Some(v) = args.cst {
        param.cst = v;
    }
    if let Some(v) = args.plate {
        param.plate = v;
    }
    if let Some(v) = args.dir {
        param.dir = v;
    }
    if let Some(v) = args.step_file {
        param.step_file = v;
    }

    if !param.bem && !param.lifting_line {
        eprintln!("proply-rs: select a design loop (set `bem` or `lifting_line`, or pass --bem / --lifting-line)");
        exit(1);
    }

    let resolution_m = (param.radius - param.hub_radius) / param.resolution as f64;

    let store: Arc<Mutex<PolarStore>> = Arc::new(Mutex::new(PolarStore::load(
        proply_rs::cache::default_cache_path().as_str(),
    )));

    let mut p = Prop::new(param.clone(), resolution_m, store.clone());
    p.n_blades = param.blades;
    p.set_plate_mode(param.plate);

    let m = Motor::new(param.motor_Kv, param.motor_no_load_current, param.motor_winding_resistance);
    let (optimum_torque, optimum_rpm) = m.get_qmax(param.motor_volts);
    let power = m.get_pmax(param.motor_volts);

    println!("\nPROPLY: Automatic propeller Design\n\n");
    println!(
        "Optimum Motor Torque {:5.3} Nm at {:5.1} RPM, power={:5.1} Watts",
        optimum_torque, optimum_rpm, power
    );
    println!("Spanwise resolution (mm) {:4.2}", resolution_m * 1000.0);
    println!("{}", param);
    let dv = optimize::dv_from_thrust(param.thrust, param.radius, param.forward_airspeed);
    println!(
        "Airspeed at propellers (hovering): {:4.2} m/s",
        param.forward_airspeed + dv
    );
    println!("\n\n");

    let mut thrust = param.thrust;
    let goal_torque = optimum_torque * 1.5;
    let (mut q, mut t) = if param.lifting_line {
        p.lift_line_design(optimum_rpm, thrust, param.ar)
    } else {
        p.full_optimize(optimum_rpm, thrust)
    };
    println!("Total Thrust: {:5.2}, Torque: {:5.3}", t, q);
    if param.auto {
        while q > goal_torque {
            thrust *= 0.95 * goal_torque / q;
            // Re-run with the same design loop that produced the initial
            // design (lifting-line or BEM).
            (q, t) = if param.lifting_line {
                p.lift_line_design(optimum_rpm, thrust, param.ar)
            } else {
                p.full_optimize(optimum_rpm, thrust)
            };
            println!("Total Thrust: {:5.2} (N), Torque: {:5.2} (Nm)", t, q);
        }
    }

    let step_filename = if param.step_file.is_empty() {
        format!("{}/{}.step", param.dir, param.name)
    } else {
        param.step_file
    };
    match step_out::write_prop(&mut p, param.n) {
        Ok(text) => {
            step_out::write_step_file(&step_filename, &text).unwrap_or_else(|e| {
                eprintln!("cannot write {}: {}", step_filename, e);
                exit(1);
            });
            println!("Wrote {}", step_filename);
        }
        Err(e) => {
            eprintln!("proply-rs: {}", e);
            exit(1);
        }
    }

    // Persist any newly simulated polars.
    store.lock().unwrap().save();
}
