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

struct Args {
    param: String,
    n: usize,
    bem: bool,
    auto: bool,
    resolution: usize,
    dir: String,
    step_file: String,
    plate: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        param: "prop_design.json".into(),
        n: 40,
        bem: false,
        auto: false,
        resolution: 40,
        dir: ".".into(),
        step_file: String::new(),
        plate: false,
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
            "--n" => a.n = value()?.parse().map_err(|_| "bad --n".to_string())?,
            "--bem" => a.bem = true,
            "--auto" => a.auto = true,
            "--naca" => {} // the only foil family supported by this port
            "--resolution" => a.resolution = value()?.parse().map_err(|_| "bad --resolution".to_string())?,
            "--dir" => a.dir = value()?,
            "--step-file" => a.step_file = value()?,
            "--plate" => a.plate = true, // testing: analytic flat-plate polars
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

    if !args.bem {
        eprintln!("proply-rs: the --bem design loop is required to produce a propeller");
        exit(1);
    }

    let param = match DesignParameters::from_file(&args.param) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("proply-rs: {}", e);
            exit(1);
        }
    };

    let resolution_m = (param.radius - param.hub_radius) / args.resolution as f64;

    let store: Arc<Mutex<PolarStore>> = Arc::new(Mutex::new(PolarStore::load(
        proply_rs::cache::default_cache_path().as_str(),
    )));

    let mut p = Prop::new(param.clone(), resolution_m, store.clone());
    p.n_blades = param.blades;
    p.set_plate_mode(args.plate);

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
    let (mut q, mut t) = p.full_optimize(optimum_rpm, thrust);
    println!("Total Thrust: {:5.2}, Torque: {:5.3}", t, q);
    if args.auto {
        while q > goal_torque {
            thrust *= 0.95 * goal_torque / q;
            (q, t) = p.full_optimize(optimum_rpm, thrust);
            println!("Total Thrust: {:5.2} (N), Torque: {:5.2} (Nm)", t, q);
        }
    }

    let step_filename = if args.step_file.is_empty() {
        format!("{}/{}.step", args.dir, param.name)
    } else {
        args.step_file
    };
    match step_out::write_prop(&mut p, args.n) {
        Ok(text) => {
            std::fs::write(&step_filename, text)
                .unwrap_or_else(|e| {
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
