//! Write just the hub solid as a STEP file, for debugging FreeCAD imports
//! without running the full design loop.
//!
//! Usage: cargo run --release -p proply-rs --example hub_step -- <param.json> [out.step]

use std::process::exit;

use proply_rs::design_parameters::DesignParameters;
use proply_rs::step_out;

fn main() {
    let mut args = std::env::args().skip(1);
    let param_file = args.next().unwrap_or_else(|| {
        eprintln!("usage: hub_step <param.json> [out.step]");
        exit(2);
    });
    let out = args.next().unwrap_or_else(|| "hub.step".into());
    let param = match DesignParameters::from_file(&param_file) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("hub_step: {}", e);
            exit(1);
        }
    };
    match step_out::hub_only_step(&param) {
        Ok(text) => {
            step_out::write_step_file(&out, &text).unwrap_or_else(|e| {
                eprintln!("cannot write {}: {}", out, e);
                exit(1);
            });
            println!("Wrote {}", out);
        }
        Err(e) => {
            eprintln!("hub_step: {}", e);
            exit(1);
        }
    }
}
