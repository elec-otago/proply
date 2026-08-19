// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Benchmark the serial vs parallel alpha sweep on the canonical workload
//! (NACA 0012 at Re = 1e6, 80-point sweep from -20° to 20°), and report
//! the per-point parity between the two.
//!
//! Usage: cargo run --release -p rust-foil --example sweep_bench

use std::time::Instant;

use rust_foil::XFoil;

fn main() {
    let mut engine = XFoil::new();
    engine.set_show_output(false);
    engine.naca(12);
    engine.set_reynolds(1.0e6);
    engine.set_max_iter(150);

    // Serial (on a fresh engine clone, as the design loop would start).
    let mut ser = engine.clone();
    let t0 = Instant::now();
    let serial = ser.aseq(-20.0, 20.0, 80);
    let t_ser = t0.elapsed();

    // Parallel.
    let t0 = Instant::now();
    let par = engine.aseq_par(-20.0, 20.0, 80);
    let t_par = t0.elapsed();

    println!("threads:  {}", rayon::current_num_threads());
    println!("serial:   {:8.3?}", t_ser);
    println!("parallel: {:8.3?}", t_par);
    println!("speedup:  {:.2}x", t_ser.as_secs_f64() / t_par.as_secs_f64());

    // Parity summary between the two sweeps.
    assert_eq!(serial.len(), par.len());
    let mut max_cl = 0.0f64;
    let mut max_cd = 0.0f64;
    let mut conv_diff = 0usize;
    for (s, p) in serial.iter().zip(par.iter()) {
        max_cl = max_cl.max((s.1 - p.1).abs());
        max_cd = max_cd.max((s.2 - p.2).abs());
        if s.5 != p.5 {
            conv_diff += 1;
        }
    }
    println!("max |dCL| = {:.5}, max |dCD| = {:.6}, conv-flag diffs = {}", max_cl, max_cd, conv_diff);
}
