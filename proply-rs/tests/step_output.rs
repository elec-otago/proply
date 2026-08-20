// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Round-trip test for the STEP writer: build a small synthetic propeller,
//! write the STEP text, and re-read it with step-io's own parser.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use proply_rs::blade_element::BladeElement;
use proply_rs::cache::PolarStore;
use proply_rs::design_parameters::DesignParameters;
use proply_rs::foil::Naca4;
use proply_rs::prop::Prop;
use proply_rs::step_out;

fn synthetic_prop() -> Prop {
    let param = DesignParameters::default(); // radius 0.0625, hub 0.005, 2 blades
    let store: Arc<Mutex<PolarStore>> =
        Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json")));
    let mut p = Prop::new(param, 0.002, store.clone());
    p.n_blades = 2;
    // Four stations hub -> tip with simple foils (no design loop needed for
    // the geometry).
    for r in [0.006, 0.02, 0.04, 0.06] {
        let foil = Rc::new(RefCell::new(Naca4::new(0.012, 0.12, 0.06, 0.4)));
        let be = BladeElement::new(r, 0.002, foil, 0.4, 10000.0, 1.0, store.clone());
        p.blade_elements.push(be);
    }
    p
}

#[test]
fn step_round_trip() {
    let mut p = synthetic_prop();
    let text = step_out::write_prop(&mut p, 12).expect("write_prop");
    // The file is a valid Part 21 document.
    assert!(text.starts_with("ISO-10303-21;"), "wrong header");
    assert!(text.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"), "no B-rep");

    // Round-trip through step-io's own reader.
    let (model, report) = step_io::read(text.as_bytes()).expect("step_io read");
    let _ = model;
    assert!(
        report.dropped.is_empty(),
        "round-trip dropped {} entities",
        report.dropped.len()
    );

    // Sanity: the writer emits both a blade and a hub solid plus the
    // assembly placements (blade + hub instances).
    assert!(text.contains("'blade'"), "blade part missing");
    assert!(text.contains("'hub'"), "hub part missing");
    assert!(text.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), "no assembly");
    assert!(text.contains("MANIFOLD_SOLID_BREP"), "no manifold solid");
}

#[test]
fn step_writer_errors_with_single_station() {
    let mut p = synthetic_prop();
    p.blade_elements.truncate(1);
    let res = step_out::write_prop(&mut p, 12);
    assert!(res.is_err(), "single station should be rejected");
}

#[test]
fn stations_are_vertically_centred_on_the_hub() {
    // A twisted station must be centred on z = 0 (the hub's mid-plane).
    // The twist rotation pivots on the chord line, which leaves a twisted
    // section riding to one side of z = 0 — the blade root then protrudes
    // above the hub instead of being centred on it.
    let store: Arc<Mutex<PolarStore>> =
        Arc::new(Mutex::new(PolarStore::load("/nonexistent/cache.json")));
    for twist in [0.0, 0.4, -0.4] {
        let foil = Rc::new(RefCell::new(Naca4::new(0.012, 0.12, 0.06, 0.4)));
        let be = BladeElement::new(0.006, 0.002, foil, twist, 10000.0, 1.0, store.clone());
        let (lower, upper) = be.get_foil_points(24, 0.0);
        let mut zmin = f64::INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        for p in lower.iter().chain(upper.iter()) {
            zmin = zmin.min(p[2]);
            zmax = zmax.max(p[2]);
        }
        assert!(
            (0.5 * (zmin + zmax)).abs() < 1.0e-12,
            "twist {}: station centre at z = {}",
            twist,
            0.5 * (zmin + zmax)
        );
        assert!(zmax - zmin > 1.0e-3, "twist {}: degenerate station", twist);
    }
}

#[test]
fn hub_step_round_trip() {
    let param = DesignParameters::default();
    let text = step_out::hub_only_step(&param).expect("hub_only_step");

    // The file is a valid Part 21 document.
    assert!(text.starts_with("ISO-10303-21;"), "wrong header");
    assert!(text.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"), "no B-rep");

    // Round-trip through step-io's own reader.
    let (model, report) = step_io::read(text.as_bytes()).expect("step_io read");
    let _ = model;
    assert!(
        report.dropped.is_empty(),
        "round-trip dropped {} entities",
        report.dropped.len()
    );

    // The hub is a single closed shell of four faces: outer cylinder side,
    // bore side, and the two annular caps.
    assert!(text.contains("'hub'"), "hub part missing");
    assert!(text.contains("MANIFOLD_SOLID_BREP"), "no manifold solid");
    assert_eq!(text.matches("ADVANCED_FACE").count(), 4, "expected 4 faces");
    assert_eq!(text.matches("CIRCLE").count(), 8, "expected 8 arc edges");

    // The annular caps are ring faces: FACE_BOUND (inner) entities in
    // addition to FACE_OUTER_BOUND.
    assert!(text.contains("FACE_OUTER_BOUND"), "no outer bound");
    assert!(text.contains("FACE_BOUND"), "no inner (bore) bound");
}

#[test]
fn hub_rejects_bore_larger_than_hub() {
    let param = DesignParameters {
        center_hole: 0.006, // 6 mm bore vs 5.05 mm hub radius
        ..Default::default()
    };
    let res = step_out::hub_only_step(&param);
    assert!(res.is_err(), "bore larger than the hub should be rejected");
}
