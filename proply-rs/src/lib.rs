// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! proply-rs: propeller design in Rust.
//!
//! A port of the Python `proply` propeller design package.  Blade element
//! momentum theory is used to size each radial blade element, and the
//! `rust-foil` crate (a Rust port of XFOIL) supplies airfoil polars.
//! The final propeller is written out as a STEP (AP242) file.

pub mod arad;
pub mod blade_element;
pub mod cache;
pub mod design_parameters;
pub mod foil;
pub mod lift_line;
pub mod motor;
pub mod nurbs;
pub mod optimize;
pub mod pchip;
pub mod polyfit;
pub mod prop;
pub mod simulator;
pub mod smooth;
pub mod step_out;
pub mod units;
pub mod yaml_out;
