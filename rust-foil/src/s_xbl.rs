// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Shared BL routines (port of `s_xbl.f90`).

use crate::state::{IVX, Xfoil};

/// Sets the BL Newton system line number corresponding to each BL station.
pub fn iblsys(xf: &mut Xfoil) {
    let mut iv = 0;
    for is in 0..2 {
        for ibl in 2..=xf.nbl[is] as usize {
            iv += 1;
            xf.isys[is][ibl] = iv as i32;
        }
    }

    xf.nsys = iv;
    if xf.nsys > 2 * IVX {
        panic!("IBLSYS: BL system array overflow.");
    }
}
