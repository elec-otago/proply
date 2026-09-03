// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Triangle-mesh (PLY) export of a designed propeller.
//!
//! The STEP writer ([`crate::step_out`]) serialises the design as a NURBS
//! B-rep that only a CAD kernel can re-mesh; for lightweight headless
//! previews (`render-step`) the same in-memory geometry — the per-station
//! foil outlines wrapped onto the rotor cylinders, plus the hub — is
//! exported here as a plain triangle soup in one ASCII PLY document
//! (assembly frame: the propeller axis is Z, blades rotated around it).
//!
//! The geometry deliberately mirrors [`crate::step_out::write_prop`] and
//! [`BladeElement::get_foil_points`] exactly (metres, not the STEP's
//! millimetres), so the mesh and the STEP always describe the same blade.
//! Winding order is not guaranteed to be outward — the renderer shades
//! two-sided — but every blade and the hub are closed shells.

use crate::prop::Prop;

/// One triangle: three vertex positions (metres, assembly frame).
type Tri = ([f64; 3], [f64; 3], [f64; 3]);

/// The propeller as an ASCII PLY document: vertices + triangle faces.
pub fn ply_prop(prop: &mut Prop, n_points: usize) -> String {
    let mut tris: Vec<Tri> = Vec::new();
    blade_tris(prop, n_points, &mut tris);
    hub_tris(&prop.param, &mut tris);
    to_ply(&tris)
}

/// Triangles of every blade (rotated copies of the single blade solid at
/// azimuth 0), in the assembly frame.
fn blade_tris(prop: &mut Prop, n_points: usize, out: &mut Vec<Tri>) {
    if prop.n_blades == 0 || prop.blade_elements.len() < 2 {
        return;
    }
    // One blade at azimuth 0, built the way step_out does: per-station
    // (lower, upper) outlines from the element's foil, wrapped on the
    // cylinder of the station radius with the scimitar offset applied.
    let radii: Vec<f64> = prop.blade_elements.iter().map(|be| be.r).collect();
    let scims: Vec<f64> = radii.iter().map(|&r| prop.get_scimitar_offset(r)).collect();
    let stations: Vec<(Vec<[f64; 3]>, Vec<[f64; 3]>)> = prop
        .blade_elements
        .iter()
        .zip(scims.iter())
        .map(|(be, s)| be.get_foil_points(n_points, *s))
        .collect();
    let mut single: Vec<Tri> = Vec::new();
    loft_tris(&stations, &mut single);

    // The blade solid lives with its LE on the +X ray (scimitar 0); copy
    // it n_blades times around the axis.
    for k in 0..prop.n_blades {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / prop.n_blades as f64;
        let (s, c) = angle.sin_cos();
        for (a, b, d) in &single {
            let rot = |p: &[f64; 3]| [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
            out.push((rot(a), rot(b), rot(d)));
        }
    }
}

/// Skin one blade solid between its station outlines: the upper and lower
/// loft surfaces between consecutive stations, the thin trailing-edge cap,
/// and fan caps over the root and tip station profiles.
fn loft_tris(stations: &[(Vec<[f64; 3]>, Vec<[f64; 3]>)], out: &mut Vec<Tri>) {
    let n = stations[0].0.len();
    let m = stations.len();
    let quad = |out: &mut Vec<Tri>, a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| {
        out.push((a, b, c));
        out.push((a, c, d));
    };
    for i in 0..m - 1 {
        let (l0, u0) = &stations[i];
        let (l1, u1) = &stations[i + 1];
        // Upper and lower loft surfaces (profiles run LE -> TE on both).
        for j in 0..n - 1 {
            quad(out, u0[j], u0[j + 1], u1[j + 1], u1[j]);
            quad(out, l0[j], l1[j], l1[j + 1], l0[j + 1]);
        }
        // Trailing-edge cap (the two profiles' TE points coincide in
        // chordwise order; connect upper TE to lower TE along the span).
        quad(
            out,
            u0[n - 1],
            u1[n - 1],
            l1[n - 1],
            l0[n - 1],
        );
    }
    // Root and tip caps: fan the closed profile outline (upper LE..TE,
    // then lower TE..LE) around its centroid.
    for idx in [0usize, m - 1] {
        let (l, u) = &stations[idx];
        let mut loop_pts: Vec<[f64; 3]> = Vec::with_capacity(2 * n);
        loop_pts.extend_from_slice(u); // LE .. TE on the upper surface
        for k in (0..n).rev() {
            loop_pts.push(l[k]); // TE .. LE on the lower surface
        }
        let mut ctr = [0.0; 3];
        for p in &loop_pts {
            for d in 0..3 {
                ctr[d] += p[d];
            }
        }
        for d in 0..3 {
            ctr[d] /= loop_pts.len() as f64;
        }
        for k in 0..loop_pts.len() {
            let a = loop_pts[k];
            let b = loop_pts[(k + 1) % loop_pts.len()];
            out.push((a, ctr, b));
        }
    }
}

/// The hub: a closed cylinder (outer radius = hub radius, height =
/// `hub_depth`, centred on z = 0) — the same solid `step_out` writes
/// (minus the mounting-bore hole, invisible at render scale).
fn hub_tris(param: &crate::design_parameters::DesignParameters, out: &mut Vec<Tri>) {
    let r = param.hub_radius + 0.05e-3;
    let (z0, z1) = (-param.hub_depth / 2.0, param.hub_depth / 2.0);
    if r <= 1e-9 || z1 <= z0 {
        return;
    }
    const SEG: usize = 48;
    let ring = |z: f64| -> Vec<[f64; 3]> {
        (0..SEG)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / SEG as f64;
                [r * a.cos(), r * a.sin(), z]
            })
            .collect()
    };
    let bottom = ring(z0);
    let top = ring(z1);
    let quad = |out: &mut Vec<Tri>, a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| {
        out.push((a, b, c));
        out.push((a, c, d));
    };
    for i in 0..SEG {
        let j = (i + 1) % SEG;
        quad(out, bottom[i], top[i], top[j], bottom[j]);
    }
    for (z, ring_pts) in [(z0, &bottom), (z1, &top)] {
        let _ = z;
        let ctr = [0.0, 0.0, ring_pts[0][2]];
        for i in 0..SEG {
            let a = ring_pts[i];
            let b = ring_pts[(i + 1) % SEG];
            out.push((a, b, ctr));
        }
    }
}

/// Serialise the triangle soup as ASCII PLY.
fn to_ply(tris: &[Tri]) -> String {
    let mut verts: Vec<[f64; 3]> = Vec::with_capacity(tris.len() * 3);
    let mut faces: Vec<[u32; 3]> = Vec::with_capacity(tris.len());
    for (a, b, c) in tris {
        let i = verts.len() as u32;
        verts.push(*a);
        verts.push(*b);
        verts.push(*c);
        faces.push([i, i + 1, i + 2]);
    }
    let mut s = String::with_capacity(verts.len() * 40 + faces.len() * 24 + 64);
    s.push_str("ply\nformat ascii 1.0\n");
    s.push_str(&format!("element vertex {}\n", verts.len()));
    s.push_str("property float x\nproperty float y\nproperty float z\n");
    s.push_str(&format!("element face {}\n", faces.len()));
    s.push_str("property list uchar int vertex_indices\nend_header\n");
    for v in &verts {
        s.push_str(&format!("{:.6} {:.6} {:.6}\n", v[0], v[1], v[2]));
    }
    for f in &faces {
        s.push_str(&format!("3 {} {} {}\n", f[0], f[1], f[2]));
    }
    s
}
