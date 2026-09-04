// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Triangle-mesh (PLY) export of a designed propeller.
//!
//! The STEP writer ([`crate::step_out`]) serialises the design as a NURBS
//! B-rep that only a CAD kernel can re-mesh; for lightweight headless
//! previews (`render-step`) the same in-memory geometry — the per-station
//! foil outlines wrapped onto the rotor cylinders, plus the hub — is
//! exported here as an indexed triangle mesh in one ASCII PLY document
//! (assembly frame: the propeller axis is Z, blades rotated around it).
//!
//! The geometry deliberately mirrors [`crate::step_out::write_prop`] and
//! [`BladeElement::get_foil_points`] exactly (metres, not the STEP's
//! millimetres), so the mesh and the STEP always describe the same blade.
//! Winding order is not guaranteed to be outward — the renderer shades
//! two-sided — but every blade and the hub are closed shells.
//!
//! Every vertex carries texture coordinates and a part marker: blades
//! (part 1) have `u` = chordwise fraction (0 = LE .. 1 = TE) and `v` =
//! spanwise fraction (0 = hub .. 1 = tip); the hub (part 0) is a plain
//! cylinder matching the design (radius `hub_radius`, height `hub_depth`,
//! with the `center_hole` mounting bore through it) with `u` = azimuth
//! fraction and `v` = height fraction.  `render-step` uses the two parts
//! to shade the hub with a material and pattern distinct from the blades,
//! and overlays zebra/grid patterns on the blade whose distortion exposes
//! station-to-station undulation.

use crate::prop::Prop;

/// A vertex: position (metres, assembly frame) + texture coordinates +
/// part (0 = hub, 1 = blade).
#[derive(Clone, Copy)]
struct V {
    pos: [f64; 3],
    u: f64,
    v: f64,
    part: u8,
}

/// One station's two profile outlines (lower, upper): chordwise points
/// running leading edge -> trailing edge in metres.
type Outline = (Vec<[f64; 3]>, Vec<[f64; 3]>);

/// Indexed mesh under construction.
struct Mesh {
    verts: Vec<V>,
    tris: Vec<[usize; 3]>,
}

impl Mesh {
    fn new() -> Self {
        Mesh {
            verts: Vec::new(),
            tris: Vec::new(),
        }
    }
    fn push(&mut self, v: V) -> usize {
        self.verts.push(v);
        self.verts.len() - 1
    }
    fn tri(&mut self, a: usize, b: usize, c: usize) {
        self.tris.push([a, b, c]);
    }
    fn quad(&mut self, a: usize, b: usize, c: usize, d: usize) {
        self.tris.push([a, b, c]);
        self.tris.push([a, c, d]);
    }
    fn to_ply(&self) -> String {
        let mut s = String::with_capacity(self.verts.len() * 56 + self.tris.len() * 24 + 96);
        s.push_str("ply\nformat ascii 1.0\n");
        s.push_str(&format!("element vertex {}\n", self.verts.len()));
        s.push_str("property float x\nproperty float y\nproperty float z\n");
        s.push_str("property float u\nproperty float v\nproperty int part\n");
        s.push_str(&format!("element face {}\n", self.tris.len()));
        s.push_str("property list uchar int vertex_indices\nend_header\n");
        for v in &self.verts {
            s.push_str(&format!(
                "{:.6} {:.6} {:.6} {:.6} {:.6} {}\n",
                v.pos[0], v.pos[1], v.pos[2], v.u, v.v, v.part
            ));
        }
        for t in &self.tris {
            s.push_str(&format!("3 {} {} {}\n", t[0], t[1], t[2]));
        }
        s
    }
}

/// The propeller as an ASCII PLY document: indexed vertices + faces.
pub fn ply_prop(prop: &mut Prop, n_points: usize) -> String {
    let mut mesh = Mesh::new();
    let n_st = prop.blade_elements.len();
    // Per-station scimitar offsets (the interpolator needs &mut Prop) and
    // the station outlines, shared by the blade solid.
    let radii: Vec<f64> = prop.blade_elements.iter().map(|be| be.r).collect();
    let scims: Vec<f64> = radii.iter().map(|&r| prop.get_scimitar_offset(r)).collect();
    let outlines: Vec<Outline> = prop
        .blade_elements
        .iter()
        .zip(scims.iter())
        .map(|(be, s)| be.get_foil_points(n_points, *s))
        .collect();
    if n_st >= 2 && prop.n_blades > 0 {
        // Build one blade at azimuth 0, then rotate copies into the mesh.
        for k in 0..prop.n_blades {
            let angle = 2.0 * std::f64::consts::PI * k as f64 / prop.n_blades as f64;
            let (sin, cos) = angle.sin_cos();
            let rot = |p: &[f64; 3]| [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1], p[2]];
            blade_into(&outlines, n_st, n_points, rot, &mut mesh);
        }
    }
    hub_into(&prop.param, &mut mesh);
    mesh.to_ply()
}

/// Skin one blade solid between its station outlines into `mesh`, applying
/// `rot` to every vertex (azimuth 0 = identity).  Outline point `(i, j)`
/// on side `side` (0 = upper, 1 = lower) gets texture coordinates
/// `u = j/(n-1)` (chordwise, LE -> TE) and `v = i/(n_st-1)` (spanwise).
fn blade_into<F: Fn(&[f64; 3]) -> [f64; 3]>(
    outlines: &[Outline],
    n_st: usize,
    n: usize,
    rot: F,
    mesh: &mut Mesh,
) {
    let v_of = |i: usize| i as f64 / (n_st - 1) as f64;
    let u_of = |j: usize| j as f64 / (n - 1) as f64;
    let push = |mesh: &mut Mesh, i: usize, j: usize, side: u8| -> usize {
        let (l, u) = &outlines[i];
        let p = if side == 0 { u[j] } else { l[j] };
        mesh.push(V {
            pos: rot(&p),
            u: u_of(j),
            v: v_of(i),
            part: 1,
        })
    };
    let quad = |mesh: &mut Mesh, a: (usize, usize, u8), b: (usize, usize, u8), c: (usize, usize, u8), d: (usize, usize, u8)| {
        let (ai, aj, as_) = a;
        let (bi, bj, bs) = b;
        let (ci, cj, cs) = c;
        let (di, dj, ds) = d;
        let ia = push(mesh, ai, aj, as_);
        let ib = push(mesh, bi, bj, bs);
        let ic = push(mesh, ci, cj, cs);
        let id = push(mesh, di, dj, ds);
        mesh.quad(ia, ib, ic, id);
    };
    for i in 0..n_st - 1 {
        for j in 0..n - 1 {
            // Upper loft then lower loft (both outlines run LE -> TE).
            quad(mesh, (i, j, 0), (i, j + 1, 0), (i + 1, j + 1, 0), (i + 1, j, 0));
            quad(mesh, (i, j, 1), (i + 1, j, 1), (i + 1, j + 1, 1), (i, j + 1, 1));
        }
        // Trailing-edge cap: upper TE to lower TE along the span.
        quad(
            mesh,
            (i, n - 1, 0),
            (i + 1, n - 1, 0),
            (i + 1, n - 1, 1),
            (i, n - 1, 1),
        );
    }
    // Root and tip caps: fan the closed outline loop (upper LE..TE, then
    // lower TE..LE) around its centroid.
    for i in [0usize, n_st - 1] {
        let mut loop_pts: Vec<(usize, u8)> = Vec::with_capacity(2 * n);
        for j in 0..n {
            loop_pts.push((j, 0));
        }
        for j in (0..n).rev() {
            loop_pts.push((j, 1));
        }
        let mut ctr = [0.0; 3];
        for (j, side) in &loop_pts {
            let (l, u) = &outlines[i];
            let p = if *side == 0 { u[*j] } else { l[*j] };
            for (c, q) in ctr.iter_mut().zip(p.iter()) {
                *c += q;
            }
        }
        for c in ctr.iter_mut() {
            *c /= loop_pts.len() as f64;
        }
        let cc = mesh.push(V {
            pos: rot(&ctr),
            u: 0.5,
            v: v_of(i),
            part: 1,
        });
        for k in 0..loop_pts.len() {
            let (ja, sa) = loop_pts[k];
            let (jb, sb) = loop_pts[(k + 1) % loop_pts.len()];
            let a = push(mesh, i, ja, sa);
            let b = push(mesh, i, jb, sb);
            mesh.tri(a, cc, b);
        }
    }
}

/// The hub: a plain cylinder matching the design data exactly — outer
/// radius `hub_radius`, height `hub_depth`, centred on z = 0, with the
/// mounting-bore hole through it (`center_hole` in the design JSON, or
/// half the hub radius when absent) — the same solid `step_out` writes.
/// Part 0, with its own texture coordinates (u = azimuth fraction around
/// the axis, v = height fraction) so render-step can shade the hub with a
/// material and pattern distinct from the blades.
fn hub_into(param: &crate::design_parameters::DesignParameters, mesh: &mut Mesh) {
    let r = param.hub_radius;
    let (z0, z1) = (-param.hub_depth / 2.0, param.hub_depth / 2.0);
    if r <= 1e-9 || z1 <= z0 {
        return;
    }
    const SEG: usize = 64;
    let ring = |mesh: &mut Mesh, radius: f64, z: f64, v: f64, part: u8| -> Vec<usize> {
        (0..SEG)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / SEG as f64;
                mesh.push(V {
                    pos: [radius * a.cos(), radius * a.sin(), z],
                    u: i as f64 / SEG as f64,
                    v,
                    part,
                })
            })
            .collect()
    };
    // Outer cylindrical side.
    let bottom = ring(mesh, r, z0, 0.0, 0);
    let top = ring(mesh, r, z1, 1.0, 0);
    for i in 0..SEG {
        let j = (i + 1) % SEG;
        mesh.quad(bottom[i], top[i], top[j], bottom[j]);
    }
    // The mounting bore from the design data: a genuine hole only (a bore
    // of essentially the full hub radius would leave nothing of the hub).
    // The bore's wall is its own part (2) so the renderer shades it as
    // the dark cavity the hole is, instead of the lit hub material.
    let bore = param.center_hole();
    let has_bore = bore > 0.0 && bore < 0.99 * r;
    let bore_wall_b = has_bore.then(|| ring(mesh, bore, z0, 0.0, 2));
    let bore_wall_t = has_bore.then(|| ring(mesh, bore, z1, 1.0, 2));
    if let (Some(bb), Some(bt)) = (&bore_wall_b, &bore_wall_t) {
        for i in 0..SEG {
            let j = (i + 1) % SEG;
            mesh.quad(bb[i], bt[i], bt[j], bb[j]);
        }
    }
    // End caps: annular covers from the outer rim down to the bore (a few
    // concentric ring fans for robustness; the innermost ring is the bore
    // circle itself).  Without a bore the caps are filled discs, covered
    // by a centre fan as before.
    let bore_frac = bore / r;
    let mids: Vec<f64> = [0.72, 0.45, 0.2]
        .into_iter()
        .filter(|&f| f > bore_frac + 1.0e-6)
        .collect();
    for (z, outer, bore_cap) in [
        (z0, &bottom, has_bore.then(|| ring(mesh, bore, z0, 0.5, 0))),
        (z1, &top, has_bore.then(|| ring(mesh, bore, z1, 0.5, 0))),
    ] {
        let mut rings: Vec<Vec<usize>> = vec![outer.clone()];
        for &f in &mids {
            rings.push(ring(mesh, r * f, z, 0.5, 0));
        }
        if let Some(bc) = bore_cap {
            rings.push(bc);
        }
        for w in rings.windows(2) {
            for i in 0..SEG {
                let j = (i + 1) % SEG;
                mesh.quad(w[0][i], w[0][j], w[1][j], w[1][i]);
            }
        }
        if !has_bore {
            let cc = mesh.push(V {
                pos: [0.0, 0.0, z],
                u: 0.5,
                v: 0.5,
                part: 0,
            });
            let inner = &rings[rings.len() - 1];
            for i in 0..SEG {
                let a = inner[i];
                let b = inner[(i + 1) % SEG];
                mesh.tri(a, b, cc);
            }
        }
    }
}
