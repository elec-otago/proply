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
//! Every vertex carries texture coordinates: `u` is the chordwise
//! fraction (0 = leading edge .. 1 = trailing edge), `v` the spanwise
//! fraction (0 = hub .. 1 = tip) along the blade, and `part` is 0 for the
//! hub and 1 for the blades.  `render-step` overlays zebra / grid
//! patterns on these coordinates, whose distortion exposes the
//! station-to-station surface undulation plain shading hides.

use crate::prop::Prop;

/// One station's two profile outlines (lower, upper): chordwise points
/// running leading edge -> trailing edge in metres.
type Outline = (Vec<[f64; 3]>, Vec<[f64; 3]>);

/// A vertex: position (metres, assembly frame) + texture coordinates +
/// part (0 = hub, 1 = blade).
#[derive(Clone, Copy)]
struct V {
    pos: [f64; 3],
    u: f64,
    v: f64,
    part: u8,
}

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
    // the station outlines, used by both the blade solid and the hub sizing.
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
    let (root_chord, root_thick) = root_section_extent(&outlines);
    hub_into(&prop.param, root_chord, root_thick, &mut mesh);
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

/// The hub: a closed cylinder (outer radius = hub radius, height =
/// `hub_depth`, centred on z = 0) — the same solid `step_out` writes
/// (minus the mounting-bore hole, invisible at render scale).
/// Chord (straight-line distance between the root's leading and trailing
/// points) and thickness (z-extent) of the innermost blade section, used
/// to size the hub spinner that hides the roots.
fn root_section_extent(outlines: &[Outline]) -> (f64, f64) {
    if outlines.is_empty() {
        return (0.0, 0.0);
    }
    let (lower, upper) = &outlines[0];
    let n = lower.len();
    let chord = if n >= 2 {
        let d = [
            lower[0][0] - lower[n - 1][0],
            lower[0][1] - lower[n - 1][1],
            lower[0][2] - lower[n - 1][2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    } else {
        0.0
    };
    let mut zmin = f64::INFINITY;
    let mut zmax = f64::NEG_INFINITY;
    let mut any = false;
    for p in lower.iter().chain(upper.iter()) {
        zmin = zmin.min(p[2]);
        zmax = zmax.max(p[2]);
        any = true;
    }
    let thick = if any { zmax - zmin } else { 0.0 };
    (chord, thick)
}

/// The hub rendered as a rounded (barrel) spinner rather than the true
/// slim cylinder: a real hub visually contains the blade roots, and the
/// design's hub radius equals the blade-root radius, so a bare cylinder
/// at that radius leaves the wrapped blade roots visible around and
/// behind it.  Size the spinner from the actual root section so it
/// occludes the roots without ballooning: max radius a little past the
/// root's wrapped chord, half-height covering the root thickness (and the
/// true hub depth), with a rounded (superellipse) profile and closed
/// ends.
fn hub_into(param: &crate::design_parameters::DesignParameters, root_chord: f64, root_thick: f64, mesh: &mut Mesh) {
    let base_r = param.hub_radius + 0.05e-3;
    if base_r <= 1e-9 {
        return;
    }
    // Enough to cover the root section (wrapped on the base_r cylinder,
    // so a slightly larger body of revolution hides it from outside).
    let r_max = base_r + 0.55 * root_chord.max(root_thick * 2.0);
    let hz = (0.5 * param.hub_depth).max(0.6 * root_thick).max(0.002);
    const SEG: usize = 56;
    const RINGS: usize = 12;
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for k in 0..=RINGS {
        let t = -1.0 + 2.0 * k as f64 / RINGS as f64; // -1 .. 1
        let z = t * hz;
        // Rounded profile: bulge to r_max mid-height, taper to base_r at
        // the ends (a superellipse keeps the shoulders full enough to be a
        // hub, not a lens).
        let width = (1.0 - t.abs().powf(2.5)).max(0.0);
        let r = base_r + (r_max - base_r) * width.powf(0.5);
        let ring: Vec<usize> = (0..SEG)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / SEG as f64;
                mesh.push(V {
                    pos: [r * a.cos(), r * a.sin(), z],
                    u: 0.0,
                    v: 0.0,
                    part: 0,
                })
            })
            .collect();
        rings.push(ring);
    }
    for w in rings.windows(2) {
        for i in 0..SEG {
            let j = (i + 1) % SEG;
            let (a, b, c, d) = (w[0][i], w[0][j], w[1][j], w[1][i]);
            mesh.quad(a, b, c, d);
        }
    }
    // Close the two ends with fan caps at the (tiny) end rings.
    for ring in [&rings[0], &rings[RINGS]] {
        let z = mesh.verts[ring[0]].pos[2];
        let cc = mesh.push(V {
            pos: [0.0, 0.0, z],
            u: 0.0,
            v: 0.0,
            part: 0,
        });
        for i in 0..SEG {
            let a = ring[i];
            let b = ring[(i + 1) % SEG];
            mesh.tri(a, b, cc);
        }
    }
}
