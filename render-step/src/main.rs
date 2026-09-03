// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! `render-step`: render a proply-rs propeller design to a PNG, headlessly.
//!
//! The propeller geometry lives in the STEP file's companion triangle mesh
//! (a `.ply` written beside every `.step` by the design pipeline), so this
//! tool takes the STEP file, loads the sibling mesh, and rasterises it
//! with a small software z-buffer — no OpenGL, no CAD kernel, deterministic
//! output on any machine.  The camera looks at the rotor from ~45° off its
//! axis and frames the whole propeller from its bounding box with a margin.
//!
//! Usage:
//! ```text
//! render-step --step build/out/honda_gx35.step [--png out.png] [--size WxH]
//! ```
//! The default output name replaces the `.step` extension with `.png`.

use std::env;
use std::fs::File;
use std::io::BufWriter;

const TILT_DEG: f64 = 45.0; // camera elevation off the rotor axis
const VFOV_DEG: f64 = 32.0; // vertical field of view
const MARGIN: f64 = 0.98; // model half-diagonal -> frame half-height
const SS: usize = 2; // supersampling factor per axis

fn main() {
    let mut step = None;
    let mut png = None;
    let mut w = 1920usize;
    let mut h = 1080usize;
    let mut it = env::args().skip(1);
    while let Some(a) = it.next() {
        let mut val = || it.next().unwrap_or_default();
        match a.as_str() {
            "--step" => step = Some(val()),
            "--png" => png = Some(val()),
            "--size" => {
                let s = val();
                let mut p = s.split('x');
                w = p.next().and_then(|v| v.parse().ok()).unwrap_or(w);
                h = p.next().and_then(|v| v.parse().ok()).unwrap_or(h);
            }
            _ => {
                eprintln!("render-step: unknown argument {a}");
                std::process::exit(2);
            }
        }
    }
    let step = match step {
        Some(s) => s,
        None => {
            eprintln!("render-step: usage: render-step --step <file.step> [--png out.png]");
            std::process::exit(2);
        }
    };
    let png = png.unwrap_or_else(|| {
        let p = std::path::Path::new(&step);
        p.with_extension("png").to_string_lossy().into_owned()
    });

    let ply = {
        let p = std::path::Path::new(&step);
        p.with_extension("ply")
    };
    let (verts, faces) = match load_ply(&ply.to_string_lossy()) {
        Ok(v) => v,
        Err(e) => {
            eprintln!(
                "render-step: cannot load mesh {} (re-run the design with --mesh-file): {}",
                ply.to_string_lossy(),
                e
            );
            std::process::exit(1);
        }
    };
    if faces.is_empty() {
        eprintln!("render-step: mesh {} has no faces", ply.to_string_lossy());
        std::process::exit(1);
    }

    // Frame the camera from the model's bounding box.
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for v in &verts {
        for d in 0..3 {
            lo[d] = lo[d].min(v[d]);
            hi[d] = hi[d].max(v[d]);
        }
    }
    let center = [
        0.5 * (lo[0] + hi[0]),
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
    ];
    let diag = ((0..3).map(|d| (hi[d] - lo[d]).powi(2)).sum::<f64>()).sqrt();

    let (w2, h2) = (w * SS, h * SS);
    let (zbuf, frame) = render(
        &verts,
        &faces,
        center,
        diag,
        w2,
        h2,
    );

    // Box-filter the supersampled frame down to the requested size.
    let out = downsample(&frame, w2, h2, w, h);
    write_png(&png, w, h, &out).unwrap_or_else(|e| {
        eprintln!("render-step: cannot write {}: {}", png, e);
        std::process::exit(1);
    });
    let _ = zbuf;
    eprintln!("render-step: wrote {} ({}x{})", png, w, h);
}

/// Software z-buffer render of the triangle mesh: perspective camera at
/// ~45° off the rotor axis, fitted to the model with a margin, two-sided
/// Lambert shading.  Returns the depth buffer (unused) and the RGB frame.
fn render(
    verts: &[[f64; 3]],
    faces: &[[u32; 3]],
    center: [f64; 3],
    diag: f64,
    w: usize,
    h: usize,
) -> (Vec<f32>, Vec<u8>) {
    let tilt = TILT_DEG.to_radians();
    let s = tilt.sin();
    // Camera-to-model direction: 45° off +Z, in the X-Z diagonal plane.
    let vd = [
        s / std::f64::consts::SQRT_2,
        s / std::f64::consts::SQRT_2,
        tilt.cos(),
    ];
    let fov = VFOV_DEG.to_radians();
    let r = 0.5 * diag;
    let dist = MARGIN * r / (fov / 2.0).tan();
    let eye = [
        center[0] - vd[0] * dist,
        center[1] - vd[1] * dist,
        center[2] - vd[2] * dist,
    ];
    // Camera basis: forward = vd; image up = world +Z projected.
    let zc = vd; // forward
    let up = [0.0, 0.0, 1.0];
    let mut xc = cross(&up, &zc);
    normalize(&mut xc);
    let yc = cross(&zc, &xc);
    let f = (h as f64 / 2.0) / (fov / 2.0).tan();

    let mut zbuf = vec![f32::INFINITY; w * h];
    let mut frame = vec![240u8; w * h * 3]; // background

    let light = {
        let mut l = [0.55, 0.42, 0.72];
        normalize(&mut l);
        l
    };
    let to_cam = |v: &[f64; 3]| -> [f64; 3] {
        [
            dot(v, &xc) - dot(&eye, &xc),
            dot(v, &yc) - dot(&eye, &yc),
            dot(v, &zc) - dot(&eye, &zc),
        ]
    };
    let mut cam: Vec<[f64; 3]> = verts.iter().map(to_cam).collect();
    let _ = &mut cam;

    for face in faces {
        let (a, b, c) = (cam[face[0] as usize], cam[face[1] as usize], cam[face[2] as usize]);
        // Backface cull (view space: +z toward the camera).
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        // No backface culling: the PLY's winding is not guaranteed outward
        // (the renderer shades two-sided instead).
        let _ = cross(&ab, &ac);
        // World-space normal for shading.
        let wa = verts[face[0] as usize];
        let wb = verts[face[1] as usize];
        let wc = verts[face[2] as usize];
        let mut wn = cross(
            &[wb[0] - wa[0], wb[1] - wa[1], wb[2] - wa[2]],
            &[wc[0] - wa[0], wc[1] - wa[1], wc[2] - wa[2]],
        );
        normalize(&mut wn);
        if dot(&wn, &[eye[0] - wa[0], eye[1] - wa[1], eye[2] - wa[2]]) < 0.0 {
            for w in wn.iter_mut() {
                *w = -*w;
            }
        }
        let shade = (0.35 + 0.65 * dot(&wn, &light).max(0.0)).clamp(0.0, 1.0);
        let col = [
            (0.72 * 255.0 * shade) as u8,
            (0.78 * 255.0 * shade) as u8,
            (0.83 * 255.0 * shade) as u8,
        ];

        // Rasterise with edge functions.
        let x0 = ((a[0] / a[2] * f + w as f64 / 2.0).floor().max(0.0)) as usize;
        // (bounds computed per triangle below)
        let proj = |p: &[f64; 3]| -> (f64, f64, f64) {
            (p[0] / p[2] * f + w as f64 / 2.0, h as f64 / 2.0 - p[1] / p[2] * f, p[2])
        };
        let (ax, ay, az) = proj(&a);
        let (bx, by, _bz) = proj(&b);
        let (cx, cy, _cz) = proj(&c);
        let _ = (x0, az);
        let min_x = (ax.min(bx).min(cx).floor().max(0.0)) as usize;
        let max_x = (ax.max(bx).max(cx).ceil().min(w as f64 - 1.0)) as usize;
        let min_y = (ay.min(by).min(cy).floor().max(0.0)) as usize;
        let max_y = (ay.max(by).max(cy).ceil().min(h as f64 - 1.0)) as usize;
        let area = ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)).abs();
        if area < 1e-12 {
            continue;
        }
        for py in min_y..=max_y {
            let row = py * w;
            for px in min_x..=max_x {
                let x = px as f64 + 0.5;
                let y = py as f64 + 0.5;
                let w0 = ((bx - ax) * (y - ay) - (by - ay) * (x - ax)) / area;
                let w1 = ((cx - bx) * (y - by) - (cy - by) * (x - bx)) / area;
                let w2 = 1.0 - w0 - w1;
                if w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9 {
                    continue;
                }
                let z = w0 * az + w1 * _bz + w2 * _cz;
                let i = row + px;
                if z < zbuf[i] as f64 {
                    zbuf[i] = z as f32;
                    let j = i * 3;
                    frame[j] = col[0];
                    frame[j + 1] = col[1];
                    frame[j + 2] = col[2];
                }
            }
        }
    }
    (zbuf, frame)
}

fn downsample(src: &[u8], sw: usize, sh: usize, dw: usize, dh: usize) -> Vec<u8> {
    let (sx, sy) = (sw / dw, sh / dh);
    let mut out = vec![0u8; dw * dh * 3];
    for y in 0..dh {
        for x in 0..dw {
            let mut acc = [0u64; 3];
            for yy in 0..sy {
                for xx in 0..sx {
                    let j = ((y * sy + yy) * sw + x * sx + xx) * 3;
                    for c in 0..3 {
                        acc[c] += src[j + c] as u64;
                    }
                }
            }
            let n = (sx * sy) as u64;
            let j = (y * dw + x) * 3;
            for c in 0..3 {
                out[j + c] = (acc[c] / n) as u8;
            }
        }
    }
    out
}

fn write_png(path: &str, w: usize, h: usize, rgb: &[u8]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let wtr = BufWriter::new(file);
    let mut enc = png::Encoder::new(wtr, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    writer.write_image_data(rgb)?;
    Ok(())
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn normalize(v: &mut [f64; 3]) {
    let n = dot(v, v).sqrt().max(1e-12);
    for w in v.iter_mut() {
        *w /= n;
    }
}

/// Minimal ASCII-PLY reader (vertices + triangle faces).
type Mesh = (Vec<[f64; 3]>, Vec<[u32; 3]>);

fn load_ply(path: &str) -> Result<Mesh, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut nv = 0usize;
    let mut nf = 0usize;
    let mut in_header = true;
    let mut verts: Vec<[f64; 3]> = Vec::new();
    let mut faces: Vec<[u32; 3]> = Vec::new();
    for raw in text.lines() {
        let line = raw.trim();
        if in_header {
            if line == "end_header" {
                in_header = false;
            } else if let Some(rest) = line.strip_prefix("element vertex ") {
                nv = rest.parse().map_err(|_| "bad vertex count".to_string())?;
            } else if let Some(rest) = line.strip_prefix("element face ") {
                nf = rest.parse().map_err(|_| "bad face count".to_string())?;
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        if verts.len() < nv {
            let x: f64 = it.next().ok_or("bad vertex")?.parse().map_err(|_| "bad x")?;
            let y: f64 = it.next().ok_or("bad vertex")?.parse().map_err(|_| "bad y")?;
            let z: f64 = it.next().ok_or("bad vertex")?.parse().map_err(|_| "bad z")?;
            verts.push([x, y, z]);
        } else if faces.len() < nf {
            let k: usize = it.next().ok_or("bad face")?.parse().map_err(|_| "bad count")?;
            let a: u32 = it.next().ok_or("bad face")?.parse().map_err(|_| "bad a")?;
            let b: u32 = it.next().ok_or("bad face")?.parse().map_err(|_| "bad b")?;
            let c: u32 = it.next().ok_or("bad face")?.parse().map_err(|_| "bad c")?;
            if k != 3 {
                return Err(format!("only triangle faces are supported (got {k})"));
            }
            faces.push([a, b, c]);
        } else {
            break;
        }
    }
    if verts.len() != nv || faces.len() != nf {
        return Err(format!(
            "short mesh: {}/{} vertices, {}/{} faces",
            verts.len(),
            nv,
            faces.len(),
            nf
        ));
    }
    Ok((verts, faces))
}
