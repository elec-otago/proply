//! STEP (AP242) output via the step-io crate.
//!
//! The propeller is written as an assembly: one blade part (a single
//! watertight NURBS solid: upper loft, lower loft, trailing-edge cap and
//! the two end caps) placed `n_blades` times around the hub, plus a hub
//! part with the mounting-hole through-bore.
//!
//! Blade faces: 1 = upper loft (same_sense=false), 2 = lower loft, 3 = TE
//! cap, 4 = hub end cap (same_sense=false), 5 = tip end cap.  Shared edges
//! are traversed in opposite directions by the two adjacent faces, so the
//! shell is a consistently oriented manifold.

use step_io::build::{CurveInput, FaceBoundInput, Frame, SurfaceInput, Vertex};
use step_io::generated::model::EdgeCurveId;

use crate::nurbs::{
    common_knots, common_params, interpolate_with_knots, line, loft, polyline, ruled, NurbsCurve,
    NurbsSurface,
};
use crate::prop::Prop;

const MM: f64 = 1000.0;

type Pt = [f64; 3];

fn scale(p: Pt, s: f64) -> Pt {
    [p[0] * s, p[1] * s, p[2] * s]
}

fn err<E: std::fmt::Display>(e: E) -> String {
    format!("step: {}", e)
}

/// Convert a non-rational curve to the step-io input type.
fn to_step_curve(c: &NurbsCurve) -> step_io::build::NurbsCurve {
    step_io::build::NurbsCurve {
        degree: c.degree,
        control_points: c.control_points.clone(),
        weights: vec![1.0; c.control_points.len()],
        knots: c.knots.clone(),
    }
}

/// Convert a non-rational surface to the step-io input type.
fn to_step_surface(s: &NurbsSurface) -> step_io::build::NurbsSurface {
    let n_u = s.control_points.len();
    let n_v = s.control_points[0].len();
    step_io::build::NurbsSurface {
        degree_u: s.degree_u,
        degree_v: s.degree_v,
        control_points: s.control_points.clone(),
        weights: vec![vec![1.0; n_v]; n_u],
        knots_u: s.knots_u.clone(),
        knots_v: s.knots_v.clone(),
    }
}

/// Build the STEP text for the propeller described by `prop`.
pub fn write_prop(prop: &mut Prop, n_points: usize) -> Result<String, String> {
    let mut b = step_io::StepBuilder::new().map_err(err)?;
    b.header(&step_io::build::HeaderInput {
        file_name: Some(format!("{}.step", prop.param.name)),
        originating_system: Some("proply-rs".into()),
        ..Default::default()
    });

    let assembly = b.part(&prop.param.name).map_err(err)?;
    let blade_part = b.part("blade").map_err(err)?;
    let hub_part = b.part("hub").map_err(err)?;

    // --- Blade solid -------------------------------------------------------
    let radii: Vec<f64> = prop.blade_elements.iter().map(|be| be.r).collect();
    let scim: Vec<f64> = radii
        .iter()
        .map(|r| prop.get_scimitar_offset(*r))
        .collect();
    let stations: Vec<(Vec<Pt>, Vec<Pt>)> = prop
        .blade_elements
        .iter()
        .zip(scim.iter())
        .map(|(be, s)| {
            let (lower, upper) = be.get_foil_points(n_points, *s);
            (
                lower.iter().map(|p| scale(*p, MM)).collect(),
                upper.iter().map(|p| scale(*p, MM)).collect(),
            )
        })
        .collect();
    let n_st = stations.len();
    if n_st < 2 {
        return Err("need at least two blade stations for a STEP solid".into());
    }

    // Common parametrization across every profile (upper and lower).
    let mut all_profiles: Vec<Vec<Pt>> = Vec::with_capacity(2 * n_st);
    for (lower, upper) in &stations {
        all_profiles.push(lower.clone());
        all_profiles.push(upper.clone());
    }
    let u = common_params(&all_profiles);
    let knots = common_knots(n_points, &u);

    let interp = |pts: &[Pt]| -> NurbsCurve {
        NurbsCurve {
            degree: 3,
            control_points: interpolate_with_knots(pts, &u, &knots),
            knots: knots.clone(),
        }
    };

    let lower_curves: Vec<NurbsCurve> = stations.iter().map(|(l, _)| interp(l)).collect();
    let upper_curves: Vec<NurbsCurve> = stations.iter().map(|(_, u)| interp(u)).collect();

    // Shared boundary curves.
    let le_points: Vec<Pt> = stations.iter().map(|(l, _)| l[0]).collect(); // LE_j = lower[0] == upper[0]
    let te_upper_points: Vec<Pt> = stations.iter().map(|(_, u)| u[n_points - 1]).collect();
    let te_lower_points: Vec<Pt> = stations.iter().map(|(l, _)| l[n_points - 1]).collect();

    let le_curve = polyline(&le_points);
    let te_upper_curve = polyline(&te_upper_points);
    let te_lower_curve = polyline(&te_lower_points);
    let tip_te_seg = line(te_upper_points[n_st - 1], te_lower_points[n_st - 1]);
    let hub_te_seg = line(te_lower_points[0], te_upper_points[0]);

    // The five surfaces.
    let upper_loft = loft(&upper_curves);
    let lower_loft = loft(&lower_curves);
    let te_cap = ruled(&te_upper_curve, &te_lower_curve);
    let hub_cap = ruled(&lower_curves[0], &upper_curves[0]);
    let tip_cap = ruled(&lower_curves[n_st - 1], &upper_curves[n_st - 1]);

    // Vertices.
    let v_le: Vec<Vertex> = le_points
        .iter()
        .map(|p| b.vertex(*p).map_err(err))
        .collect::<Result<_, String>>()?;
    let v_te_u: Vec<Vertex> = te_upper_points
        .iter()
        .map(|p| b.vertex(*p).map_err(err))
        .collect::<Result<_, String>>()?;
    let v_te_l: Vec<Vertex> = te_lower_points
        .iter()
        .map(|p| b.vertex(*p).map_err(err))
        .collect::<Result<_, String>>()?;

    // Edges (authored directions).
    let e_le = b.edge(v_le[0], v_le[n_st - 1], CurveInput::Nurbs(to_step_curve(&le_curve))).map_err(err)?;
    let e_te_u = b
        .edge(v_te_u[0], v_te_u[n_st - 1], CurveInput::Nurbs(to_step_curve(&te_upper_curve)))
        .map_err(err)?;
    let e_te_l = b
        .edge(v_te_l[0], v_te_l[n_st - 1], CurveInput::Nurbs(to_step_curve(&te_lower_curve)))
        .map_err(err)?;
    let e_tip_te = b
        .edge(v_te_u[n_st - 1], v_te_l[n_st - 1], CurveInput::Nurbs(to_step_curve(&tip_te_seg)))
        .map_err(err)?;
    let e_hub_te = b
        .edge(v_te_l[0], v_te_u[0], CurveInput::Nurbs(to_step_curve(&hub_te_seg)))
        .map_err(err)?;
    let e_hub_u = b
        .edge(v_le[0], v_te_u[0], CurveInput::Nurbs(to_step_curve(&upper_curves[0].clone())))
        .map_err(err)?;
    let e_hub_l = b
        .edge(v_le[0], v_te_l[0], CurveInput::Nurbs(to_step_curve(&lower_curves[0].clone())))
        .map_err(err)?;
    let e_tip_u = b
        .edge(v_le[n_st - 1], v_te_u[n_st - 1], CurveInput::Nurbs(to_step_curve(&upper_curves[n_st - 1].clone())))
        .map_err(err)?;
    let e_tip_l = b
        .edge(v_le[n_st - 1], v_te_l[n_st - 1], CurveInput::Nurbs(to_step_curve(&lower_curves[n_st - 1].clone())))
        .map_err(err)?;

    // Faces.  Loop directions are CCW when viewed from the exterior.
    let face_upper = b
        .face(
            SurfaceInput::Nurbs(to_step_surface(&upper_loft)),
            false, // natural surface normal points into the solid
            vec![FaceBoundInput::outer(vec![
                (e_le, true),      // hub -> tip
                (e_tip_u, true),   // LE -> TE at the tip
                (e_te_u, false),   // tip -> hub
                (e_hub_u, false),  // TE -> LE at the hub
            ])],
        )
        .map_err(err)?;
    let face_lower = b
        .face(
            SurfaceInput::Nurbs(to_step_surface(&lower_loft)),
            true,
            vec![FaceBoundInput::outer(vec![
                (e_le, false),     // tip -> hub
                (e_hub_l, true),   // LE -> TE at the hub
                (e_te_l, true),    // hub -> tip
                (e_tip_l, false),  // TE -> LE at the tip
            ])],
        )
        .map_err(err)?;
    let face_te = b
        .face(
            SurfaceInput::Nurbs(to_step_surface(&te_cap)),
            true,
            vec![FaceBoundInput::outer(vec![
                (e_te_u, true),    // hub -> tip along the upper TE
                (e_tip_te, true),  // upper -> lower at the tip
                (e_te_l, false),   // tip -> hub along the lower TE
                (e_hub_te, true),  // lower -> upper at the hub
            ])],
        )
        .map_err(err)?;
    let face_hub_cap = b
        .face(
            SurfaceInput::Nurbs(to_step_surface(&hub_cap)),
            false, // natural normal is +r; outward is -r
            vec![FaceBoundInput::outer(vec![
                (e_hub_u, true),    // LE -> TE along the upper profile
                (e_hub_te, false),  // upper -> lower
                (e_hub_l, false),   // TE -> LE along the lower profile
            ])],
        )
        .map_err(err)?;
    let face_tip_cap = b
        .face(
            SurfaceInput::Nurbs(to_step_surface(&tip_cap)),
            true,
            vec![FaceBoundInput::outer(vec![
                (e_tip_l, true),    // LE -> TE along the lower profile
                (e_tip_te, true),   // lower -> upper
                (e_tip_u, false),   // TE -> LE along the upper profile
            ])],
        )
        .map_err(err)?;

    let _blade_solid = b
        .solid(blade_part, "blade", vec![face_upper, face_lower, face_te, face_hub_cap, face_tip_cap])
        .map_err(err)?;

    // --- Hub solid (with the mounting-hole through-bore) --------------------
    build_hub(&mut b, hub_part, &prop.param)?;

    // --- Assembly -----------------------------------------------------------
    let z = [0.0, 0.0, 1.0];
    let x = [1.0, 0.0, 0.0];
    b.place(assembly, hub_part, Frame { origin: [0.0, 0.0, 0.0], axis: z, ref_dir: x })
        .map_err(err)?;
    for k in 0..prop.n_blades {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / prop.n_blades as f64;
        b.place(
            assembly,
            blade_part,
            Frame { origin: [0.0, 0.0, 0.0], axis: z, ref_dir: [angle.cos(), angle.sin(), 0.0] },
        )
        .map_err(err)?;
    }

    b.finish().map_err(|e| format!("step finish: {}", e))
}

/// Build the hub solid: a cylinder (outer side plus two planar caps).
///
/// The Python prop's mounting-hole through-bore is not included: step-io's
/// multi-shell output (rings with inner bounds, or BREP_WITH_VOIDS) is not
/// read back correctly by FreeCAD's OCCT importer (every face arrives as an
/// open shell).  A plain cylinder imports reliably; the bore is a trivial
/// pocket operation in any CAD package.
///
/// The circles are SINGLE full-circle edges shared by the side face and the
/// caps; adjacent faces must reference the same edge entities or the shell
/// does not close.
fn build_hub(
    b: &mut step_io::StepBuilder,
    hub_part: step_io::build::Part,
    param: &crate::design_parameters::DesignParameters,
) -> Result<(), String> {
    let hub_r = (param.hub_radius + 0.05e-3) * MM;
    let hh = param.hub_depth * MM;
    let z0 = -hh / 2.0;
    let z1 = hh / 2.0;

    // One seam vertex per circle level (angle 0).
    let a_bot = b.vertex([hub_r, 0.0, z0]).map_err(err)?;
    let a_top = b.vertex([hub_r, 0.0, z1]).map_err(err)?;

    let circle_frame = |z: f64| Frame {
        origin: [0.0, 0.0, z],
        axis: [0.0, 0.0, 1.0],
        ref_dir: [1.0, 0.0, 0.0],
    };
    let outer_bot_full = b
        .edge(a_bot, a_bot, CurveInput::Circle(circle_frame(z0), hub_r))
        .map_err(err)?;
    let outer_top_full = b
        .edge(a_top, a_top, CurveInput::Circle(circle_frame(z1), hub_r))
        .map_err(err)?;
    let outer_seam = b.edge(a_bot, a_top, CurveInput::Line).map_err(err)?;

    let side = b
        .face(
            SurfaceInput::Cylinder(
                Frame { origin: [0.0, 0.0, z0], axis: [0.0, 0.0, 1.0], ref_dir: [1.0, 0.0, 0.0] },
                hub_r,
            ),
            true,
            vec![FaceBoundInput::outer(vec![
                (outer_bot_full, true),
                (outer_seam, true),
                (outer_top_full, false),
                (outer_seam, false),
            ])],
        )
        .map_err(err)?;

    let cap_bot = b
        .face(
            SurfaceInput::Plane(Frame { origin: [0.0, 0.0, z0], axis: [0.0, 0.0, -1.0], ref_dir: [1.0, 0.0, 0.0] }),
            true,
            vec![FaceBoundInput::outer(vec![(outer_bot_full, false)])],
        )
        .map_err(err)?;

    let cap_top = b
        .face(
            SurfaceInput::Plane(Frame { origin: [0.0, 0.0, z1], axis: [0.0, 0.0, 1.0], ref_dir: [1.0, 0.0, 0.0] }),
            true,
            vec![FaceBoundInput::outer(vec![(outer_top_full, true)])],
        )
        .map_err(err)?;

    b.solid(hub_part, "hub", vec![side, cap_bot, cap_top])
        .map_err(err)?;
    Ok(())
}

/// Write a STEP containing only the hub (for debugging FreeCAD imports).
pub fn hub_only_step(param: &crate::design_parameters::DesignParameters) -> Result<String, String> {
    let mut b = step_io::StepBuilder::new().map_err(err)?;
    b.header(&step_io::build::HeaderInput {
        file_name: Some("hub.step".into()),
        originating_system: Some("proply-rs".into()),
        ..Default::default()
    });
    let hub_part = b.part("hub").map_err(err)?;
    build_hub(&mut b, hub_part, param)?;
    b.finish().map_err(|e| format!("step finish: {}", e))
}

/// Type-check helper (never used): asserts the edge id type used in face
/// bounds matches the builder's return type.
#[allow(dead_code)]
fn _assert_edge_type(_: EdgeCurveId) {}
