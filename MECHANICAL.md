# Mechanical blade-thickness sizing

How the mechanical stiffness (beam-deflection) calculation works and how
it chooses the airfoil thickness, for the `mech_thickness` design option.

Historically proply sizes the section thickness with a *geometric* power
law: a `p = 0.3` curve running from `hub_depth` at the root down to a
tenth of it at the tip.  That law says nothing about how hard the blade
is being pushed — a 1 N toy prop and a 100 N industrial prop with the
same `hub_depth` get the same sections — and it makes the *hub's*
thickness decide the *blade's* airfoil sections.

The mechanical law replaces that with a stiffness calculation: the blade
is treated as a cantilever beam anchored at the hub, the thrust load is
integrated along the blade, the deflection it would cause is estimated,
and the thickness is chosen so the blade does not deform too much.  The
hub thickness is deliberately not involved.

## The mechanical model

### 1. The load: thrust per unit span

A converged design hands back the thrust of every radial station.  These
element loads are *annular* — the thrust of the whole disk annulus at
that radius.  Dividing by the panel width `dr` and the blade count `B`
gives the load per unit span carried by **one blade**:

```text
q(r) = dT(r) / dr / B        [N/m]
```

`q(r)` acts in the z direction: normal to the rotor disk, which is
exactly the direction thrust tugs the blade out of plane.

### 2. The bending moment

Moment at radius `r` is the integral of the load outboard of `r`,
integrated from the tip inwards (the tip carries nothing, the root
carries everything):

```text
M(r) = ∫_r^R q(s) (s − r) ds        [N·m]
```

### 3. The beam: deflection from the moment — the twist

The blade bends as an Euler–Bernoulli beam fixed at the hub.  The
cross-section is approximated as a solid rectangle of chord `c` and
thickness `t`, and it is **twisted by the local twist `θ(r)`** (the
section's pitch about the spanwise axis).  The deflection is in the
z direction, so the stiffness is the section's second moment about the
rotor-disk chord axis — for the rotated rectangle:

```text
I(r) = c³(r)·t(r)/12 · sin²θ  +  c(r)·t³(r)/12 · cos²θ        [m⁴]
     = c(r)·t(r)/12 · ( (t·cos θ)² + (c·sin θ)² )
```

The two terms are the z-projections of the section's two dimensions,
squared and summed: `t·cos θ` is the **z-component of the thickness**
(flatwise bending), and `c·sin θ` is the **z-component of the chord**
(chordwise bending).  Because the chord is typically much longer than the
thickness, the chord's projection dominates at high twist — the twisted
section's resistance to the z load "tends to the chord", i.e. to its
total z-extent `c·sin θ + t·cos θ`.  This is the "take the twist into
account, as the chord can make the beam stiffer" part of the TODO, and it
is why the untwisted `I = c t³/12` overestimates the thickness a pitched
blade needs.  The curvature follows from the moment and the material's
elastic modulus `E`:

```text
w''(r) = M(r) / (E · I(r)),     w(r_hub) = 0,  w'(r_hub) = 0
```

### 4. Sizing: thickness for a constant-curvature bend

The thickness is laid out so the blade bends with **constant
curvature** — a smooth arc with no curvature concentration anywhere.
Fixing `w'' = κ` and closing the tip deflection at the allowed value
`δ` over the beam length `L = R − r_hub` (an arc of curvature `κ` over
length `L` drops `δ = κ L² / 2`, so `κ = 2δ / L²`), and substituting
`w''` and the twist-aware `I` into the curvature equation, the section
thickness satisfies the cubic

```text
t(r)³·cos²θ  +  t(r)·c(r)²·sin²θ  =  6·M(r)·L² / (E·c(r)·δ)   (=: K)
```

with

```text
L = R − r_hub         (the beam length)
δ = deflection_fraction · R      (the allowed tip deflection)
```

For an untwisted section this collapses to the closed form
`t = (K)^(1/3)`; at high twist the linear (chord-projection) term takes
over, `t ≈ K / (c²·sin²θ)`, and the thickness is governed by the chord,
not the foil thickness.  This is what the "deflection caused by the
thrust integrated along the blade" buys you: the moment distribution
picks the *shape* (thickest where the moment peaks, tapering toward the
tip), the deflection limit picks the *scale*, and the twist decides how
much of the section's stiff chord lies in the load direction.

### 5. The floor

At the tip `M → 0`, so the sizing would taper the section to a
knife edge.  The sized thickness is floored at a minimum fraction of the
**local chord** (default `thickness_floor: 0.06`, the thinnest ARA-D
table).  The floor also keeps the polars inside the foil families'
support.  Wherever the floor binds, the blade is stiffer than the sizing
requires, so the realized deflection is *at most* the allowed one.

## How the calculation impacts the foil thickness

- **The thickness follows the load, not the hub geometry.**  The root
  section is sized by the bending moment the thrust creates there; a
  harder-pushing prop gets a thicker blade, an easy one a thinner blade,
  and `hub_depth` plays no part.
- **The response to load is gentle.**  The cubic root damps everything:
  doubling the thrust (`M` doubles) adds only `2^(1/3) ≈ 1.26`× the
  thickness; a stiffer material (`E` doubled) needs `2^(−1/3) ≈ 0.79`×;
  a looser deflection limit (`δ` doubled, e.g. 5% → 10% of R) also gives
  `0.79`×.  The blade does not chase every load change.
- **Chords stiffen the blade — and twist aims them at the load.**  The
  flatwise term `c t³ cos²θ / 12` means a wide-chord section needs less
  thickness (`c` doubled → absolute thickness × `2^(−1/3)`).  The
  twist-aware term `c³ t sin²θ / 12` goes further: wherever the section
  is pitched, its *chord* projects onto the load direction and carries
  the bending.  On the example below the strongly twisted root sections
  are sized to the floor — the wide, pitched root is already stiff in z,
  so the calculation leaves it at the minimum `t/c`.  (A section's total
  z-extent is `c·sin θ + t·cos θ`; at high twist it tends to the chord.)
- **The tip section is a floor decision.**  The tip carries almost no
  moment, so its thickness is whatever the floor is — the law never
  sizes a useful tip from the load alone.  Raise `thickness_floor` for a
  chunkier tip.
- **It replaces the power law's numbers.**  On the web-default design
  (68 mm radius, three blades, 3 GPa modulus, 5% of R allowed): with the
  twist included, root `t/c ≈ 0.08` (versus 0.71 by the geometric law
  and 0.28 by the untwisted mechanical law), most outboard sections ride
  the `thickness_floor`, and the sized blade is predicted to deflect
  3.00 mm — inside the 3.4 mm allowed — where a blade built to the
  geometric law is predicted to deflect 4.9 mm.
- **The section `t/c` flows into the aerodynamics.**  The sized absolute
  thickness is divided by the design's own chord (the law is stored as a
  radius → `t/c` curve, so the sized section scales exactly with the
  final chord), the station foils are built at that `t/c` — NACA
  4-series and CST set it directly, ARA-D interpolates its tables — and
  the airfoil polars are re-simulated for the new sections.  The design
  is then re-run on the sized law, so the reported torque, thrust and
  geometry all belong to the mechanically sized blade.

## Using it

In the design JSON:

```json
{
    "mech_thickness": true,
    "modulus": "3 GPa",
    "deflection_fraction": 0.05,
    "thickness_floor": 0.06
}
```

or on the command line: `--mech-thickness`, `--modulus <P>`,
`--deflection-fraction <F>`, `--thickness-floor <F>`, and in the web
demo via the **Mechanical thickness (beam sizing)** checkbox in the
Propeller Specifications tab.  The **Elastic modulus (GPa)** field on
the same tab sets `modulus` — raise it for a stronger material and the
mechanical law sizes much thinner foils (the thickness goes as
`E^(−1/3)`): on the example below, root `t/c` drops from ≈ 0.08 at the
nylon/ABS default (3 GPa) to the `thickness_floor` (0.06) everywhere
with carbon fibre (100 GPa).

| key | meaning | unit | default |
| --- | --- | --- | --- |
| `mech_thickness` | use the mechanical law instead of the geometric power law | bool | `false` |
| `modulus` | elastic modulus `E` of the blade material (moulded nylon/ABS ≈ 3 GPa, carbon ≈ 70–150 GPa) | Pa (`"3 GPa"`, `"3000 MPa"`) | `3e9` |
| `deflection_fraction` | allowed tip deflection as a fraction of the prop radius `R` | — | `0.05` |
| `thickness_floor` | minimum section thickness as a fraction of the local chord | — | `0.06` |

The pipeline runs the design once (with the geometric law) to get the
converged station loads, sizes the thickness law from them, and re-runs
the design on the sized law — the reported operating point matches the
mechanically sized blade.  The design summary records what happened:

```yaml
design:
  thickness_law: mechanical
  tip_deflection_mm: 3.002042
```

## Example: the web-default design

Converged station data and the resulting sections (BEM, plate polars,
E = 3 GPa, `deflection_fraction: 0.05`, `thickness_floor: 0.06`), with
the local twist that enters the sizing:

```text
  r (m)    chord (mm)   twist (°)   t (mm)    t/c
 0.0060      8.0         25.4       0.65     0.081
 0.0129      7.9         18.8       0.47     0.060  (floor)
 0.0198      6.9         15.8       0.42     0.060  (floor)
 0.0267      5.4         15.2       0.32     0.060  (floor)
 0.0336      3.7         15.7       0.22     0.060  (floor)
 0.0404      2.4         16.2       0.20     0.084
 0.0473      1.8         15.9       0.23     0.129
 0.0542      1.7         14.3       0.27     0.165
 0.0611      1.7         10.7       0.21     0.123
 0.0680      1.8          5.0       0.11     0.060  (floor)
```

The root is the thickest section (highest moment), but the twist does
most of the work there: at 25° pitch a 8.0 mm chord presents ≈ 3.4 mm of
itself to the z load (its z-extent is `c·sin θ + t·cos θ ≈ 3.5 mm`,
versus `t·cos θ ≈ 0.6 mm` of thickness projection alone), so the sizing
needs only a thin section.  The strongly twisted inboard part of the
blade rides the `thickness_floor`; the mild-twist outboard section
thickens towards mid-span where the chord tapers.  Predicted tip
deflection: 3.00 mm of the 3.40 mm allowed.  A lifting-line design of
the same prop (thinner, more lightly twisted sections) sizes on the same
law, with the floor binding across most of the blade.

## Assumptions and limitations

- **Rectangular section approximation.**  The inertia
  `I = c³t/12·sin²θ + ct³/12·cos²θ` treats the section as a solid
  rectangle twisted by the local pitch; real airfoil sections are
  cambered, tapered in thickness and carry the twist in their geometry,
  so the stiffness estimate is approximate (the TODO asks for an
  *approximated* deflection).  The `t³` and `sin²θ` terms dominate, which
  is what the sizing is sensitive to.
- **Static, single-material beam.**  No fatigue, no centrifugal or
  torsional coupling, no shear deflection; `E` is a single number.  The
  twist enters as the pure section rotation about the spanwise axis;
  the bending/twist coupling of a slender pre-twisted beam is not
  modelled.
- **The load is the design's own.**  The law is sized on the converged
  operating point's thrust distribution and the design is re-run once on
  the sized law — a single refinement, not a fixed-point iteration to
  full aerodynamic–structural convergence (the polars do shift a little
  with thickness; the cubic keeps the effect small).
- **Uniform curvature is a design choice.**  Any thickness *shape*
  scaled to meet the tip-deflection limit satisfies the constraint; the
  constant-curvature layout was picked because it concentrates no
  curvature anywhere.  The floor then takes over where the shape would
  fall below it.
- **Opt-in.**  The geometric power law remains the default so existing
  designs are unchanged; enable the mechanical law with `mech_thickness`.

See [CHANGES.md](CHANGES.md) (2026-08-30, "Mechanical airfoil-thickness
law") for the implementation notes.
