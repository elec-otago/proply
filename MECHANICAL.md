# Mechanical blade-thickness sizing

How the mechanical strength calculation works and how it chooses the
airfoil thickness, for the `mech_thickness` design option.

Historically proply sizes the section thickness with a *geometric* power
law: a `p = 0.3` curve running from `hub_depth` at the root down to a
tenth of it at the tip.  That law says nothing about how hard the blade
is being pushed — a 1 N toy prop and a 100 N industrial prop with the
same `hub_depth` get the same sections — and it makes the *hub's*
thickness decide the *blade's* airfoil sections.

The mechanical law replaces that with a strength calculation: the blade
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

### 3. The beam: deflection from the moment

The blade bends as an Euler–Bernoulli beam fixed at the hub.  The
cross-section is approximated as a solid rectangle of width `c` (the
chord) and height `t` (the airfoil thickness), so the second moment of
area about the flapping axis is

```text
I(r) = c(r) · t(r)³ / 12             [m⁴]
```

and the curvature follows from the moment and the material's elastic
modulus `E`:

```text
w''(r) = M(r) / (E · I(r)),     w(r_hub) = 0,  w'(r_hub) = 0
```

The **chord enters the stiffness** (`I ∝ c · t³`) — this is the "the
chord can make the beam stiffer" part of the TODO.  A wide section is
stiff even when thin.  The twist enters through the *design's* chord
distribution: both design loops build their chords with the twist
(the geometric chord cap divides by `cos(twist)`, and the lifting-line
loop optimises the chord against the twisted flow), so the chord curve
that stiffens the blade already reflects the twist.

### 4. Sizing: thickness for a constant-curvature bend

The thickness is laid out so the blade bends with **constant
curvature** — a smooth arc with no curvature concentration anywhere.
Fixing `w'' = κ` and closing the tip deflection at the allowed value
`δ` over the beam length `L = R − r_hub` (an arc of curvature `κ` over
length `L` drops `δ = κ L² / 2`, so `κ = 2δ / L²`), and substituting
`w''` and `I` into the curvature equation, the thickness comes out
closed form:

```text
t(r) = ( 6 · M(r) · L² / (E · c(r) · δ) )^(1/3)
```

with

```text
L = R − r_hub         (the beam length)
δ = deflection_fraction · R      (the allowed tip deflection)
```

This is what the "deflection caused by the thrust integrated along the
blade" buys you: the moment distribution picks the *shape* (thickest
where the moment peaks, tapering toward the tip), and the deflection
limit picks the *scale*.

### 5. The floor

At the tip `M → 0`, so the closed form would taper the section to a
knife edge.  The sized thickness is floored at a minimum fraction of the
**local chord** (default `thickness_floor: 0.06`, the thinnest ARA-D
table).  The floor also keeps the polars inside the foil families'
support.  Wherever the floor binds, the blade is stiffer than the closed
form, so the realized deflection is *at most* the allowed one.

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
- **Chords stiffen the blade.**  `I ∝ c t³` means a wide-chord section
  needs less thickness: with `c` doubled the sized absolute thickness
  drops by `2^(−1/3)`, and the `t/c` *fraction* by `2^(−4/3)`.  On
  lifting-line designs, where the mid-span chord bulges, the floor binds
  over that whole bulge (see the example below) — the calculation
  explicitly says "this blade is already stiff there".
- **The tip section is a floor decision.**  The tip carries almost no
  moment, so its thickness is whatever the floor is — the law never
  sizes a useful tip from the load alone.  Raise `thickness_floor` for a
  chunkier tip.
- **It replaces the power law's numbers.**  On the web-default design
  (68 mm radius, three blades, 3 GPa modulus, 5% of R allowed): root
  `t/c ≈ 0.28` (≈ 1.8 mm) versus 0.71 by the geometric law, and the
  sized blade is predicted to deflect 3.38 mm — inside the 3.4 mm
  allowed — where a blade built to the geometric law is predicted to
  deflect 4.9 mm.  By the beam model the old law is over-thick at the
  root and under-stiff overall.
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
`E^(−1/3)`): on the example below, root `t/c` drops from 0.28 at the
nylon/ABS default (3 GPa) to 0.09 with carbon fibre (100 GPa), with the
mid- and outer sections landing on the `thickness_floor`.

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
  tip_deflection_mm: 3.379012
```

## Example: the web-default design

Converged station data and the resulting sections (BEM, plate polars,
E = 3 GPa, `deflection_fraction: 0.05`, `thickness_floor: 0.06`):

```text
  r (m)    chord (mm)   thrust (N)   t (mm)    t/c
 0.0060      6.6          0.006      1.84     0.278
 0.0129      7.0          0.034      1.46     0.208
 0.0198      7.0          0.078      1.17     0.166
 0.0267      6.6          0.126      0.96     0.146
 0.0336      5.8          0.167      0.81     0.141
 0.0404      4.8          0.191      0.72     0.151
 0.0473      4.0          0.190      0.69     0.175
 0.0542      3.4          0.168      0.66     0.195
 0.0611      3.1          0.144      0.46     0.148
 0.0680      3.1          0.114      0.19     0.060  (floor)
```

Thickness peaks at the root (highest moment), decays towards the tip and
lands on the floor at the tip where the moment vanishes.  Predicted tip
deflection: 3.38 mm of the 3.40 mm allowed.  A lifting-line design of
the same prop widens the mid-span chord enough that the floor binds
across `r ≈ 0.03–0.04 m` — the wide chord is already stiff, so the
calculation leaves the section at the minimum (`t/c = 0.06`).

## Assumptions and limitations

- **Rectangular section approximation.**  `I = c t³ / 12` treats the
  section as a solid rectangle; real airfoil sections are cambered and
  tapered in thickness, so the stiffness estimate is approximate (the
  TODO asks for an *approximated* deflection).  The `t³` term dominates,
  which is what the sizing is sensitive to.
- **Static, single-material beam.**  No fatigue, no centrifugal or
  torsional coupling, no shear deflection; `E` is a single number.
- **The load is the design's own.**  The law is sized on the converged
  operating point's thrust distribution and the design is re-run once on
  the sized law — a single refinement, not a fixed-point iteration to
  full aerodynamic–structural convergence (the polars do shift a little
  with thickness; the cubic root keeps the effect small).
- **Uniform curvature is a design choice.**  Any thickness *shape*
  scaled to meet the tip-deflection limit satisfies the constraint; the
  constant-curvature layout was picked because it concentrates no
  curvature anywhere.  The floor then takes over where the shape would
  fall below it.
- **Opt-in.**  The geometric power law remains the default so existing
  designs are unchanged; enable the mechanical law with `mech_thickness`.

See [CHANGES.md](CHANGES.md) (2026-08-30, "Mechanical airfoil-thickness
law") for the implementation notes.
