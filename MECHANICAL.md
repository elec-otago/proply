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

### 3. The beam: deflection from the moment — the real section, the twist and the camber

The blade bends as an Euler–Bernoulli beam fixed at the hub, with the
section's **real** bending inertia (not the enclosing rectangle) about
the rotor-disk chord axis.  The section of chord `c` and thickness `t`
is twisted by the local twist `θ(r)` (its pitch about the spanwise
axis), so the z-bending inertia of the *rotated* section is its two
principal moments, each weighted by the twist — the z-projections of the
section's two dimensions:

```text
I(r) = i_flat(r)·c(r)·t³(r)/12 · cos²θ  +  i_edge(r)·c³(r)·t(r)/12 · sin²θ
```

`t·cos θ` is the **z-component of the thickness** (flatwise bending) and
`c·sin θ` the **z-component of the chord** (chordwise bending); the
shape factors `i_flat` and `i_edge` scale the two moments from the
rectangle to the **actual foil shape** ([`foil::SectionShape`], computed
from the section's shape points):

- a real section is roughly *half* as stiff as its rectangle even
  symmetric (`i ≈ 0.47` flatwise, `≈ 0.44` edgewise for the NACA
  4-series), so the law sizes ~29% thicker sections than the rectangle
  model;
- **camber raises `i_flat`** — a curved (cambered) section is stiffer
  than a flat one, because its mean line carries area away from the
  bending axis (`i_flat ≈ 0.62` at 2% camber, `≈ 1.35` at 6% camber, at
  12–15% thickness), while `i_edge` barely moves.  This is the "a curved
  blade is stiffer than a flat one" part of the question, quantified.

Because the chord is typically much longer than the thickness, the
chord's projection dominates at high twist — the twisted section's
resistance to the z load "tends to the chord", i.e. to its total
z-extent `c·sin θ + t·cos θ`.  This is the "take the twist into
account, as the chord can make the beam stiffer" part of the TODO.  The
curvature follows from the moment and the material's elastic modulus
`E`:

```text
w''(r) = M(r) / (E · I(r)),     w(r_hub) = 0,  w'(r_hub) = 0
```

### 4. Sizing: thickness for a constant-curvature bend

The thickness is laid out so the blade bends with **constant
curvature** — a smooth arc with no curvature concentration anywhere.
Fixing `w'' = κ` and closing the tip deflection at the allowed value
`δ` over the beam length `L = R − r_hub` (an arc of curvature `κ` over
length `L` drops `δ = κ L² / 2`, so `κ = 2δ / L²`), and substituting
`w''` and the twist- and shape-aware `I` into the curvature equation,
the section thickness satisfies the cubic

```text
i_flat·cos²θ·t(r)³  +  i_edge·sin²θ·c(r)²·t(r)  =  6·M(r)·L² / (E·c(r)·δ)   (=: K)
```

with

```text
L = R − r_hub         (the beam length)
δ = deflection_fraction · R      (the allowed tip deflection)
```

For an untwisted section this collapses to the closed form
`t = (K / (i_flat·cos²θ))^(1/3)` (the shape factor included); at high
twist the linear (chord-projection) term takes over,
`t ≈ K / (i_edge·c²·sin²θ)`, and the thickness is governed by the chord,
not the foil thickness.  This is what the "deflection caused by the
thrust integrated along the blade" buys you: the moment distribution
picks the *shape* (thickest where the moment peaks, tapering toward the
tip), the deflection limit picks the *scale*, the twist decides how much
of the section's stiff chord lies in the load direction, and the camber
decides how stiff the section itself is.

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
  flatwise term means a wide-chord section needs less thickness (`c`
  doubled → absolute thickness × `2^(−1/3)`), and wherever the section
  is pitched its *chord* projects onto the load direction and carries
  the bending (`i_edge·c³ t·sin²θ / 12` — a section's total z-extent is
  `c·sin θ + t·cos θ`; at high twist it tends to the chord).
- **Cambered (curved) sections are stiffer than flat ones.**  The law
  sizes with the real section's inertia about the chord line
  (`i_flat`), and camber raises it: 2% camber at 12% thickness gives
  `i_flat ≈ 0.62` (≈ 10% thinner than the symmetric section), 6% camber
  at 15% gives `≈ 1.35` (≈ 30% thinner).  On the example below, adding
  6% camber drops the whole inboard section to the `thickness_floor`.
- **The tip section is a floor decision.**  The tip carries almost no
  moment, so its thickness is whatever the floor is — the law never
  sizes a useful tip from the load alone.  Raise `thickness_floor` for a
  chunkier tip.
- **It replaces the power law's numbers.**  On the web-default design
  (68 mm radius, three blades, 3 GPa modulus, 5% of R allowed), sized
  with the real NACA sections (zero camber): root `t/c ≈ 0.17` (versus
  0.71 by the geometric law and 0.08 by the rectangle-based mechanical
  law — a real symmetric section is only ≈ 0.47× as stiff as its
  rectangle), tapering to `≈ 0.09–0.23` outboard and the floor at the
  tip; the sized blade closes the tip deflection at 3.38 mm of the
  3.40 mm allowed.  (A blade built to the geometric power law deflects
  only ≈ 2.5 mm under the same beam model — its root `t/c ≈ 0.71`
  overshoots the stiffness budget: it buys far more stiffness than the
  deflection limit needs, at roughly four times the root thickness.)
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
demo via the **Mechanical** tab: the **Mechanical thickness (beam
sizing)** checkbox enables the law, and **Elastic modulus (GPa)**,
**Deflection limit** and **Minimum thickness** set `modulus`,
`deflection_fraction` and `thickness_floor`.  Raise the modulus for a
stronger material and the mechanical law sizes much thinner foils (the
thickness goes as `E^(−1/3)`): on the example below, root `t/c` drops
from ≈ 0.17 at the nylon/ABS default (3 GPa) to the `thickness_floor`
(0.06) with carbon fibre (100 GPa).

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
E = 3 GPa, `deflection_fraction: 0.05`, `thickness_floor: 0.06`), with
the local twist that enters the sizing:

```text
  r (m)    chord (mm)   twist (°)   t (mm)    t/c
 0.0060     11.5         26.6       1.97     0.171
 0.0129     13.6         15.3       1.61     0.118
 0.0198     14.3         11.4       1.35     0.095
 0.0267     13.5         11.1       1.23     0.092
 0.0336     11.3         11.6       1.21     0.107
 0.0404      8.6         11.2       1.22     0.142
 0.0473      6.3          9.3       1.20     0.191
 0.0542      4.9          6.4       1.12     0.229
 0.0611      4.4          3.8       0.76     0.172
 0.0680      4.4          4.3       0.26     0.060  (floor)
```

Sized with the **real NACA sections** (`i_flat ≈ 0.47`, camber 0): the
sections are ≈ 29% thicker than the rectangle model everywhere, the root
— the highest moment, still carrying the twist and camber geometry — is
`t/c ≈ 0.17` (versus 0.71 by the geometric law and 0.08 by the
rectangle-based law), and the deflection closes at 3.38 mm of the
3.40 mm allowed.  Adding 6% camber (`--camber 0.06`) shows the stiffness
of the curved section: `i_flat` jumps to ≈ 1.35 and the whole inboard
blade sizes down onto the `thickness_floor` (root `t/c` 0.06, predicted
deflection 2.77 mm).  A lifting-line design of the same prop sizes on
the same law.

## Assumptions and limitations

- **Section model.**  The inertia uses the real section shape from the
  foil's shape points — the actual airfoil (camber included) scaled by
  the twist — with `i_flat`/`i_edge` relative to the enclosing
  rectangle, so the camber-stiffness is captured (a curved section is
  genuinely stiffer in the model).  The beam itself is still the
  Euler–Bernoulli approximation the TODO asks for: no shear deflection,
  no bending/twist coupling of the slender pre-twisted beam, and the
  section is treated as solid (its internal web-like distribution is
  not modelled).
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
