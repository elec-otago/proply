# TODO for the web based version

Document each change once done in CHANGES.md. Then remove from TODO.md

* Design a mechanism for airfoil thickness. The basic idea should be mechanical, treating the blade as a beam and approximating deflection. Then calculating the deflection caused by the thrust integrated along the blade (in the z direction). The thickness should be chosen to keeping the blade shape from deforming too much. Also take the twist into account as the chord can make the beam stiffer. The hun thickness should not be involved in deciding this.
