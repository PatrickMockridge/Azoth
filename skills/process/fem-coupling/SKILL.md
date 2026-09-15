# Fem coupling

Link a NeqSim process simulation and engineering documents to a finite-element model of the solid: heat conduction through a layered wall, transient cooldown, species diffusion in porous rock, and the thermal and pressure stress that follow. Merges P&ID, STID, datasheet, insulation-specification and inspection inputs into a traceable design basis, converts a flashed NeqSim fluid into a film coefficient, a Biot and Fourier number and a mesh and time-step target, solves the layered one-dimensional problem with a dependency-free finite-element solver verified against the closed-form resistance, generates a structured Gmsh mesh in two dimensions or swept into three (revolved pipe or vessel wall, extruded plate or block) and a runnable scikit-fem or FEniCSx case, screens which backend is defensible (scikit-fem, FEniCSx, SfePy, MFEM, OpenSeesPy, PyNite), renders the mesh and the solved field off-screen with PyVista including surface, cut-plane and clipped three-dimensional views, gates the study on discretisation, mesh independence, energy balance and boundary placement, and reduces the field to the U-value, U-multiplier, hot-spot factor and no-touch time a one-dimensional NeqSim model consumes. USE WHEN: a task needs a temperature or stress field inside a solid that a one-dimensional heat-transfer coefficient cannot produce - a local insulation defect, a support or clamp short-circuit, a buried or non-radial soil path, a nozzle or wall discontinuity, a cooldown or thermal-shock transient, diffusion through a porous medium - when a three-dimensional geometry or a rendered field is needed, or when an existing thermal or stress finite-element report must be qualified before its numbers are trusted.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P11.

## When to Use

Link a NeqSim process simulation and engineering documents to a finite-element model of the solid: heat conduction through a layered wall, transient cooldown, species diffusion in porous rock, and the thermal and pressure stress that follow. Merges P&ID, STID, datasheet, insulation-specification and inspection inputs into a traceable design basis, converts a flashed NeqSim fluid into a film coefficient, a Biot and Fourier number and a mesh and time-step target, solves the layered one-dimensional problem with a dependency-free finite-element solver verified against the closed-form resistance, generates a structured Gmsh mesh in two dimensions or swept into three (revolved pipe or vessel wall, extruded plate or block) and a runnable scikit-fem or FEniCSx case, screens which backend is defensible (scikit-fem, FEniCSx, SfePy, MFEM, OpenSeesPy, PyNite), renders the mesh and the solved field off-screen with PyVista including surface, cut-plane and clipped three-dimensional views, gates the study on discretisation, mesh independence, energy balance and boundary placement, and reduces the field to the U-value, U-multiplier, hot-spot factor and no-touch time a one-dimensional NeqSim model consumes. USE WHEN: a task needs a temperature or stress field inside a solid that a one-dimensional heat-transfer coefficient cannot produce - a local insulation defect, a support or clamp short-circuit, a buried or non-radial soil path, a nozzle or wall discontinuity, a cooldown or thermal-shock transient, diffusion through a porous medium - when a three-dimensional geometry or a rendered field is needed, or when an existing thermal or stress finite-element report must be qualified before its numbers are trusted.

## Inputs

To be declared by the calculation that backs this skill.

## Outputs

To be declared by the calculation that backs this skill.

## How a calculation runs

No azoth calculation runs yet. The skill exists so that the task, and its place in
the roadmap, is named rather than lost.

## Python usage pattern

None — the calculation is not ported.

## The keycard

Not applicable until the calculation is ported.

## Validation checklist

- [ ] Confirm the task actually needs this calculation.
- [ ] Use NeqSim's validated engine for real work until the tranche lands.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A placeholder used as a design result | The calculation is not ported | Use NeqSim's validated engine |

## Limitations

Azoth has not ported this calculation; nothing here produces a number.

## Related Azoth functionality

None yet. Azoth backs this at tranche P11. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
