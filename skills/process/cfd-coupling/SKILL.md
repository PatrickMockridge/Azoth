# Cfd coupling

Link a NeqSim process simulation and engineering documents to a CFD study, single-phase or multiphase. Merges P&ID, STID, datasheet and plant-data inputs into a traceable design basis, converts a flashed NeqSim fluid into CFD boundary conditions, takes both phases and the interfacial tension from a multiphase flash and screens which multiphase model is defensible, writes and runs a complete OpenFOAM case (steady single-phase RANS or transient volume of fluid) for arbitrary geometry, reads the solved fields back, gates the study on wall treatment / mesh independence / turbulence model, and converts local-versus-bulk results into enhancement factors for one-dimensional models. Tonal-noise requests are first gated on source topology, internal geometry, synchronized spectra, event conditions, acoustic terminations, and structural boundaries; steady RANS is never presented as tonal-source diagnosis. USE WHEN: a task needs local flow detail a one-dimensional model cannot generate - velocity or shear peaks at bends, welds, restrictions, tees, headers or tube bundles, flow maldistribution across a bundle or manifold, stratified or slug two-phase behaviour in a line, a pressure-drop check on real geometry, a CFD/FEM report qualification, or a fail-closed readiness assessment for aeroacoustic or flow-induced tonal noise.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P11.

## When to Use

Link a NeqSim process simulation and engineering documents to a CFD study, single-phase or multiphase. Merges P&ID, STID, datasheet and plant-data inputs into a traceable design basis, converts a flashed NeqSim fluid into CFD boundary conditions, takes both phases and the interfacial tension from a multiphase flash and screens which multiphase model is defensible, writes and runs a complete OpenFOAM case (steady single-phase RANS or transient volume of fluid) for arbitrary geometry, reads the solved fields back, gates the study on wall treatment / mesh independence / turbulence model, and converts local-versus-bulk results into enhancement factors for one-dimensional models. Tonal-noise requests are first gated on source topology, internal geometry, synchronized spectra, event conditions, acoustic terminations, and structural boundaries; steady RANS is never presented as tonal-source diagnosis. USE WHEN: a task needs local flow detail a one-dimensional model cannot generate - velocity or shear peaks at bends, welds, restrictions, tees, headers or tube bundles, flow maldistribution across a bundle or manifold, stratified or slug two-phase behaviour in a line, a pressure-drop check on real geometry, a CFD/FEM report qualification, or a fail-closed readiness assessment for aeroacoustic or flow-induced tonal noise.

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
