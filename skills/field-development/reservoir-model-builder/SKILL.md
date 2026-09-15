# Reservoir model builder

Set up a screening-level reservoir model from whatever data exists, on a data-maturity ladder from a single public headline volume up to a full static-model parameter set, and refine it as data arrives. Enforces a data-first workflow: six modelling ingredients (geometry, petrophysics, fluid, SCAL, contacts, volumes) each have a ranked source ladder, every rung above the one in use must carry a recorded outcome of used/blocked/absent before the model may be built, a blocked source is reported as an access request rather than an acquisition programme, and any downgrade is surfaced with what it gave up so recovery factor can be decomposed into displacement times sweep. USE WHEN: a task needs a reservoir model for a field where only open data is available (for example an NCS field on public resource pages), needs a best-guess structural model when there is NO seismic, log or contact data at all (assumed play-typical trap style, layered stratigraphy, fluid contact and culmination solved to honour published volumes, structure-aware well placement, and a full assumption register), needs volumetrics from area/net pay/porosity/Sw, needs hydrostatic pressure and geothermal temperature defaults from depth, needs a recovery factor and drive mechanism from analogues, needs a well count and productivity index from permeability, needs to turn an ALREADY-BUILT static model (retrieved grid dimensions plus PORO/PERMX/NTG/SATNUM/EQLNUM cell arrays) into validated OPM Flow GRID and PROPS include files while catching the failure modes that produce a deck which runs and is wrong (array length mismatch, zero-based region arrays, fractions stored in percent, PERMZ left equal to PERMX, undefined sentinels) and reconciling the rebuilt in-place volume against a reported P90/P50/P10, or needs a NeqSim SimpleReservoir/WellFlow specification with a provenance trail and a ranked data-acquisition plan.

This is an `advisory` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. No azoth backing is planned — this is advisory or data, not a calculation.

## When to Use

Set up a screening-level reservoir model from whatever data exists, on a data-maturity ladder from a single public headline volume up to a full static-model parameter set, and refine it as data arrives. Enforces a data-first workflow: six modelling ingredients (geometry, petrophysics, fluid, SCAL, contacts, volumes) each have a ranked source ladder, every rung above the one in use must carry a recorded outcome of used/blocked/absent before the model may be built, a blocked source is reported as an access request rather than an acquisition programme, and any downgrade is surfaced with what it gave up so recovery factor can be decomposed into displacement times sweep. USE WHEN: a task needs a reservoir model for a field where only open data is available (for example an NCS field on public resource pages), needs a best-guess structural model when there is NO seismic, log or contact data at all (assumed play-typical trap style, layered stratigraphy, fluid contact and culmination solved to honour published volumes, structure-aware well placement, and a full assumption register), needs volumetrics from area/net pay/porosity/Sw, needs hydrostatic pressure and geothermal temperature defaults from depth, needs a recovery factor and drive mechanism from analogues, needs a well count and productivity index from permeability, needs to turn an ALREADY-BUILT static model (retrieved grid dimensions plus PORO/PERMX/NTG/SATNUM/EQLNUM cell arrays) into validated OPM Flow GRID and PROPS include files while catching the failure modes that produce a deck which runs and is wrong (array length mismatch, zero-based region arrays, fractions stored in percent, PERMZ left equal to PERMX, undefined sentinels) and reconciling the rebuilt in-place volume against a reported P90/P50/P10, or needs a NeqSim SimpleReservoir/WellFlow specification with a provenance trail and a ranked data-acquisition plan.

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

None yet. No azoth backing is planned — this is advisory or data, not a calculation. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
