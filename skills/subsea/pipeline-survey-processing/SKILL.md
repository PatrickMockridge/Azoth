# Pipeline survey processing

Educational as-built pipeline survey processing that turns raw survey rows (KP, depth to top of pipe, seabed depth, easting/northing or latitude/longitude) into a cleaned, sign-normalised, resolution-filtered pipeline profile with flagged erroneous points, free-span and cover/burial candidates, a repeat-survey change comparison, a traceable processing log, and an elevation profile a NeqSim pipe model can consume. USE WHEN: a task needs a public, screening-level pipeline profile built from survey or inspection data before detailed DNV-RP-F105 free-span, DNV-RP-F109 on-bottom-stability, DNV-RP-F114 pipe-soil, or NeqSim thermal-hydraulic analysis.

This is an `data-retrieval` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. No azoth backing is planned — this is advisory or data, not a calculation.

## When to Use

Educational as-built pipeline survey processing that turns raw survey rows (KP, depth to top of pipe, seabed depth, easting/northing or latitude/longitude) into a cleaned, sign-normalised, resolution-filtered pipeline profile with flagged erroneous points, free-span and cover/burial candidates, a repeat-survey change comparison, a traceable processing log, and an elevation profile a NeqSim pipe model can consume. USE WHEN: a task needs a public, screening-level pipeline profile built from survey or inspection data before detailed DNV-RP-F105 free-span, DNV-RP-F109 on-bottom-stability, DNV-RP-F114 pipe-soil, or NeqSim thermal-hydraulic analysis.

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
