# Compressor antisurge recycle

Set up anti-surge recycle control for a centrifugal compressor in NeqSim, including compressor-chart generation, steady-state AntiSurgeRecycleCalculator use, dynamic AntiSurgeController PI control, and CompressorAntiSurgeApplication topology binding for hot/cold recycle valves and speed runback. USE WHEN: a task needs to protect a NeqSim compressor from surge with a recycle (spill-back) loop and a compressor performance chart is either supplied or must be generated.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. The specification puts **automation** beyond P12, so no tranche backs it.

## When to Use

Set up anti-surge recycle control for a centrifugal compressor in NeqSim, including compressor-chart generation, steady-state AntiSurgeRecycleCalculator use, dynamic AntiSurgeController PI control, and CompressorAntiSurgeApplication topology binding for hot/cold recycle valves and speed runback. USE WHEN: a task needs to protect a NeqSim compressor from surge with a recycle (spill-back) loop and a compressor performance chart is either supplied or must be generated.

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

None yet. The specification puts **automation** beyond P12, so no tranche backs it. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
