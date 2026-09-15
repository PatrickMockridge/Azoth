# Pvt regression characterization factor

Public weighted multi-target regression of a split/characterization factor against measured PVT and separator data (saturation pressure, GOR, stock-tank-oil density, formation volume factor). USE WHEN: a task must calibrate one heavy-end characterization factor so a fluid model reproduces several measured PVT/separator quantities at once, with per-target weights and residual reporting, before rigorous NeqSim EOS regression.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P4.

## When to Use

Public weighted multi-target regression of a split/characterization factor against measured PVT and separator data (saturation pressure, GOR, stock-tank-oil density, formation volume factor). USE WHEN: a task must calibrate one heavy-end characterization factor so a fluid model reproduces several measured PVT/separator quantities at once, with per-target weights and residual reporting, before rigorous NeqSim EOS regression.

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

None yet. Azoth backs this at tranche P4. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
