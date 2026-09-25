# Utility balance screening

Educational utility-balance screening for instrument air demand, fuel gas Wobbe index, and cooling water duty/flow with a simple capacity margin roll-up (NORSOK U-001 / ISA-7.0.01 style). USE WHEN: a task needs a public, screening-level estimate of instrument air demand, cooling water flow, fuel gas Wobbe index compliance, and utility capacity utilisation before detailed utility-system design.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. The `standards/` tree is **beyond P12** in the specification's order, so no tranche backs it.

## When to Use

Educational utility-balance screening for instrument air demand, fuel gas Wobbe index, and cooling water duty/flow with a simple capacity margin roll-up (NORSOK U-001 / ISA-7.0.01 style). USE WHEN: a task needs a public, screening-level estimate of instrument air demand, cooling water flow, fuel gas Wobbe index compliance, and utility capacity utilisation before detailed utility-system design.

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

None yet. The `standards/` tree is **beyond P12** in the specification's order, so no tranche backs it. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
