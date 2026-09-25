# Separator modelling

Educational separator screening indicators for gas load, residence time, and capacity warnings. USE WHEN: a task needs public separator capacity screening without proprietary design methods.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. azoth's `process.separator` is an **equilibrium flash split**; the gas-load, residence-time and capacity indicators this screening wants are not ported.

## When to Use

Educational separator screening indicators for gas load, residence time, and capacity warnings. USE WHEN: a task needs public separator capacity screening without proprietary design methods.

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

None yet. azoth's `process.separator` is an **equilibrium flash split**; the gas-load, residence-time and capacity indicators this screening wants are not ported.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
