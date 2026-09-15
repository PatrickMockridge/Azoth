# E300 fluid io

Read, write, and add water to Eclipse E300 fluid files for NeqSim. USE WHEN: a task needs to load an E300 file into a NeqSim fluid, write a NeqSim fluid to E300, or add water to a fluid or E300 file with the public PVTsim water parameters (volume shift 0.084004, parachor 10.0, kij 0.5).

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P4.

## When to Use

Read, write, and add water to Eclipse E300 fluid files for NeqSim. USE WHEN: a task needs to load an E300 file into a NeqSim fluid, write a NeqSim fluid to E300, or add water to a fluid or E300 file with the public PVTsim water parameters (volume shift 0.084004, parachor 10.0, kij 0.5).

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
