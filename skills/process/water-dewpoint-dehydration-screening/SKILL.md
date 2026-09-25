# Water dewpoint dehydration screening

Educational gas water-content and dehydration screening using the public GPSA Bukacek saturated-water-content correlation. USE WHEN: a task needs a public, screening-level estimate of saturated water content in natural gas, a check against a sales-gas water spec, or a decision on whether a stated stream composition is already DEHYDRATED or merely cooled and knocked out (the saturation test), before detailed dehydration design or before asserting any hydrate/ice finding.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. This is a property correlation, which is tranche P1's tree, and it is not among the ones P1 ported.

## When to Use

Educational gas water-content and dehydration screening using the public GPSA Bukacek saturated-water-content correlation. USE WHEN: a task needs a public, screening-level estimate of saturated water content in natural gas, a check against a sales-gas water spec, or a decision on whether a stated stream composition is already DEHYDRATED or merely cooled and knocked out (the saturation test), before detailed dehydration design or before asserting any hydrate/ice finding.

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

None yet. This is a property correlation, which is tranche P1's tree, and it is not among the ones P1 ported. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
