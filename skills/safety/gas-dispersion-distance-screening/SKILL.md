# Gas dispersion distance screening

Educational Gaussian point-source gas-dispersion screening that estimates the downwind distance to a target concentration (LFL or toxic ppm) using the public Pasquill-Gifford / Briggs rural model. USE WHEN: a task needs a public, screening-level distance-to-LFL or distance-to-ppm estimate for a continuous gas release before detailed dispersion or consequence analysis.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P11.

## When to Use

Educational Gaussian point-source gas-dispersion screening that estimates the downwind distance to a target concentration (LFL or toxic ppm) using the public Pasquill-Gifford / Briggs rural model. USE WHEN: a task needs a public, screening-level distance-to-LFL or distance-to-ppm estimate for a continuous gas release before detailed dispersion or consequence analysis.

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
