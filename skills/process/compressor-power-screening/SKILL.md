# Compressor power screening

Educational centrifugal-compressor power screening using the public polytropic-head equation (API 617 / API 619 style). USE WHEN: a task needs a public, screening-level estimate of polytropic head, discharge temperature, and gas power for a single compression stage before detailed compressor selection.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. azoth's `process.compressor` is the **isentropic step over an efficiency**; the polytropic-head route this screening uses is not ported.

## When to Use

Educational centrifugal-compressor power screening using the public polytropic-head equation (API 617 / API 619 style). USE WHEN: a task needs a public, screening-level estimate of polytropic head, discharge temperature, and gas power for a single compression stage before detailed compressor selection.

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

None yet. azoth's `process.compressor` is the **isentropic step over an efficiency**; the polytropic-head route this screening uses is not ported.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
