# Noise screening

Standards-based gas-valve and restriction noise screening at a stated receiver distance using either a current measured A-weighted level or a conservative pressure-drop energy model. USE WHEN: a task needs noise triage, receiver/workplace assessment, or routing to detailed IEC 60534-8-3 prediction while keeping acoustic-induced-vibration assessment separate.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. The specification puts **mechanical design** beyond P12, so no tranche backs it.

## When to Use

Standards-based gas-valve and restriction noise screening at a stated receiver distance using either a current measured A-weighted level or a conservative pressure-drop energy model. USE WHEN: a task needs noise triage, receiver/workplace assessment, or routing to detailed IEC 60534-8-3 prediction while keeping acoustic-induced-vibration assessment separate.

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

None yet. The specification puts **mechanical design** beyond P12, so no tranche backs it. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
