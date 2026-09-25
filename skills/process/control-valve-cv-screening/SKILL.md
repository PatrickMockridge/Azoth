# Control valve cv screening

Educational control-valve sizing screening that estimates the required flow coefficient (Kv/Cv) and flags choked flow per public IEC 60534-2-1 / ISA-75.01 equations for liquid and gas service. USE WHEN: a task needs a public, screening-level required Cv and a choked-flow flag before detailed control-valve selection.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. azoth's `hydraulics.control_valve_cv` is the **liquid, non-choked** relation and gives the flow *from* a coefficient; the required-coefficient direction and the IEC 60534-2-1 choked-flow criterion are not ported.

## When to Use

Educational control-valve sizing screening that estimates the required flow coefficient (Kv/Cv) and flags choked flow per public IEC 60534-2-1 / ISA-75.01 equations for liquid and gas service. USE WHEN: a task needs a public, screening-level required Cv and a choked-flow flag before detailed control-valve selection.

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

None yet. azoth's `hydraulics.control_valve_cv` is the **liquid, non-choked** relation and gives the flow *from* a coefficient; the required-coefficient direction and the IEC 60534-2-1 choked-flow criterion are not ported.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
