# Two phase flow regime screening

Educational two-phase flow-regime screening that classifies a horizontal gas-liquid flow pattern from superficial velocities using a simplified public Mandhane-style map and flags slug risk. USE WHEN: a task needs a public, screening-level flow-regime flag (stratified, slug, annular, bubble) and slug-risk triage to pair with slug-flow and flow-induced-vibration work.

This is an `advisory` skill: azoth will not compute it, and the map below is the content.
**No tranche backs it**: the map is NeqSim's `fluidmechanics/`, which azoth does not port, and azoth's own `hydraulics.*` is where a flow-regime calculation belongs if one is written.

## When to Use

Educational two-phase flow-regime screening that classifies a horizontal gas-liquid flow pattern from superficial velocities using a simplified public Mandhane-style map and flags slug risk. USE WHEN: a task needs a public, screening-level flow-regime flag (stratified, slug, annular, bubble) and slug-risk triage to pair with slug-flow and flow-induced-vibration work.

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

None, and none is planned under this name. NeqSim's `fluidmechanics/` is not the port source, so the backing would be azoth's own `hydraulics.*`. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
