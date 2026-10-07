# Gas scrubber sizing screening

Educational gas-scrubber sizing screening using the public Souders-Brown / K-factor relation for vertical separators with a mist-eliminator gas-load check. USE WHEN: a task needs a public, screening-level estimate of the Souders-Brown velocity, required vessel diameter, velocity utilisation, and mist-eliminator load for a vertical gas scrubber before detailed separator design.

This is a `screening` skill and, for the sizing it wants, still a placeholder: azoth's
`process.gas_scrubber` is an **equilibrium flash split** and, given a vessel's diameter and a
design K, reports that vessel's Souders-Brown **capacity utilisation** — but the required
diameter, the K-factor velocity itself and the mist-eliminator load this screening is for are
still not ported.

## When to Use

Educational gas-scrubber sizing screening using the public Souders-Brown / K-factor relation for vertical separators with a mist-eliminator gas-load check. USE WHEN: a task needs a public, screening-level estimate of the Souders-Brown velocity, required vessel diameter, velocity utilisation, and mist-eliminator load for a vertical gas scrubber before detailed separator design.

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

Azoth has not ported the sizing; nothing here produces a required diameter or a mist-eliminator
load. The utilisation is the one of the four numbers this skill wants that the library answers,
and only when the vessel's diameter and its design K are both stated.

## Related Azoth functionality

`process.gas_scrubber`'s `capacity_utilization`, which is the utilisation this skill's checklist
asks for. It is an **equilibrium flash split** otherwise, and the required diameter, the K-factor
velocity and the mist-eliminator load are not ported.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
