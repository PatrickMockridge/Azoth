# Surf cooldown screening

Educational SURF flowline/riser cooldown and no-touch-time screening placeholder with public assumptions. USE WHEN: a task needs a quick, public estimate of how long an insulated subsea flowline stays above its hydrate formation temperature after shutdown, and should be directed to validated NeqSim methods for real cooldown and hydrate calculations.

This is a `screening` skill with a mixed backing. **Its hydrate half is ported and the cooldown is not**: the temperature it compares against is `eos.hydrate_formation_temperature`, and the cooldown estimate is `pvtsimulation/flowassurance/SurfCooldownAnalyzer`, which **no tranche backs** - it is long-tail work, carried under Tier 4 of the port roadmap rather than behind a P-number.

## When to Use

Educational SURF flowline/riser cooldown and no-touch-time screening placeholder with public assumptions. USE WHEN: a task needs a quick, public estimate of how long an insulated subsea flowline stays above its hydrate formation temperature after shutdown, and should be directed to validated NeqSim methods for real cooldown and hydrate calculations.

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

`eos.hydrate_formation_temperature` is the hydrate temperature this compares against, and it is ported. The cooldown itself is not: NeqSim's is `pvtsimulation/flowassurance/SurfCooldownAnalyzer`, composed of `PipelineCooldownCalculator` and that hydrate temperature, and neither cooldown calculator is ported. They sit in `pvtsimulation/flowassurance/`, but the tranche that directory is mapped to covers the hydrate, wax, asphaltene and scale families and not these - a cooldown is a thermal transient, not a thermodynamic equilibrium, which is why it is long-tail rather than tranche work.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
