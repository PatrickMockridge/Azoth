# Production network routing

Educational production-network routing screening that routes wells through manifolds and flowlines/risers to a host and rolls up an arrival pressure, and screens multiwell flow regulated by a facility inlet/separator pressure. USE WHEN: a task needs a public, screening-level inflow rate per well, aggregated manifold rates, a platform arrival-pressure roll-up, or pressure-regulated multiwell flow from an inlet/separator and reservoir pressures before detailed NeqSim inflow and multiphase-hydraulics design.

This is an `data-retrieval` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. No azoth backing is planned — this is advisory or data, not a calculation.

## When to Use

Educational production-network routing screening that routes wells through manifolds and flowlines/risers to a host and rolls up an arrival pressure, and screens multiwell flow regulated by a facility inlet/separator pressure. USE WHEN: a task needs a public, screening-level inflow rate per well, aggregated manifold rates, a platform arrival-pressure roll-up, or pressure-regulated multiwell flow from an inlet/separator and reservoir pressures before detailed NeqSim inflow and multiphase-hydraulics design.

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

None yet. No azoth backing is planned — this is advisory or data, not a calculation. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
