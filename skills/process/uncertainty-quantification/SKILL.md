# Uncertainty quantification

Monte Carlo uncertainty quantification, tornado sensitivity and global sensitivity analysis for NeqSim tasks. Supplies inverse-CDF marginals (uniform, triangular, normal, log-normal fitted from a stated P10/P90, deterministic), unit-hypercube samplers (seeded pseudo-random, Latin hypercube, Halton), a technical/economic model split that caches the expensive flowsheet stage so price and cost parameters never trigger a re-solve, linear-interpolation percentiles in the ascending p10<=p50<=p90 convention the task gate enforces, a scale-invariant split-half convergence check, a swing-ranked one-at-a-time tornado, optional SALib Saltelli/Sobol' and Morris backends that expose the interaction effects a tornado is blind to, an optional chaospy polynomial-chaos surrogate that reaches the same statistics in one to two orders of magnitude fewer model evaluations, and emission of the uncertainty block that the task report generator and CI gate consume. USE WHEN: a task must report P10/P50/P90, a tornado diagram or a probability of a negative outcome, a Monte Carlo loop wraps an expensive NeqSim simulation, parameters must be ranked before a sensitivity budget is committed, an interaction between uncertain inputs is suspected, or an existing uncertainty block must be audited for sample count, convergence and percentile convention.

This is a `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. The specification puts `statistics/` **beyond P12**, so no tranche backs it.

## When to Use

Monte Carlo uncertainty quantification, tornado sensitivity and global sensitivity analysis for NeqSim tasks. Supplies inverse-CDF marginals (uniform, triangular, normal, log-normal fitted from a stated P10/P90, deterministic), unit-hypercube samplers (seeded pseudo-random, Latin hypercube, Halton), a technical/economic model split that caches the expensive flowsheet stage so price and cost parameters never trigger a re-solve, linear-interpolation percentiles in the ascending p10<=p50<=p90 convention the task gate enforces, a scale-invariant split-half convergence check, a swing-ranked one-at-a-time tornado, optional SALib Saltelli/Sobol' and Morris backends that expose the interaction effects a tornado is blind to, an optional chaospy polynomial-chaos surrogate that reaches the same statistics in one to two orders of magnitude fewer model evaluations, and emission of the uncertainty block that the task report generator and CI gate consume. USE WHEN: a task must report P10/P50/P90, a tornado diagram or a probability of a negative outcome, a Monte Carlo loop wraps an expensive NeqSim simulation, parameters must be ranked before a sensitivity budget is committed, an interaction between uncertain inputs is suspected, or an existing uncertainty block must be audited for sample count, convergence and percentile convention.

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

None yet. The specification puts `statistics/` **beyond P12**, so no tranche backs it. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
