# Benchmark reference data

Turn the mandatory benchmark-validation step of a NeqSim task into a reproducible, source-traceable comparison. Supplies a registry of independent reference sources with authority tiers, validated ranges and stated uncertainties (IAPWS-95, IAPWS-IF97, Span-Wagner CO2, Setzmann-Wagner methane, Span nitrogen, Bucker-Wagner ethane, Lemmon propane, GERG-2008, CoolProp HEOS, NIST WebBook), an offline anchor table of published critical points, triple points, boiling points and one ambient liquid density that runs with no optional dependency and no network, an optional CoolProp backend for reference values at any state, and a comparison layer that grades PASS/WARN/FAIL, rejects a reference that does not outrank the model basis, records whether the deviation is inside the reference's own uncertainty, enforces the three-point minimum, and emits the exact benchmark_validation block that the task report generator and CI gate consume. USE WHEN: a task must validate NeqSim output against independent reference data, a benchmark notebook is being written, a benchmark_validation block must be produced for results.json, a reported deviation must be traced to a citable source, or an existing benchmark claim must be checked for independence and resolution before it is trusted.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P4.

## When to Use

Turn the mandatory benchmark-validation step of a NeqSim task into a reproducible, source-traceable comparison. Supplies a registry of independent reference sources with authority tiers, validated ranges and stated uncertainties (IAPWS-95, IAPWS-IF97, Span-Wagner CO2, Setzmann-Wagner methane, Span nitrogen, Bucker-Wagner ethane, Lemmon propane, GERG-2008, CoolProp HEOS, NIST WebBook), an offline anchor table of published critical points, triple points, boiling points and one ambient liquid density that runs with no optional dependency and no network, an optional CoolProp backend for reference values at any state, and a comparison layer that grades PASS/WARN/FAIL, rejects a reference that does not outrank the model basis, records whether the deviation is inside the reference's own uncertainty, enforces the three-point minimum, and emits the exact benchmark_validation block that the task report generator and CI gate consume. USE WHEN: a task must validate NeqSim output against independent reference data, a benchmark notebook is being written, a benchmark_validation block must be produced for results.json, a reported deviation must be traced to a citable source, or an existing benchmark claim must be checked for independence and resolution before it is trusted.

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

None yet. Azoth backs this at tranche P4. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
