# Reference fluid synthetic generation

Public helpers to generate representative or synthetic fluid cases from a common reference fluid by adjusting a split/characterization factor, match that factor to measured PVT/separator data, blend well/fluid compositions into a field composition by molar-rate allocation, and — when there is NO PVT report and possibly no sample at all — build a declared best-guess fluid basis from fundamentals and named analogues: reservoir temperature from a provincial geothermal gradient, fluid type from GOR/degAPI bands, GOR and stock-tank gravity interpolated on an analogue depth trend with low/base/high cases, a seed light-ends + C7+ cut composition ready for a NeqSim EOS characterization, a per-parameter provenance and assumption register, and a ranked data-acquisition plan. USE WHEN: a task must calibrate a heavy-end split factor against measurements, produce field-level or per-case representative fluids from a reference model, combine several wells/fluids into one allocated field fluid, or establish a fluid for a discovery or prospect that has no laboratory PVT, before rigorous NeqSim characterization.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P4.

## When to Use

Public helpers to generate representative or synthetic fluid cases from a common reference fluid by adjusting a split/characterization factor, match that factor to measured PVT/separator data, blend well/fluid compositions into a field composition by molar-rate allocation, and — when there is NO PVT report and possibly no sample at all — build a declared best-guess fluid basis from fundamentals and named analogues: reservoir temperature from a provincial geothermal gradient, fluid type from GOR/degAPI bands, GOR and stock-tank gravity interpolated on an analogue depth trend with low/base/high cases, a seed light-ends + C7+ cut composition ready for a NeqSim EOS characterization, a per-parameter provenance and assumption register, and a ranked data-acquisition plan. USE WHEN: a task must calibrate a heavy-end split factor against measurements, produce field-level or per-case representative fluids from a reference model, combine several wells/fluids into one allocated field fluid, or establish a fluid for a discovery or prospect that has no laboratory PVT, before rigorous NeqSim characterization.

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
