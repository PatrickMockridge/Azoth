# Near well and injectivity

Derive what the rock will give and take, and hand it to NeqSim: productivity and injectivity indices, their evolution as saturation fronts develop, and the SCAL basis behind them. Standardises on OPM Flow as the reservoir simulator, pyscal for relative permeability and resdata for output, and covers converting a NeqSim compositional fluid into a black-oil (PVTO/PVDG) or gas-condensate (VAPOIL/PVTG/PVDO) PVT table that OPM Flow will actually accept. USE WHEN: a productivity or injectivity index is about to be assumed, injectors must be checked against a voidage requirement, productivity decay through the bubble point matters, a gas condensate needs retrograde dropout and condensate banking represented, a reservoir model must be sized backwards from a mandated production profile, a NeqSim fluid must become a PVTO/PVDG/PVTG/PVTW deck section, or an Eclipse-format reservoir model must be built and run.

This is an `screening` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. Azoth backs this at tranche P6.

## When to Use

Derive what the rock will give and take, and hand it to NeqSim: productivity and injectivity indices, their evolution as saturation fronts develop, and the SCAL basis behind them. Standardises on OPM Flow as the reservoir simulator, pyscal for relative permeability and resdata for output, and covers converting a NeqSim compositional fluid into a black-oil (PVTO/PVDG) or gas-condensate (VAPOIL/PVTG/PVDO) PVT table that OPM Flow will actually accept. USE WHEN: a productivity or injectivity index is about to be assumed, injectors must be checked against a voidage requirement, productivity decay through the bubble point matters, a gas condensate needs retrograde dropout and condensate banking represented, a reservoir model must be sized backwards from a mandated production profile, a NeqSim fluid must become a PVTO/PVDG/PVTG/PVTW deck section, or an Eclipse-format reservoir model must be built and run.

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

None yet. Azoth backs this at tranche P6. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
