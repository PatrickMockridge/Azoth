# Gas dispersion distance screening

Educational Gaussian point-source gas-dispersion screening that estimates the downwind distance to a target concentration (LFL or toxic ppm) using the public Pasquill-Gifford / Briggs rural model. USE WHEN: a task needs a public, screening-level distance-to-LFL or distance-to-ppm estimate for a continuous gas release before detailed dispersion or consequence analysis.

This is a `screening` skill: it produces no azoth number. The library carries the
calculation as `process.gas_dispersion` and `process.gaussian_plume` - azoth-id-ok: named to
disclaim, because this page runs neither - and the skill is not promoted to back them until a
measurement says so.

## When to Use

Educational Gaussian point-source gas-dispersion screening that estimates the downwind distance to a target concentration (LFL or toxic ppm) using the public Pasquill-Gifford / Briggs rural model. USE WHEN: a task needs a public, screening-level distance-to-LFL or distance-to-ppm estimate for a continuous gas release before detailed dispersion or consequence analysis.

## Inputs

The source strength, the effective release height, the wind speed, the Pasquill class, the
terrain category, the receptor's position and the target concentration. The library's
declaration of these is `specs/models/process/gaussian_plume.toml`; the screening page does
not restate it.

## Outputs

The two dispersion coefficients, the concentration at the receptor, the ground centreline
concentration and the downwind distance to the target. The library declares the same five.

## How a calculation runs

No azoth calculation runs from this skill. The library's screening is the process id named
above, and it is not promoted as this skill's backing until a measurement compares them.

## Python usage pattern

None from this skill. The library's own entry points are named below; azoth-id-ok: named to
disclaim, since this page does not call them.

## The keycard

Not applicable: the screening takes its own scalar inputs, and the library's id takes the
released composition through the databank.

## Validation checklist

- [ ] Confirm the task actually needs this calculation.
- [ ] Use the library's process id for a number, and NeqSim's validated engine where the
  dense-gas or toxic endpoint is the question.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A placeholder used as a design result | The calculation is not run from this page | Use the library's process id or NeqSim |

## Limitations

This page produces no number. Its public model is the passive Gaussian plume; the library's
screening id also carries the dense-gas and toxic endpoints. azoth-id-ok: named to disclaim.

## Related Azoth functionality

`process.gas_dispersion` carries the ambient-flash screening with both arms - azoth-id-ok:
named to disclaim, because this skill runs none of them. The plume arithmetic alone is
`process.gaussian_plume` - azoth-id-ok: disclaimed - and `process.probit` is the dose-response
- azoth-id-ok: disclaimed - which this page also runs not. The validated engine today is
NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
