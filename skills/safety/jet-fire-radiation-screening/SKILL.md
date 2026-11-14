# Jet fire radiation screening

Educational jet-fire thermal-radiation screening that estimates the radiative flux at a distance, or the distance to a target flux, using the public single-point-source flame radiation model. USE WHEN: a task needs a public, screening-level jet-fire radiation flux or a distance to a personnel/equipment thermal limit before detailed fire-consequence analysis.

This is a `screening` skill: it produces no azoth number. The library carries the
calculation as `process.jet_fire` - azoth-id-ok: named to disclaim, because this page runs
it not - and the skill is not promoted to back it until a measurement says so.

## When to Use

Educational jet-fire thermal-radiation screening that estimates the radiative flux at a distance, or the distance to a target flux, using the public single-point-source flame radiation model. USE WHEN: a task needs a public, screening-level jet-fire radiation flux or a distance to a personnel/equipment thermal limit before detailed fire-consequence analysis.

## Inputs

The fuel mass flow, the lower heating value, the radiative fraction, the ambient humidity,
a distance or a target flux. The library's declaration of these is
`specs/models/process/jet_fire.toml`; the screening page does not restate it.

## Outputs

The total radiated power, the transmissivity and incident flux at a distance, and the
distance to a target flux. The library declares the same four.

## How a calculation runs

No azoth calculation runs from this skill. Its screening is `process.jet_fire` - azoth-id-ok: disclaimed, because this page runs it not - and it is not promoted as this skill's backing until a measurement compares them.

## Python usage pattern

None from this skill. The library's own entry point is named below; azoth-id-ok: named to
disclaim, since this page does not call it.

## The keycard

Not applicable: the screening takes its own scalar inputs.

## Validation checklist

- [ ] Confirm the task actually needs this calculation.
- [ ] Use the library's process id for a number, and NeqSim's validated engine where the
  dense-gas or fireball consequence is the question.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A placeholder used as a design result | The calculation is not run from this page | Use the library's process id or NeqSim |

## Limitations

This page produces no number. Its public model is the single-point-source flame radiation;
the library's id also carries the fixed-point inverse and the transmissivity it shares with
`process.bleve` - azoth-id-ok: named to disclaim, because this page runs it not.

## Related Azoth functionality

`process.jet_fire` is the point-source radiation - azoth-id-ok: named to disclaim, because
this skill runs it not. `process.bleve` shares its transmissivity - azoth-id-ok: disclaimed
- `process.pool_fire` is the pool-fire view factor - azoth-id-ok: disclaimed - and
`process.vce` is the explosion overpressure - azoth-id-ok: disclaimed. The validated engine
today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
