# Olga multiphase simulator

Run the OLGA transient multiphase flow simulator from Python: locate the installation, generate its PVT table and hydrate equilibrium curve from a NeqSim fluid, validate a genkey case with a rule check, launch the batch engine with the right flags, decode the engine exit code, and read .tpl trend and .ppl profile results. USE WHEN: a task must execute an OLGA case, batch or sweep OLGA runs, build a two- or three-phase OLGA PVT table or a HYDRATECURVE for HYDRATECHECK, diagnose an OLGA failure or licence error, or benchmark OLGA against the NeqSim pipeline models (TwoFluidPipe, PipeBeggsAndBrills) so that a multiphase flow result can be quoted with a known accuracy.

This is a `data-retrieval` skill: it drives OLGA and reads what OLGA writes.
**No tranche backs the simulator**: OLGA is a commercial tool and not NeqSim physics, so azoth will not compute what this runs. The one part of it that is azoth's is the hydrate equilibrium curve, `eos.hydrate_formation_temperature`.

## When to Use

Run the OLGA transient multiphase flow simulator from Python: locate the installation, generate its PVT table and hydrate equilibrium curve from a NeqSim fluid, validate a genkey case with a rule check, launch the batch engine with the right flags, decode the engine exit code, and read .tpl trend and .ppl profile results. USE WHEN: a task must execute an OLGA case, batch or sweep OLGA runs, build a two- or three-phase OLGA PVT table or a HYDRATECURVE for HYDRATECHECK, diagnose an OLGA failure or licence error, or benchmark OLGA against the NeqSim pipeline models (TwoFluidPipe, PipeBeggsAndBrills) so that a multiphase flow result can be quoted with a known accuracy.

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

`eos.hydrate_formation_temperature` computes the hydrate equilibrium temperature a `HYDRATECURVE` is built from. The simulator, the case files and the trend and profile readers are OLGA's and no tranche will back them. The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
