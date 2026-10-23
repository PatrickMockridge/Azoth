# Sand erosion screening

Sand-erosion and remaining-wall-life screening with azoth. USE WHEN: a task needs an estimate
of the sand-erosion rate, the cumulative wall loss and the remaining wall life on a
solids-bearing line before detailed DNV RP O501 work.

This is an `azoth` skill. `eos.erosion_prediction` computes the API RP 14E erosional velocity,
the DNV RP O501 sand-erosion rate at the chosen fitting and material, the cumulative loss over
the design life, the remaining wall and the risk band.

## When to Use

- When a produced-sand rate has to be turned into a wall-loss estimate or a remaining life.
- When the fitting matters: an elbow, tee, blind tee, reducer, choke or weld each carry their
  own geometry factor, impact angle and eroded area.
- When a material change is being compared: the erosion constant, the velocity exponent and
  the wall density differ by material.

## Inputs

- `mixture_density`, `mixture_velocity`, `mixture_viscosity` — the flowing mixture.
- `pipe_diameter`, `wall_thickness`, `corrosion_allowance`, `design_life` — the line.
- `sand_rate` — the produced sand rate, in kg/s.
- `sand_particle_diameter`, `sand_particle_density` — the solids.
- `geometry`, `bend_radius_ratio` — the fitting and, for an elbow, its radius over diameter.
- `pipe_material` — the wall.
- `api_c_factor` — the API RP 14E C factor, 100 by default.

## Outputs

- `erosional_velocity`, `velocity_ratio` — the API RP 14E check.
- `erosion_rate`, `cumulative_erosion`, `remaining_wall_thickness`.
- `within_api_limits`, `within_erosion_limits`, `risk_level` — `low`, `medium`, `high` or
  `critical`.

## How a calculation runs

The erosional velocity is `C / sqrt(rho)` in the class's lb/ft³ and ft/s, converted back to
m/s. The DNV rate is `K F(alpha) Vp^n (dp/dp_ref)^0.2 GF (rho_p/rho_ref)`, where the particle
velocity is the mixture velocity scaled by a Stokes-number ratio, `K` and `n` come from the
material, and the geometry factor, impact angle and eroded area come from the fitting. A
straight pipe takes the class's own empirical branch instead. A rate above 0.1, 0.5 or
1.0 mm/yr - or a velocity ratio above 0.8, 1.0 or 1.5 - sets the risk band.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
r = azoth.eos.erosion_prediction(
    mixture_density=q(100.0, "kg/m**3"), mixture_velocity=q(10.0, "m/s"),
    pipe_diameter=q(0.1524, "m"), wall_thickness=q(0.010, "m"),
    sand_rate=q(50.0 / 86400.0, "kg/s"), sand_particle_diameter=q(0.00025, "m"),
    sand_particle_density=q(2650.0, "kg/m**3"), mixture_viscosity=q(0.001, "Pa*s"),
    geometry="elbow", bend_radius_ratio=1.5, pipe_material="carbon_steel",
    api_c_factor=100.0, corrosion_allowance=q(0.003, "m"),
    design_life=q(25 * 365.25 * 86400.0, "s"),
)
print(r.risk_level, r.remaining_wall_thickness.to("mm"))
```

## The keycard

The material constants are the calculation's own, not databank fields; a different erosion
constant is a different correlation rather than a component override, so this screen does not
read a keycard.

## Validation checklist

- [ ] `sand_rate` is in kg/s and the particle diameter in metres - the class's own units are
      kg/day and mm, and the conversion is the caller's.
- [ ] The wall thickness and corrosion allowance are the design values, not the nominal pipe.
- [ ] The geometry is the worst fitting on the line, not the straight run.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A rate a thousand times too high | The sand rate was passed in kg/day | This API takes kg/s |
| A rate a thousand times too low | The particle diameter was passed in mm | This API takes metres |
| `low` risk with a thin wall | The band reads the rate and the velocity ratio, not the wall | Read `remaining_wall_thickness` and `within_erosion_limits` |

## Limitations

A screening correlation, not a CFD erosion model: the geometry factors, impact angles and
eroded areas are the class's, and the straight-pipe branch is its own empirical reduction.
Particle tracking, sand-size distributions and multiphase flow-pattern effects are not
modelled.

## Related Azoth functionality

`eos.erosion_prediction` — the model page under `docs/src/eos/`; the line's hydraulic state it
is fed from is `eos.pipe_beggs_and_brills` or the caller's own.

## References

- `docs/src/eos/erosion_prediction.md`.
- API RP 14E and DNV-RP-O501.
- NeqSim: https://github.com/equinor/neqsim — `ErosionPredictionCalculator`.
