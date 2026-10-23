# Surf cooldown screening

Subsea flowline cooldown and no-touch-time screening with azoth. USE WHEN: a task needs to
estimate how long an insulated flowline stays above its hydrate formation temperature after
shutdown, before a transient thermal study.

This is an `azoth` skill. `eos.surf_cooldown` computes the whole screen: the flash the fluid
density and specific heat come from, the hydrate equilibrium temperature at the operating
pressure, the lumped cooldown, and the crossing of the hydrate temperature plus a margin.

## When to Use

- When a shutdown leaves a wet line cooling towards the seabed and the question is how long
  before hydrates can form.
- When an insulation or U-value choice has to be ranked before a detailed thermal transient.
- When a requirement is already stated and the verdict is what matters: does the line meet it?

## Inputs

- `components`, `z` — the fluid by name and its mole fractions; **water must be one of them**
  unless the screen is deliberately reporting no hydrate risk.
- `operating_pressure`, `flash_temperature` — the state the fluid and the hydrate curve are
  read at.
- `internal_diameter`, `wall_thickness`, `insulation_thickness`, `coating_thickness`,
  `insulation_conductivity` — the pipe.
- `seabed_temperature` — the ambient.
- `overall_u_value` or `external_htc` — a supplied U, or the layer resistances.
- `hydrate_margin`, `required_no_touch_time`, `total_time`, `time_step` — the margin, the
  requirement and the march's resolution.
- `require_water` — set it when a dry fluid must be refused rather than reported as no risk.

## Outputs

- `verdict` — `ok`, `marginal`, `critical` or `no_hydrate_risk`.
- `no_touch_time` — the crossing time; absent when the class's answer is positive infinity.
- `hydrate_equilibrium_temperature`, `water_present`, `water_mole_fraction`.
- `fluid_density`, `fluid_specific_heat`, `overall_u_value`, `time_constant` — what the march
  actually used, so the answer can be checked.

## How a calculation runs

The fluid is flashed at the operating temperature and pressure, and the phases' total mass over
their Peneloux-corrected volume is the density; their mole-weighted molar `cp` over the total
mass is the specific heat. The hydrate equilibrium temperature comes from
`eos.hydrate_formation_temperature` at the operating pressure. `eos.pipeline_cooldown` then
marches one lumped thermal mass (fluid, steel and insulation) by explicit Euler, and the
no-touch time is where the fluid crosses the hydrate temperature plus the margin. The verdict
compares that time with the requirement at its three-quarter point, or with the absolute twelve
and six hours when no requirement is given.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
r = azoth.eos.surf_cooldown(
    ["methane", "ethane", "propane", "water"], [0.82, 0.06, 0.03, 0.09],
    operating_pressure=q(120.0, "bar"), flash_temperature=q(65.0, "degC"),
    internal_diameter=q(0.254, "m"), wall_thickness=q(0.0159, "m"),
    insulation_thickness=q(0.060, "m"), coating_thickness=q(0.006, "m"),
    insulation_conductivity=q(0.22, "W/(m*K)"), seabed_temperature=q(4.0, "degC"),
    overall_u_value=q(2.5, "W/(m**2*K)"), required_no_touch_time=q(8.0, "h"),
)
print(r.verdict, r.no_touch_time.to("h"))
```

## The keycard

The hydrate's guest parameters are the databank's own columns, and a keycard overrides a
component's constants the same way it overrides any other. The steel conductivity, insulation
density and specific heat are the class's constants and are not inputs.

## Validation checklist

- [ ] The fluid carries water, or `water_present` is knowingly false.
- [ ] The operating pressure and the insulation's conductivity are the design values.
- [ ] A `no_hydrate_risk` verdict on a wet line has been read with `water_present` in view, not
      as a proof that the line is safe.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `no_hydrate_risk` on a wet line | The hydrate solve failed, which the class reports as no risk | Read `water_present`, and set `require_water` in a design workflow |
| A no-touch time that ignores the insulation | `overall_u_value` was supplied and the layer inputs were left unused | Leave the U unset to run the layer resistances, or supply the design U |
| A verdict with no requirement | The absolute twelve- and six-hour bands decided | State `required_no_touch_time` |

## Limitations

A lumped single thermal mass: no axial temperature profile, no soil or burial model, and no
distributed insulation defects. The hydrate equilibrium is thermodynamic and inhibitor-free;
the screen says when the line can reach the hydrate curve, not whether hydrates nucleate.

## Related Azoth functionality

`eos.pipeline_cooldown` — the march this screen drives; `eos.hydrate_formation_temperature` —
the equilibrium it compares against; and `azoth-hydrate-margin-check` for the margin alone.

## References

- `docs/src/eos/surf_cooldown.md` and `docs/src/eos/pipeline_cooldown.md`.
- NeqSim: https://github.com/equinor/neqsim — `SurfCooldownAnalyzer`.
