# Hydrate margin check

Educational hydrate operating-margin screening placeholder with public assumptions. USE WHEN: a task needs a quick check of whether an operating point keeps a safe temperature margin above a hydrate equilibrium temperature and should be directed to validated NeqSim methods for real calculations.

This is an `azoth` skill. It computes the hydrate equilibrium temperature with the ported
`eos.hydrate_formation_temperature` and takes the margin from it; what it does not compute is
anything the operating point itself needs — a flowing temperature, a heat-transfer path or an
inhibitor's effect.

## When to Use

- When an operating point needs a temperature margin above the hydrate curve.
- When the question is which *structure* forms, because the margin differs between them.
- When a produced-water or wet-gas line is being screened before a full flow-assurance study.

## Inputs

- `components` — the substances by name, and **water must be one of them**.
- `P` — the pressure the margin is measured at.
- `z` — the overall mole fractions, summing to one.
- `eos` — the cubic, `srk` or `pr`.

## Outputs

- `temperature` — the hydrate equilibrium temperature at that pressure.
- `structure` — `structure_i` or `structure_ii`, whichever the fill is more stable in.
- `iterations` and `residual` — the bisection's diagnostics.

## How a calculation runs

The hydrate's fugacity is written as `phi_w(T, P)` against the fluid's own water fugacity,
and the equilibrium temperature is where the two meet: `eos.hydrate_formation_temperature`
bisects on temperature, with the structure chosen by comparing the two structures'
*exponents* so that the fluid's water fugacity cancels exactly.

**The margin is subtraction, not physics.** `T_operating - T_hydrate(P)` is the whole of it,
and its two halves are only as good as the operating temperature and the pressure it is
quoted at.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
components = ["methane", "ethane", "propane", "water"]
z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]

equilibrium = azoth.eos.hydrate_formation_temperature(components, P=q(100.0, "bar"), z=z, eos="srk")
print(equilibrium.structure)  # structure_ii
margin = q(285.0, "K") - equilibrium.temperature
print(margin.to("delta_degC"))  # the subcooling, if the margin is negative
```

## The keycard

The guest parameters — the Langmuir constants, the cavity radii, the Lennard-Jones pair — are
the databank's own columns, resolved from the names, and a keycard overrides a component's
constants the same way it overrides any other.

## Validation checklist

- [ ] Water is in the feed; a fluid without it has no formation temperature and is refused.
- [ ] The pressure the margin is quoted at is the pressure at the *point of interest*.
- [ ] Read `structure`: a margin for the wrong structure is a margin for a hydrate that is
      not the one forming.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `InvalidInputError: no water` | The feed has no water component | Add it; a hydrate is water's |
| A margin quoted without a pressure | The equilibrium temperature is a function of both | State the pressure it belongs to |
| A negative margin read as safe | The sign | Subcooled means *inside* the hydrate region |

## Limitations

An inhibitor, a kinetic subcooling allowance and a salinity correction are not modelled. The
equilibrium is thermodynamic: it says where a hydrate can form, not how long it takes, and a
line sitting inside the region with no nucleation is a real outcome this cannot see.

## Related Azoth functionality

`eos.hydrate_formation_temperature` — the model page under `docs/src/eos/`. Its inverse is
`eos.hydrate_formation_pressure`, and how much forms once inside the region is
`eos.hydrate_fraction`.

## References

- `docs/src/eos/hydrate_formation_temperature.md`, and the guest data's provenance in
  `databank/manifest.toml`.
- NeqSim: https://github.com/equinor/neqsim — `HydrateFormationTemperatureFlash`.
