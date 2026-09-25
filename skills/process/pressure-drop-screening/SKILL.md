# Pressure-drop screening

Single-phase line pressure drop from the Darcy-Weisbach relation, with the friction
factor implicit or explicit, and the pressure gradient that follows from dividing by
the line's length.

## When to Use

- When a task needs a screening pressure drop for a single-phase line, or a pressure
  gradient to compare against a line-sizing guideline.

## Inputs

- Density, velocity, diameter and viscosity for the Reynolds number.
- A friction factor, length, diameter, density and velocity for the drop.
- For a line with a fluid and a state rather than a density, `process.pipe`, which
  takes a composition, a pressure, a temperature, a length, a diameter and a
  roughness.

## Outputs

- `re` and `regime`; `dp`; and from `process.pipe`, `pressure_drop` and the outlet
  state.

## How a calculation runs

`reynolds_number` gives the regime; the friction factor is implicit (Colebrook-White)
or explicit (Swamee-Jain, Haaland); `darcy_weisbach` gives the drop. The gradient the
guidelines screen against is `dp / L` — that division is arithmetic on the answer, not
a calculation of its own, and no guideline is shipped.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
re = azoth.hydraulics.reynolds_number(
    q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
)
re.re, re.regime  # 149401.2, FlowRegime.TURBULENT

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
r.dp  # 22455.0 Pa
```

## The keycard

The fitting coefficients come from `data/fittings/crane_k_factors.csv`, which ships as
estimated placeholders. Supply real coefficients via a keycard's `fittings` section
before any design use.

## Validation checklist

- [ ] The line is single-phase. The Darcy-Weisbach form takes one density, so applying
      it to a flashing line is a modelling error the arithmetic cannot see.
- [ ] A transitional Reynolds number carries `TRANSITIONAL_FLOW`, and the friction
      factor is indeterminate there.
- [ ] Read `warnings`: omitting viscosity skips the regime check rather than failing.
- [ ] Screen the gradient against a guideline, not the total drop.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `RANGE_CHECK_SKIPPED` | Viscosity omitted | Pass `mu`, or accept the skipped regime check |
| A wrong pressure drop | Placeholder fitting coefficients used | Supply real `fittings` via a keycard |
| The regime wrong | A transitional `re` read as turbulent | Read `regime` |
| A two-phase line screened as one | The relation assumes a single density | Use `process.pipe`, which flashes |

## Limitations

The drop assumes one phase at a constant density; `process.pipe` is the line model
that flashes, and it is a model rather than a calc because a composition vector has
nowhere in the scalar registry to go. The guideline itself — NORSOK P-002, GPSA — is
not shipped, and the shipped fitting coefficients are `estimated_dummy`.

## Related Azoth functionality

`hydraulics.reynolds_number`, `hydraulics.friction_factor_colebrook`,
`hydraulics.friction_factor_swamee_jain`, `hydraulics.friction_factor_haaland`,
`hydraulics.darcy_weisbach`, `hydraulics.crane_k_factors` — pages under
`docs/src/hydraulics/` — and `process.pipe` for a line carrying a fluid.

## References

- `docs/src/hydraulics/darcy_weisbach.md`.
- `docs/src/process/pipe.md`.
