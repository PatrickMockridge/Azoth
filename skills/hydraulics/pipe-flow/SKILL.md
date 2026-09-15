# Pipe flow

From the Reynolds number through the friction factor to the Darcy-Weisbach pressure
drop, with fitting losses available separately. This is the hydraulics slice a pipe
sizing or pressure-drop task starts from.

## When to Use

- When a task needs a Reynolds number, a friction factor, a straight-pipe pressure
  drop, or fitting resistance coefficients.

## Inputs

- Density, velocity, diameter, viscosity; a friction factor and pipe geometry for the
  pressure drop; a list of fitting ids and a `f_t` for the coefficients.

## Outputs

- `re` and `regime`, the friction factor `f`, `dp`, and the fitting `K` factors.

## How a calculation runs

`reynolds_number` feeds every other calc; the friction factor is implicit
(Colebrook-White) or explicit (Swamee-Jain, Haaland). Pipe *with* fittings is a
composition done by the `azoth pipe` CLI, not a calc of its own — the two losses use
different methods and adding them is a modelling decision the caller should see.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
re = azoth.hydraulics.reynolds_number(
    q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
)
r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
r.dp  # 22455 Pa, with RANGE_CHECK_SKIPPED because viscosity was omitted
```

## The keycard

The fitting coefficients come from `data/fittings/crane_k_factors.csv`, which ships as
estimated placeholders. Supply real coefficients via a keycard's `fittings` section
before any design use.

## Validation checklist

- [ ] Omit viscosity and the regime check is skipped — read `warnings`, or pass `mu`.
- [ ] A transitional Reynolds number carries `TRANSITIONAL_FLOW`; the friction factor
      is indeterminate there.
- [ ] Treat the shipped `crane_k_factors.csv` as placeholder, not data.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `RANGE_CHECK_SKIPPED` | Viscosity omitted | Pass `mu`, or accept the skipped check |
| A wrong pressure drop | Placeholder fitting coefficients used | Supply real `fittings` via a keycard |
| The regime wrong | A transitional `re` read as turbulent | Read `regime` |

## Limitations

The `crane_k_factors.csv` values are estimated dummies, wrong by up to a factor of two
while looking reasonable. Choked-flow throat area (`choked_flow_area`) is the
isentropic basis, not relief-valve sizing to a standard.

## Related Azoth functionality

`hydraulics.reynolds_number`, `hydraulics.friction_factor_colebrook`,
`hydraulics.friction_factor_swamee_jain`, `hydraulics.friction_factor_haaland`,
`hydraulics.crane_k_factors`, `hydraulics.darcy_weisbach` — pages under
`docs/src/hydraulics/`.

## References

- `docs/src/hydraulics/darcy_weisbach.md`.
