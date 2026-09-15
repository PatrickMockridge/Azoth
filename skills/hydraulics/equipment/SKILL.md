# Pumps, valves, orifices and choked flow

The equipment-side hydraulics: shaft power, control-valve and orifice flow, and the
choked-flow throat area a relief path needs.

## When to Use

- When a task needs a pump's shaft power, a valve or orifice flow, or a choked-flow
  throat area.

## Inputs

- Pump: density, flow, delivered head, efficiency. Valve: `Cv`, pressure drop,
  specific gravity. Orifice: bore, pressure drop, density, discharge coefficient.
  Choked: mass flow, upstream pressure/density, isentropic exponent.

## Outputs

- `power` (pump), `q` (valve/orifice), `a` (choked throat) — each a quantity.

## How a calculation runs

These are correlations: `pump_power` from flow, head and efficiency;
`control_valve_cv` from the US-convention `Cv`; `orifice_flow` from a supplied
discharge coefficient `Cd`; `choked_flow_area` from the isentropic critical-flow
relation.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
power = azoth.hydraulics.pump_power(q(998.0, "kg/m**3"), q(0.01, "m**3/s"), q(30.0, "m"), 0.75)
power.power  # shaft power, watt
```

## The keycard

`orifice_flow`'s discharge coefficient `Cd` may be omitted when a keycard supplies
`coefficients.hydraulics.orifice_flow.Cd`. An explicit argument always wins; a
missing coefficient with no card is an error, not a default.

## Validation checklist

- [ ] The pump's `eta` is in `(0, 1]`; `H` is delivered head, not including the
      pump's own losses.
- [ ] The orifice `Cd` is supplied, or the keycard supplies it.
- [ ] `choked_flow_area` assumes the flow is choked; that cannot be checked here.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A wrong power | `eta` outside `(0, 1]`, or losses double-counted in `H` | Read the spec's assumption |
| `InvalidInputError` | `orifice_flow` without a `Cd` | Pass `Cd`, or the keycard |
| The wrong throat area | Unchoked flow sized as choked | Confirm choked flow first |

## Limitations

`choked_flow_area` is the isentropic basis, not relief-valve sizing to a standard: the
API 520 de-rating coefficients are the caller's to apply. `control_valve_cv` does not
reproduce IEC 60534's coefficient tables.

## Related Azoth functionality

`hydraulics.pump_power`, `hydraulics.control_valve_cv`, `hydraulics.orifice_flow`,
`hydraulics.choked_flow_area` — pages under `docs/src/hydraulics/`.

## References

- `docs/src/hydraulics/choked_flow_area.md`.
