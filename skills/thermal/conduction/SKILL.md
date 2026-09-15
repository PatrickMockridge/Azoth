# Steady conduction

The first heat-transfer namespace: steady conduction through a plane wall. Its
existence is the point — it proved the calc pipeline is domain-agnostic rather than
shaped around pipe flow.

## When to Use

- When a task needs the steady heat flow through a plane wall from its conductivity,
  area, temperature difference and thickness.

## Inputs

- Thermal conductivity `k`, area `A`, a temperature *difference* `dT`, and thickness
  `L`.

## Outputs

- `q`, the heat flow, following the sign of `dT`.

## How a calculation runs

`q = k * A * dT / L`, the steady-state conduction relation. `dT` is a difference, not
an absolute temperature: pass `Q(30, "delta_degC")` or `Q(30, "K")`.

## Python usage pattern

```python
import azoth

r = azoth.thermal.conduction_plane_wall(
    azoth.ureg.Quantity(0.5, "W/(m*K)"),
    azoth.ureg.Quantity(2.0, "m**2"),
    azoth.ureg.Quantity(30.0, "K"),
    azoth.ureg.Quantity(0.1, "m"),
)
r.q  # heat flow, watt
```

## The keycard

None: `conduction_plane_wall` takes its four quantities directly and reads no
coefficient.

## Validation checklist

- [ ] `dT` is a difference — `delta_degC` or kelvin — not an absolute temperature.
- [ ] `k`, `A` and `L` are positive.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A wrong sign or magnitude | `dT` given as an absolute temperature | Pass a difference |
| `OutOfRangeError` | A non-positive `k`, `A` or `L` | Check the inputs |

## Limitations

One correlation, one geometry: steady conduction through a plane wall. Transient and
multi-layer cases are not implemented.

## Related Azoth functionality

`thermal.conduction_plane_wall` — page under `docs/src/thermal/`.

## References

- `docs/src/thermal/conduction_plane_wall.md`.
