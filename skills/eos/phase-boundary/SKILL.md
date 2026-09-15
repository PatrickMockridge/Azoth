# Phase boundaries

Where a mixture first gives off vapour (bubble), first condenses (dew), saturates
(pure component), or reaches its critical point. These are the models that delimit
the two-phase region.

## When to Use

- When a task needs a bubble-point or dew-point pressure at a temperature.
- When a pure component's saturation pressure, or a mixture's critical point, is
  wanted.

## Inputs

- A `Mixture` and a liquid (`x`) or vapour (`y`) composition; or `Tc`, `Pc`, `omega`
  and `T` for a pure saturation pressure.

## Outputs

- `pressure` (bubble/dew), `p_sat` (pure), `tc`/`pc`/`vc` (critical), plus the
  incipient composition and K-values.

## How a calculation runs

These are *models*: each fixes a procedure — an outer iteration holding one phase and
solving for the pressure at which the other first appears. A mixture above its
critical condition has neither a bubble nor a dew point, and that is reported as a
real state, not a solver failure.

## Python usage pattern

```python
import azoth

methane = azoth.eos.component("methane")
r = azoth.eos.pure_saturation(
    Tc=methane.Tc, Pc=methane.Pc, omega=methane.omega, T=azoth.ureg.Quantity(120.0, "K")
)
r.p_sat  # the saturation pressure, below Tc
```

## The keycard

A pure saturation takes `Tc`, `Pc`, `omega` directly — the caller's values, or the
databank's via `component(name)`. A keycard overrides them per parameter.

## Validation checklist

- [ ] `T` is below `Tc` for a saturation pressure; above it there is none.
- [ ] A bubble takes the liquid composition `x`; a dew takes the vapour `y`.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `OutOfRangeError` | `T >= Tc` for a pure saturation | Saturation needs `T < Tc` |
| The wrong boundary | `x` passed where `y` belongs | Bubble uses `x`, dew uses `y` |
| A critical point read as a split | The model does not decide the feed | Read the spec's assumptions |

## Limitations

A bubble/dew model takes the held composition as given — whether that phase is stable
is not asked. A pure component's bubble point *is* its saturation pressure, so
`bubble_pressure` refuses a one-component mixture.

## Related Azoth functionality

`eos.bubble_pressure`, `eos.dew_pressure`, `eos.pure_saturation`, `eos.critical_point`
— model pages under `docs/src/eos/`.

## References

- `docs/src/eos/critical_point.md` — the Heidemann-Khalil construction and why the
  naive `dP/dV = d2P/dV2 = 0` route is wrong.
