# Pump hydraulics screening

Hydraulic and shaft power for a centrifugal pump from the head it develops, the flow it
moves and its efficiency; the affinity laws and the NPSH margin are built from this
against a vendor curve.

## When to Use

- When a task needs a screening shaft power for a pump, or the hydraulic power a duty
  requires before a vendor selection.

## Inputs

- Density, volumetric flow, head in metres of the pumped fluid, and the overall
  efficiency: hydraulic power delivered over shaft power supplied.

## Outputs

- `power`, the shaft power the machine must be supplied with.

## How a calculation runs

`pump_power` is the whole calculation: `rho * g * q * H / eta`. The affinity laws are
ratios of this result at two speeds or diameters, and the NPSH margin is
`(p_suction - p_vapour) / (rho * g)` plus the elevation and line losses — both are
arithmetic over the answer and over a state, not calculations of their own. The vapour
pressure is `eos.pure_saturation`; the required NPSH is the vendor's curve, which
nothing here supplies.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
r = azoth.hydraulics.pump_power(q(998.0, "kg/m**3"), q(0.01, "m**3/s"), q(30.0, "m"), 0.75)
r.power  # 3914.81468 W
```

## The keycard

Not applicable: this is water-density arithmetic over four numbers a caller supplies.
If the pumped fluid's density comes from a fluid rather than a measurement, it is the
density that needs a keycard, not this calculation.

## Validation checklist

- [ ] `eta` is in `(0, 1]`: zero is singular and above one would beat the reversible
      machine, and both raise.
- [ ] `H` is head in metres of the *pumped fluid*, not of water. A head quoted in
      another fluid's metres is a factor with no symptom.
- [ ] The vapour pressure in an NPSH margin is at the *pump suction* temperature.
- [ ] Treat the result as a screening estimate: a real pump's efficiency is a point on
      a curve, and this takes one number for it.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| Power too high or low | A density the fluid does not have | Take `rho` from the state, not from water |
| A singular error | `eta = 0` | Supply a real efficiency |
| A negative head | The sign convention reversed | A negative head is a turbine, which is a different calculation |

## Limitations

The efficiency is a single number rather than a curve, so the result is a screening
estimate and not a vendor prediction. NPSH-available is assembled by the caller;
NPSH-required is the vendor's. Choked flow and relief sizing are not this calculation.

## Related Azoth functionality

`hydraulics.pump_power` (`docs/src/hydraulics/pump_power.md`), and `eos.pure_saturation`
for the vapour pressure an NPSH margin needs.

## References

- `docs/src/hydraulics/pump_power.md`.
