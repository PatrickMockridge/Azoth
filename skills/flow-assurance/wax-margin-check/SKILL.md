# Wax margin check

Educational wax appearance margin screening with public assumptions. USE WHEN: a task needs a quick check of the temperature at which wax first appears in an oil and the margin an operating temperature keeps above it, and should be directed to validated NeqSim methods for real calculations.

This is an `azoth` skill. The wax amount is ported (`eos.tp_multiflash_wax`) and the appearance
temperature is found from it; what it does not model is the deposit — how much of the wax that
forms sticks to a wall.

## When to Use

- When an oil needs a wax appearance temperature (WAT) and a margin above it.
- When the question is how much wax is stable at a given temperature, not only whether any is.
- When a cut's own properties have to be built before either, because the databank has no row
  for a pseudo-component.

## Inputs

- `components` — the substances by name. **At least one must be a wax former**, which the
  databank states in its `waxformer` column: the n-alkanes from n-hexane upward.
- `T`, `P` and `z` for the flash; `eos`, the cubic the fluid and each cut's reference liquid
  are built on.
- For a cut: `molar_mass` and `density`, through `eos.tbp_fraction_properties`.

## Outputs

- `eos.tp_multiflash_wax` — `wax_fraction`, `phase_count`, `beta`, `x`, `converged`.
- `eos.tbp_fraction_properties` — `tc`, `pc`, `boiling_temperature`, `acentric_factor`,
  `attraction_exponent`.
- `eos.wax_solid_fugacity` — `fugacity_coefficient`, the solid's own.

## How a calculation runs

The wax phase is a solid whose coefficient is `phi_liquid exp(fusion)` over the component's
heat of fusion and triple point — a **function of the state and the component alone**, so it
is one vector for the whole solve rather than something rebuilt per trial. The amount comes
from `eos.tp_multiflash_wax`, which is `eos.tp_multiflash` with that one more phase in the set.

**The WAT is a root, not a model.** There is no appearance-temperature calculation: wax is
either stable at a temperature or it is not, so the appearance temperature is where
`wax_fraction` reaches zero and is found by bisection on the flash. The example below does it
in seventeen steps.

**`converged` is false where nothing forms**, because the three-phase set has a phase with
nowhere to go and the answer returned is the two-phase flash's — which is the same state. That
is the signal a WAT search bisects on.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
components = ["methane", "n-heptane", "nc14", "nc20"]
z = [0.7, 0.1, 0.1, 0.1]

formed = azoth.eos.tp_multiflash_wax(components, T=q(265.0, "K"), P=q(5.0, "bar"), z=z, eos="srk")
print(formed.wax_fraction)  # 0.0707

absent = azoth.eos.tp_multiflash_wax(components, T=q(275.0, "K"), P=q(5.0, "bar"), z=z, eos="srk")
print(absent.wax_fraction, absent.converged)  # 0.0, False
```

## The keycard

The wax flag and the melt data (`heatoffusion`, `triplepointtemperature`) are the databank's
own columns. A pseudo-component the databank has no row for takes its properties from
`eos.tbp_fraction_properties` and its melt data from the caller — which is what NeqSim's
`addTBPWax` does.

## Validation checklist

- [ ] At least one component is a wax former, or no wax can form and the answer is trivially
      zero.
- [ ] A cut carries a heat of fusion and a triple point; a zero is refused rather than read as
      a substance that does not melt.
- [ ] `converged` is read: it says which of the two flashes answered.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `wax_fraction` always zero | No wax former in the feed | Check the `waxformer` column |
| A WAT below the operating range | Bisecting on `wax_fraction` alone | Bisect on `converged` as well |
| A cut with a placeholder melt | `nc15` and up carry `6000` J/mol at `273.0` K | Take the melt data from the assay |

## Limitations

An inhibitor, a pour-point depressant, an ageing deposit and an energy balance are not
modelled: `wax_fraction` is how much wax is *stable*, not how much deposits. **A fluid of the
higher alkanes should precipitate far above 265 K and this one does not**, because
`nc15`–`nc24` carry a placeholder melt in the databank; that is the table rather than the
model, and a real assay's melt data is what fixes it.

## Related Azoth functionality

`eos.tp_multiflash_wax`, `eos.wax_solid_fugacity`, `eos.tbp_fraction_properties`,
`eos.tp_multiflash` — model pages under `docs/src/eos/`.

## References

- `docs/src/eos/tp_multiflash_wax.md` and `docs/src/eos/tbp_fraction_properties.md`.
- NeqSim: https://github.com/equinor/neqsim — `TPmultiflashWAX`, `ComponentWax.fugcoef2`.
