# Hydrate screening

Educational hydrate formation screening placeholder with public assumptions. USE WHEN: a task needs a quick, public screen of whether a wet stream forms a hydrate, in which structure, at what temperature or pressure, and how much of it forms once inside the region, and should be directed to validated NeqSim methods for real calculations.

This is an `azoth` skill. All three questions are ported: the equilibrium temperature
(`eos.hydrate_formation_temperature`), its inverse (`eos.hydrate_formation_pressure`) and the
amount (`eos.hydrate_fraction`).

## When to Use

- When a wet stream needs a first pass: does it form a hydrate at all at this condition.
- When the operating pressure is the free variable and the equilibrium temperature is not.
- When the question is how much water leaves the fluid once the region is entered.

## Inputs

- `components` — the substances by name, with water among them.
- `P` and `z` for the temperature; `T` and `z` for the pressure; `T`, `P` and `z` for the
  amount.
- `eos` — the cubic, `srk` or `pr`.

## Outputs

- `eos.hydrate_formation_temperature` — `temperature`, `structure`.
- `eos.hydrate_formation_pressure` — `pressure`, `structure`.
- `eos.hydrate_fraction` — `beta`, `structure`, and `balance_error`.

## How a calculation runs

All three sit on one equality: the water's fugacity in the fluid equals the hydrate's, with
the structure chosen by comparing the two structures' exponents rather than their coefficients
so that the fluid's own water fugacity cancels. The temperature and the pressure are the same
solve on different axes; `eos.hydrate_fraction` adds the amount, taking the cage composition
from **both** cavities and asserting the material balance at the answer — NeqSim's own
distributes the guests by the small cage alone and leaves a phase whose mole fractions sum to
`1.1201`, so the two do not agree and this is the one that closes.

**`balance_error` is the check to read.** It is the feed's recovery summed over the phases
minus one, and a fraction reported beside a large one is a fraction of a state that does not
exist.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
components = ["methane", "ethane", "propane", "water"]
z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]

print(azoth.eos.hydrate_formation_pressure(components, T=q(285.0, "K"), z=z, eos="srk").pressure)
formed = azoth.eos.hydrate_fraction(components, T=q(283.15, "K"), P=q(100.0, "bar"), z=z, eos="srk")
print(formed.beta, formed.structure, formed.balance_error)
```

## The keycard

The guest constants and the cavity data are the databank's own columns, resolved from the
names. A keycard overrides a component's own constants, not the cavity geometry.

## Validation checklist

- [ ] Water is one of the `components`; a dry stream has no hydrate and is refused.
- [ ] `balance_error` is small — read it before the fraction.
- [ ] The structure is reported, not assumed: I and II have different curves.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `InvalidInputError: no water` | The feed is dry | A hydrate is water's; there is no answer |
| A fraction read without `balance_error` | The state may not close | Read both |
| The wrong curve | Assuming structure I | Read `structure` |

## Limitations

An inhibitor, a salt correction and a kinetic induction time are not modelled. `beta` is
thermodynamic: it says how much hydrate is stable, not how much has formed.

## Related Azoth functionality

`eos.hydrate_formation_temperature`, `eos.hydrate_formation_pressure`,
`eos.hydrate_fraction` — model pages under `docs/src/eos/`.

## References

- `docs/src/eos/hydrate_fraction.md`, which records the divergence from NeqSim's composition
  and the measurement behind it.
- NeqSim: https://github.com/equinor/neqsim — `HydrateFormationPressureFlash`,
  `HydrateFormationTemperatureFlash`, `TPHydrateFlash`.
