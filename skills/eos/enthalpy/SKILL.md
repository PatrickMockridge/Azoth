# Enthalpy, entropy and heat capacity

The energy side of the EOS namespace: ideal-gas heat capacity from a polynomial, the
heat of vaporisation, liquid heat capacity, and the absolute molar enthalpy and
entropy of a mixture.

## When to Use

- When a task needs an ideal-gas `cp` at a temperature, a heat of vaporisation, or a
  mixture's molar enthalpy/entropy at a state.

## Inputs

- The five `cp` polynomial coefficients (`cp_a`..`cp_e`) and a temperature; for
  enthalpy/entropy, a `Mixture`, an ideal-gas model, `T`, `P`, `z` and the
  compressibility root.

## Outputs

- `cp` (a heat capacity), `h` and `s` (molar enthalpy and entropy), each a quantity.

## How a calculation runs

The coefficients carry the powers of temperature in their units, so `cp_a` is a heat
capacity and `cp_b` is one per kelvin. `molar_enthalpy_entropy` assembles
`H = H_ig(T_ref) + ∫Cp dT + H_dep`, where the ideal-gas model carries the datum and
nothing checks it.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
cp = azoth.eos.ideal_gas_cp(
    q(30.0, "J/(mol*K)"),
    q(0.0, "J/(mol*K**2)"),
    q(0.0, "J/(mol*K**3)"),
    q(0.0, "J/(mol*K**4)"),
    q(0.0, "J/(mol*K**5)"),
    T=q(300.0, "K"),
)
cp.cp  # 30 J/(mol*K), with a constant polynomial
```

## The keycard

A keycard can supply the `cp_a`..`cp_e` coefficients for a component; the databank
carries them for every substance it ships. A substance a keycard adds needs its own,
or it has no enthalpy path.

## Validation checklist

- [ ] The ideal-gas datum is the same across calls being compared; subtracting two
      enthalpies from different reference values gives a plausible number, not an error.
- [ ] `T` lies inside the range the coefficients were fitted over — it is not checked.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A plausible wrong difference | Two enthalpies from different reference values | Use the same datum |
| A negative `cp` | The polynomial evaluated outside its fitted range | The result carries the warning |
| The wrong root | `compressibility` for the wrong phase | Use the root for the phase wanted |

## Limitations

A polynomial evaluated outside its fitted range turns over; the result carries
`OUT_OF_VALID_RANGE` — the case the check catches is `cp` going negative, not the more
likely one where the value is a few per cent wrong.

## Related Azoth functionality

`eos.ideal_gas_cp`, `eos.molar_enthalpy_entropy`, `eos.liquid_heat_capacity`,
`eos.heat_of_vaporization` — pages under `docs/src/eos/`.

## References

- `docs/src/eos/molar_enthalpy_entropy.md`.
