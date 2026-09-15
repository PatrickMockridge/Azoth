# The cubic equation of state

A cubic EOS decomposes into constitutive coefficients — the alpha-function
coefficient, the reduced attraction parameters, the compressibility factor, the
fugacity coefficient — each with its own source and worked example. This skill is how
an agent drives that chain for Peng-Robinson, Soave-Redlich-Kwong and Redlich-Kwong.

## When to Use

- When a task needs a compressibility factor, fugacity coefficient, molar volume or
  mass density from a cubic.
- When a cubic's alpha correlation, Peneloux volume shift or binary mixing must be
  applied.

## Inputs

- Acentric factor (`omega`, a plain float) and reduced temperature/pressure for the
  alpha function; reduced `a` and `b` for the factor and departure.

## Outputs

- `kappa`, `a_reduced`/`b_reduced`, `z`, `ln_phi`, `v`, `rho` — each a typed field on
  a result object.

## How a calculation runs

Everything is in reduced variables, so the calc is dimensionless end to end and
carries no unit. The coefficients come first, deliberately: a transposed digit in the
`omega**2` coefficient produces an EOS that still runs and is slightly wrong
everywhere, and no check on `z` can see it, because `z` is computed *from* the
coefficient.

## Python usage pattern

```python
import azoth

kappa = azoth.eos.pr_kappa(0.0115)  # methane's acentric factor
ab = azoth.eos.pr_alpha_ab(kappa.kappa, Tr=2.0, Pr=1.0)
z = azoth.eos.pr_z_factor(ab.a_reduced, ab.b_reduced)
```

## The keycard

A keycard can override `omega` per component, and name a model variant (`cubic_eos` /
`peng_robinson` / `peng_robinson` alpha / `classical_kij` mixing). A variant the build
does not implement is refused when the card is loaded.

## Validation checklist

- [ ] `omega` is a plain float, not a quantity — the acentric factor is dimensionless.
- [ ] The alpha correlation matches the cubic; each has its own coefficient calc.
- [ ] `z` is read as the admissible roots, not any real root of the polynomial.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A slightly-wrong answer everywhere | A coefficient digit transposed | Compare the coefficient calc's worked example |
| `UnitMismatchError` | `omega` passed as a quantity | Pass a plain float |
| The middle root used | It is a root, not a state | Use the smallest/largest admissible roots |

## Limitations

The middle root of the cubic is deliberately not returned: it is a root of the
polynomial and not a state the EOS describes. A `kappa` below zero (helium) returns
with `OUT_OF_VALID_RANGE`, not an error — a value out of range is still a value.

## Related Azoth functionality

`eos.pr_kappa`, `eos.srk_kappa`, `eos.rk_alpha_ab`, `eos.pr78_kappa`, `eos.twu_kappa`,
`eos.prsv_kappa`, `eos.pr_alpha_ab`, `eos.srk_alpha_ab`, `eos.pr_z_factor`,
`eos.srk_z_factor`, `eos.pr_departure`, `eos.srk_departure`, `eos.rk_departure`,
`eos.pr_molar_volume`, `eos.pr_mass_density`, `eos.pr_peneloux_shift`,
`eos.srk_peneloux_shift`, `eos.vdw1f_mix_binary`, `eos.rachford_rice_binary` — pages
under `docs/src/eos/`.

## References

- `docs/src/architecture/specification.md`, part one — the cubic family (P2).
