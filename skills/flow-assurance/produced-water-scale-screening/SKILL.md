# Produced water scale screening

Public produced-water brine builder and screening-level scale evaluation. USE WHEN: a task needs to turn an ion analysis, preset, or TDS value into a NeqSim-ready electrolyte ion mapping and a quick scale/mixing-incompatibility screening (BaSO4, SrSO4, CaSO4, CaCO3), and should be directed to validated NeqSim checkScalePotential methods for design-grade work.

This is an `azoth` skill. The saturation ratio and the amount precipitated are ported
(`eos.scale_saturation_ratio`, `eos.salt_precipitation`) and the activity coefficients come
from the ported Pitzer phase; **the ion mapping is the caller's** — a brine is named with the
databank's own species names, and no analysis-to-species conversion is done here.

## When to Use

- When a produced water needs a saturation ratio against a mineral's solubility product.
- When two waters are being checked for mixing incompatibility: the ratio of the mixture, not
  of either alone, is the answer.
- When the question is how much mineral comes out, not only whether any can.

## Inputs

- `components` — the brine by name. `water` must be one of them, and every ion the mineral is
  built from must be too; the databank's species names are `Na+`, `Cl-`, `Ca++`, `Ba++`,
  `SO4--`, `CO3--`, `HCO3-` and so on.
- `z` — the overall mole fractions over the brine as named.
- `salt` — the mineral, a row of `COMPSALT.csv`: `NaCl`, `BaSO4`, `CaCO3`, `FeS` and the rest.
- `T` and `P`.

## Outputs

- `eos.scale_saturation_ratio` — `saturation_ratio`, `ion_activity_product`,
  `solubility_product`.
- `eos.salt_precipitation` — `precipitated_moles`, `initial_saturation_ratio`,
  `final_saturation_ratio`, `extent_of_maximum`.

## How a calculation runs

The ratio is the ion activity product over the solubility product, `IAP/Ksp`, both taken in
the brine's own activity coefficients: the ions come from `eos.pitzer_phase`, which returns
`gamma` for every species and the water's activity. **Above one the brine is supersaturated;
at one it is at equilibrium; below one nothing can form.**

`eos.salt_precipitation` turns that into an amount by moving the ions down until the ratio
reaches one, re-solving the Pitzer coefficients at each extent — so the bracket is `[0,` the
extent at which an ion runs out`]` and `extent_of_maximum` says which end the answer was at.
One means an ion was exhausted; less means the ratio crossed one first.

**A mineral the brine cannot form is refused, not reported as a zero**, because a zero would be
an answer about a mineral that is not there.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
brine = ["water", "Na+", "Cl-", "CO3--", "HCO3-"]
z = [
    0.802568218298555,
    0.0963081861958266,
    0.0963081861958266,
    0.00321027287319422,
    0.00160513643659711,
]

pitzer = azoth.eos.pitzer_phase(brine, T=q(298.15, "K"), x=z)
ratio = azoth.eos.scale_saturation_ratio(
    "NaCl",
    x1=z[1],
    x2=z[2],
    x_water=z[0],
    gamma1=pitzer.gamma[1],
    gamma2=pitzer.gamma[2],
    water_activity=pitzer.water_activity,
    T=q(298.15, "K"),
    P=q(10.0, "bar"),
)
print(ratio.saturation_ratio)  # 1.2721: a quarter over halite's product

taken = azoth.eos.salt_precipitation(brine, salt="NaCl", T=q(298.15, "K"), P=q(10.0, "bar"), z=z)
print(taken.precipitated_moles, taken.extent_of_maximum)
```

## The keycard

The ions and their charges are the databank's rows; the Pitzer parameters — the `beta` terms,
`theta` and the temperature coefficients — are the vendored dataset's. A species pair the
dataset has no `theta` for is refused rather than run with a zero.

## Validation checklist

- [ ] The brine names every ion the mineral is built from, and water.
- [ ] `T` and `P` are the ones at the point of interest: `Ksp` is a function of both.
- [ ] A mixing-incompatibility question is asked of the *mixture's* composition.
- [ ] `extent_of_maximum` is read: one means an ion ran out.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `InvalidInputError: is built from` | The brine does not carry one of the mineral's ions | Name it, or ask about a mineral it can form |
| A ratio computed with `gamma = 1` | The ideal-solution assumption | Use `eos.pitzer_phase`'s `gamma` |
| A ratio read as an amount | A ratio is not a mass | `eos.salt_precipitation` gives the amount |

## Limitations

The multi-mineral case — several solids competing for the same ions — is not here: the
complementarity loop is written and parked, because a list of mineral names has no shape in
the spec system. One mineral at a time is what this answers. Kinetics, an inhibitor and the
scale that adheres are all outside it.

## Related Azoth functionality

`eos.scale_saturation_ratio`, `eos.salt_precipitation`, `eos.pitzer_phase` — model pages under
`docs/src/eos/`.

## References

- `docs/src/eos/scale_saturation_ratio.md` and `docs/src/eos/salt_precipitation.md`.
- NeqSim: https://github.com/equinor/neqsim — `CheckScalePotential`, `MultiSaltPrecipitation`,
  `CalcSaltSatauration`.
