# Units and the vocabulary

Units cross azoth's API as `pint` quantities. A bare number where a length is
expected is a `UnitMismatchError`, not a silent thousand-fold error. The canonical
unit for every spec input is compiled from one vocabulary table, so a unit added
anywhere else is a unit that disagrees.

## When to Use

- When a caller passes quantities to a calculation and must get the unit right.
- When a spec input declares a unit and a caller wonders what it may pass.
- When a result carries a unit and must be converted before use.

## Inputs

- Quantities built with `azoth.ureg.Quantity`.
- Plain floats only for genuinely dimensionless inputs — an acentric factor, a
  Reynolds number, a friction factor.

## Outputs

- The result's fields as `pint` quantities in their declared units.

## How a calculation runs

Both kernels extract the SI base magnitude before doing arithmetic, so the two run
the same operations on the same numbers rather than each trusting a units library to
arrive at them by a different route. `azoth.ureg` is the pint registry, and
`azoth.ureg.Quantity` builds a quantity. `azoth.Q` is a type alias for
annotations, not a constructor.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity
q(1000.0, "Pa").to("kPa")  # 1.0 kilopascal
azoth.ureg.Quantity(1.5, "m/s")  # the same constructor, shorter

# A bare float where a length is expected is refused, not misread:
# azoth.hydraulics.darcy_weisbach(0.02, 100.0, ...)  -> UnitMismatchError
```

## The keycard

A keycard coefficient must declare a unit, checked against the spec's declared unit
for that input. A value in the wrong unit is a factor with no symptom — the unit
declaration is what catches it.

## Validation checklist

- [ ] A dimensioned input is a quantity; a dimensionless one is a plain float.
- [ ] The unit is the one the calc page declares — not merely "a pressure".
- [ ] A conversion is explicit, not assumed: `q(...).to("kPa")`.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `UnitMismatchError` | A bare float where a quantity is expected | Build `q(value, "unit")` |
| A plausible wrong answer | The wrong unit accepted and passed through | Check the unit against the spec |
| A factor of a thousand between kernels | A conversion written by hand | The vocabulary is the one table; no second copy |

## Limitations

The vocabulary is closed: `specs/vocabulary/vocabulary.toml` is the one hand-written
source, and the schema's unit enum, the Rust table and the Python map are compiled
from it. A unit outside it is refused at load, not guessed at.

## Related Azoth functionality

- `specs/vocabulary/vocabulary.toml` — the 24-unit vocabulary.
- `lean/Azoth/Vocabulary.lean`, `lean/Azoth/Dim.lean` — the proved unit theorems.
- `docs/src/calculus/dimensions.md` and `docs/src/calculus/vocabulary.md`.

## References

- `docs/src/architecture/specification.md`, part three — the vocabulary rule.
