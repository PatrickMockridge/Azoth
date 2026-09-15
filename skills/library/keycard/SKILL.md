# Author a keycard

A keycard is one TOML file that overrides or extends what azoth ships: a component's
critical constants, a binary interaction parameter, a coefficient, a fluid table, a
fitting, or a named model variant. The card is **where responsibility sits** — the
library implements, the engineer decides.

## When to Use

- When a caller's data differs from the shipped databank and must be supplied.
- When a component, `kij`, coefficient, fitting or model variant is to be declared.
- When a value that a calculation reads must carry its own provenance, not the
  library's.

## Inputs

- A TOML document with `schema_version = 2` and the sections the library reads:
  `keyholder`, `components`, `kij`, `coefficients`, `models`, `fittings`, `fluids`.

## Outputs

- A `Keycard` value — immutable, and the whole of what the file said.

## How a calculation runs

`azoth.keycard.load(path)` reads a file and returns a value; `use(mapping)` builds one
from an already-parsed document. Neither stores anything, and there is no
process-wide card. A call that should read a card is handed one; a call handed none
reads what the library ships.

## Python usage pattern

```python
import azoth

card = azoth.keycard.use(
    {
        "schema_version": 2,
        "keyholder": {"name": "example"},
        "components": {"methane": {"omega": {"value": 0.0115, "unit": "dimensionless"}}},
    }
)

azoth.eos.component("methane", card=card)  # the card's omega, the databank's Tc and Pc
```

A name already in the databank is overridden parameter by parameter; a new name must
be complete. A parameter or unit the card's vocabulary does not know is refused when
the card is loaded, not when it is finally used.

## The keycard

The card is data and cannot be anything else. Every section is named choices from
closed vocabularies plus numbers, so nothing in a card can make the library do
something it does not already implement.

## Validation checklist

- [ ] `schema_version` is `2`; a stale version is refused rather than guessed at.
- [ ] Every coefficient declares a unit, checked against the spec's declared unit.
- [ ] A new component supplies everything a cubic reads — `Tc`, `Pc`, `omega`.
- [ ] A model variant names only the vocabularies this build implements.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A value silently ignored | A section or parameter the format does not define | Refused on load — read the error |
| A wrong answer that looks fine | A coefficient in the wrong unit | Declare the unit; it is checked |
| Two answers for one call | A card assumed global | Pass `card=` explicitly |
| A partial component accepted | Missing `Tc`/`Pc`/`omega` | Complete it or it is refused |

## Limitations

The card is **not shared and not stored**: there is no `current()` and no `clear()`.
A caller with two datasets passes a different card to each call rather than running
two processes. The shipped data carries a `verify_status` column; a user's keycard
rows deliberately do not — the engineer who supplied the value knows whether they
trust it.

## Related Azoth functionality

- `keycard.example.toml` — the template.
- `python/src/azoth/keycard.py` — the reader.
- `specs/schema/keycard.schema.json` and `tools/check_user_data.py` — the schema and
  its checker.
- `docs/src/calculus/capability.md` — the keycard as a capability.

## References

- `docs/src/architecture/specification.md`, part two — the keycard.
