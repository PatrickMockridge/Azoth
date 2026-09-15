# Run a calculation with azoth

azoth computes equations of state, flashes, hydraulics and heat transfer, each
calculation implemented twice — a Python reference and a Rust extension bound by
PyO3 — and held to the same answer. This skill is how an agent calls it: the import,
the units, the result contract, and the keycard.

## When to Use

- When a task needs a number azoth computes — a pressure drop, a flash split, a
  saturation pressure, a compressibility factor.
- When an agent must read a result's warnings before using it.
- When a caller needs to override or extend the shipped data with a keycard.

Do not use this skill to add a new calculation — that is the port discipline, not a
call.

## Inputs

- Domain functions under `azoth.hydraulics.*`, `azoth.eos.*` and `azoth.thermal.*`.
- Quantities built with `azoth.ureg.Quantity`; plain floats only
  for genuinely dimensionless inputs such as an acentric factor, a Reynolds number
  or a friction factor.
- A keycard, passed as `card=...` where one is wanted.

## Outputs

- A result object with typed fields (`r.dp`, `r.beta`, `r.phase`, ...).
- `r.warnings` — every range check that fired or was skipped.
- `r.is_clean` — `False` the moment any warning is present.

## How a calculation runs

Two kernels run the same arithmetic and are compared by tolerance. Which one answers
is chosen at import: `azoth.backends()` lists what is available,
`azoth.active_backend()` what is in use, and `azoth.use_backend(...)` selects. Units
cross the API as `pint` quantities, so a bare number where a length is expected
raises `UnitMismatchError` rather than being silently misread.

## Python usage pattern

```python
import azoth

q = azoth.ureg.Quantity

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
r.dp  # 22455.0 pascal
r.warnings  # (Warning(RANGE_CHECK_SKIPPED, ...),)
r.is_clean  # False

fluid = azoth.eos.from_names(["methane", "n-butane"])
flash = azoth.eos.pt_flash(fluid, T=q(330.0, "K"), P=q(2.5e6, "Pa"), z=[0.6, 0.4])
flash.beta, flash.phase
```

## The keycard

`azoth.keycard.load(path)` reads a TOML file and returns a value; it stores nothing
and there is no process-wide card. A call that should read a card is handed one, and
a call handed none reads what the library ships:

```python
card = azoth.keycard.load("keycard.example.toml")
methane = azoth.eos.component("methane", card=card)
```

An explicit argument always beats the card. A name or unit the card's vocabulary does
not know is refused when the card is loaded, not when it is finally used.

## Validation checklist

- [ ] Check `r.warnings`, not just the answer: a value out of its validated range is
      still returned, with a warning.
- [ ] Pass a quantity where the spec declares a unit; a bare float is a
      `UnitMismatchError`.
- [ ] Read `r.is_clean` when a result must carry no caveat at all.
- [ ] Pass the keycard explicitly; do not assume a card loaded elsewhere is in force.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `UnitMismatchError` | A bare float where a quantity is expected | Build `q(value, "unit")` |
| A plausible wrong answer | A quantity in the wrong unit | Check the unit against the calc page |
| A result used as validated | `r.warnings` never read | Read warnings, or `r.is_clean` |
| Two answers for one call | A card assumed global | Pass `card=` explicitly |
| A pressure drop that looks fine | The shipped `crane_k_factors.csv` is placeholder data | Supply your own fittings via a keycard |

## Limitations

azoth is early: the implemented set is the table in `README.md`. The shipped fitting
coefficients under `data/fittings/crane_k_factors.csv` are estimated placeholders,
not engineering data — a pressure drop from them can be wrong by a factor of two and
look reasonable. See "Not for design work yet" in `README.md`.

## Related Azoth functionality

The ids below are namespaced by domain and are what a caller reaches for; each has a
generated page under `docs/src/`:

- `hydraulics.darcy_weisbach`, `hydraulics.reynolds_number` — pipe flow.
- `eos.pr_kappa`, `eos.pr_alpha_ab`, `eos.pr_z_factor` — the Peng-Robinson cubic.
- `eos.pt_flash`, `eos.bubble_pressure`, `eos.dew_pressure` — phase equilibria.
- `thermal.conduction_plane_wall` — steady conduction.

The CLI `azoth pipe` composes several hydraulics calcs for one command.

## References

- `README.md` — quick start and the implemented list.
- `docs/src/architecture/specification.md` — the normative statement.
- `keycard.example.toml` — the keycard template.
- `NOTICE` — what ships and where it came from.
