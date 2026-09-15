# Flashes and stability

A flash is the deliverable of the EOS namespace: how much of a feed is vapour, what
the two phases are made of, and the fugacity coefficients that make it so. This skill
is how an agent runs a temperature-pressure flash and its companions.

## When to Use

- When a feed's phase split at a temperature and pressure is wanted.
- When an enthalpy (`ph_flash`) or entropy (`ps_flash`) specification drives the split.
- When the question is whether a feed is stable as a single phase (`stability_test`).

## Inputs

- A `Mixture` (via `azoth.eos.from_names`), a temperature/pressure, and a composition
  `z` that sums to one.

## Outputs

- `beta` (the vapour fraction), `x` and `y` (the two compositions), `phase`, and the
  `ln_phi` columns.

## How a calculation runs

A flash is a *model*, not a calculation: its spec fixes a procedure (successive
substitution) rather than an equation. Successive substitution finds *a* stationary
point; whether the feed was stable is a different question, answered by
`stability_test` (Michelsen's tangent-plane criterion).

## Python usage pattern

```python
import azoth

fluid = azoth.eos.from_names(["methane", "n-butane"])
r = azoth.eos.pt_flash(
    fluid, T=azoth.ureg.Quantity(330.0, "K"), P=azoth.ureg.Quantity(2.5e6, "Pa"), z=[0.6, 0.4]
)
r.beta, r.phase  # 0.842…, two_phase
```

## The keycard

The mixture's `kij` comes from the databank (or a keycard's `kij` section, which wins
per pair). A keycard's `models` section names the cubic variant the flash runs on.

## Validation checklist

- [ ] `z` sums to one; it is checked rather than renormalised.
- [ ] Read `phase`, not `beta`: `beta` is `None` where there is no vapour fraction.
- [ ] A flash is not a stability proof — run `stability_test` for that.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `InvalidInputError` | `z` does not sum to one | It is refused, not silently corrected |
| A wrong split | Reading `beta` when the feed is single-phase | Read `phase` |
| A false single-phase verdict | The flash converged to `x = y = z` | Run `stability_test` |

## Limitations

There is no stability test inside the flash: a converged split is a stationary point,
not a proof the feed was unstable. A negative flash (`beta` outside `[0, 1]`) is a
real reading, reported with `OUT_OF_VALID_RANGE`.

## Related Azoth functionality

`eos.pt_flash`, `eos.ph_flash`, `eos.ps_flash`, `eos.stability_test` — model pages
under `docs/src/eos/`.

## References

- `docs/src/eos/pt_flash.md` and `docs/src/eos/stability_test.md`.
