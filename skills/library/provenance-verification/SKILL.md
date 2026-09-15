# Provenance and verification

azoth's answers are worth only what their validation is worth, and validation that
cannot be read cannot be checked. This skill is how a caller or agent reads what a
result rests on: where a value came from, whether it is `verified`, `unverified` or
an `estimated_dummy`, and what a warning means.

## When to Use

- When a result must be attributed to a source before it is trusted.
- When a caller needs to know whether a value is a placeholder.
- When deciding whether a divergence from NeqSim is a finding or a failure.

## Inputs

- A result's `warnings` and `is_clean`.
- A shipped data row's `verify_status` column.
- `NOTICE` and `databank/manifest.toml`.

## Outputs

- A statement of what a value rests on, and what has and has not been confirmed.

## How a calculation runs

Three distinct things are often conflated, and azoth keeps them apart:

- **Warnings are not errors.** A value out of its validated range is still returned,
  with a warning. What the library must never do is return it silently.
- **A check that could not run is not a check that passed.** `RANGE_CHECK_SKIPPED`
  says "never checked", not "fine".
- **`verified` means the code is what it claims to be**, not that the correlation is
  right for your fluid, roughness or Reynolds number.

## Python usage pattern

```python
import azoth

r = azoth.hydraulics.darcy_weisbach(
    0.02,
    azoth.ureg.Quantity(100.0, "m"),
    azoth.ureg.Quantity(0.1, "m"),
    azoth.ureg.Quantity(998.0, "kg/m**3"),
    azoth.ureg.Quantity(1.5, "m/s"),
)
r.is_clean  # False
[w.code for w in r.warnings]  # [RANGE_CHECK_SKIPPED]
```

## The keycard

The shipped data carries a `verify_status` column because it is the library's
statement about itself. A user's keycard rows do not — the engineer who supplied the
value knows whether they trust it, and a status field would only invite them to
assert something nobody can check.

## Validation checklist

- [ ] Check `warnings` before using a result, not after.
- [ ] Read `NOTICE` for what ships and where it came from.
- [ ] Treat a divergence from NeqSim as a finding, never a silent correction.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A placeholder used as data | The `estimated_dummy` status not read | Check `verify_status` |
| A result used as validated | `warnings` ignored | Read warnings, or `is_clean` |
| Confidence manufactured | A status field nobody can check | Record provenance, not status |

## Limitations

None of this proves the correlation is right for a specific fluid or operating range.
A verified artifact means the code is what it claims to be; the sources are cited in
each spec, and the papers are the caller's to read.

## Related Azoth functionality

- `NOTICE` — what ships and where it came from.
- `databank/manifest.toml` — the column-by-column vendoring record.
- `provenance.json` — the per-build record (commit and SHA-256 of every spec,
  implementation, test and data file).

## References

- `docs/src/architecture/specification.md` — verification and the port's three rules.
- `README.md`, "Verifying a result".
