# Validation cases

A validation case is an **external** check on a calculation: a worked example
from a standard, a textbook, a published paper, or an independent calculation,
recorded as data and run by CI.

It is deliberately separate from the tests in a spec. Those are generated from
the spec itself and check that the code does what the spec says. A validation
case checks whether the *spec was right* - it comes from outside the registry, so
it can disagree with it.

## The format

```json
{
  "id": "crane_tp410_example_3_5",
  "calc": "hydraulics.darcy_weisbach",
  "source": {
    "attribution": "Crane TP-410 Example 3-5",
    "verification": "source_needed",
    "notes": "Why this is or is not confirmed."
  },
  "derivation": "How the expected value was obtained.",
  "inputs": { "f": 0.02, "L": 100.0, "D": 0.1, "rho": 998.0, "v": 1.5 },
  "expected": { "dp": 22455.0 },
  "tolerance": 1.0e-9
}
```

`calc` must be a calc id from the registry. Inputs are named exactly as the
spec's `inputs` block and are in the units that block declares, so a validation
case needs no unit handling of its own - the runner converts using the spec.

`source.verification` is one of:

| Value | Meaning |
|---|---|
| `verified` | A person checked the attribution against the primary source |
| `unverified` | The attribution is believed but unconfirmed |
| `source_needed` | The attribution could not be established, or the source could not be consulted |

A `source_needed` case may still be run - the *arithmetic* is ours and can be
checked - but the runner reports it, and the case never claims to validate
against a source nobody has read.

## Running them

```bash
pytest python/tests/validation
```

Every case is also run in CI. A mismatch fails the build.

## What a passing case does and does not mean

A passing case means the implementation agrees with the recorded numbers. It does
not mean those numbers are right: that depends entirely on `source.verification`,
and on whether the case's inputs are inside the range the calc is validated over.

Read the case's `notes` before treating it as evidence of anything.

## Adding one

1. Put the numbers in `<source>/<case>.json` under `validation/`.
2. Say where they came from, honestly. If you found the example in a secondary
   reference, say so rather than naming the standard it claims to derive from.
3. Derive the expected values yourself where you can, and say how.
4. Never transcribe a table or worked example from a copyrighted standard. Cite
   the equation and work the arithmetic - the arithmetic is not copyrightable
   and is the part worth checking anyway.
