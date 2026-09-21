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

`calc` names an id from **either** registry - a calculation or a model. Inputs are
named exactly as the spec's `inputs` block and are in the units that block
declares, so a validation case needs no unit handling of its own - the runner
converts using the spec.

A model's case is worth having for a different reason than a calc's. A calc's spec
asserts an *equation*, and a published worked example can check that the equation
was transcribed correctly. A model's spec asserts a *procedure* - a scheme, a
tolerance, an iteration cap - none of which any source states, so the only check on
it that does not come from the spec itself is an external one on where its answer
lands.

Two things follow. `expected` may be a number, a vector (a composition, or one
K-value per component), an enum member as its spec spelling, or `null` for an
output that is genuinely absent - `eos.pt_flash` reports no vapour fraction when a
feed has no two-phase solution, and a case can assert that. And an input that is an
*object* rather than a number needs an entry in the runner's `ARGUMENT_BUILDERS`,
which today holds one: `eos.pt_flash` takes a `Mixture`, so its case names its
components and the runner resolves them through the databank.

## NeqSim is the ground truth, and `neqsim/` is how to reach it

This library is a port of NeqSim. Where a published source and NeqSim
disagree about what a calculation does, NeqSim is what azoth is trying to be, so a
disagreement with NeqSim is the more interesting failure.

`validation/neqsim/` holds a small Java driver that builds a system in NeqSim and
prints a flash's answer in the fields azoth reports. It is not built by CI - it needs
a JVM and NeqSim's jar - so the numbers it prints are recorded in the cases beside
it, with the fluid and the state written down so anyone can run it again. `FlashTp`
carries the two commands at the top.

The jar is not vendored. `databank/sources/neqsim/` carries NeqSim's *data files*,
because those are what the databank is compiled from and they are small; the jar is
20 MB of compiled Java that nothing here links against.

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
