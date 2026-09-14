# The keycard

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

**A keycard is how you add to azoth, or override what it ships, without writing Rust
or Python.** It is one YAML file. `keycard.example.yaml` at the repository root is the
template; copy it, fill in what you are entitled to use, and:

```python
import azoth

azoth.keycard.load("keycard.yaml")
azoth.eos.component("methane")      # now your values, not the databank's
```

Check it first — the same file, before it reaches a calculation:

```bash
python tools/check_user_data.py keycard.yaml
```

## The sections

| Section | Overrides | Read by |
|---|---|---|
| `keyholder` | — | nothing in azoth; it records who is asserting the right to use the values |
| `components` | `data/components/components.csv` | `eos.component`, `eos.from_names` |
| `kij` | `data/components/kij.csv` | every mixture's mixing rule |
| `coefficients` | — | a calculation's named argument, as a default |
| `models` | — | `eos.from_model` |
| `fittings` | `data/fittings/crane_k_factors.csv` | `tools/gen_user_data.py` — then a rebuild |
| `fluids` | `data/fluids/` | `tools/gen_user_data.py` — then a rebuild |

The last two are different in kind from the rest, and the difference is worth knowing
before you write a keycard. `keycard.yaml` is a **source**, not a runtime input: the four
sections above them are read when a keycard is loaded, and `fittings` and `fluids` are
compiled into the shipped data files by `tools/gen_user_data.py` instead. That is what
lets the Rust core embed those files with `include_str!` while Python reads the same
bytes from disk, which is the property `python/tests/test_data_agreement.py` rests on.
So a fitting or a fluid takes effect after a `gen_user_data.py` run and a rebuild — and
neither is visible to `azoth pipe`, which is a Rust binary and never sees a
Python-loaded keycard.

Everything is by name. A component you supply with a name the databank already has
replaces the parameters you list and keeps the rest — correcting one value does not
mean restating the others, and does not silently lose them. A name the databank does
not have adds a substance, and then all three of `Tc`, `Pc` and `omega` are required:
a partial component is refused rather than completed from something similar, which
would be inventing data.

## It is data, and it cannot be anything else

Every section is named choices from closed vocabularies plus numbers. There is nothing
in a keycard that this library has not already implemented, and nothing a keycard can
make it do.

That is enforced where the keycard is **loaded**, not where it is used. A model
declaring a cubic this build does not run is refused by name, rather than being stored
and silently evaluated as Peng-Robinson — which would be a wrong answer with no
symptom, and therefore the worst kind.

The vocabularies are narrower than they will be. `shape` admits one cubic and `alpha`
one alpha function, because one is what both implementations run and a test holds the
vocabulary to them together. Widening either is a change that happens *alongside* the
implementation, never before it.

## Coefficients are a default, and an explicit argument always wins

```python
azoth.hydraulics.orifice_flow(..., Cd=0.61)   # this call uses 0.61; the keycard is not read
azoth.hydraulics.orifice_flow(...)            # this call uses the keycard's value
```

The direction matters. A caller who writes a value in their own call is stating it for
that call, and overriding it from a file they loaded an hour ago would be the worst
kind of surprise.

Most coefficients need no keycard at all — pass them. The section exists for the case
where the same value is used everywhere and belongs in one place. See
[Copyright and licensed data](./copyright.md) for why so many of them are arguments in
the first place.

**A coefficient must declare a unit**, and it is checked against the spec's declared
unit for that input. A value in the wrong unit is a factor with no symptom: entirely
ordinary-looking, and wrong.

## One keycard at a time

`load()` sets *the* keycard, process-wide. There is one rather than a stack, because
the alternative is a library whose answers depend on call order — the same calculation
returning two results in one process because someone loaded a keycard between them.
`azoth.keycard.clear()` goes back to what ships.

Two datasets are two processes, which is also the only way to compare them honestly.

## What a keycard is not

**It is not a place to record provenance.** A row may carry a `citation` and nothing
requires one: no status field, no check, no rule that two fields agree. That is
deliberate rather than an omission — nobody can check whether a person read a standard,
so a required field would be a form to fill in rather than a fact, and a form teaches
people to fill it in.

What the library owes you instead is disclosure: [What ships](./data.md) says what the
library carries and where it came from. The rest — whether a value is right, and whether
you may use it — is engineering judgement, and it is yours.

**And that is the point of the file, not a gap in it.** The card is the capability
declaration: it records what you are entitled to compute with, which is why the
accountability sits with you rather than with the library. [S6](./spec.md#s6-the-library-implements-the-engineer-decides)
is where that is stated — *the library implements; the engineer decides*. A row is right
or wrong in the databank or in your card, and no check in this library can make it
otherwise. If you are adding to azoth, that section is also the reason not to build one.

**It is not shared.** A keycard is per-process and per-caller. Passing one around
between threads that expect different data is not supported, and the answer is one
keycard per process.
