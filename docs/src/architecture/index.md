# Architecture

azoth is a chemical-engineering calculation and data library: the thermodynamics,
hydraulics and heat transfer of a process, implemented twice and held to each other.
It is a port of Equinor's [NeqSim](https://github.com/equinor/neqsim), and what it
carries is what NeqSim carries, named as NeqSim names it.

Everything on this page is stated once, here. A later page that repeats it is wrong,
and a later page that contradicts it is a defect.

## One document, many readings

```
specs/**/*.toml          one file per calculation, model, case, and the unit vocabulary
specs/schema/*.json      the contract each of those is held to
        │
        ├─► tools/gen_registry.py   ─► crates/azoth-*/src/spec_gen.rs
        │                              python/src/azoth/_registry_gen.py
        ├─► tools/gen_models.py     ─► crates/azoth-eos/src/model_gen.rs
        │                              python/src/azoth/_models_gen.py
        ├─► tools/gen_vocabulary.py ─► crates/azoth-core/src/unit_vocab_gen.rs
        │                              python/src/azoth/core/_units_gen.py
        │                              specs/schema/unit.schema.json
        │                              lean/Azoth/Vocabulary.lean, Gate.lean
        ├─► tools/gen_docs.py       ─► docs/src/**, including docs/src/SUMMARY.md
        ├─► tools/gen_stub.py       ─► python/src/azoth/_core.pyi
        └─► tools/provenance.py     ─► provenance.json

databank/sources/neqsim/*.csv ─► tools/gen_databank.py ─► data/components/*.csv
keycard.toml                  ─► tools/gen_user_data.py ─► data/fittings/*.csv
                                                           data/fluids/*.csv
```

**Every generated file on the right comes from a file on the left, and CI regenerates
all of them and fails on a diff.** So a specification is not a document that is
supposed to match the code; it is the thing the code was made from. `tools/` also
holds the checkers - `spec_lint`, `prose_lint`, `check_links`, `check_user_data`,
`check_manifest`, `check_lean_axioms` - which generate nothing and are what fail the
build.

## What is true, and where it is written

| what is true | where it is |
|---|---|
| the format of a calculation, a model and a case | `specs/schema/{calc,model,case}.schema.json` |
| which units a spec may name, and their dimensions | `specs/vocabulary/vocabulary.toml` |
| the dimension group, and when two dimensions are equal | `lean/Azoth/Dim.lean` |
| each calculation's equation, domain and cases | `specs/calcs/**/*.toml` |
| what a range check does at run time | `crates/azoth-core/src/spec.rs` |
| what is taken from NeqSim, column by column | `databank/manifest.toml` |
| the licence obligations for the vendored data | `NOTICE`, `LICENSE-CC-BY-4.0` |
| what the library promises and does not | [the specification](./specification.md) |

`tools/check_lean_axioms.py` reads `lean/Azoth/Axioms.lean` and refuses any axiom the
development was not allowed to use, so the Lean half is checked rather than trusted.

## Two implementations, one answer

Python is the reference implementation; Rust is the second, and the two are bound by
PyO3. Units cross the public API as `pint` quantities on one side and `uom` quantities
on the other, and both extract the SI base magnitude before doing arithmetic - so the
two run the same operations on the same numbers rather than each trusting a units
library to arrive at them by a different route.

Every case in every spec runs through both, and the comparison is **by tolerance and
not bit-equality**: `log10` and `sqrt` are not correctly-rounded in general and `libm`
differs between platforms. Where a procedure iterates, the iteration counts are
required to match. Bit-equality is asserted separately and only for one platform and
build.

Rust is not here for speed. For a single call the PyO3 boundary costs more than the
arithmetic saves. Its value is being an independent implementation of the same
document - which is also why the data files are embedded on both sides with
`include_str!` and a test compares the embedded bytes against the file on disk rather
than comparing parsed values: two files can parse into the same numbers.

## The port

azoth is a port of NeqSim, and the rules a port follows - what it names, what it may
claim, and what accepts it - are [Part one of the specification](./specification.md).
The record of what was taken from which upstream file, column by column, is
`databank/manifest.toml`; `tools/check_manifest.py` refuses a `not-ported` row that does
not name the NeqSim class that would close it, so the file is a backlog rather than a
dumping ground.

## The rest of this section

- [The specification](./specification.md) is the normative statement: what azoth is
  for, what it owes a reader, and the policy decisions behind the shape of it.
- [Spec files](./spec-files.md) is the format of everything under `specs/`.
- [The calculus of thermodynamic dimensionality](../calculus/index.md) is the formal
  layer: the types, the vocabulary, and what a keycard is as a capability.
