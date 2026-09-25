# Architecture

azoth is a chemical-engineering calculation and data library: the thermodynamics,
hydraulics and heat transfer of a process, in Python for the ecosystem and Rust for the
engine. It is a port of Equinor's [NeqSim](https://github.com/equinor/neqsim), and what it
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
        ├─► tools/provenance.py     ─► provenance.json
        └─► tools/gen_model_inputs.py ─► crates/azoth-process/src/model_inputs_gen.rs

specs/unit_ops/**/*.toml ─► tools/gen_palette.py ─► crates/azoth-process/src/palette_gen.rs

databank/sources/neqsim/*.csv ─► tools/gen_databank.py ─► data/components/*.csv
keycard.toml                  ─► tools/gen_user_data.py ─► data/fittings/*.csv
                                                           data/fluids/*.csv
```

**Two of those read a specification nothing else compiles from.** The process palette is
declared in `specs/unit_ops/`, one file per entry, and it is read *at run time* by
`load_palette` rather than generated — because the schema is the `UnitOpSpec` type and the checker
and a front-end have to read the same declaration. `gen_palette.py` therefore embeds those files as
text (a browser has no filesystem) and `gen_model_inputs.py` compiles the one thing the model
generator drops: the inputs' names, kinds and optionals, which a form needs and a bound does not
carry.

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

## Two implementations, mirrored

Every calculation is written twice by hand, from one spec: a Rust kernel and a Python
kernel, implementing the same declared arithmetic and never generated from each other.
Neither is a wrapper over the other and neither is authoritative — they are **mirrors**, and
the library is the pair. The spec names both, in its `implementations` block, and every
calculation in `specs/` names both.

**What differs between them is what each is good for, not which one computes.** Rust is safe
at compile time, expresses the process calculus natively in its type system, and composes
into very large — or many concurrent — simulations that run in minutes where Python would
take days; it is also the half that runs without Python at all. Python is the
ecosystem-facing half: basic calculations, notebooks, the data-science and AI-SDK tooling
around them, and the agentic layer. PyO3 binds the two into one library, and units cross the
public API as `pint` quantities on one side and `uom` quantities on the other, both
extracting the SI base magnitude before doing arithmetic — so the two run the same operations
on the same numbers.

Because both exist, every case in every spec runs through both, and the comparison is
**by tolerance and not bit-equality**: `log10` and `sqrt` are not correctly-rounded in
general and `libm` differs between platforms. Where a procedure iterates, the iteration
counts are required to match. Bit-equality is asserted separately and only for one
platform and build. The data files are embedded on both sides with `include_str!`, and a
test compares the embedded bytes against the file on disk rather than comparing parsed
values, because two files can parse into the same numbers.

**That comparison is a test, not the reason there are two.** Each half is wanted for what it
is good at — Python for the ecosystem the library has to live in, Rust for compile-time
safety and for composing large or concurrent work — so writing the arithmetic twice is what
that costs, and the comparison is what the cost buys.

What the comparison proves is bounded, and the bound is worth stating because it is easy to
read past. Two kernels written independently against one declaration can drift apart, which
is what the test catches. They were still both written from the same declaration, so their
agreement says nothing about whether the declaration is right. Whether the spec was right is
what the validation cases ask, and `validation/README.md` says so in its own words.

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
