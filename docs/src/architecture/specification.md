# The specification

**An opinionated port of [NeqSim](https://github.com/equinor/neqsim) to Rust, with every
calculation mirrored in Python.**

This page is normative: where a docstring, a comment or a habit disagrees with it, this
page wins and the other thing is a bug. It has three parts — the port, the keycard, and
the calculus.

---

# Part one: the port

## Until the port is complete there is no other priority

Everything below is subordinate to finishing the port. There is no feature work, no
polish, no refactor and no discussion that outranks it. When the port is complete this
section is deleted and the priority moves to whatever is next.

## What the port is

NeqSim 3.20.0 is 1,285,392 lines of Java across 3,371 files. The part being ported is its
**calculation layer**:

| Package | Files | Lines |
|---|---|---|
| `thermo/` — phase models, components, systems, mixing rules | 370 | 146,563 |
| `thermodynamicoperations/` — flashes, saturation ops, envelopes | 128 | 53,035 |
| `physicalproperties/` — transport and interfacial properties | 108 | 17,353 |
| **Total** | **~606** | **~217,000** |

Inside that: **~70 phase models** in seven families (cubic, activity, associating CPA and
SAFT, reference/Helmholtz, electrolyte, solid/hydrate, bases), **~105 flash operations**
and ~60 saturation operations, **~60 physical-property methods**, 25 alpha functions and
15 mixing rules, over a **~90-field component model**.

azoth has **43 ids** — 34 calculations and 9 models.

## The ledger

`databank/manifest.toml` records every column of every vendored NeqSim table, and for each
column that is not carried across, the NeqSim class that would close it, under a reason
prefix from a closed vocabulary: `not-ported`, `not-yet`, `not-a-value`,
`empty-upstream`, `licence`, `superseded-by`.

`tools/check_manifest.py` validates it and prints the tally. That tally is the port's
progress bar, and it is the checker's to print rather than this page's to copy: run

```bash
python tools/check_manifest.py
```

## How a calculation is ported

**The spec is a data sheet.** A file under `specs/` declares and does not explain; the
format is specified in [Spec files](./spec-files.md). It names the inputs and outputs with
their units, the valid range, the assumptions that are *not* checked, and the cases.

**Two kernels are hand-written, and only two.** One Rust file and one Python file. They
are separate implementations of the same declared arithmetic, and they are compared case
by case in the same process with the iteration counts required to match. This is the
property the port exists to preserve, and it is the one thing that must not be generated.

**Everything that names them is generated.** A new calculation currently costs eleven hand
edits across Rust and Python on top of its two kernels — result structs, transport
classes, registration lines, `__all__` entries, dispatch wrappers. Every one of those is a
function of the spec's `inputs:`, `outputs:` and the id. They are being moved into a
generator, which reduces a ported calculation to **four hand-written files**: the spec,
the Rust kernel, the Python kernel, and the test.

**A spec declares every result field, in both directions.** A field the spec does not name
cannot be emitted by a generator, and a declared output with no field is a promise the
code cannot keep. `python/tests/test_registry_contract.py` holds both directions.

## Verification

**A port names three things**: the paper for the method, the NeqSim class for the port, and
its own notes for what it changed. **A port never upgrades a verification status** —
reading someone's Java is not reading the paper, so everything ported ships `unverified`.

**A port is accepted on azoth's own tests, never on its provenance.** That NeqSim
implements something is evidence it can be implemented, not evidence it is right. Its
`CriticalPointFlash` is correct by inspection and validated nowhere.

**NeqSim runs as a differential oracle, and it does not gate the build.** The checkout at
`/tmp/neqsim-check/neqsim` is built (Java 21, `target/neqsim-3.20.0.jar`) and ships 29
paired input/output flash cases plus literature benchmarks. A divergence from NeqSim is a
finding, never a failure: the oracle is a second opinion, not a source of truth.

## The order

Each tranche closes a set of manifest columns, and a tranche is not done until its columns
are vendored and reachable.

| | Tranche | NeqSim |
|---|---|---|
| **P0** | The port engine: generate the glue, and the oracle | — |
| P1 | Property models — viscosity, conductivity, diffusivity, surface tension, liquid Cp, Antoine, heat of vaporisation | `physicalproperties/`, `interfaceproperties/` |
| P2 | The cubic family and the 25 alpha functions | `thermo/phase/`, `attractiveeosterm/` |
| P3 | Mixing rules | `thermo/mixingrule/` |
| P4 | Reference EOS — GERG-2004/2008, IAPWS-95, Leachman, Span-Wagner, EOS-CG, BWR | `thermo/phase/` |
| P5 | Activity models — NRTL, UNIFAC, UNIQUAC, Wilson, Van Laar | `thermo/phase/` |
| P6 | The full flash set and the phase envelope | `flashops/`, `saturationops/` |
| P7 | Associating — CPA, UMR-CPA, PC-SAFT, SAFT-VR-Mie | `thermo/phase/` |
| P8 | Electrolytes | `thermo/phase/` |
| P9 | Solids and flow assurance — hydrates, wax, asphaltene, scale, freezing | `flashops/`, `pvtsimulation/` |
| P10 | Reactions | `chemicalreactions/` |
| P11 | Unit operations | `process/equipment/` |
| P12 | Flowsheets | — |

Beyond P12: mechanical design, safety, cost, electrical, automation, the MCP server,
`standards/`, `statistics/`, `fluidmechanics/` and `pvtsimulation/`. They are last because
they are engineering deliverables rather than thermodynamics — a statement about **order**,
not about whether they belong.

---

# Part two: the keycard

## The keycard is where responsibility sits

**A keycard declares which data its holder is entitled to use, and therefore which physics
the library runs for them.** It is not a configuration file. The shipped baseline is
NeqSim's; a user's card is theirs; and the difference between the two is who is
accountable for a value being right.

This is why the card is the layer that matters. It determines:

- **the equation of state**, through `models` — which cubic, which alpha function, which
  mixing rule;
- **the chemical data**, through `components` and `kij` — critical constants, acentric
  factors, binary interaction parameters;
- **the equipment data**, through `fittings` — the Crane coefficients a pressure drop is
  computed from;
- **the property tables**, through `fluids`.

Change a card and every answer changes. A result resting on a value nobody was licensed to
supply is a result nobody should have produced, and no schema can discharge that.

## What a keycard is, formally

One TOML file. The format is `specs/schema/keycard.schema.json`; the readers are
`python/src/azoth/keycard.py` and `azoth_eos::card`. Its sections are:

| Section | Stage | Overrides |
|---|---|---|
| `keyholder` | runtime | nothing — it records whose card this is |
| `components` | runtime | `data/components/components.csv` |
| `kij` | runtime | `data/components/kij.csv` |
| `coefficients` | runtime | a calculation's named argument, as a default |
| `models` | runtime | the declared cubic variants |
| `fittings` | compiled | `data/fittings/crane_k_factors.csv` |
| `fluids` | compiled | `data/fluids/` |

**The two stages are part of the format.** A `runtime` section is read by the loaded
keycard the moment a calculation is called. A `compiled` section is consumed by a
generator that writes a shipped data file, so a change to one takes effect after a
regeneration **and a rebuild**. Each section carries an `x-azoth-stage` annotation in the
schema, and the compiler derives its work list by reading them.

## What a keycard cannot be

**It is data, and it cannot be anything else.** Every section is named choices from closed
vocabularies plus numbers. Nothing in a keycard can make this library do something it does
not already implement, and a name outside a vocabulary is refused when the card is
*loaded* rather than when it is finally used — a model that silently fell back to
Peng-Robinson would be a wrong answer with no symptom.

**It is not a place to record provenance.** A row may carry a `citation`; nothing requires
one, nothing validates one, and no field records a status. An earlier version required a
`verification.status` on every spec and a machine-fetchable `source_ref` beside it, with a
rule that the two agree. Nobody can check whether a person read a standard, so the field
was a form to fill in rather than a fact — and a form teaches people to fill it in, which
manufactures confidence rather than producing it.

**It is not shared, and it is not stored.** `azoth.keycard.load` reads a file and
returns what it says; the calls that should read it are handed it, and a call handed
none reads the data this library ships. A library whose answers depend on call order
is a library that returns two results for one calculation, and a card in force is
exactly that — so there is no card in force, no `current()` to read one and no
`clear()` to undo one. A caller with two datasets passes a different card to each
call rather than running two processes.

## What the library owes instead of a status field

Disclosure, and it is specific:

- `NOTICE` says what ships and where it came from.
- Every calculation's page names its source and states, in plain words, what has and has
  not been confirmed.
- Shipped data carries a `verify_status` column — `estimated_dummy`, `unverified` or
  `verified` — derived at generation time from whether a row's citation says the value is a
  placeholder. It drives the warning a user gets when a result rests on the Crane
  coefficients: seven values that are not from any standard and can make a pressure drop
  wrong by a factor of two while looking entirely reasonable.

The library's own data carries that column because it is the library's statement about
itself. A user's rows do not, because the engineer who supplied a value knows whether they
trust it, and a status field would only invite them to assert something nobody can check.
The two are not the same kind of thing and compressing them into one rule loses something
real.

## The rules a keycard follows

- **Everything the library reads, a keycard must be able to express.** A parameter the
  databank carries and the card cannot is a value a user cannot supply.
- **A name already in the databank is overridden parameter by parameter.** Correcting one
  value does not mean restating the others and does not silently lose them.
- **A new name must be complete.** A partial component is refused rather than completed
  from a similar substance, which would be inventing data.
- **An explicit argument always wins.** A caller who writes a value in their own call is
  stating it for that call; overriding it from a file loaded an hour ago is the worst kind
  of surprise.
- **A coefficient must declare a unit**, checked against the spec's declared unit for that
  input. A value in the wrong unit is a factor with no symptom.

---

# Part three: the calculus

**The calculus is normative for the types, and this page is normative for the
policy.** They do not overlap, and neither overrides the other. Where a signature,
a dimension or a channel declaration disagrees with
[The calculus of thermodynamic dimensionality](../calculus/index.md), the calculus
wins and the other thing is a bug. Where a rule about what is ported, how a
calculation is declared, or where responsibility for data sits disagrees with this
page, this page wins. A type is not a policy, and the two are not two answers to
one question.

Four rules from it are stated here because they constrain everything else on this
page:

- **One vocabulary table.** `specs/vocabulary/vocabulary.toml` is the one
  hand-written source of the canonical units a spec may declare, and the schema's
  unit enum, the Rust conversion table and the Python map are compiled from it. A
  unit added anywhere else is a unit that disagrees.

- **The calculus is hand-written; the vocabulary is data.** The same split
  [How a calculation is ported](#how-a-calculation-is-ported) makes for kernels,
  extended by one file: the Lean development is written by hand, and everything
  that merely names a unit is generated from the table.

- **No conversion factor is written down anywhere.** A factor is a number `uom`
  and `pint` each already know. A third copy is a copy that can disagree with
  both, and it is why `mm` once put a factor of a thousand between the two
  implementations while every test passed.

- **Pi and rho are a language, not a spec field.** The calculus at
  [Processes and channels](../calculus/process.md) states a unit operation as a
  process on typed, directional channels, and that statement is what tranche P11
  will build the tier against. **Nothing declares a port today**: the unit-operation
  tier was deleted rather than repaired, so there is no model carrying a port, no
  lint rule holding one against its inputs and outputs, and no interpreter for the
  calculus. A page describing a port is a specification, and it is marked as one.
