# Specification

This page says what azoth **is**. It is normative: where another page, a docstring, a
comment or a habit disagrees with it, this page wins and the other thing is a bug. It is
deliberately short on hedging and long on decisions, because a document that describes
what a project would like to be is a document nobody can act on.

If you are adding a calculation, read [S9](#s9-what-a-contribution-costs) and stop
there. If you are deciding whether something belongs in azoth at all, read
[S8](#s8-the-scope). If you are wondering why it is written in Rust, read
[S2](#s2-why-rust-and-not-java).

---

## S1. What azoth is

**An opinionated port of NeqSim to Rust, with every calculation mirrored in Python.**

Both words are load-bearing.

*Port*: the algorithms come from [NeqSim](https://github.com/equinor/neqsim), Equinor's
open-source Java process simulator, under Apache-2.0. They are carried across rather
than re-derived, because a second derivation of the Soave alpha function is not
knowledge, it is a second chance to make a mistake. What that obliges is attribution,
which the repository discharges once in `NOTICE` rather than per file — see
[azoth and NeqSim](./comparison/neqsim.md).

*Opinionated*: where azoth's structure differs from NeqSim's, that is a decision that
can be defended, and this page defends it. The differences are not accidents of a
different language. They are:

| | NeqSim | azoth |
|---|---|---|
| Truth lives in | the Java source | `specs/**/*.yaml`, which generates the code and the docs |
| Nothing is registered by | — | anyone; the id is the address (S5) |
| Data ships as | a bundled databank | the same databank, vendored; the keycard is the user's |
| Correctness is claimed by | a test suite | two implementations that must agree case by case |
| The unit of composition is | an object graph of units | a document (S4) — designed, not built |

## S2. Why Rust, and not Java

NeqSim is a large, working, thirty-year-old Java codebase. Writing azoth in Rust
throws that away and pays for it. The case for doing it anyway, and the price:

**For:**

- **`pip install azoth` does not need a JVM.** NeqSim on Python means a JDK *and*
  JPype, and a user who wanted to compute a Reynolds number has installed a Java
  runtime to do it. That is the single largest practical difference and it is why the
  Python mirror exists at all.
- **One wheel per platform, across Python versions.** The extension is built against
  `abi3-py312`, so a single artifact serves every later Python rather than a matrix.
- **A small distribution.** NeqSim ships a shaded JAR plus a runtime. azoth's core is a
  few megabytes, because it vendors data rather than a platform.
- **Invariants the compiler holds, not conventions it hopes for.** `unsafe_code =
  "forbid"` is a lint setting, not a policy. Physical quantities cross the API as `uom`
  types, so passing a bare float where a pressure belongs does not compile. `match` is
  exhaustive and there are no nulls. In Java each of these is a thing a careful
  programmer does; here each is a thing a careless one cannot avoid.
- **No warm-up and no collector.** A CLI that starts in milliseconds is a CLI that can
  be used in a shell pipeline.
- **One static binary** for the command line, with no runtime to ship alongside it.

**Against, stated because a specification that only lists advantages is advertising:**

- **Rust is harder to contribute to**, and the pool of chemical engineers who write it
  is smaller than the pool who write Python or Java. This is a real cost and it is paid
  on every contribution.
- **The borrow checker is genuinely bad at a mutable object graph of mutually
  referencing units** — a stream that knows its source, a unit that owns its streams,
  a flowsheet that owns its units and is reachable from all of them. Fighting that is
  how you get `Rc<RefCell<...>>` everywhere and a codebase nobody can reason about.
  **This is why a flowsheet is a document rather than an object graph** (S4), and it is
  the clearest example of a language shaping the design instead of merely hosting it.
- **Compile times** are real, and a full `cargo test --workspace` is not instant.

## S3. A core in the middle

`azoth-core` holds the vocabulary every other crate needs and nothing domain-specific:
units, errors, warnings, the spec runtime, range checks, and the solvers. **Every
domain crate depends on it and on no sibling.** `azoth-eos` does not know that
hydraulics exists, and cannot come to know without an edit to its `Cargo.toml` that a
reviewer will see.

```
                          azoth-core
          units · errors · warnings · spec runtime · solvers
                                │
        ┌──────────────┬────────┴────────┬──────────────┐
    azoth-eos     azoth-thermal   azoth-hydraulics   azoth-cli
        │
   azoth-process
        └──────────────┴────────┬────────┴──────────────┘
                          azoth-python
```

The rule does work rather than decorating a diagram:

- **A new domain is a new crate**, not a fourth branch inside an existing module. The
  boundary is in the build system, so it is checked by the compiler rather than by
  review.
- **Nothing can reach across**, so a change in one domain cannot silently alter
  another's behaviour. Cross-domain composition lives above the domains, in
  `azoth-python` today and in `azoth-process` for the process layer (S4).
- **`azoth-core` stays small.** It may not grow a dependency on a domain, and it may
  not grow a calculation.

**`azoth-process` is the one crate that depends on a sibling, and it is a deliberate
exception rather than a lapse.** A unit operation is a flash call plus arithmetic -
that is the whole of the Pareto argument for the scope of this port, and it is read
from NeqSim's source rather than asserted. A process layer that could not call the
flashes would not be a process layer, so `azoth-process` depends on `azoth-eos`. It
belongs to the tier above the domains rather than inside one, and the diagram above
shows it that way: it is drawn *under* `azoth-eos` because that is what it consumes,
and it joins the composition tier that `azoth-python` and `azoth-cli` already occupy
rather than sitting beside the domains.

`azoth-test-support` (shared test fixtures) is depended on by tests only and is not part
of the runtime layering.

## S4. Composable, not monolithic

Four levels, each composing the one below rather than replacing it. Every level is
useful on its own, and the API at each level is the one below plus structure.

1. **A calculation** — one equation, scalars in, a result out. `darcy_weisbach`.
2. **A model** — a *procedure*, or a computation over vectors. A PT flash iterates; a
   critical point searches. These are not calcs because a model's spec is where a
   procedure over **vectors** lives: a composition vector and a matrix of interaction
   parameters have nowhere to go in the calc registry, whose inputs are scalars. Where
   either kind of spec iterates, it names the scheme, the convergence rule and the
   tolerance, because two implementations that run different loops diverge.
3. **A unit operation** — a transformation of streams. A separator takes one stream and
   gives two; a compressor takes a stream and a duty.
4. **A flowsheet** — units, connections, and a solver over them.

**A flowsheet will be a document, and it is designed rather than built.** There is no
flowsheet schema, no format and no code for one in this tree; [Roadmap](./roadmap.md)
records it as not started. The design is YAML: units, their connections, and the solver.
Both implementations would execute the same document, which is what would make a
flowsheet diffable, hashable, reviewable in a pull request and reproducible from a file
rather than reconstructible from a script.

This is a deliberate refusal, and S2 is where the reason lives: the obvious design — a
mutable graph of units holding references to each other — is exactly the shape Rust is
worst at, and building it would produce a codebase held together by interior mutability.
A document has none of that. The imperative API in Python builds a document; it does not
replace one.

The cost is real and is recorded in [Roadmap](./roadmap.md): a flowsheet has no single
retraceable worked example the way a calculation does. What replaces it is a
hand-computable case plus conservation checks that hold at every answer.

## S5. It ships data, and the keycard extends it

**These two things are the entire registry. There is no registration step.**

**What ships: the vendored NeqSim databank.** 173 components and 516 binary interaction
parameters, generated into `data/components/` from NeqSim's `COMP.csv` and `INTER.csv`
by `tools/gen_databank.py`. Each component row carries its citation verbatim — NeqSim
v3.20.0, Equinor and NTNU, Apache-2.0 — and the interaction-parameter table carries none,
because `COMP.csv` has no citation column to copy; that attribution lives in `NOTICE`.
It is vendored rather than depended upon, because NeqSim is Apache-2.0,
because the data is the output of work at Equinor and NTNU rather than something
invented, and because a calculation library that ships no components is a calculator.

**What a user adds: the keycard.** One YAML file. Its sections are `keyholder`,
`components`, `kij`, `fluids`, `fittings`, `coefficients` and `models`. Everything in it
overrides or extends what ships, by name. A keycard is how a user adds to azoth
**without writing Rust or Python** — a component parameter, a fluid property table, a
fitting coefficient, a discharge coefficient and a cubic-EOS variant are all data.

**No contributor registers anything.** A calculation is found by its id, and the id is
its address:

```
"hydraulics.darcy_weisbach"
   → python/src/azoth/hydraulics/reference/darcy_weisbach.py
   → crates/azoth-hydraulics/src/darcy_weisbach.rs
```

The module path and the function name **follow from the id**. Adding a calculation adds
files and requires no registration step. What enumeration exists — the generated
registries, and the list of what is implemented in this book — is a build artefact a
generator emits, not a list a person maintains. Deleting registration rather than
generating it is deliberate: a generated list is still a list that can be wrong.

**The burden and the responsibility sit with the data.** If a value is wrong it is wrong
in the vendored databank or in a keycard. Equinor and NTNU are accountable for the
first; the engineer holding the keycard is accountable for the second, and for their
right to use it. Not the library, and not a schema.

## S6. Provenance is the engineer's job, not the library's

A keycard **may** carry a citation. Nothing requires one, nothing validates one, and no
field records a status.

This is the position rather than an omission, and it was arrived at by removing the
opposite. An earlier version required a `verification.status` on every spec and a
machine-fetchable `source_ref` beside it, with a rule that the two agree. Nobody can
check whether a person read a standard. So the field was a form to fill in rather than a
fact, and a form teaches people to fill it in — which is worse than having no field at
all, because it manufactures confidence.

**What the library owes instead is disclosure, and it is specific:**

- `NOTICE` says what ships and where it came from.
- Every calculation's page names its source and states, in plain words, what has and
  has not been confirmed — see [S7](#s7-every-rust-calculation-is-mirrored-in-python)
  on how that is checked, and the `## Notes` section on any calculation page for what
  it looks like.

The rest is engineering judgement, which is a professional responsibility and not one a
YAML linter can discharge.

**One asymmetry is deliberate, and is stated here so that it is not mistaken for a
leftover.** The data *this repository ships* carries a `verify_status` column; a
keycard's rows do not. The two are not the same kind of thing:

- A keycard is the user's, and the library does not ask. The engineer who supplied a
  value knows whether they trust it, and a status field would only invite them to
  assert something nobody can check.
- The shipped data is the *library's own* statement about itself, and it is what makes
  disclosure concrete rather than a promise. It is derived at generation time from
  whether a row's citation says the value is a placeholder, and it drives the warning a
  user gets when a result rests on the Crane coefficients — seven values that are not
  from any standard and can make a pressure drop wrong by a factor of two while looking
  entirely reasonable.

Compressing those two into one rule either way loses something real. Asking users for a
status produces forms; dropping the shipped column produces a library that ships
placeholders and says nothing.

## S7. Every Rust calculation is mirrored in Python

Not a binding to a black box, and not a second-class path. **Both implementations are
real, both are readable, and they are compared against each other case by case.** Where
a procedure iterates, the **iteration counts are required to match**, because two
implementations that agree on an answer reached by different paths have not been shown
to agree on anything durable.

Three reasons, in order of weight:

1. **Two independent implementations are the one property NeqSim cannot have.** A
   single implementation can be wrong in a way its own tests encode. Two, written from
   one specification, disagree at the point where one of them is wrong.
2. **conda and Jupyter with no toolchain.** A user who wants azoth in a notebook should
   not need a Rust compiler, and the Python path is the one they will actually read.
3. **The Python path is the readable one.** `numpy`-free, dependency-light, written to
   be read beside the equation it implements.

The rule has a consequence worth stating: **there is no third lane.** There is no
"Python-only, for convenience" and no "Rust-only, for speed". A calculation that exists
in one language does not exist.

## S8. The scope

**This is a statement of where effort goes, not a gate.** It exists so that "should we
port this?" has an answer that does not have to be re-argued each time.

**In scope:**

- The cubic equation-of-state core: Peng-Robinson and its variants, mixing rules,
  departures, the Helmholtz energy and what follows from it.
- Flashes: TP, PH, PS. Stability and the critical point.
- Properties: density, enthalpy, entropy, heat capacity, viscosity, thermal
  conductivity.
- The Pareto set of unit operations: mixer, splitter, separator, valve, heater, cooler,
  compressor, pump, expander, heat exchanger.
- Flowsheets, including recycle.
- Reports, the command line, the agent surface.
- Data and the keycard.

**Out of scope**, recorded with reasons in [Roadmap](./roadmap.md): distillation
columns, transient pipeline simulation, reactors, networks, power systems, CPA and SAFT
and GERG and electrolytes, mechanical design, cost estimation, safety systems, PVT
laboratory analysis, hydrates, wax and asphaltene.

The boundary is drawn on a Pareto argument rather than on difficulty: **nearly every
NeqSim unit operation is a flash call plus arithmetic**, and the surrounding thousands
of lines are performance charts, entrainment models, geometry sizing and mechanical
design. What is in scope is what makes a flowsheet run. What is out is what makes one
acceptable to a detailed-design review, and that is a different product.

## S9. What a contribution costs

| To add | You write | You edit |
|---|---|---|
| A calculation | the spec, one Rust file, one Python file | the glue that names them, below |
| A unit operation | the same, over streams | the same |
| A component, `kij` pair, coefficient or model variant | a keycard | **no code at all** |
| A fluid table or a fitting coefficient | a keycard | re-run `tools/gen_user_data.py`, then rebuild |
| A flowsheet | nothing yet — it is designed, not built | — |

**No registry decides what exists.** There is no dispatch table and no registration
call, because a calculation's id *is* its address. What is **not** free is the glue that
attaches the two implementations to their names in each language — the PyO3 wrapper, the
transport class, the namespace's `__all__` and the two batch arms — which is typed by
hand, and which `python/tests/test_registration_completeness.py` and the batch coverage
test fail loudly on when one is missed, naming the file to go and edit.
[CONTRIBUTING.md](https://github.com/PatrickMockridge/Azoth/blob/main/CONTRIBUTING.md)
has the current list, measured rather than remembered.

What is *not* free, and is not meant to be:

- **A worked example, and it must be retraceable by hand.** This is the one requirement
  the linter kept when it dropped the rest, because it is the cheapest check in the
  registry and it is what pins a number to something a reader can follow. Skipping it is
  allowed; skipping it silently is not.
- **Both implementations**, per S7.
- **Units-safe signatures and typed errors.** A bare float where a pressure belongs does
  not compile; a range violation is a warning on the result, never a silent clamp and
  never an exception.

Everything else — the docs page, the registry entry, the range checks, the test cases,
the type stubs — is generated from the spec. The generator emits shape; **you write
arithmetic.**
