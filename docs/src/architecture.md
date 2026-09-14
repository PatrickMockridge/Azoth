# How azoth is put together

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

[The specification](./spec.md) is the normative page: where this one disagrees
with it, this one is wrong. What this page is for is orientation — what the
pieces are, how they nest, and where a new piece of your own would go.

## A core in the middle

`azoth-core` holds the vocabulary every other crate needs and nothing
domain-specific: units, errors, warnings, the spec runtime, range checks, and
the solvers. It contains no engineering calculation, and that is a rule rather
than a description of how it happens to look today.

**Every domain crate depends on it and on no sibling.** `azoth-eos` does not
know that hydraulics exists, and cannot come to know without an edit to its
`Cargo.toml` that a reviewer will see. So a new domain is a **new crate**, not a
fourth branch inside an existing module — the boundary lives in the build
system, where the compiler checks it rather than a reviewer having to.

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

`azoth-process` is the one crate that depends on a sibling, and it is a
deliberate exception rather than a lapse. A unit operation **is** a flash call
plus arithmetic; a process layer that could not call the flashes would not be a
process layer. It belongs to the tier above the domains, and the diagram draws it
that way. `azoth-test-support` is depended on by tests only and is not part of
the runtime layering at all.

## Four levels, each composing the one below

Every level is useful on its own, and the API at each level is the one below plus
structure. Nothing is replaced on the way up.

| Level | What it is | An example |
|---|---|---|
| **A calculation** | one equation, scalars in, a result out | [`hydraulics.darcy_weisbach`](./hydraulics/darcy_weisbach.md) |
| **A model** | a *procedure*, or a computation over vectors | [`eos.pt_flash`](./eos/pt_flash.md) |
| **A unit operation** | a transformation of streams | [`process.separator`](./process/separator.md) |
| **A flowsheet** | units, connections, and a solver over them | — |

The line between the first two is where the arguments are. A calculation's spec
fixes an **equation**; a model's fixes a **loop** — a scheme, a tolerance, an
iteration cap — because two implementations that run different loops diverge, and
one that agrees only to within its own tolerance has not been shown to agree at
all. A model's inputs are therefore vectors and matrices where a calculation's
are scalars, which is the other reason a composition vector has nowhere in the
calc registry to go.

**The fourth level is specified and not built.** A flowsheet is a *document* —
YAML naming units, their connections and the solver — rather than an object graph
of units holding references to each other, and the reasoning is in
[specification, S4](./spec.md). What is true today is that nothing executes such
a document: the ladder currently stops at unit operations, and
[Roadmap](./roadmap.md) has the programme rather than this page.

## Where your thing goes

Most changes are smaller than they look, and the commonest mistake is writing
code for something that is data.

| You want to | You write | You edit |
|---|---|---|
| Change a component's `Tc`, `Pc` or `omega` | a row in a keycard | nothing |
| Add a binary interaction parameter | a row in a keycard | nothing |
| Add or correct a fluid property table | a section in a keycard | nothing |
| Supply a fitting's equivalent length | a row in a keycard | nothing |
| Set a coefficient once instead of passing it everywhere | a section in a keycard | nothing |
| Name a cubic variant and its substances once | a keycard `models:` entry | nothing |
| Add an equation | the spec, one Rust file, one Python file | nothing |
| Add a unit operation | the same, over streams | nothing |
| Plug in your own source of fluid properties | a `PropertyProvider` | nothing |
| Add a whole new domain | the spec files and a new crate | its `Cargo.toml` and the workspace members |

**Nothing is registered.** A calculation is found by its id, and the id *is* its
address: `hydraulics.darcy_weisbach` names
`azoth/hydraulics/reference/darcy_weisbach.py` and
`azoth-hydraulics/src/darcy_weisbach.rs`, by convention. Adding one adds files
and edits no list, because there is no list. What enumeration exists — for
`azoth describe` and for the list of what is implemented in this book — is a
build artefact a generator emits, not something a person maintains. Deleting
registration rather than generating it was deliberate:
[specification, S5](./spec.md) says why.

The last row is the only one that touches the build system, and that is the
point of the crate rule rather than an inconvenience.

## The extension philosophy, stated as a refusal

**There are three ways to extend azoth, and there is deliberately no fourth.**

1. **Data, through the keycard.** A component parameter, an interaction
   parameter, a fluid table, a fitting coefficient, a coefficient a calculation
   takes, a named cubic variant. Everything here overrides or extends what ships,
   by name, and none of it is code — a keycard cannot make the library do
   anything it could not already do. [The keycard](./keycard.md) is the page.

2. **Code, by contributing.** A spec plus one Rust file and one Python file.
   Compile-time, generator-driven, and held to the same tests as everything else.

3. **A fluid property source, on the Python side.** `PropertyProvider` is a
   two-method protocol — a density and a viscosity — for code that turns a fluid
   name into properties. It exists because that is the one place where the right
   answer depends on data this library cannot ship.

**What there is not, and will not be, is a runtime plugin registry.** No
`register()` call, no dynamic loading, no trait a user implements to add a
calculation at run time. That is a refusal with reasons rather than an item
nobody has got to:

- **It would create a third lane, which [S7](./spec.md) forbids.** Every
  calculation here exists twice, once in Rust and once in Python, and the two are
  compared against each other case by case. A plugin is written in one language.
  Admitting it would mean admitting calculations that exist in one
  implementation, which is exactly the class of thing the two-implementation rule
  exists to exclude — and the rule earns its keep precisely where it is
  inconvenient.
- **It would be a second source of truth for something that already has one.**
  The id is the address. A registry is a list that can disagree with the tree it
  describes, and a list that can be wrong is a failure this project is organised
  against.

The consequence is worth stating plainly rather than leaving to be discovered: if
what you want to add is an equation, you are writing Rust and Python, not
configuring something. If what you want to add is data, you never touch either.

## What holds it together

Two rules do most of the work, and both are enforced rather than documented and
hoped for.

**The specs are the source of truth.** `specs/calcs/**/*.yaml` and
`specs/models/**/*.yaml` define each calculation: equation, source, valid range,
assumptions, worked example, tests. The registries, the range checks, the type
stub and every page under the namespaces are generated from them, and CI
regenerates and fails on any difference. The generator emits **shape**; the
arithmetic is written by hand and held to the spec by tests.

**Each calculation is data, and the data is yours to hold.** If a value is wrong
it is wrong in the databank azoth vendors or in a keycard. Equinor and NTNU are
accountable for the first — [azoth and NeqSim](./comparison/neqsim.md) says what
that means and [What ships](./data.md) says what is in it. The engineer holding
the keycard is accountable for the second, and for their right to use it.
[Specification, S6](./spec.md) is where that position is argued, including why no
field records a verification status.

## See also

- [The specification](./spec.md) — the normative version of everything above.
- [The keycard](./keycard.md) — extending it without writing code.
- [The batch API](./batch.md) — one calculation over arrays, and what it gives up.
- [Roadmap](./roadmap.md) — what is next, and what is deliberately not.
