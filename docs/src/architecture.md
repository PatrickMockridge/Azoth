# How azoth is put together

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

[The specification](./spec.md) is the normative page: where this one disagrees
with it, this one is wrong. What this page is for is orientation — what the
pieces are, how they nest, and where a new piece of your own would go.

## A core in the middle

`azoth-core` holds the vocabulary the other crates share and nothing
domain-specific: units, errors, warnings, results, range checks, the spec runtime
and the solvers. The crate's own doc comment states the rule — it
"deliberately contains no engineering calculations".

**A domain crate depends on `azoth-core` and on no sibling**, with one exception.
`azoth-process` also depends on `azoth-eos`, because a unit operation is a flash
call plus arithmetic, and a process layer that could not call the flashes would
not be a process layer. So a new domain is a **new crate** rather than a fourth
branch inside an existing module.

The domain crates are `azoth-eos`, `azoth-hydraulics`, `azoth-thermal` and
`azoth-process`. The other two crates are not domains, and they do depend on
siblings — `azoth-cli` and `azoth-python` are entry points, and an entry point has
to reach everything it exposes:

| Crate | Depends on | What it is |
|---|---|---|
| `azoth-core` | — | units, errors, warnings, results, range checks, the spec runtime, the solvers |
| `azoth-eos` | core | equations of state, flashes, mixture properties |
| `azoth-hydraulics` | core | pipe flow, fittings, valves, pumps |
| `azoth-thermal` | core | conduction |
| `azoth-process` | core, **eos** | unit operations — the one domain crate that takes a sibling |
| `azoth-cli` | core, hydraulics | the `azoth` binary |
| `azoth-python` | core, eos, hydraulics, thermal, process | the compiled extension |
| `azoth-test-support` | core | a dev-dependency of the namespace crates only |

## Four levels

| Level | What it is | An example |
|---|---|---|
| **A calculation** | one equation, scalars in, a result out | [`hydraulics.darcy_weisbach`](./hydraulics/darcy_weisbach.md) |
| **A model** | a *procedure*, or a computation over vectors | [`eos.pt_flash`](./eos/pt_flash.md) |
| **A unit operation** | a transformation of streams | [`process.separator`](./process/separator.md) |
| **A flowsheet** | units, connections, and a solver over them | — |

The first three are built. **The fourth is not.** S4 argues that a flowsheet
should be a *document* — YAML naming units, their connections and the solver —
rather than an object graph of units holding references to each other, but that
is a design in a specification rather than a schema, a format or any code:
`specs/` holds `calcs/` and `models/` and nothing else, and nothing in the tree
executes such a document. [Roadmap](./roadmap.md) records it as not started.

The difference between the first two levels is narrower than it looks, and worth
stating precisely because the names invite the wrong answer. A spec under
`specs/models/` whose `kind` is `procedure` names an `algorithm`: a scheme, a
convergence rule, a tolerance and an iteration cap. That is not, however,
something a calculation cannot have — a calc spec may name a `solver` block when
its equation is implicit, and two of them do
([`hydraulics.friction_factor_colebrook`](./hydraulics/friction_factor_colebrook.md)
and [`eos.pr_z_factor`](./eos/pr_z_factor.md)). What a model adds is a
**procedure over vectors**: a composition vector and a matrix of interaction
parameters have nowhere to go in the calc registry, whose inputs are scalars.

## Where your thing goes

Most changes are smaller than they look, and the commonest mistake is writing
code for something that is data.

| You want to | You write | You edit |
|---|---|---|
| Change a component's `Tc`, `Pc` or `omega` | a row in a keycard | nothing |
| Add a binary interaction parameter | a row in a keycard | nothing |
| Set a coefficient once instead of passing it everywhere | a section in a keycard | nothing |
| Name a cubic variant and its substances once | a keycard `models:` entry | nothing |
| Add or correct a fluid property table | a section in a keycard | re-run `tools/gen_user_data.py`, rebuild |
| Supply a fitting's equivalent length | a row in a keycard | re-run `tools/gen_user_data.py`, rebuild |
| Plug in your own source of fluid properties | a `PropertyProvider` | nothing |
| Add an equation | the spec, one Rust file, one Python file | nothing |
| Add a unit operation | the same, over streams | nothing |
| Add a whole new domain | the spec files and a new crate | its `Cargo.toml` and the workspace members |

The two rows that need a rebuild are the ones that feed the **data files** rather
than the running library. A keycard's `components`, `kij`, `coefficients` and
`models` sections are read at run time. Its `fluids` and `fittings` sections are
not: `keycard.yaml` is a *source* that `tools/gen_user_data.py` compiles into
`data/fluids/*.csv` and `data/fittings/*.csv`, and the Rust side embeds those with
`include_str!` while the Python side locates them on disk. That is what makes the
two implementations read byte-identical data, and it is also why a fluid or a
fitting cannot take effect without a rebuild.

**Nothing is registered.** A calculation is found by its id, and the id *is* its
address: `hydraulics.darcy_weisbach` names
`azoth/hydraulics/reference/darcy_weisbach.py` and
`azoth-hydraulics/src/darcy_weisbach.rs`, by convention. Adding one adds files and
edits no list by hand. The generated registries do change — `spec_gen.rs`,
`_registry_gen.py` and `_models_gen.py` are rewritten by the generators — but a
person does not maintain them.

The last row is the only one that touches the build system, and that is the point
of the crate rule rather than an inconvenience.

## How it extends

There are three ways in, and they are the ones the specification sets out rather
than a taxonomy invented here.

**Data, through the keycard.** A component parameter, an interaction parameter, a
fluid table, a fitting coefficient, a coefficient a calculation takes, a named
cubic variant — everything overrides or extends what ships, by name, and none of
it is code. What a keycard cannot do is make the library do something it could not
already do. [The keycard](./keycard.md) is the page;
[specification, S5](./spec.md) is the rule.

**Code, by contributing.** A spec plus one Rust file and one Python file. Both
implementations are held to the same tests, and the cross-language comparison runs
every case through each. [Specification, S9](./spec.md) states the whole contract.

**A fluid property source, on the Python side.** `PropertyProvider` is a
protocol with a `name` property, a `density()` and a `dynamic_viscosity()`, for
code that turns a fluid name into properties. It exists because that is the one
place where the right answer depends on data this library cannot ship.

What is *not* here is a runtime plugin registry: there is no `register()` call, no
dynamic loading, and no trait or protocol a user implements to add a calculation
at run time. The consequence is worth stating plainly rather than leaving to be
discovered — if what you want to add is an equation, you are writing Rust and
Python, not configuring something; if what you want to add is data, you never
touch either.

## What holds it together

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
