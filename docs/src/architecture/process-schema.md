# The process schema

This page states the schema a GUI reads and writes: the palette of unit operations,
the flowsheet that wires them, and the rules the checker holds a flowsheet to. It is
the concrete form of [the process calculus](../calculus/process.md), which is normative
for the types — a channel is a field record with a polarity, conservation is linearity,
and feedback is restriction. Where a format here disagrees with that page, the calculus
wins.

The schema lives in Rust. `crates/azoth-process` owns the types, the checker and the
first kernels; the format is TOML, deserialised by the same serde types that define the
contract. There is no JSON Schema and no Python linter for this layer — the Rust type is
the schema, and validation is `azoth_process::validate`.

## The material stream

The palette adopts one record for the arrows, the fields the balance lemma names:

| field | dimension | meaning |
|---|---|---|
| `n` | `molar_flow` | molar flow, mol/s |
| `z` | `dimensionless` (vector) | composition, one entry per component |
| `P` | `pressure` | pressure |
| `T` | `thermodynamic_temperature` | temperature |
| `h` | `molar_energy` | molar enthalpy, J/mol |

A port declares this record inline; the flowsheet checker enforces that two connected
ports carry the same record, field by field, comparing dimensions by exponent tuple
rather than by name.

## A unit-operation spec

One file under `specs/unit_ops/` declares one palette entry — its id, where its port
shape comes from, the parameters a user fills in, and its ports:

```toml
id = "unit_ops.pump"
name = "Pump"
[source]
standard = "NeqSim process/equipment/pump/Pump.java"
[parameters.outlet_pressure]
unit = "Pa"
description = "the pressure the pump raises the stream to"
[[ports]]
name = "inlet"
direction = "in"
[ports.fields]
n = { dimension = "molar_flow" }
z = { dimension = "dimensionless", shape = "vector" }
P = { dimension = "pressure" }
T = { dimension = "thermodynamic_temperature" }
h = { dimension = "molar_energy" }
[[ports]]
name = "outlet"
direction = "out"
[ports.fields]
# the same record
```

A port's `multiplicity` is `one` by default and `many` for a mixer's inlets or a
splitter's outlets. A field's `shape` is `scalar` by default and `vector` for
composition.

A parameter declares `required = true` when a kernel cannot run without it — the fact a
form's asterisk and the checker's `MissingParameter` both need, and one the declaration
carried no word for until P12's close-out. It is held to two sources that agree: the
`[inputs]` block of `specs/models/process/<id>.toml`, which marks an input `optional = true`
when a caller may leave it out, and the shim in
`crates/azoth-process/src/executor/dispatch.rs`, which reads a required parameter with
`Parameters::si`/`number`/`flag`/`text`/`vector` and an optional one with the `optional_*`
family. `tests/palette.rs` re-checks the agreement on every run.

The palette is grouped by family under `specs/unit_ops/`: two-port machines, separators,
mixer/splitter, the heat exchanger, columns, reactors and utility units. Two kinds of
NeqSim equipment are deliberately **not** palette entries: the solver/control blocks
(`Adjuster`, `SetPoint`, `Calculator`, `Recycle`, `VirtualStream`) are P12 machinery, and
the domain modules (renewables, reservoir) are not unit operations.

## A flowsheet spec

A flowsheet is a **self-contained simulation**: it declares its inputs — what the user
specifies — and its products, which the run calculates. Boundary streams are named and
instances wired by named connections:

```toml
id = "flowsheets.demo"
name = "Feed, pump, heat, separate, recycle"
products = ["vapour_product"]
[[inputs]]
name = "feed_1"
components = ["methane", "n-butane"]
n = 1.0
z = [0.9, 0.1]
P = 5.0e5
T = 300.0
[[instances]]
id = "p1"
unit = "unit_ops.pump"
[instances.parameters]
outlet_pressure = 2.0e6
[[connections]]
from = "feed_1"
to = "p1.inlet"
[[recycles]]
stream = "recycle_1"
from = "sep1.liquid"
to = "mix1.feed"
```

An input's fields are the material stream's own — `n`, `z`, `P` and `T` from the record
above — so no unit is written and the field's dimension supplies it, exactly as an
instance parameter's bare number takes its unit from the palette entry. **`h` is the one
field of the record a user does not write**: it is a state function of `(T, P, z)`, so the
kernel derives it, and a written `h` is refused rather than ignored. **`products` is names
only**, because an output is calculated and has nothing to state.

`products` is written *before* `[[inputs]]`, and that order is not cosmetic: TOML has no
way to reopen a key after an array of tables, so a value written after one belongs to the
last element of it. For the same reason `[layout]` is written *after* every array of tables:
it is a table, and the writer emits its tables in field order.

**`[layout]` is a view and not physics.** `[layout.instances]`, `[layout.inputs]` and
`[layout.products]` each map a name to an `[x, y]` pair, and nothing in the run reads one: the
checker ignores it, the executor never sees it, and a flowsheet written by hand carries none —
the projection derives a deterministic layout from the connection graph and uses a stored
position only where one exists. It is in the document so that one artifact is the whole editor
state, which is what makes the round trip lossless for a figure drawn rather than only for a
document typed. **A table and not a pair of fields on `[[instances]]`**, because `products` is
names only — an output is calculated and has nothing to state — so a per-instance position
could not cover the boundary, and one table covers all three kinds uniformly.

`from` is a feed name or `instance.port` (an outlet); `to` is a product name or
`instance.port` (an inlet). A recycle is a connection that closes a loop, its stream
named as the tear.

**An outlet declared `multiplicity = "many"` is addressed by position**: `from =
"split1.products[0]"` is the first of the streams that port returns, and each one is bound
under that same path, so a connection writes what the run produced. A `one` outlet has no
position and an inlet has none either — a `many` *inlet* takes several connections to the
same port name, which is how a mixer is wired — and the checker refuses an index written
where it does not belong rather than dropping it.

A `[[recycles]]` entry also carries the tear's convergence, all of it optional and every
default `Recycle`'s own: `flow_tolerance`, `composition_tolerance`, `temperature_tolerance`
and `pressure_tolerance` (each `1e-2`), `max_iterations` (`10`), `minimum_flow` (`1e-20`
kg/hr) and `acceleration_method` (`direct_substitution`; `wegstein` is carried, `broyden`
is refused by name with the measurement that closes it). **`minimum_flow` is a switch and
not a tolerance**: a tear below it is deactivated outright, its residuals *declared* zero
rather than measured, and its loop stops at one pass.

## The checker

`azoth_process::validate` runs the calculus's rules:

- every instance names a palette unit op, and **exactly** its declared parameters — every one
  it cannot run without given, none it does not declare, and **each value of the kind its
  declaration states** (the kind is not in the palette, so it is read from the model's own input
  declaration — the same table a form chooses a widget from);
- a connection joins an outlet (or feed) to an inlet (or product);
- a Port-to-Port connection joins dimension-compatible field records;
- an input's record could be a stream — one mole fraction per substance, and a fluid that
  names at least one;
- **linearity** — a `one` port is consumed/produced exactly once, a `many` port at
  least once, and every feed/product is used exactly once;
- no name is used for two things a reader has to tell apart — an instance and a feed, an instance
  and a product, or a feed and a product;
- **a connection and a tear may each be declared once**: a connection is identified by its two
  endpoints and a tear by its name, because those are what the executor resolves them by, so a
  second entry under the same identity feeds one port the same stream twice. Measured — a
  duplicated connection takes a mixer's outlet from 1 mol/s to 2 and a duplicated tear from 1.5 to
  2, and neither `OverfedPort` nor the feed rule fires on a `many` inlet, so both documents used to
  pass and compute something nobody asked for;
- a tear's `acceleration_method` is one this port can apply, read from the class's own two refusals
  rather than from a second list of names;
- every loop in the instance graph passes through a declared recycle.

The last rule about an input is a *shape* rule and no more: the checker compares dimensions
by exponent tuple and resolves no names, so whether a named substance exists is the
databank's answer and is refused at the run.

`azoth check --flowsheet <file>` runs it from the command line; a cargo test walks
`specs/unit_ops/` and `specs/flowsheets/` and holds every shipped file to it.

## The wire projection

A front-end does not read the document; it reads a projection of it, and the projection is stated
here because two of its rules are what make an editor mechanical rather than a translation.

```text
a node    instance:sep1          input:feed_1          product:vapour_product
a handle  sep1.liquid            s1.products[0]        feed_1
an edge   e0 … over the connections then the recycles, in document order
```

**A node's id is `{role}:{name}`.** Instance-versus-boundary is disjoint by the name rule above;
boundary-versus-boundary is not, which is why the prefix cannot be dropped — a bare name would
merge an input and a product of the same name into one node.

**A handle's id is the stream's own path** — `feed_1`, `sep1.liquid`, `s1.products[0]` — so the
handle an edge attaches to, the session's key for the same stream and the edge's `data.path` are
one string rather than three spellings reconciled in three places. A feed's handle is the feed's
name and a product's is the product's, which makes every edge the same shape (`source` → `target`)
and needs no special case at the boundary. A `one` port draws one handle; a `many` **inlet** draws
one however many connections reach it, which is how a mixer is wired; a `many` **outlet** draws one
per stream it returns.

**How many streams a `many` outlet returns is a value and not a declaration** — the checker says so
itself — so the projection draws the union of what the document addresses and what the instance's
vector parameter says, and an entry with two vector parameters would be ambiguous and leave the
wiring alone to decide.

The graph is emitted in **xyflow's** node/edge shape — `id`, `type`, `position`, `data`, `source`,
`target`, `sourceHandle`, `targetHandle` — so the editor's model is the schema's. A `[[connections]]`
entry has no id in the schema, so its position in the document is its identity and the edge id is
that position.

**The projection never refuses a document.** A canvas has to draw a broken flowsheet, because
drawing it is how a user fixes it: an instance naming a unit op the palette does not carry gets a
node with no ports, and an edge naming an instance nobody declared keeps the id it implies — which
the checker's `UnknownInstance` explains, and which is reported rather than repaired.

## The kernels

A kernel is a unit operation's arithmetic, a pure function of its inlets. **Twenty-seven of the
29 palette entries carry a registered `process.*` id** — a spec, a kernel, a Python reference, a
case set and a NeqSim capture — composing the calculations in `azoth-eos` rather than adding new
physics. Two do not: `unit_ops.simple_absorber` is refused on measured evidence, and
`unit_ops.rate_based_packed_column` is a second physics carried by the distillation workstream
that added it beside the column.

**The flowsheet executor is built**, and it is `crates/azoth-process`'s `executor` module: the
`unit_ops.*` → kernel dispatch table, the execution order, the tear with `Recycle`'s own four
tolerances, a session whose every value has a stable path, a JSON codec and structured
diagnostics. It runs `specs/flowsheets/demo.toml` and converges its recycle, held to a NeqSim
`ProcessSystem` capture, and the document runs from three places — the library, `azoth run
--flowsheet` and `azoth.process.run_flowsheet` — with no argument but the file. **One entry it
refuses by name**: `unit_ops.packed_column` has a kernel
and a registered id, and its declaration describes the packing rather than the column, so the
executor says so rather than running a machine the declaration does not describe. What is still
owed is named in `ROADMAP.md`.
