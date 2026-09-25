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
last element of it.

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
  it cannot run without given, and none it does not declare;
- a connection joins an outlet (or feed) to an inlet (or product);
- a Port-to-Port connection joins dimension-compatible field records;
- an input's record could be a stream — one mole fraction per substance, and a fluid that
  names at least one;
- **linearity** — a `one` port is consumed/produced exactly once, a `many` port at
  least once, and every feed/product is used exactly once;
- every loop in the instance graph passes through a declared recycle.

The last rule about an input is a *shape* rule and no more: the checker compares dimensions
by exponent tuple and resolves no names, so whether a named substance exists is the
databank's answer and is refused at the run.

`azoth check --flowsheet <file>` runs it from the command line; a cargo test walks
`specs/unit_ops/` and `specs/flowsheets/` and holds every shipped file to it.

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
