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

The palette is grouped by family under `specs/unit_ops/`: two-port machines, separators,
mixer/splitter, the heat exchanger, columns, reactors and utility units. Two kinds of
NeqSim equipment are deliberately **not** palette entries: the solver/control blocks
(`Adjuster`, `SetPoint`, `Calculator`, `Recycle`, `VirtualStream`) are P12 machinery, and
the domain modules (renewables, reservoir) are not unit operations.

## A flowsheet spec

A flowsheet names boundary streams and wires instances by named connections:

```toml
id = "flowsheets.demo"
feeds = ["feed_1"]
products = ["vapour_product"]
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

`from` is a feed name or `instance.port` (an outlet); `to` is a product name or
`instance.port` (an inlet). A recycle is a connection that closes a loop, its stream
named as the tear.

## The checker

`azoth_process::validate` runs the calculus's rules:

- every instance names a palette unit op, and only its declared parameters;
- a connection joins an outlet (or feed) to an inlet (or product);
- a Port-to-Port connection joins dimension-compatible field records;
- **linearity** — a `one` port is consumed/produced exactly once, a `many` port at
  least once, and every feed/product is used exactly once;
- every loop in the instance graph passes through a declared recycle.

`azoth check --flowsheet <file>` runs it from the command line; a cargo test walks
`specs/unit_ops/` and `specs/flowsheets/` and holds every shipped file to it.

## The kernels

A kernel is a unit operation's arithmetic, a pure function of its inlets. Fourteen live in
`crates/azoth-process/src/kernels/` — splitter, mixer, separator, throttling valve, heat
exchanger, pump, heater, cooler, filter, compressor, expander, pipe, manifold and gas scrubber —
composing the calculations in `azoth-eos` rather than adding new physics. The flowsheet
**executor** — turning a flowsheet into something that runs — is tranche P12, and is not built:
this page's subject is the schema the
executor and a GUI will both read and write.
