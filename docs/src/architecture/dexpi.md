# DEXPI exchange

DEXPI (Data Exchange in the Process Industry) is the process industry's vendor-neutral
format for the two drawings that surround a simulation: the P&ID and the PFD/BFD.
`azoth-dexpi` writes a flowsheet, and its run, as a DEXPI 2.0 **Process** model, and reads
one back as a scaffold. This page states what the exchange carries, what it costs, and what
it refuses.

The exchange is the file half of the interoperation surface; [the middleware](middleware.md)
is the live half. Neither is a second copy of the other: the middleware drives a document
that is already loaded, while DEXPI is how a plant-engineering tool hands one over and
receives one back.

## What this does not claim

A file that passes the official schema and these rules is **not** a DEXPI EV certificate,
and it is **not** a named-CAE round trip. Every export report states
`NOT_A_DEXPI_EV_CERTIFICATE` and `QUALIFICATION_REQUIRED` for exactly that reason. Importing
into the precise CAE product and version a project uses, exporting again, and reviewing the
differences is a separate, accountable step.

## The release the files name

The Process profile imports the Core and Process information models at
`https://data.dexpi.org/models/2.0.0/Core.xml` and `.../Process.xml`; the Plant profile
imports Core and `.../Plant.xml`. DEXPI 2.0.1 publishes a byte-identical XML schema, so one
digest pins both releases:

```
sha256 0bd6568b7c52fd59ac4ef6d8547ab2b2aa890c009a3a1bfbc5d7cc18c53b36e6
       crates/azoth-dexpi/schema/DEXPI_XML_Schema.xsd
```

The schema's own header says `version="1.0.0"`; it is the 2.0 release, and the digest — not
that attribute — is its identity. The file is CC BY 4.0, unmodified, attributed in `NOTICE`;
the crate around it is AGPL-3.0-or-later.

DEXPI 1.4 — the last P&ID release serialised as Proteus XML — and the native 2.0 Plant
(P&ID) model are **not ported** in this tranche. The classes that would close them are named
at the end of this page.

## The metamodel

A DEXPI 2.0 document is a `Model` with `Import` elements and one tree of `Object`s. An
object carries a `type` — a name reference such as `Process/Process.Pumping` — and three
kinds of child: `Components` (composition), `Data` (a typed property) and `References` (ids
or names). Elements are unqualified; there is no target namespace.

`Data property` is an XSD `nameReference`, not a free string: a bare name may carry no dot,
and a dot is legal only after a `prefix/`. That is why azoth's extension properties are
spelled with a slash after the vendor prefix, not a dot.

The reader refuses a document type declaration, caps the document at 64 MiB and nesting at
128, and treats an element it does not know as an error rather than a subtree walked past.

## The mapping

The mapping is data, not a match arm. `specs/dexpi/process_units.toml` states one row per
palette entry, `tools/gen_dexpi_mapping.py` compiles it, and the generator refuses a palette
entry with no row or a row naming no entry — so adding a unit operation cannot leave the
exchange silently incomplete.

A row's disposition is what the exchange costs:

- **faithful** — the DEXPI Process type carries what the unit does;
- **lossy** — the reviewed type is the nearest counterpart and the row's `note` says what is
  lost; an export that includes it requires `--allow-lossy-mappings`;
- **refused** — no reviewed Process type exists, and the export fails naming the unit.
  Gaining a mapping is a review, not a flag; `--map <unit>=<type>` is the explicit override,
  and it is recorded in the report.

The Process mapping covers **31 palette entries**: 11 faithful, 10 lossy and 10 refused.

The reviewed Process types are NeqSim's `Dexpi20ProcessModelWriter.type` set — `Source`,
`Sink`, `TransportingFluids`, `Compressing`, `Pumping`, `Distilling`, `SeparatingByGravity`,
`Cooling`, `ExchangingThermalEnergy`, `HeatingInFurnace`, `StoringFluids`,
`RegulatingFlow`, `MixingSimple` and `SplittingMaterial`. A type outside that set has not
been reviewed against the DEXPI Process model, and this exchange does not invent one: the
Process model has no reacting, flaring, ejecting or manifold step, so reactors, the flare,
the ejector, the manifold, the component splitter and the two four-port rotating machines
are refused by name.

The refusals are not a judgement about DEXPI. The Plant model has reactor and ejector
classes, so those units gain a mapping when the Plant profile lands.

## Quantities and units

The Process profile's state is mass flow in kg/h, absolute pressure in bar, and temperature
in °C — the references NeqSim's writer emits and the official Core physical-quantity model
names. Import converts through one explicit table:

| quantity | reference | SI |
|---|---|---|
| mass flow | `...MassFlowRateUnit.KilogramPerHour` | kg/s = value / 3600 |
| pressure | `...PressureAbsoluteUnit.Bar` | Pa = value × 10⁵ |
| temperature | `...TemperatureUnit.DegreeCelsius` | K = value + 273.15 |

An unreviewed reference is refused by name. Gauge pressure is refused rather than given an
atmosphere, Fahrenheit and Rankine are refused, and a non-finite magnitude is refused —
never written as `NaN` or `INF`.

## Identity

Structural XML ids are counters — `ProcessStep1`, `Stream1`, `Port1` — so uniqueness is by
construction and the `xsd:ID` pattern is met. The azoth identity travels in the standard
`Identifier` and `Label` fields, which are free strings and preserve a name such as
`10-VA-001` exactly.

The `Model` name is a DEXPI `name` and may carry no dot, so a flowsheet id such as
`flowsheets.demo` is sanitised for it; with `--with-azoth-extensions` the exact id and
display name travel on the engineering model.

Extension properties are how azoth round-trips what the standard model cannot carry. They
are schema-valid and a consumer that does not know them ignores them:

| property | on | carries |
|---|---|---|
| `Azoth/UnitOp` | process step | the palette id, recovering the many-to-one mappings |
| `Azoth/InstanceId` | process step | the azoth instance id |
| `Azoth/Port` | material port | the palette port name |
| `Azoth/Index` | material port | the position of a `many` port's stream |
| `Azoth/From`, `Azoth/To` | stream | the azoth endpoints |
| `Azoth/Recycle` | process step | the tear name of a projected recycle |
| `Azoth/FlowsheetId`, `Azoth/FlowsheetName` | engineering model | the document's identity |

## Export

The flowsheet is checked and run before it is written: an error diagnostic refuses the
export, because a DEXPI file is a handover artifact and a document that does not validate is
not one to hand over. A run that does not converge is refused unless `--allow-unconverged`
is given, in which case the numbers are published with a warning. A stream the run did not
produce carries no quantity; nothing is defaulted.

One process step is written per boundary and per instance, one `MaterialPort` per stream
endpoint, and one `Process.Stream` per material edge — `connections` and `recycles` alike.
DEXPI has no recycle: with `--project-recycles` a tear becomes a `TransportingFluids` step
with two streams, which is NeqSim's own projection; without it the export refuses.

## Import

DEXPI carries no composition, so mass flow cannot become molar flow and a feed record cannot
be written. An import therefore produces a **scaffold**: instances, connections, products and
recycles, every feed left undeclared, and no unit carrying parameters. `azoth check` on the
result reports `unknown_feed` and `missing_parameter` for each site; that verdict is the
statement that the document is not runnable yet.

Port names are the lossy part, because a DEXPI material port is connection-scoped rather
than named the way a palette port is. The extension properties recover them exactly. Without
them the importer assigns by direction and document order against the palette and reports
`port_assignment_assumed` — a heat exchanger's hot and cold sides are the dangerous case, and
a strict import refuses rather than guesses.

Where several azoth units share one Process type — `Distilling` covers every column,
`SeparatingByGravity` covers the separator family — the importer chooses a documented default
and reports `ambiguous_step_type`; the extension property removes the ambiguity. A type with
no reviewed azoth unit is kept as `dexpi.unmapped`, which the checker refuses, rather than
substituted with a plausible one.

## Validation

Three layers, three different claims:

1. **Semantic rules** (`azoth dexpi validate`): the exact imports, unique ids, resolvable
   references, the reviewed type set, one Process model with at least one step and stream,
   port/stream reciprocity, and finite unit-qualified quantities.
2. **The official schema** (`--xsd`): `xmllint --nonet --noout --schema` against the pinned
   XSD. A missing `xmllint` is a failure, never a silent pass — a check that did not run must
   not be reported as one that did.
3. **An independent parser**: the Python test loads every export with pyDEXPI in strict mode.

**Import has no upstream oracle.** NeqSim has no native DEXPI 2.0 Process reader, so the
writer is the oracle for export and the reader is held to the DEXPI model, the semantic
rules, and export→import→export round-trip. The export oracle is
`validation/neqsim/Dexpi20ProcessProbe.java` and its committed capture.

## The Plant model

A Plant exchange is a P&ID. It carries no process state - no mass flow, pressure or
temperature - and what it does carry is the part a PFD loses: the tag of every item, the
nozzle each stream leaves by, and the line it belongs to.

One tagged plant item is written per mapped instance, with one `Nozzle` per palette port
(each owning a `PipingNode`). Each material edge becomes a `PipingNetworkSegment` with a
`Pipe` joining two nodes. A flowsheet that states `line` on a connection is grouped into one
`PipingNetworkSystem` per line, whose `LineNumber` is that line; a connection without one
goes under `AzothProcessConnections`, and import turns that placeholder back into no line.
Feeds and products are `FlowIn`/`FlowOut` off-page connectors inside the segment that uses
them.

**A valve is not a tagged item.** The DEXPI Plant model places a valve as a piping
component inside a segment, and a tagged item is a different class. This writer places no
piping components yet, so `unit_ops.throttling_valve` is refused on write and read-only;
placing valves needs the segment's `Items` component and is the next step. pyDEXPI, loading
a valve as a tagged item, rejects the file - which is how this was found.

The Plant mapping writes **9 of 31 palette entries** and reads 23.

Reading is wider than writing by design: NeqSim's reviewed Plant writer emits a narrow set,
while a real P&ID carries every valve, column and reactor class. A class this port reads but
does not write maps to its azoth unit on import, one class to one unit; a class in neither
list is kept as `dexpi.unmapped`, which the checker refuses.

Import is the same scaffold policy as Process: instances, connections and line numbers, with
every feed undeclared and no parameter invented. A loop with no marker is found by the
importer, torn at the lexicographically first edge in it, and reported as `recycle_inferred`;
an azoth export writes the real tear in `Azoth/Recycle`.

## Proteus

DEXPI 1.x serialises the same drawing as Proteus XML: a `PlantModel` root, `ComponentClass`
attributes instead of `type` references, `Connection` elements with `FromID`/`ToID`, and
`GenericAttributes`. `crates/azoth-dexpi/src/proteus.rs` **projects** that document into the
DEXPI 2.0 Plant metamodel, and the Plant importer reads the result - so the topology rules,
port assignment, line numbers and recycle inference are shared rather than written twice.

What the projection drops is named: positions and presentation are graphics; inline fittings
this port has no unit for (a reducer, a blind flange, a strainer, an orifice plate) are
joined through with an `inline_component_skipped` diagnostic; connections ending at a
component with no azoth unit are counted in a `topology_gap` note. A valve in a segment is
kept as `unit_ops.throttling_valve`.

The official DEXPI 1.3 sample is not vendored - its redistribution terms are not clear enough
to copy - so the tests use a fixture written for them, and the sample is read by hand during
development.

## Not ported

**A valve cannot be written yet.** Placing a `GlobeValve` in a segment's `Items` component
and splitting the segment's pipe around its nodes is the step that closes it; until then
`unit_ops.throttling_valve` is refused on Plant write.

Instrumentation, control loops and SIS/HIPPS semantics need a document concept azoth does
not have - there is no instrument model - and are a tranche of their own, not a serialiser.
The Plant profile reads the instrumentation classes so a real P&ID validates, and reports
the objects as unplaced. Graphics need the diagram representation groups;
`DexpiLayoutEngine` (3,572 lines) is the class that would close automatic routing, and it is
not ported. Multi-area `ProcessModel` packages, the conformance-assessment machinery and the
revision-impact reports are NeqSim's `Dexpi20ProcessModelPackage*` family and are not ported
either: one file, one exchange, until a caller needs a package.

*Enforcement: check — `tools/gen_dexpi_mapping.py --check` refuses a palette entry with no mapping row or a row naming no entry, and the vendored schema's SHA-256; `crates/azoth-dexpi/tests/` holds the mapping, units, refusals, round trips and the two NeqSim oracles (Process and Plant); `crates/azoth-dexpi/tests/process_export.rs` and `plant.rs` validate every export with `xmllint` against the pinned XSD and fail rather than skip when the tool is absent; `python/tests/test_dexpi.py` loads both exports with pyDEXPI in strict mode; and `docs/claims.toml` holds the mapping counts to the tables.*
