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
statement that the document is not runnable yet. When the file carries graphics, the import
also writes the drawing's arrangement into `[layout]`, so the scaffold opens where the drawing
drew it; `--no-layout` drops it.

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
   port/stream reciprocity, and finite unit-qualified quantities. Graphics are held to their
   structure, not the unit profile: a `Core/Diagram.*` or `Plant/Diagram.*` class this port
   does not render is a warning, and a graphical `Represents` that points at nothing is a
   warning rather than a broken reference.
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

**A valve is a piping component, not a tagged item.** The DEXPI Plant model places a valve
inside a segment's `Items`, with its own two `PipingNode`s, and its `Connections` chain a
pipe to the valve and a pipe from it. The writer contracts the two flowsheet connections
that reach a valve's ports into one segment, places the `GlobeValve` in `Items`, and splits
the pipes around it; a chain of valves becomes one segment with one item each. The reader
takes the segment's pipes as its edges and the segment's `Items` as instances, so the valve
comes back as `unit_ops.throttling_valve`. pyDEXPI, loading a valve as a tagged item,
rejects the file - which is how this was found.

The Plant mapping writes **10 of 31 palette entries** and reads 23.

Reading is wider than writing by design: NeqSim's reviewed Plant writer emits a narrow set,
while a real P&ID carries every valve, column and reactor class. A class this port reads but
does not write maps to its azoth unit on import, one class to one unit; a class in neither
list is kept as `dexpi.unmapped`, which the checker refuses.

Import is the same scaffold policy as Process: instances, connections and line numbers, with
every feed undeclared and no parameter invented. A loop with no marker is found by the
importer, torn at the lexicographically first edge in it, and reported as `recycle_inferred`;
an azoth export writes the real tear in `Azoth/Recycle`.

### Instrumentation and control

The flowsheet's `[[instruments]]` and `[[signals]]` (documented in
[the process schema](process-schema.md#the-control-layer)) exchange as the DEXPI 2.0
instrumentation objects. One `ProcessInstrumentationFunction` is written per instrument, with
the tag split into the model's own parts - category, modifier and number - and with the
`Location` as a reference to `Plant/Enumerations.LocationClassification.*`. The safety class
and voting architecture go in the model's `SafetyRelevanceClass` and
`VotingSystemRepresentation`; the tag, service, measured and actuated endpoints, range, safe
state and intervals travel in `Azoth/…` properties, because the model has no field for them
here. The instrument references the plant item it watches or commands through
`SensingLocation`/`ActuatingLocation`.

A signal is a `SignalConveyingFunction` under its source instrument, with `Source` and
`Target` references to the two instruments and a `SignalConveyingType` reference for the
carrier. Import reads the whole layer back, including the extensions when they are present.

**Not ported in this layer**: `ProcessSignalGeneratingFunction` and `ActuatingFunction` (the
model's sensor and final-element decomposition, which the instrument record does not split),
`InstrumentationLoopFunction` grouping (a loop is derivable from `loop_number`), alarm and
trip setpoints, and instrumentation graphics. The SIL and voting recorded are exchange
evidence, not SIL verification.

The **Proteus** reader carries the same layer: a `ProcessInstrumentationFunction` becomes an
instrument with its tag, location and the item its `ProcessSignalGeneratingFunction` is
located in, and a `SignalConveyingFunction` between two instruments becomes the signal. A
signal ending at an actuator is reported as `instrumentation_signal_skipped` rather than
invented as an instrument.

## Graphics

A drawing is **presentation, not process state**. DEXPI 2.0 carries it in the `Core/Diagram`
package: a `Diagram` holding `RepresentationGroup`s, each linked by `Represents` to the
conceptual object it draws, each holding a `Static` of `PolyLine`/`Polygon`/`Text` primitives.
`crates/azoth-dexpi/src/drawing/` reads that model onto a neutral `Drawing` - groups of
polylines, polygons, ellipses, arcs and text, with the representing object's id - and writes
the reviewed subset back.

The reader is wider than the writer. It resolves `ShapeCatalogue`/`ShapeUsage`, nested groups,
labels, connector lines and node positions, and understands the geometry a vendor file
carries. A drawing class or a shape it cannot render is a `DEXPI_GRAPHICS_UNSUPPORTED`
warning, never a refusal: a CAD file's furniture is not a defect in its process model. Both
`Core/Diagram.*` and `Plant/Diagram.*` are the graphics namespace for that reason, and an
unresolved `Represents` is a warning rather than a broken reference.

**Coordinates.** DEXPI and CAD are y-up; the editor's `[layout]` is y-down. The reader flips
once (around the diagram's `MaxY`, or the Proteus `Extent/Max`) and normalises to a zero
origin; the writer flips back. The arrangement therefore survives an exchange **up to a
translation**, which is what a layout view claims. An import with graphics places each
instance, feed and product in `[layout]`; `--no-layout` drops it.

**What an export draws.** `--with-drawing` adds a PFD to a Process export and a simplified
P&ID to a Plant export. A unit is a box with its name and palette or Plant class, a stream is
an orthogonal polyline with an arrow and its label (and, when a run wrote it, kg/h), a
boundary is a pentagon, a recycle is a diamond with dashed streams, and a P&ID adds line
numbers, instrument bubbles and dashed signals. This is deliberately **not a symbol library**:
the primitives are generic Core shapes, and a CAD tool replaces them with its own ISO symbols.
Routing is one orthogonal step, with a return lane for a backward edge; `DexpiLayoutEngine`
and its 3,572 lines are not ported, and the drawing places what the canvas places.

### The Drawing workspace

`azoth-dexpi`'s `view::inspect` answers one question about a file: **what would azoth
model?** It reads the document and its drawing, runs the same import the CLI runs, and marks
each conceptual object `modelled`, `lossy`, `ambiguous`, `unmapped`, `topology` or
`control`; a drawn group with no conceptual object behind it is `drawing_only`. The statuses
come from the same mapping tables the importers read - there is no second mapping. The answer
also carries the drawing, the scaffold TOML the file would become, and the positions the
import wrote into `[layout]`.

The editor's **Drawing workspace** is a read-only viewer over that answer: it imports a file,
draws it, lists what would be modelled and what would not, and offers "Open as flowsheet"
(the scaffold into the editor) and "Export PFD"/"Export P&ID" (the open document back out).
It is **not a CAD editor**: no symbol is moved or authored, and the figure is not saved back
as a drawing. The same call is on `/dexpi` over `azoth serve` and in the wasm module, so the
viewer needs no server; a file with no graphics gets a schematic derived from its topology,
labelled as such.

## Proteus

DEXPI 1.x serialises the same drawing as Proteus XML: a `PlantModel` root, `ComponentClass`
attributes instead of `type` references, `Connection` elements with `FromID`/`ToID`, and
`GenericAttributes`. `crates/azoth-dexpi/src/proteus.rs` **projects** that document into the
DEXPI 2.0 Plant metamodel, and the Plant importer reads the result - so the topology rules,
port assignment, line numbers and recycle inference are shared rather than written twice.

The projection itself carries no graphics; the same parsed tree feeds
`drawing/proteus.rs`, which reads the `ShapeCatalogue`, each component's `Position`/`Scale`,
the segment `CenterLine`s and the `Label` texts into the neutral drawing. Its groups use the
projected object ids, so a drawn item and the imported instance share one key. Inline fittings
this port has no unit for (a reducer, a blind flange, a strainer, an orifice plate) are joined
through in the topology with an `inline_component_skipped` diagnostic and kept as
drawing-only groups; connections ending at a component with no azoth unit are counted in a
`topology_gap` note. A valve in a segment is kept as `unit_ops.throttling_valve`.

The official DEXPI 1.3 sample is not vendored - its redistribution terms are not clear enough
to copy - so the tests use a fixture written for them, and the sample is read by hand during
development.

## Not ported

The instrumentation layer exchanges; the pieces named in that section are not ported yet.
A mixer or splitter is a junction of pipes rather than a tagged item and is still refused
on Plant write; a valve is the one piping component this writer places.
The middleware's command model has no instrument edit - a document's instruments are written
by hand or by an exchange, not by a form - and the editor's graph does not draw them.

Graphics exchange as generic primitives, but **the viewer is not a CAD editor**: it does not
move, create or author a symbol, and it does not route a stream. `DexpiLayoutEngine` (3,572
lines) is the class that would close automatic routing, and it is not ported; an export
places the canvas's arrangement with one orthogonal step. Multi-area `ProcessModel` packages,
the conformance-assessment machinery and the revision-impact reports are NeqSim's
`Dexpi20ProcessModelPackage*` family and are not ported either: one file, one exchange, until
a caller needs a package.

*Enforcement: check — `tools/gen_dexpi_mapping.py --check` refuses a palette entry with no mapping row or a row naming no entry, and the vendored schema's SHA-256; `crates/azoth-dexpi/tests/` holds the mapping, units, refusals, round trips, the drawing read/write and placement tests, and the two NeqSim oracles (Process and Plant); `crates/azoth-dexpi/tests/drawing_export.rs` draws a PFD and a P&ID from the shipped fixtures and checks every `Represents`; `process_export.rs` and `plant.rs` validate every export with `xmllint` against the pinned XSD and fail rather than skip when the tool is absent; `python/tests/test_dexpi.py` loads the process, plant and drawing exports with pyDEXPI in strict mode; and `docs/claims.toml` holds the mapping counts to the tables.*
