# Exchanging a DEXPI drawing

`azoth dexpi` writes a flowsheet as a DEXPI 2.0 PFD or P&ID, reads one back as a simulation
scaffold, and draws both directions with generic `Core/Diagram` primitives. The format and
its limits are stated in [the DEXPI exchange](../architecture/dexpi.md); this page is the
walkthrough.

## A PFD from a simulation

```sh
azoth dexpi export \
  --flowsheet specs/flowsheets/demo.toml \
  --profile process \
  --out demo.dexpi.xml \
  --allow-lossy-mappings --project-recycles --with-azoth-extensions --with-drawing
```

The Process exchange runs the flowsheet first and writes each stream's mass flow (kg/h),
absolute pressure (bar) and temperature (°C). `--with-drawing` adds the PFD: one box per
step and boundary, a labelled polyline per stream, and a diamond for a projected recycle.
A run that does not converge refuses unless `--allow-unconverged` publishes its numbers with
a warning; `--topology-only` skips the run and draws no values.

## A simplified P&ID

```sh
azoth dexpi export \
  --flowsheet plant.toml \
  --profile plant \
  --out plant.dexpi.xml \
  --allow-lossy-mappings --with-azoth-extensions --with-drawing
```

The Plant exchange carries tags, nozzles, line numbers and the instrumentation layer. The
drawing adds items, off-page connectors, one segment per stream, instrument bubbles and
dashed signals. It is a **simplified** P&ID: one segment per stream, and a throttling valve
is placed inside the segment's `Items` with the pipes split around it; other inline fittings
are not written.

Both profiles use generic primitives, not ISO 10628 symbols. The file is the handoff: a CAD
tool replaces the symbols and keeps the topology, tags and values.

## Reading a P&ID

```sh
azoth dexpi import --input plant.dexpi.xml --out scaffold.toml
azoth check --flowsheet scaffold.toml
```

The import detects the profile from the file's imports, or takes `--profile plant` explicitly;
a Proteus `PlantModel` (DEXPI 1.x) is accepted on the same command. What comes back is a
**scaffold**: instances, connections, products, line numbers and instruments are there; the
feeds are undeclared and no unit carries parameters, because DEXPI has no composition.

`azoth check` is the statement of what remains: `unknown_feed` for each feed to declare and
`missing_parameter` for each unit. The drawing's arrangement is written into `[layout]`, so
the scaffold opens in the editor where the drawing placed it. `--no-layout` drops it.

The import report names what it lost: `recycle_inferred` for a loop with no tear marker,
`port_assignment_assumed` where a DEXPI port is connection-scoped, `ambiguous_step_type`
where several azoth units share one DEXPI type, and `unsupported_plant_class` for an item
kept as `dexpi.unmapped`. `--strict` turns every warning into a refusal.

## From Python

```python
from azoth import dexpi

xml, report = dexpi.export_process(
    open("specs/flowsheets/demo.toml").read(),
    "specs/unit_ops",
    allow_lossy_mappings=True,
    project_recycles=True,
    with_azoth_extensions=True,
    drawing=True,
)
toml, imported = dexpi.import_process(xml, "specs/unit_ops")
print(imported["layout_positions"])
```

`export_plant(..., drawing=True)` does the same for a P&ID, and `import_plant` /`import_proteus`
take `layout=False` to drop the arrangement. Every call returns the text and its report as
JSON, the same two documents the CLI prints.

## The Drawing workspace

The editor's **Drawing** tab is the same view as `azoth dexpi view`, in the browser. Import a
DEXPI 2.0 P&ID/PFD or a Proteus file; the workspace draws it, lists every object's coverage -
modelled, lossy, ambiguous, unmapped, topology, control, drawing-only - and says what the file
would become: instances, connections, instruments, and the check errors still owed.

- **Open as flowsheet** loads the scaffold into the editor, where `azoth check`'s diagnostics
  mark the feeds and parameters still to declare.
- **Export PFD** and **Export P&ID** write the *open simulation* back as a DEXPI file with the
  drawing; a refusal lists its diagnostics in the panel.
- **Download SVG** writes the figure as shown, for a document or a slide.

It is a viewer, not a CAD: it does not move or author a symbol. The arrangement is carried by
`[layout]`; hand the XML to CAD for symbols and revisions.

`azoth dexpi view --input FILE --json` is the headless form, and `dexpi.view(xml)` the
Python one.

## What the exchange does not promise

A file that passes the schema and these rules is not a DEXPI EV certificate and not a named-CAE
round trip; every export report says so. The drawing is azoth's own generic representation,
not a reproduction of a vendor's symbols, and an export regenerates it from `[layout]` rather
than preserving the original CAD geometry. The viewer reads a drawing; it never edits one.
