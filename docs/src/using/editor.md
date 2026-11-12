# The editor

The flowsheet editor is a process flowsheet in a browser: a palette of unit operations, a
drawing, a window for the selected object, and a workbook of every stream. The calculations run
in the page — the engine is compiled to WebAssembly — so nothing is installed and no case is
uploaded.

**[Open the editor →](https://patrickmockridge.github.io/azoth/)**

This page is part of the user manual; [Using azoth](./index.md) is its start.

![The editor showing the demo flowsheet](../images/editor.png)

## First five minutes

![The welcome card over the demo flowsheet](../images/welcome.png)

1. The demo opens on its own: two feeds meeting at a mixer, then a pump, a heater, a separator and a liquid recycle. Each mixer feed has its own inlet point.
2. Press **Solve**. The status changes from `stale` to `solved`, and each connection shows the
   temperature and pressure the run reached.
3. Click a unit operation to open its Inspector. **Design** holds the editable parameters, and
   double-clicking the node opens the same Design fields in the node window itself. **Conditions**
   holds the feed or product stream, **Worksheet** the streams on its ports, **Results** whatever
   the model reached besides its streams, and **Provenance** the model, hashes and data authority
   behind it.
4. Open **Workbook** at the bottom for every stream in one table. **Messages** holds the model
   check's findings and any run error; clicking an entry selects the object it is about.
5. **Case → Demo cases** opens other worked cases.

## The drawing

- Drag from an outlet handle to an inlet handle to connect two objects. Every connection point
  says `IN` or `OUT`; a two-stream exchanger labels `hot` and `cold` on each side. **One line uses
  one point**: a mixer has a separate inlet for every feed, and a splitter a separate outlet for
  every product. Connectors are drawn as right-angle routes with rounded corners.
- Double-click a node for its window: a unit operation's Design fields are editable in place,
  followed by its status, streams, results and messages. The button in the window hands the other
  sheets to the Inspector.
- Click an object to select it; drag a unit operation to move it.
- Delete or Backspace, or **Remove** in the object's window, removes the selected object. A removed
  feed takes its connections with it; **Undo** restores it.
- A connection can be turned into a **tear** — a recycle — from its Connections sheet. A tear is
  drawn as a moving line and is iterated until its residual is small.
- Connectors are routed around the node plates; a recycle takes its own return lane below the
  drawing, so it never runs back over a unit or along its own outbound line.

## Feeds, products and units

A **feed** is four numbers and a fluid: molar flow, pressure, temperature and composition. Use
the component picker rather than spelling substance names; the databank is what decides whether a
name exists.

A **product** is a name only: its stream is calculated.

The unit set at the top right changes what is *shown* — SI or field units, including °C and °F.
It never changes the document, and solving has no opinion about it.

## Appearance

The top bar's **Appearance** control chooses a layout and a colour, and remembers the pair for
this browser. The layouts are **Macintosh**, **Windows XP** and **Cyberpunk**; each can be
**light**, **neutral** or **dark**. The default is Windows XP light.

The layouts are not colour schemes over one screen. **Windows XP** is a single MDI frame: a Luna
title bar, a menu row, a toolbar and Explorer panes around a document window. **Macintosh** is a
grey desktop: the menu bar is at the top and the navigator, drawing, inspector, workbook and
status are separate windows — drag a title bar to move one, drag its bottom-right corner to
resize it, and use **View → Reset window positions** to put them back. **Cyberpunk** is a
full-bleed drawing with translucent HUD rails; the switches at the right of the top strip hide
the navigator, inspector, workbook or command deck, and **Commands** in the strip — or
**Ctrl/Cmd+K** anywhere in that layout — opens a command palette over the same commands as the
toolbar and menus. The workbook rail's top edge resizes it.

![The Macintosh desktop](../images/editor-mac.png)

![The Cyberpunk HUD in its dark colour](../images/editor-cyberpunk.png)

The profiles are cosmetic: they change nothing in the document, the stream values or the unit
set, and the report keeps its fixed paper palette when printed. The choice is a reader's, like
the unit set, and is not written into a case.

## What the numbers mean

A result has two readings, and the split is deliberate.

**The human reading** is the **Results** sheet and the report's cover: a plain label, the value
and unit in the unit set in force, and the library's own sentence about how the number is
obtained — often the equation behind it, so an engineer can check the answer by hand. Only fields
the run published are drawn. `false` and `0` are values; a field the design did not configure is
not a dash row, and a parameter that was never set is not in the report.

**The machine record** is the collapsed *Machine-readable fields* block on the same sheet, the
**Provenance** tab, and the report's appendix. It keeps the model's own field names in the
library's units, the parameters the run was given, the model id, the NeqSim class it was ported
from, the data authority (the standard keycard this build ships, derived from NeqSim values) and
the spec and kernel hashes. **Download run record (JSON)** on the report writes the same record —
the document, the raw session, the provenance and the diagnostics — so another tool or agent can
reproduce the run and compare it field by field.

**Verification is a record field, not a badge on the answer.** `verified`,
`partially_verified`, `source_needed` and `unverified` describe how far the *model* is checked
against evidence outside its own examples. `partially_verified` means the kernels reproduce the
examples pinned in the model's own specification and no independent validation case is recorded;
it is the normal state for this port. The status is printed as the library spells it, once, in
Provenance and the report's appendix, with the sentence that defines it. **Caveats** on a result
sheet are the calculation's own warnings about a number it did return — a bound that could not be
checked, a correlation at the edge of its range.

The port is incomplete, and what is missing is recorded in the [road map](../agentic/roadmap.md)
and at the front of this book. **These are screening values, not design deliverables.**

## Saving and printing

**Save** writes the case under its name as a TOML document; **Open** reads one back. The header's
case-name field sets that name. The workbook exports to CSV,
and **Report** prints or saves a case file with the drawing, the streams, the results, the
caveats and the source record. The report also downloads the run as JSON, which is the
machine-readable form of the same run.

The **Calculations** switch opens a separate workspace for the library's standalone calculations;
the **View** menu carries the same two commands. See [calculation sheets](./calculations.md) for
the workspace and [reading a calculation page](./reading-a-calculation.md) for the reference page
a sheet links to.
