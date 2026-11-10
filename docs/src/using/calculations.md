# Calculation sheets

The editor's **Calculations** workspace runs one registered calculation at a time. A sheet holds
the calculation's inputs and its last answer; a **calculation book** holds several sheets and is
saved as one TOML file. The flowsheet stays loaded when you switch between the two workspaces.

## Choose a calculation

The picker on the left lists every calculation the library registers, grouped by namespace
(hydraulics, heat transfer, equations of state, characterization and reactions), with a short
**Common** list at the top and a search box. Each entry carries a one-line statement of what the
calculation is for, so the list can be read rather than guessed at. A calculation whose inputs or
outputs are tables is listed but cannot be added yet; its row says so rather than disappearing.

Adding a calculation creates a sheet and selects it. A book can hold as many sheets as the case
needs, and the tabs above the sheet move between them.

## Read the sheet

A sheet is meant to be read top to bottom, before a number is typed:

- **About this calculation** states the purpose, links to the full reference page, and holds the
  equation and the source references.
- **Inputs** is one labelled line per field. The explanation is the visible label; the symbol and
  unit are its metadata; the checked limits are stated in the unit set in force.
- **When this applies** lists the conditions the calculation does not check at runtime. They are
  the caller's responsibility, and a result inside every checked limit can still be meaningless
  if they do not hold.
- **Results** shows every declared output with its meaning and unit, even before the first run.

## Fill and calculate

Every input is a field with the unit the calculation declares. The unit set in the title strip
changes what the fields show and accept; the sheet itself stores SI, exactly as a flowsheet's
streams do, so switching units never changes the stored number.

- A required input is marked with `*` until it is set.
- An optional input can be left empty, or cleared back to empty with **clear**.
- **Load worked example** fills the calculation's own published example. It is the fastest way to
  see what a sheet expects.
- **Calculate** runs the calculation. Editing an input marks the answer **stale**; nothing is
  recomputed until Calculate is pressed, because a flash or a rating is not free.

## Read the answer

Each result carries its value, unit and meaning. Below it:

- **Warnings** are conditions that do not stop the calculation but change how much the number
  should be trusted — a transitional flow regime, a result outside the validated range.
- A red sentence is the calculation's own refusal, naming the input that caused it.
- **Model check** gives the model id, its verification status, the spec and code hashes, and the
  test count. It links to the implementation notes on the calculation's manual page.

The printable **Report** shows the purpose, the inputs and outputs with their meanings, the
assumptions, references, warnings and model check for every sheet, one sheet per page. **Save**
writes the book; **New** and **Open** replace it, with a confirmation when there are unsaved
edits. Undo and redo cover the edits in the current book.

## What is not here yet

The calculations with table-shaped inputs or outputs — fitting lists and the zone and shell-and-
tube design records — have no sheet editor in this version. They remain available through the
library and the command line. The calculation sheets are screening tools: the port is incomplete,
and every shipped model is marked unverified or partially verified, so a sheet is not a design
deliverable. [Reading a calculation page](./reading-a-calculation.md) explains the reference
sections in full.
