# Using azoth

The editor runs in the browser: the engine is compiled to WebAssembly, so nothing is installed
and no case is uploaded. It has two workspaces, and the one you are in never throws the other
away.

## The flowsheet

Build a process from **feeds**, **unit operations** and **products**, joined by streams. The
palette adds a unit operation; the property window fills in its design; **Solve** runs the
document. Start with [the editor](./editor.md).

## The calculation book

The **Calculations** switch opens a separate book of standalone calculations. Each sheet says what
the calculation is for, what every input means, what it will return, and when the answer applies.
A book holds several sheets, is saved in TOML, and prints as a report. Start with
[calculation sheets](./calculations.md).

## Reading the reference

Every calculation has a page in this book, generated from the same specification the code is
generated from. [Reading a calculation page](./reading-a-calculation.md) explains the sections,
from the purpose summary and input meanings through the checked limits, the conditions the
calculation does not check, and the worked example.

The sidebar's calculation reference is the long form of what a sheet links to. A sheet's
**Read the full manual** link opens the page for its calculation.

## Install and call the library

The library is a Rust crate with a Python package and a command line. The
[front page](../index.md) has the install route and the batch API; [the middleware](../architecture/middleware.md)
describes the doors a front end uses.

[Python scripting and Jupyter](./python.md) is the route for a script or a notebook:
`import azoth`, the `pint` quantities the API takes, the worked notebooks under
[`notebooks/`](https://github.com/PatrickMockridge/azoth/tree/main/notebooks), the batch API,
and the provenance every result carries.

## Before you rely on a number

The port is incomplete. Every shipped model is marked `unverified` or `partially_verified`, and
the pages state what was not carried. **These are screening values, not design deliverables.**
