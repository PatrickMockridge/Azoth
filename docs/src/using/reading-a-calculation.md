# Reading a calculation page

Every calculation in the reference is generated from its specification. The order is deliberate,
and a sheet in the editor links here when the short form is not enough. The
[pump-power page](../hydraulics/pump_power.md) is a good worked example.

## What this is for

The paragraph under the title is the calculation's purpose: what it computes and the engineering
decision it supports. It is the first thing to read, because a number from the wrong calculation
is worse than no number.

## Equation and source

**Equation** shows the relation in LaTeX and in the form the library evaluates. **Source** names
the standard, paper or textbook, and — where it matters — what about the attribution is not
confirmed.

## Notes

**Notes** is implementation detail: where this port departs from its source, which branches are
not carried, and known deviations. It is audit material for judging a number, not a condition the
caller can satisfy. The browser does not show it; the sheet links here for it.

## Inputs and outputs

The tables give each field's name, unit and meaning. The editor shows the same meaning beside the
control, so a sheet and this page cannot describe a field differently.

## Checked limits

**Valid range** lists the bounds the calculation enforces. Each bound says what happens when it is
violated — the run is refused, or it returns a value with a warning — and why the bound exists.
The editor shows the same limits beside each input, in the unit set you chose, and quotes the
bound's own sentence when a value is outside it.

## Conditions for a valid result

These are conditions the calculation does **not** check at runtime: phase, steady state, the
meaning of a property, the applicability of a correlation. They are the caller's responsibility.
A result that is inside every checked limit can still be meaningless if these do not hold.

## Worked example

A fully specified case with the arithmetic written out. **Load worked example** in the editor
fills the same inputs, so the number on the sheet can be retraced by hand.

## Tests and how far this is checked

The page lists the tests that exercise the calculation and the provenance record: the verification
status, the hashes of the specification and the kernels, and whether any external validation case
names it. `partially_verified` is the normal case for this port; it means the kernels reproduce
the examples pinned in the calculation's own specification and no independent validation case is
recorded. The editor prints the status as the library spells it, with this sentence beside it, in
the sheet's **Model check** — not as a badge on the answer, which the run's own calculated/out-of-
date state describes.

## Implementation

The Python and Rust entry points, and any data files the calculation reads. A result carries the
same provenance at runtime, so a number can be traced back to the page that explains it.
