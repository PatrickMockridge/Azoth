# chemeng

Open, validated, citable chemical engineering calculations.

Every calculation in this book ships with its equation, where the equation came
from, the range in which it is validated, its assumptions, a worked example, and
the tests that exercise it. Nothing here is prose written alongside code: the
pages under [Hydraulics](./hydraulics/index.md) are **generated from the same
specification files the code is generated from**, so they cannot drift from it.
CI regenerates them and fails on any difference.

## The two ideas this library is built around

**Warnings are not errors.** A value outside the range in which a correlation was
validated is still a value. Refusing to return it would be less useful than
returning it with a warning. What the library must never do is return it
*silently*.

**A check that could not run is not a check that passed.** When an optional input
is missing, the range check that depends on it reports `RANGE_CHECK_SKIPPED`
rather than quietly succeeding. "Checked and fine" and "never checked" are
different states, and the API keeps them different.

Both are visible in every result: each carries its warnings, and each can be
asked whether it is clean.

## Reading a calculation page

Each page has the same shape, and the order is deliberate:

| Section | Why it is there |
|---|---|
| Equation | In LaTeX for a reader, and in the form the library evaluates, so you can check one against the other |
| Source | Where the equation came from. Marked `TODO: source needed` when it is not confirmed |
| Verification | Whether a person has checked the citation. **Read this before trusting a number** |
| Inputs and outputs | Names, units, and what they mean |
| Valid range | The bounds, what happens when one is violated, and why each exists |
| Assumptions | What is **not** checked at runtime, and is therefore your responsibility |
| Worked example | A fully specified case you can reproduce by hand |
| Tests | What is actually exercised, including what is deliberately skipped and why |

## What is implemented

The current slice is the hydraulics kernel through Darcy-Weisbach pressure drop:

- [`hydraulics.reynolds_number`](./hydraulics/reynolds_number.md)
- [`hydraulics.friction_factor_colebrook`](./hydraulics/friction_factor_colebrook.md)
- [`hydraulics.friction_factor_swamee_jain`](./hydraulics/friction_factor_swamee_jain.md)
- [`hydraulics.crane_k_factors`](./hydraulics/crane_k_factors.md)
- [`hydraulics.darcy_weisbach`](./hydraulics/darcy_weisbach.md)

Orifice, control valve, relief valve and pump calculations are not implemented.

## How the pieces fit together

Each calculation is independent, and the composition is done by the caller. The
`chemeng pipe` command performs the one shown here, and reports the two pressure
drop contributions separately rather than only their sum - they come from
different methods, and seeing which one dominates is part of judging the answer.

```mermaid
graph LR
    F["fluid properties<br/>density, viscosity"] --> V["velocity<br/>from flow and bore"]
    V --> RE["reynolds_number<br/>Re and regime"]
    V --> K
    R["roughness, diameter"] --> RR["relative roughness"]
    RR --> FF
    RE --> FF["friction factor<br/>colebrook or swamee-jain"]
    FF --> DW["darcy_weisbach<br/>straight pipe"]
    FF --> K["crane_k_factors<br/>fittings"]
    DW --> T["total pressure drop"]
    K --> T
```

Note the two arrows leaving the friction factor. The straight-pipe loss uses
`f * L/D`, and the fitting loss uses each fitting's equivalent length ratio
scaled by the *same* `f`. Feeding the fitting calculation with anything other
than the friction factor actually used for the pipe would make the two
contributions inconsistent with each other.

## Not for design work yet

The fitting coefficients in `data/fittings/crane_k_factors.csv` are **estimated
dummy values** - plausible magnitudes chosen so the software has something to run
against. They are not from Crane TP-410 or any other standard, and a pressure
drop computed from them can be wrong by a factor of two while looking entirely
reasonable. Every affected result carries an `ESTIMATED_DATA` warning, and the
[`crane_k_factors`](./hydraulics/crane_k_factors.md) page says so at the top.

Water and air properties are a different case: real published values, marked
`unverified` because they have not been checked against a primary formulation.

## How a calculation is added

1. Write a spec under `specs/calcs/<namespace>/<name>.yaml`.
2. Write one function in Python and one in Rust.
3. Declare the tests in the spec.

The documentation, the range checks both implementations enforce, and the test
cases both implementations run are all generated from that one file. If you can
write the spec, you have written the calc, the tests and the docs.
