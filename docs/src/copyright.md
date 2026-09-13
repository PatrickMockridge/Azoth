# Copyright, and licensed data

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

Some of what these calculations need is not copyrightable, and some of it is. This
page records every place that distinction changed what this library does, so that a
reader can see what is missing and why rather than discovering it one spec at a time —
and so that the next calculation has a worked answer instead of a puzzle.

## The rule

**Never reproduce copyrighted tables or text.** Cite the equation; work your own
arithmetic.

*Facts* are not copyrightable. That a 90° standard elbow has an equivalent length of
thirty diameters is a fact, and recording a fact with a citation is fine. A *table* is
copyrightable — transcribing Crane TP-410's tables into this repository would
redistribute the standard's expression of those facts, and that is not something this
project does. The same line runs through every entry below.

## What is blocked, and what each block costs

| What | Why it cannot ship | What it costs | How it is handled |
|---|---|---|---|
| **Crane TP-410 fitting coefficients** (the L_eq/D table) | The table is Crane's expression of the facts | `crane_k_factors`, and the fitting share of the `azoth pipe` CLI | `data/fittings/crane_k_factors.csv` ships **estimated dummy values**, recorded in its `verify_status` column, and a test fails the day the file is populated |
| **ISO 5167 discharge-coefficient equation** | A long fitted expression whose constants come from a table of experimental results | `orifice_flow` cannot compute `Cd` | `Cd` is a **caller input**. The spec says so, names the convention required, and argues why |
| **Crane TP-410 equation numbers** | The numbering cannot be confirmed without the standard | `reynolds_number`, `darcy_weisbach` and `crane_k_factors` cite Crane | The specs leave `source.equation` unstated rather than guessing a number, and say so in their `notes` — a gap that is visible rather than a number that was made up |
| **Crane TP-410 Example 3-5** | Reproducing a worked example from the standard is the thing the rule forbids | A `reference` test in `darcy_weisbach` | `status: skipped` with the reason recorded in `skip_reason`; the derived worked example covers the same arithmetic |
| **Perry's 8th ed. Eq. 6-42** | Same | A `reference` test in `darcy_weisbach` | Skipped, same way |
| **API 520 relief-valve sizing constants** | Fitted constants and de-rating factors | `relief_valve_area` — **not implemented** | Will take the constants as inputs, as `orifice_flow` takes `Cd`. Not built yet |
| **IEC 60534 control-valve Cv/Kv tables** | Copyrighted tables | Nothing is computed from them | `control_valve_cv` takes `Cv` as an input and is responsible for undoing the unit convention its `dimensionless` declaration hides — see below |

### And one thing that is *not* blocked, but is often mistaken for it

The water and air property tables in `data/fluids/` are real published values, and no
copyright question arises: they are facts, published widely, and this repository
reproduces no table — it records a handful of points with a citation. They are marked
`unverified`, which is **not** a statement that the numbers are doubtful: it says that
no person here has checked them against a primary formulation. That is the same
discipline the fitting registry uses, in a milder tier — `estimated_dummy` means "a
placeholder" while `unverified` means "real, unchecked".

The column is derived when the data file is generated, from whether the row's citation
says the value is a placeholder. There is deliberately no `verified` value: a tool
cannot check whether a person checked something, so offering the word would only invite
an assertion nobody can test. See [Specification, S6](./spec.md#s6-provenance-is-the-engineers-job-not-the-librarys).

## The three patterns that come out of this

Everything above resolves into one of three shapes. Knowing which one a blocked value
belongs to is the whole answer to "what do I do about it".

The first two are for values this library may **not** ship. The third is for values and
code it **may**, which is a different question and is easy to conflate with them.

### A single number becomes a function argument

A discharge coefficient, a pump efficiency, a friction factor and a Crane `f_T` are all
single values. They are **arguments**, and no file is involved:

```python
azoth.hydraulics.orifice_flow(..., Cd=0.61)         # your discharge coefficient
azoth.hydraulics.control_valve_cv(10.0, ...)        # your valve Cv
azoth.hydraulics.pump_power(..., eta=0.75)          # your efficiency, from the curve
azoth.hydraulics.darcy_weisbach(0.02, ...)          # your friction factor
azoth.hydraulics.crane_k_factors([...], f_t=0.018)
```

This is worth stating plainly, because it means **most of the blocked values need no
mechanism at all**. Five of the seven entries above are already fully usable with
licensed data, today, in both languages: you hold the standard, you read the number,
you pass it in. The library never needs to know it.

**One of those five is harder than the others, and `control_valve_cv` is where it
shows.** `Cd`, `eta` and `f` are genuinely dimensionless — ratios of one thing to
another. `Cv` is not: it is defined as so many gallons per minute at one pound per
square inch, so it carries the units `gpm/sqrt(psi)`. The schema has to declare it
`dimensionless` because no unit in the vocabulary can express a fractional power of a
non-SI unit, and that declaration is a fiction of necessity. Undoing it correctly is
then the whole job of the calculation: the conversion constant is named, derived from
the definitions of the gallon and the pound-force, and tested on both sides. A caller
who supplies a metric `Kv` instead is out by a factor of about 1.156, which is large,
silent, and produces an entirely ordinary-looking flow.

It is also not a workaround. Passing the coefficient is the better design on its own
terms — the same shape as `f` in Darcy-Weisbach and `f_t` in the equivalent-length
method — and it is why those specs do not have a licensing problem in the first place.

### A set of numbers needs data with provenance

Some things are not single values: they are tables you look up *by name*. A fitting
registry is the one that exists today; a control-valve table and a relief-valve
`K`-factor table would be the same shape when those calculations are built.

Those need a file. The repository ships placeholders, and a user supplies the real
values — see below.

### A licensed implementation may be ported, with attribution

The rule at the top of this page is about material that **cannot** be redistributed.
Where a licence grants the right, the answer is the opposite, and the two cases are
easy to mistake for each other because both involve taking someone else's work.

[NeqSim](https://github.com/equinor/neqsim) is Apache-2.0. That licence explicitly
permits redistribution and derivation, including into this repository's AGPL-3.0 code —
Apache-2.0 is one-way compatible with the GNU licences, and the obligation it imposes is
**attribution**, not copyleft reciprocity. So a NeqSim algorithm is not a Crane table.
The Crane table cannot be shipped at any price; a NeqSim implementation can be shipped
provided the credit and the statement of changes are made, which is what
[`NOTICE`](https://github.com/PatrickMockridge/Azoth/blob/main/NOTICE) at the repository
root is for — once, for everything reused, rather than restated per calculation.

**The licence is not the hard part. The verification is.** A port arrives looking
authoritative — it came from a maintained industrial library — and that is a stronger
pull than a table would be, because a table at least looks like data someone has to
check. NeqSim's `CriticalPointFlash` is the case that made this concrete: it implements
Heidemann & Khalil correctly as far as can be told by reading it, and it checks the
result **nowhere** — no pure-component test, no mixture test, and a silent `break` on a
`NaN`. Agreeing with it is not evidence of anything, because there is nothing there for
it to be right against.

So a port is governed by these, and they are in
[CONTRIBUTING.md](https://github.com/PatrickMockridge/Azoth/blob/main/CONTRIBUTING.md):

- It is not evidence. Reading someone's Java is not reading the paper the method came
  from, and a port is a second implementation of a method — never a second source for
  it. Nothing here records a status saying otherwise.
- It cites three things: the paper for the **method**, the implementation for the
  **port**, and its own `notes` for the **changes**.
- It is accepted on this library's tests, never on its provenance — and those tests have
  to be ones that would fail if the port were wrong.

**A databank is a third thing again, and the distinction is worth keeping.** A
compilation can carry its own copyright where the individual facts in it do not, so
"Apache-2.0" settles the code and does not by itself settle a data table. Where a
licensed table is vendored here it is recorded separately, with its own provenance, and
one risk is named rather than hidden: a permissive licence on someone's *compilation*
does not establish that every value inside it was cleanly sourced upstream. That is not
inspectable from here, and the mitigation is attribution rather than inspection.

## Supplying your own data

### The file

Please note: the keycard is reserved for **sets** of values that a calculation looks up
by name — a fluid's property table, a fitting registry, that sort of thing. A *single*
coefficient needs no file at all and is better passed as an argument: see the first of
the three patterns above.

`keycard.example.yaml` at the repository root is the template. Copy it to
`keycard.yaml`, fill in the values you are entitled to use, and check it:

```bash
python tools/check_user_data.py keycard.yaml
```

It has a `fittings:` section for the equivalent-length registry and a `fluids:` section
for property tables. A row **may** carry a `citation`; nothing requires one, and there
is no status field.

That is the deliberate position rather than an omission. The repository's own data
records a derived `verify_status` so that it can warn you when a result rests on its
placeholders — but a value *you* supply is yours, and a form asking you to assert
something no tool can check would teach you to fill it in rather than to know the
answer. **Provenance is your job, not this library's.**

`tools/check_user_data.py` enforces those rules and imports them from
`tools/spec_lint.py` rather than restating them, so a user's file and the repository's
files cannot drift apart in what "a valid citation" means. It checks **provenance, not
values**: nothing can tell whether a coefficient is right, only whether the row says
where it came from in a form a tool can fetch. That division is deliberate and is the
same one the specs use.

The template validates as it stands, and a test enforces that — a template that does
not pass its own checker teaches the wrong shape. Every row in it is therefore
`estimated_dummy`; the shape of a finished row is shown in comments, because an active
row claiming `verified` would be a citation nobody had checked.

### Do not commit what you generate

**If the values came from a standard you licensed, committing them redistributes
them** — which is the thing this whole page exists to avoid. `keycard.yaml` is
gitignored. The data files generated from it must not be committed either; the
repository's committed copies must stay the placeholders.

### Compiling it

```bash
python tools/check_user_data.py keycard.yaml     # validate, writes nothing
python tools/gen_user_data.py   keycard.yaml     # then compile
maturin develop                                     # and rebuild
```

```
keycard.yaml                    you write this, and keep it
       │
       │  tools/gen_user_data.py   validates with the checker's own rules, then writes
       ▼
data/fittings/crane_k_factors.csv  ─┐
data/fluids/<fluid>.csv            ─┼─►  embedded by Rust (include_str!)
                                   │     read at runtime by Python
       both languages, same bytes ─┘
```

The reason it generates files rather than being read directly is the Rust core. It
embeds its data at compile time with `include_str!`, which is what makes the two
implementations read byte-identical data — a property the cross-language agreement
tests rest on. A data file read at runtime from a path the user controls would end that
guarantee unless both sides read it, and giving Rust a runtime data path means
depending on a YAML parser there, which this project has deliberately avoided.

So the user's file is a **source**, compiled once into the canonical form — exactly the
relationship `specs/calcs/*.yaml` already has to the two languages. The cost is that
changing your data means re-running the tool and rebuilding, which suits a build from
source rather than a `pip install` of a wheel.

Three things about the tool are deliberate and worth knowing before you run it:

**It overwrites `data/` in place.** That is forced: `include_str!` needs the file at the
path the code names. The shipped placeholders are committed, so `git checkout -- data/`
brings them back — but nothing structural stops you committing licensed values by
accident. It prints a warning, and the file's own banner says so too, and the rest is
yours.

**A number keeps the text you wrote.** `0.0000200` stays `0.0000200` and does not become
`2e-05`. The same float, and not the same statement: the first says three significant
figures. A tool that reformatted your data would be discarding information you put
there on purpose.

**A fluid nothing can read is refused, not written.** Both implementations hardcode the
fluid list — Rust through `include_str!` and a `match`, Python through a `_BUILTINS`
dict — so `data/fluids/<name>.csv` for an unregistered `<name>` would be a file that
loads nothing and warns nobody. Generating it would produce exactly the kind of thing
this project is built to prevent: data that looks like it is in use. The tool refuses
and prints the five places to register the name; `example_fluid` in the template is
deliberately one of these, so the template shows the shape without pretending the
repository carries a third table.

**And if you do generate licensed data, two tests will fail** — the ones asserting the
shipped coefficients are placeholders. That is what they are for. Read this page before
changing them.

## Why not simply include the numbers

Because a validated calculation library is only worth what its validation is worth, and
a table transcribed without a licence is a legal problem attached to a technical claim.
This library's whole argument is that a number should arrive with its provenance
attached and its limits stated. Reproducing a standard's tables would undercut that in
the one way the project could not recover from: the numbers would look authoritative,
be unverifiable, and be somebody else's to license.

The alternative on offer is more work for the user and less convenience for everyone —
and it keeps the arithmetic, the range checks, the assumptions and the tests in the
open, which is the part that is worth sharing.
