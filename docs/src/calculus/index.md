# The calculus of thermodynamic dimensionality

**This section is normative for the types.** `specification.md` is normative for the port —
what is ported, how a calculation is declared, where responsibility for data sits.
Where the two meet, `specification.md` decides the policy and this section decides the
types, because a type is not a policy and the two are not two answers to one
question.

## What it is

azoth's units are not a convention. `mm` and `m` differ by a factor; `Pa` and
`Pa*s` are different dimensions; `W` and `J/s` are the same one. Each of those is
a statement about a group, not about a table of strings, and the group is the
first thing here.

A unit operation is not a tuple of arguments either. It is a process with an
inlet and an outlet, and "the mole balance closes" is either a consequence of how
its channels are used or it is a sentence in a specification that nothing checks.
That is the second thing.

Five layers, each a Lean module and a page here. **Two of the modules exist** —
`Azoth/Dim.lean` and `Azoth/Vocabulary.lean` — and the other three are the ones
their tranches will write:

| Layer | What it fixes | Lean | Page |
|---|---|---|---|
| Dimensions | what a unit's dimension is, and when two are equal | `Azoth/Dim.lean` ✅ | [Dimensions](./dimensions.md) |
| Vocabulary | which units a spec may declare, and what each one is | `Azoth/Vocabulary.lean` ✅ | [The vocabulary table](./vocabulary.md) |
| Barbs | what an observer of a channel can see, and therefore what equality means | `Azoth/Barb.lean` | [Barbs](./barbs.md) |
| Processes | a unit operation as a process on typed, directional channels | `Azoth/Process.lean` | [Processes](./process.md) |
| Reflection | feedback, serialisation, and the interoperation surface | `Azoth/Rho.lean` | [Reflection](./rho.md) |
| Capability | the keycard as authority a process holds rather than a global it reads | `Azoth/Capability.lean` | [The keycard](./capability.md) |

The two that exist are the ones written against a **vocabulary** — the canonical
unit strings a spec may declare. That vocabulary is data rather than proof: one
hand-written table, compiled into Rust, Python, JSON Schema and Lean. It is
specified separately, in [The vocabulary table](./vocabulary.md), because it is the
one part of this whose authority is a generator rather than a proof.

## Status: the vocabulary is proved, the rest is specification

**The vocabulary is proved.** The development is in `lean/`, and
`tools/check_lean_axioms.py` refuses a gap in any theorem it claims — twenty-eight
today: four about the dimension group, that its exponent vectors name it and read
back unchanged, and one per canonical unit checking that the table's exponents name
the dimension `lean-units` calls by that name. That is
[The vocabulary table](./vocabulary.md), and it is the whole of what this section
has checked against anything.

**Everything else here is a specification.** [Processes](./process.md),
[Barbs](./barbs.md), [Reflection](./rho.md) and [The keycard](./capability.md)
state the layers that tranches F4 and F7 will build, and **none of those layers
exists**: there is no channel type in Rust or Python, no `Kind`, no process
calculus, no card. A proof about them today would be a proof about nothing — it
could not be falsified by any test this repository can run, which is the only thing
a proof in this tree is for.

So each claim names its theorem and says which of three things it is:

| | Means |
|---|---|
| **Proved** | in `lean/Azoth/`, and the axiom gate covers it |
| **Specified** | the layer does not exist. The claim is what will make the tranche that builds it checkable, and it is proved then |
| **Characterised** | the layer exists and the claim describes it rather than guarding it. Nothing this repository can do would violate it, so it is stated and deliberately not proved |

The distinction between the last two is the whole of the judgement here, and it is
the one `specification.md` makes when it refuses a status field: what a claim is *worth*
depends on whether anything can falsify it. A specification is worth writing
before its layer exists, because it is what the layer is checked against. A
characterisation is worth writing too — it is what makes the dimension the right
notion rather than an arbitrary labelling — but proving it would add a theorem
nobody consults to a gate whose value comes from every entry in it being one a
change could break.

## The rule the rest depends on

A unit declares its **dimension** and nothing else. It does not declare how many
SI base units it is worth, because `uom` and `pint` each already know that, and a
number written down in a third place is a number that can disagree with both.
The table says what a unit *is*; the two units libraries say what it is *worth*;
and a test compares them rather than trusting a restatement. See
[The vocabulary table](./vocabulary.md).
