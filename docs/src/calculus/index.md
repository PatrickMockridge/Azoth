# The calculus of thermodynamic dimensionality

**This section is normative for the types.** `spec.md` is normative for the port —
what is ported, how a calculation is declared, where responsibility for data sits.
Where the two meet, `spec.md` decides the policy and this section decides the
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

Five layers, each a Lean module and a page here:

| Layer | What it fixes | Lean | Page |
|---|---|---|---|
| Dimensions | what a unit's dimension is, and when two are equal | `Azoth/Dim.lean` | [Dimensions](./dimensions.md) |
| Barbs | what an observer of a channel can see, and therefore what equality means | `Azoth/Barb.lean` | [Barbs](./barbs.md) |
| Processes | a unit operation as a process on typed, directional channels | `Azoth/Process.lean` | [Processes](./process.md) |
| Reflection | feedback, serialisation, and the interoperation surface | `Azoth/Rho.lean` | [Reflection](./rho.md) |
| Capability | the keycard as authority a process holds rather than a global it reads | `Azoth/Capability.lean` | [The keycard](./capability.md) |

The first layer is written against a **vocabulary** — the canonical unit strings a
spec may declare. That vocabulary is data rather than proof: one hand-written
table, compiled into Rust, Python and JSON Schema. It is specified separately, in
[The vocabulary table](./vocabulary.md), because it is the one part of this whose
authority is a generator rather than a proof.

## Status: written, not yet proved

Every page below states claims and names the theorem each will become. **The Lean
development is not in this tree yet**, so nothing here has been checked by a
machine. A reader should take each *Claim* as a specification of what has to be
proved, and not as a result — and the reason to say so plainly is the reason
`spec.md` gives for having no status field anywhere: a form that asks for
confidence manufactures it rather than producing it.

Each claim carries its own status line recording whether it is proved, and the
line is changed when the proof lands rather than when the proof is planned. The
one exception is [The vocabulary table](./vocabulary.md), which describes a
mechanism that is implemented and tested in this tree today.

## The rule the rest depends on

A unit declares its **dimension** and nothing else. It does not declare how many
SI base units it is worth, because `uom` and `pint` each already know that, and a
number written down in a third place is a number that can disagree with both.
The table says what a unit *is*; the two units libraries say what it is *worth*;
and a test compares them rather than trusting a restatement. See
[The vocabulary table](./vocabulary.md).
