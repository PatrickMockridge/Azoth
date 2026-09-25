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

Nine layers, each a page here and the Lean module it is stated against:

| Layer | What it fixes | Lean | Page |
|---|---|---|---|
| Dimensions | what a unit's dimension is, and when two are equal | `Azoth/Dim.lean` ✅ | [Dimensions](./dimensions.md) |
| The vocabulary table | which units a spec may declare, and what each one is | `Azoth/Vocabulary.lean` ✅ | [The vocabulary table](./vocabulary.md) |
| Numerical safety | which partial functions are total, and what a fractional exponent means | `Azoth/Pow.lean` ✅ | [Numerical safety](./numerics.md) |
| Raw and normalised variables | the derivative between a raw and a normalised variable | `Azoth/Normalisation.lean` ✅ | [Raw and normalised variables](./normalisation.md) |
| The sensitivity of a solution | how a solution moves with the equations that define it | `Azoth/Implicit.lean` (first order) | [The sensitivity of a solution](./implicit.md) |
| Reflection and feedback | feedback, serialisation, and the interoperation surface | `Azoth/Rho.lean` ✅ | [Reflection and feedback](./rho.md) |
| The keycard as a capability | authority a process holds rather than a global it reads | `Azoth/Capability.lean` ✅ | [The keycard as a capability](./capability.md) |
| Barbs | what an observer of a channel can see, and therefore what equality means | `Azoth/Barb.lean` (general barbs) | [Barbs](./barbs.md) |
| Processes and channels | a unit operation as a process on typed, directional channels | `Azoth/Process.lean` (the extensionality claim) | [Processes and channels](./process.md) |

Two further modules are the gate rather than a layer: `Azoth/Axioms.lean` carries the
`#print axioms` line for each general theorem, and the generated `Azoth/Gate.lean` the
one per canonical unit. Eleven modules in all, and every layer the calculus states has one.

The vocabulary is data rather than proof: one hand-written table, compiled into Rust,
Python, JSON Schema and Lean. It is specified separately, in
[The vocabulary table](./vocabulary.md), because it is the one part of this whose
authority is a generator rather than a proof.

## Status: what the gate covers, and what is still specified

**The axiom gate is what makes a proof count.** `tools/check_lean_axioms.py` refuses a
theorem resting on anything outside the three axioms Lean permits — 41 claims in
`Azoth/Axioms.lean`, six of them `Azoth.Dim`'s, and 37 in the generated
`Azoth/Gate.lean`, one per canonical unit, checking that the table's exponents name the
dimension `lean-units` calls by that name.

Most of the nine are proved. [Dimensions](./dimensions.md), [Numerical safety](./numerics.md),
[Raw and normalised variables](./normalisation.md) and [Reflection and
feedback](./rho.md) are complete, and [The keycard as a capability](./capability.md) is
proved against the grant a caller holds and passes (`crates/azoth-eos/src/card.rs`).
[The sensitivity of a solution](./implicit.md) is proved to first order.

Not all of them are. [Barbs](./barbs.md) proves a barbed bisimulation to be an equivalence —
the general machinery, before the *dimensional* barb that would give those two claims
something to be about, which is a later tranche — so they stay specifications.
[Processes and channels](./process.md) states three claims and **one of them is proved**:
`Azoth.Process.unit_op_is_extensional` is the adequacy claim, and it is proved at the level the
layer supports. The other two — one about conservation following from linearity, one about a
recycle having a fixed point — need something this development does not have, and the page
says which: the first is a lemma about values crossing channels and a barb records a channel,
and the second needs a metric to say what a contraction is. Naming what is missing is the
point; a claim weakened into something provable would be a different claim.

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
