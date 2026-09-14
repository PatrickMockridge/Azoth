# The keycard as a capability

`spec.md` says the keycard "is where responsibility sits": the shipped baseline is
NeqSim's, a user's card is theirs, and the difference between the two is who is
accountable for a value being right. This page says what that makes the card *as a
thing in the calculus* — and the answer is that it is authority a process holds,
not a file a process reads.

## A capability is a name

The encoding is three processes and one restriction:

```
(ν k)(  Card k G  |  Databank k  |  Calc c k  )
```

- `k` is a name, freshly bound, and it is not visible outside the bracket.
- `Card k G` holds the grant: a finite set `G` of the names its holder is entitled
  to.
- `Databank k` will answer a lookup, but only on a message that arrives on `k`.
- `Calc c k` is a calculation, and `k` is the only way it can reach the databank.

**Holding `k` is the authority.** There is no check to pass, no credential to
present and nothing to consult: a process that does not hold `k` cannot construct
a message on it, because `k` is not in its scope. That is the whole of what an
object capability is, and it is why the card is a *value* here rather than a
configuration.

This is not a metaphor laid over a lookup table. In the pi-calculus a name's scope
*is* its authority — a name not in scope cannot be used — so "the card grants the
right to use a component" and "the card's name is in scope at the call site" are
the same statement.

## Two claims

**Claim (determinism).** One card and one set of inputs determine exactly one
result:

```
∃! r.  Run c x r
```

*Status: **specified**, and this one is worth being precise about. The card exists
today as a module-level global in Python and as nothing at all in Rust, so the
claim is not merely unproved - **the implementation contradicts it**, in the two
ways the section below sets out. `Azoth.Capability.result_is_a_function_of_the_card`
is what the tranche that makes the card a value is checked against.*

This is the property `spec.md` states as a rule — *"It is not shared. One keycard
per process… A library whose answers depend on call order is a library that returns
two results for one calculation."* The claim is that rule, made formal: the result
is a function of the card and the inputs, and of nothing else.

**Claim (non-amplification).** No result derives a value whose origin the card
does not grant:

```
Run c x r   →   derives r d   →   d ∈ G
```

*Status: **specified**. Authority is not a value here yet, so there is nothing whose
exercise could fail to be a subset of what it holds.*

Non-amplification is the claim that a card cannot be *added to* by using it — that
the authority a calculation exercises is a subset of what it holds. It is the
capability-safety property, and it is the one that would fail if a calculation
could reach a component by some route other than the card.

**It has to be non-vacuous, and that is a requirement on the development rather
than a remark.** `derives` could be defined so that nothing derives anything, and
the claim would hold trivially and mean nothing. So the development carries a
witness in the other direction:

**Claim (non-vacuity).** There is a card `k` with grant `G`, a datum `d ∉ G`, and a
computation whose result relies on `d` — and for that card, `derives r d` is not
derivable.

*Status: **specified**. It is the witness that makes the claim above non-vacuous,
and it is proved in the same tranche - a gate that cannot fail is not a gate.*

The witness is what makes the claim above worth proving, and it is not optional. A
gate that cannot fail is not a gate — the same reason `spec.md` requires a bound to
name what it bounds, and the same reason the conversion test in
[The vocabulary table](./vocabulary.md) is checked by pointing a unit at the wrong
library rather than by reading the code.

## What the implementation does today

**Python holds no card, and the precedence the claim above needs is in force.**
`azoth.keycard.load` reads a file and returns what it says; it stores nothing, and
every call that reads a card is handed one:

```
card_a = load("a.yaml")
card_b = load("b.yaml")
r1 = calc(..., card=card_a)
r2 = calc(..., card=card_b)
```

Two datasets, two answers, one process, and neither call depends on the other having
happened. `python/tests/test_keycard_loader.py` asserts the absence directly - that
no `Keycard` is bound at module scope, and that there is no `current` or `clear` to
reach for - rather than asserting that a missing accessor raises, because a
module-level card still written and no longer read is the same defect one layer down.

So determinism is a description of the library rather than a requirement on it, for
the Python half.

**The Rust side still has no card at all.** `crates/azoth-eos/src/databank.rs` reads
only the files compiled into the binary, and a user's overlay never reaches it. So a
Rust caller and a Python caller with the same card get different physics from the
same inputs — which is the one thing the two-implementations rule exists to prevent.
The capability is what closes it: an overlay is a value passed in, and a value passed
in is the thing a compile-time embed cannot be. That is the next piece of this
tranche, and until it lands the non-amplification claim is half true at best.

## Disclosure, which is the other half

Non-amplification says a result cannot rest on data the card does not grant. It
does not say the result has to *say* what it rested on, and by `spec.md`'s own
account it must: *"What the library owes instead of a status field"* is disclosure.
So the result carries the origins it used, and a caller can see which names an
answer involved rather than inferring it from the card they handed in — which is
the difference between a capability that is enforced and a capability that is
auditable.

The shipped data already carries the `verify_status` column for the library's own
statement about itself. What the capability adds is the *user's* half: which of
their names, and which of the library's, an answer actually rested on.
