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

*Status: not yet proved. Will be
`Azoth.Capability.result_is_a_function_of_the_card`.*

This is the property `spec.md` states as a rule — *"It is not shared. One keycard
per process… A library whose answers depend on call order is a library that returns
two results for one calculation."* The claim is that rule, made formal: the result
is a function of the card and the inputs, and of nothing else.

**Claim (non-amplification).** No result derives a value whose origin the card
does not grant:

```
Run c x r   →   derives r d   →   d ∈ G
```

*Status: not yet proved. Will be `Azoth.Capability.no_amplification`.*

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

*Status: not yet proved. Will be `Azoth.Capability.amplification_is_possible`.*

Both statements are in the development, and only the second is what makes the
first worth proving. A gate that cannot fail is not a gate — the same reason
`spec.md` requires a bound to name what it bounds, and the same reason the
conversion test in [The vocabulary table](./vocabulary.md) is checked by pointing
a unit at the wrong library rather than by reading the code.

## What the implementation does today

Two facts, each stated because the claims above are about a design the code does
not yet have.

**The card is a module-level global in Python.** `azoth.keycard.load(path)` sets a
process-wide value, and every calculation reads it. So in one process:

```
load("a.yaml"); r1 = calc(...)
load("b.yaml"); r2 = calc(...)
```

gives two different answers to one calculation, decided by call order. That is
exactly the state `spec.md` names as the thing a keycard must not be, and the two
are not reconcilable by reading the sentence differently: a mutable that a caller
sets between two calls *is* an answer that depends on call order. Closing this is
what makes the determinism claim above a description of the library rather than a
requirement on it.

**The Rust side has no card at all.** `crates/azoth-eos/src/databank.rs` reads
only the files compiled into the binary, and a user's overlay never reaches it. So
a Rust caller and a Python caller with the same card get different physics from
the same inputs — which is the one thing the two-implementations rule exists to
prevent. The capability is what closes it: an overlay is a value passed in, and a
value passed in is the thing a compile-time embed cannot be.

Neither fact is a claim about the calculus; both are reasons the calculus is
stated now rather than after more calculations land on top of them.

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
