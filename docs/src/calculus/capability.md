# The keycard as a capability

`specification.md` says the keycard "is where responsibility sits": the shipped baseline is
NeqSim's, a user's card is theirs, and the difference between the two is who is
accountable for a value being right. This page says what that makes the card *as a
thing in the calculus* — and the answer is that it is authority a process holds,
not a file a process reads.

The page is the specification and `lean/Azoth/Capability.lean` is its proof; where the two
disagree, the page is right and the module is a bug. Each claim below names the theorem it
is proved by.

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

*Status: **proved**, by `Azoth.Capability.run_deterministic` and its
`∃!`-shaped form `Azoth.Capability.run_exists_unique`: one grant and one computation
determine exactly one answer, with no card in force that a call order could change. The
implementation states the same thing - `azoth.keycard.load` stores nothing, and
`azoth_eos::card::Card` is a value a caller holds and passes - which is what the section
below measures rather than asserts.*

This is the property `specification.md` states as a rule — *"It is not shared. One keycard
per process… A library whose answers depend on call order is a library that returns
two results for one calculation."* The claim is that rule, made formal: the result
is a function of the card and the inputs, and of nothing else.

**Claim (non-amplification).** No result derives a value whose origin the card
does not grant:

```
Run c x r   →   derives r d   →   d ∈ G
```

*Status: **proved**, by `Azoth.Capability.run_derives_in_grant`. The proof is one
step, and that is the point rather than a weakness: the content is the model's
`Computation.rests` - an answer rests on nothing the computation did not demand - together
with `Run` requiring the demands to lie inside the grant. A longer proof would mean the
gate had been put somewhere other than where a computation reaches the databank, which is
exactly the failure this claim is about.*

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

*Status: **proved**, and in two halves, because the witness is only a witness if both
are. `Azoth.Capability.the_gate_refuses` shows that under the witness grant there is **no**
answer at all - not a wrong one - while `Azoth.Capability.the_gate_can_succeed` shows that
with the one datum added the computation runs and its answer rests on it. The difference
between the two grants is exactly that datum, so `derives` is inhabited in the granted case
and the gate fails in the refused one. Without the second half the claim above could hold
by `derives` being empty.*

*`Azoth.Capability.witness_datum_is_outside` is gated alongside them: it is the fact that
makes the refusal about *this* datum rather than about the grant being empty.*

The witness is what makes the claim above worth proving, and it is not optional. A
gate that cannot fail is not a gate — the same reason `specification.md` requires a bound to
name what it bounds, and the same reason the conversion test in
[The vocabulary table](./vocabulary.md) is checked by pointing a unit at the wrong
library rather than by reading the code.

## What the implementation does today

**Python holds no card, and the precedence the claim above needs is in force.**
`azoth.keycard.load` reads a file and returns what it says; it stores nothing, and
every call that reads a card is handed one:

```
card_a = load("a.toml")
card_b = load("b.toml")
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

**Both halves read a card.** `azoth_eos::card::Card` parses the document — the same
TOML `azoth.keycard` reads — and produces the `azoth_eos::databank::Overlay` that
`entry`, `kij`, `names` and `mixture_of` take. So a Rust-native caller holds a card of
their own rather than being handed the values Python resolved, and
`python/tests/test_card_agreement.py` hands one card's text to both readers and compares
what each resolves it to — name by name, pair by pair.

**A Python caller's card reaches Rust's arithmetic by a different route.** `azoth.eos.
components.mixture_of` resolves the names and applies the card *in Python*, and
`azoth._rust_bridge` crosses the boundary with the numbers, so the flash that runs in
Rust runs on the card's values. Measured rather than reasoned about: a card shifting
methane's `Tc` to 300 K gives `beta = 0.2870200723305141` on the Python backend and
`0.2870200723305131` on the Rust one, against a baseline of `0.6824887179287704`.

**Two routes, one rule, and the rule is compared.** A card naming only `omega` keeps the
shipped `Tc` and `Pc`; that merge is implemented twice now, in `azoth.eos.components`
and in `azoth-eos::databank`, and `python/tests/test_data_agreement.py` asks both the
same question and compares the answers — one case per rule rather than per value.
That test's empty-card case is the guard that makes the rest readable: a comparison
that ignored its argument would pass it and fail every other one.

**The grant carries a model's argument, not only a calculation's, and it carries a
matrix.** A card's `coefficients` section is keyed by id and then by input name, and
`azoth.keycard.coefficient_value` is the single place the question is asked: an explicit
argument wins and the card is not consulted, otherwise the card supplies the value, and if
neither does the error names the entry to add. A vector or a matrix comes back as
**one quantity per entry** rather than one quantity wrapping the array, because pint needs
NumPy for an array magnitude and a coefficient is not worth a dependency; `Card`'s
`CoefficientValue` carries the same shape on the Rust side.

That path had one consumer until this tranche and has two now: `eos.uniquac_activity_
coefficients` takes its `aij` from it, and `eos.ge_uniquac_phase` takes the same matrix as
a *phase* parameter. UNIQUAC is the case that makes the grant load-bearing rather than
notional, because no upstream table carries a UNIQUAC interaction matrix - so `aij` can
arrive only as an argument or from a card, and a card that supplied it is a card whose
holder is accountable for it.

## Disclosure, which is the other half

Non-amplification says a result cannot rest on data the card does not grant. It
does not say the result has to *say* what it rested on, and by `specification.md`'s own
account it must: *"What the library owes instead of a status field"* is disclosure.
So the result carries the origins it used, and a caller can see which names an
answer involved rather than inferring it from the card they handed in — which is
the difference between a capability that is enforced and a capability that is
auditable.

The shipped data already carries the `verify_status` column for the library's own
statement about itself. What the capability adds is the *user's* half: which of
their names, and which of the library's, an answer actually rested on.
