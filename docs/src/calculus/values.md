# What a card's value may be

[The keycard as a capability](./capability.md) fixes *which* datums an answer may rest on. It says
nothing about *what a datum says*. `Azoth.Capability.Answer` carries a `value : Nat` and the three
claims there only ever touch which names the answer used — so a card that grants `methane.Tc` grants
the right to use whatever number is written under it, and today that may be `-1.9e9 K`.

This page is the layer that is missing, stated so the tranche that builds it has something to be
checked against. Its status is **specified**: the claims below are what make that tranche checkable,
and they are proved when it lands.

## The distinction the layer turns on

A value is not a thing the tree can prove. `Tc = 190.564 K` is empirical: it is a measurement
somebody made, and no theorem in `lean/Azoth/` bears on it. What a tree *can* do is fix the
**domain** a value has to lie in for the arithmetic over it to mean anything — which is what
`Azoth/Guards.lean` already does for a calculation's inputs, and what nothing does for a card's.

So the layer has one claim, and it is about admissibility rather than about truth:

**Claim (admissibility).** A datum a bank holds is inside the domain the layer declares for it, or
the lookup is refused.

*Status: **specified.** The layer does not exist. Nothing reads a bound on a component parameter
today: `specs/schema/component.schema.json` declares a parameter as `{value, unit, citation}` with
no `minimum` or `maximum`, and the two readers — `python/src/azoth/keycard.py` and
`crates/azoth-eos/src/card.rs` — check a value's *unit* and its *dimension* and never its
magnitude. The tranche that builds it adds the bound vocabulary to the schema, the enforcement to
both readers at the site the unit check already uses (`_quantity` and `converted`), and this
paragraph becomes the status it earns.*

## What the layer does not do

**It does not prove a value.** No bound is the right bound, and nothing here claims one is. A bound
is the model author's statement about the domain the implementation requires — `Pr`'s alpha is
defined where `Tc > 0` and `Tr ≥ 0`, so a card stating `Tc = -1` cannot be used — and what is
provable is that enforcing it makes the arithmetic total, never that the number is the physics.

**It does not make a datum's provenance checkable.** `docs/src/architecture/specification.md` says
it plainly: "Nobody can check whether a person read a standard", which is why the card carries no
status field. A bound is a fact about the domain, not a mark of quality, and the two must not be
read as one.

**It does not bound everything.** Where no domain statement exists — the acentric factor's band,
the ideal-gas polynomial's coefficients, the dielectric constants — the entry is an explicit
`unformalised` with its reason, in `lean/guards.toml`'s existing vocabulary. That is a recorded
decision rather than an omission, and it is the same shape as the class table the guards already
carry.

**It does not reach the databank's values.** The layer states a domain; the data states a magnitude.
Keeping them apart is what stops a measurement being dressed as a theorem, and it is why nothing in
`lean/` reads `data/`.

*Enforcement: nothing — the layer is specified, so there is nothing in the implementation for a failure to appear against yet. The tranche that builds it replaces this line with the checker it adds.*
