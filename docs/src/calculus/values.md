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

*Status: **checked**, and the claim is narrower than it looks. Each numeric parameter in
`specs/schema/component.schema.json` carries either a bound (`x-azoth-min`/`x-azoth-max` with a
`x-azoth-rationale`) or `x-azoth-unbounded` with the reason there is none, and a parameter
carrying neither is refused by `tools/gen_parameter_bounds.py` rather than skipped. The bounds
reach both readers through one generated table, and each refuses at the site its unit check
already uses — `azoth.keycard._within_bound` and `azoth_eos::card::check_bound`.*

*The layer's own measurement: **2** of the fifteen are bounded today: a critical temperature and a critical pressure, both of
which are absolute and so have zero, and below it, outside any state. The other thirteen state
why not, and each reason is measured rather than assumed — the acentric factor's band is refused
because the shipped table spans `-0.39003` to `2.8353`, so a sign bound would reject neon.*

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

*Enforcement: check — `tools/gen_parameter_bounds.py` refuses a parameter that states neither a bound nor a reason, and emits one table both readers take; `python/src/azoth/core/_bounds_gen.py` and `crates/azoth-core/src/parameter_bounds_gen.rs` are that table, enforced by `azoth.keycard._within_bound` and `azoth_eos::card::check_bound` at the site each reader already checks a value's unit.*
