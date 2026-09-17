/-
The keycard as a capability, as `docs/src/calculus/capability.md` states it.

`specification.md` says the keycard "is where responsibility sits": the shipped baseline
is NeqSim's, a user's card is theirs, and the difference between the two is who is
accountable for a value being right. This module fixes what that makes the card *as a
thing in the calculus*.

A capability is a name. `Card k G` holds a grant `G`; the databank answers a lookup only
on a message arriving on `k`; and a calculation that does not hold `k` cannot construct
one, because `k` is not in its scope. That is the whole of what an object capability is,
and it is why the card is a *value* here rather than a configuration.

What this module models is the part of that which is a statement about *results* rather
than about scope, because that is where the claims bite: a computation demands datums, the
grant decides whether it runs, and its answer rests on some of what it demanded.

Three claims, and the third is what keeps the second honest:

* **determinism** - one grant and one computation determine exactly one answer. Nothing
  here is a function of call order, so two runs cannot disagree.
* **non-amplification** - an answer rests only on datums the grant holds. It is the
  capability-safety property, and it is the one that would fail if a computation could
  reach a datum by some route other than the card.
* **non-vacuity** - there is a grant, a datum outside it, and a computation whose answer
  rests on that datum, with the gate refusing for the smaller grant. Without this the
  second claim could hold by `derives` being empty, and a gate that cannot fail is not a
  gate.

A grant is a predicate rather than a `Finset` or a `Set`: `Name` is a mutual inductive
over `Process` and carries no `DecidableEq`, and nothing here needs to decide membership.
Predicates also keep this module on `Azoth.Rho` alone, with no library beneath it.

The `rests` field of `Computation` is the model rather than a theorem: a computation's
answer does not rest on what it did not ask for. Everything above follows from it.
-/

import Azoth.Rho

namespace Azoth

namespace Capability

open Rho

/-- Containment of two predicates-as-sets.

Written out rather than taken from a library so that this module needs nothing beneath
`Azoth.Rho`; `A ⊆ B` is "everything `A` holds, `B` holds". -/
def HoldsWithin (A B : Name → Prop) : Prop :=
  ∀ ⦃d : Name⦄, A d → B d

/-- A grant: the names a card is entitled to. -/
abbrev Grant := Name → Prop

/-- What a computation answers with: the value, and the datums the answer rested on. -/
structure Answer where
  /-- The value the computation produced. -/
  value : Nat
  /-- The datums the answer rested on. -/
  used : Name → Prop

/-- A computation: the datums it demands of the card, its answer, and the fact that the
answer rests on nothing it did not demand.

`rests` is the model, not a theorem. It is what makes non-amplification a statement about
the card rather than about the computation's honesty: given it, an answer cannot rest on a
datum outside the grant, because the grant is what the demands had to be inside. -/
structure Computation where
  /-- The datums the computation asks the card for. -/
  demands : Name → Prop
  /-- The answer it gives, and what that answer rested on. -/
  answer : Answer
  /-- The answer rests on nothing the computation did not demand. -/
  rests : HoldsWithin answer.used demands

/-- `Run G c r`: the computation `c` runs under the grant `G` and answers `r`.

The whole of the gate: `c`'s demands must be inside `G`. There is no other condition, and
no route by which a computation reaches the databank except through what it demands - which
is what holding the name `k` amounts to. -/
def Run (G : Grant) (c : Computation) (r : Answer) : Prop :=
  HoldsWithin c.demands G ∧ r = c.answer

/-- `derives r d`: the answer `r` rested on the datum `d`. -/
def derives (r : Answer) (d : Name) : Prop :=
  r.used d

/-- **Claim (determinism).** One grant and one computation determine exactly one answer.

This is `specification.md`'s rule - *"It is not shared. One keycard per process… A library
whose answers depend on call order is a library that returns two results for one
calculation"* - made formal. There is no card in force here that a call order could change
and no second run to disagree with the first, so the result is a function of the grant and
the computation and of nothing else. -/
theorem run_deterministic {G : Grant} {c : Computation} {r r' : Answer}
    (h : Run G c r) (h' : Run G c r') : r = r' :=
  h.2.trans h'.2.symm

/-- **Claim (determinism), in the form the page states it.** -/
theorem run_exists_unique {G : Grant} {c : Computation} (h : HoldsWithin c.demands G) :
    ∃ r : Answer, Run G c r ∧ ∀ r' : Answer, Run G c r' → r' = r :=
  ⟨c.answer, ⟨h, rfl⟩, fun r' hr => run_deterministic hr ⟨h, rfl⟩⟩

/-- **Claim (non-amplification).** No answer derives a datum the grant does not hold.

The proof is one step, and it is worth saying why that is the point rather than a
weakness: the content is `Computation.rests` - an answer rests on nothing the computation
did not demand - and `Run` requiring the demands to be inside the grant. Given those two
the statement is immediate; a longer proof would mean the gate had been put somewhere
other than where a computation reaches the databank, which is exactly the failure the
claim is about. -/
theorem run_derives_in_grant {G : Grant} {c : Computation} {r : Answer} {d : Name}
    (h : Run G c r) (hd : derives r d) : G d :=
  h.1 (c.rests (by rw [h.2] at hd; exact hd))

/-- The grant the witness uses: one name, and not the one the computation wants. -/
def witnessGrant : Grant := fun n => n = Name.free 0

/-- A datum outside `witnessGrant`. -/
def witnessDatum : Name := Name.free 1

/-- A computation that demands two names and answers with a value resting on both.

Its demands straddle the witness grant: one is held, one is not. -/
def witnessComputation : Computation where
  demands := fun n => n = Name.free 0 ∨ n = Name.free 1
  answer := { value := 7, used := fun n => n = Name.free 0 ∨ n = Name.free 1 }
  rests := by
    intro d hd
    exact hd

/-- The witness datum really is outside the witness grant. -/
theorem witness_datum_is_outside : ¬ witnessGrant witnessDatum := by
  intro h
  exact Nat.noConfusion (Name.free.inj h)

/-- **The gate refuses.** Under the witness grant the computation does not run at all.

This is the half that says the gate is a gate: not that the answer is wrong, but that
there is no answer. -/
theorem the_gate_refuses : ¬ ∃ r : Answer, Run witnessGrant witnessComputation r := by
  rintro ⟨r, hdemands, -⟩
  exact witness_datum_is_outside (hdemands (Or.inr rfl))

/-- **Claim (non-vacuity).** The datum granted, the computation runs and its answer rests
on it.

Together with `the_gate_refuses` this is the witness `capability.md` requires: there is a
grant, a datum outside it, and a computation whose answer rests on that datum, and for the
smaller grant `derives r witnessDatum` is not derivable *because no `r` is*. So
`run_derives_in_grant` is not true by `derives` being empty.

`derives` is inhabited in the granted case and the gate fails in the refused one, and the
difference between the two grants is exactly the one datum. -/
theorem the_gate_can_succeed :
    ∃ r : Answer,
      Run (fun n => witnessGrant n ∨ n = witnessDatum) witnessComputation r ∧
        derives r witnessDatum := by
  refine ⟨witnessComputation.answer, ⟨?_, rfl⟩, Or.inr rfl⟩
  intro d hd
  exact hd.elim Or.inl (fun h => Or.inr h)

end Capability

end Azoth
