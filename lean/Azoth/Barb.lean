/-
Barbs and barbed bisimulation, as `docs/src/calculus/barbs.md` states them.

A barb `P ↓ c` says `P` can interact on channel `c` at once: it carries an `out`
or `in` prefix on `c` through parallel composition. Barbed bisimulation is the
observational equivalence that preserves barbs and reduction closure, in both
directions; this module fixes the barbs and the reduction closure, and defines
bisimulation as the union of all bisimulations - the largest one, which is the
tractable form of a coinductive greatest fixed point. The dimensional barb
`P ↓ (c, f, δ)` is the thermodynamic instantiation, a later tranche; the barb here
is the general one over channels.
-/

import Azoth.Rho

namespace Azoth

namespace Barb

open Rho

/-- `P ↓ c`: `P` offers on channel `c` through parallel composition. -/
inductive Barb : Process → Name → Prop where
  | out {a b : Name} {P : Process} : Barb (Process.out a b P) a
  | in {a : Name} {P : Process} : Barb (Process.in a P) a
  | par_left {P Q : Process} {a : Name} : Barb P a → Barb (Process.par P Q) a
  | par_right {P Q : Process} {a : Name} : Barb Q a → Barb (Process.par P Q) a

/-- The reflexive-transitive closure of the reduction. -/
inductive ReducesStar : Process → Process → Prop where
  | refl (P : Process) : ReducesStar P P
  | step {P Q R : Process} : Reduces P Q → ReducesStar Q R → ReducesStar P R

/-- The closure is reflexive. -/
theorem reduces_star_refl (P : Process) : ReducesStar P P :=
  ReducesStar.refl P

/-- A reduction step is in the closure. -/
theorem reduces_star_step {P Q : Process} (h : Reduces P Q) : ReducesStar P Q :=
  ReducesStar.step h (ReducesStar.refl Q)

/-- The closure is transitive. -/
theorem reduces_star_trans {P Q R : Process} (h1 : ReducesStar P Q) (h2 : ReducesStar Q R) :
    ReducesStar P R := by
  induction h1 generalizing R with
  | refl _ => exact h2
  | step hP _ ih => exact ReducesStar.step hP (ih h2)

/-- A relation is a bisimulation when it preserves barbs and reduction closure in
both directions. -/
def IsBisim (R : Process → Process → Prop) : Prop :=
  ∀ {P Q : Process}, R P Q →
    (∀ a : Name, Barb P a ↔ Barb Q a) ∧
    (∀ P', Reduces P P' → ∃ Q', ReducesStar Q Q' ∧ R P' Q') ∧
    (∀ Q', Reduces Q Q' → ∃ P', ReducesStar P P' ∧ R P' Q')

/-- Barbed bisimulation: two processes are bisimilar when some bisimulation relates
them - the union of all bisimulations, the largest one. -/
def BarbedBisim (P Q : Process) : Prop :=
  ∃ R : Process → Process → Prop, IsBisim R ∧ R P Q

/-- Bisimulation is reflexive: every process is bisimilar to itself. -/
theorem barbed_bisim_refl (P : Process) : BarbedBisim P P := by
  refine ⟨fun x y => x = y, ?_, rfl⟩
  intro x y h
  subst y
  constructor
  · intro a
    rfl
  constructor
  · intro x' hred
    exact ⟨x', ReducesStar.step hred (ReducesStar.refl x'), rfl⟩
  · intro y' hred
    exact ⟨y', ReducesStar.step hred (ReducesStar.refl y'), rfl⟩

end Barb

end Azoth
