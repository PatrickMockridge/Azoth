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

/-- Bisimulation is symmetric: the symmetric closure of a bisimulation is one. -/
theorem barbed_bisim_symm {P Q : Process} (h : BarbedBisim P Q) : BarbedBisim Q P := by
  rcases h with ⟨R, hisBisim, hPQ⟩
  refine ⟨fun x y => R x y ∨ R y x, ?_, Or.inr hPQ⟩
  intro x y hxy
  rcases hxy with hxy | hyx
  · rcases hisBisim hxy with ⟨hbarb, hredL, hredR⟩
    exact ⟨hbarb,
      (fun x' hx' => by
        rcases hredL x' hx' with ⟨y', hstar, hR⟩
        exact ⟨y', hstar, Or.inl hR⟩),
      (fun y' hy' => by
        rcases hredR y' hy' with ⟨x', hstar, hR⟩
        exact ⟨x', hstar, Or.inl hR⟩)⟩
  · rcases hisBisim hyx with ⟨hbarb, hredL, hredR⟩
    exact ⟨(fun a => (hbarb a).symm),
      (fun x' hx' => by
        rcases hredR x' hx' with ⟨y', hstar, hR⟩
        exact ⟨y', hstar, Or.inr hR⟩),
      (fun y' hy' => by
        rcases hredL y' hy' with ⟨x', hstar, hR⟩
        exact ⟨x', hstar, Or.inr hR⟩)⟩

/-- If a bisimulation relates `y` and `z`, a multi-step reduction of `y` is matched
by a multi-step reduction of `z`. -/
theorem bisim_matches_star {Rel : Process → Process → Prop} (hisBisim : IsBisim Rel)
    {y y' z : Process} (hyy' : ReducesStar y y') (hyz : Rel y z) :
    ∃ z', ReducesStar z z' ∧ Rel y' z' := by
  induction hyy' generalizing z with
  | refl _ => exact ⟨z, ReducesStar.refl z, hyz⟩
  | @step y y1 y' hred hstar ih =>
      rcases (hisBisim hyz).2.1 y1 hred with ⟨z1, hzz1, hy1z1⟩
      rcases ih hy1z1 with ⟨z', hz1z', hy'z'⟩
      exact ⟨z', reduces_star_trans hzz1 hz1z', hy'z'⟩

/-- The reverse: a multi-step reduction of the second process is matched by the
first. -/
theorem bisim_matches_star_right {Rel : Process → Process → Prop} (hisBisim : IsBisim Rel)
    {x y y' : Process} (hyy' : ReducesStar y y') (hxy : Rel x y) :
    ∃ x', ReducesStar x x' ∧ Rel x' y' := by
  induction hyy' generalizing x with
  | refl _ => exact ⟨x, ReducesStar.refl x, hxy⟩
  | @step y y1 y' hred hstar ih =>
      rcases (hisBisim hxy).2.2 y1 hred with ⟨x1, hxx1, hx1y1⟩
      rcases ih hx1y1 with ⟨x', hx1x', hx'y'⟩
      exact ⟨x', reduces_star_trans hxx1 hx1x', hx'y'⟩

/-- Bisimulation is transitive: the composition of two bisimulations is one. -/
theorem barbed_bisim_trans {P Q R : Process} (h1 : BarbedBisim P Q) (h2 : BarbedBisim Q R) :
    BarbedBisim P R := by
  rcases h1 with ⟨R1, hisBisim1, hPQ⟩
  rcases h2 with ⟨R2, hisBisim2, hQR⟩
  let S := fun x z => ∃ y, R1 x y ∧ R2 y z
  have hS : IsBisim S := by
    intro x z hxz
    rcases hxz with ⟨y, hxy, hyz⟩
    rcases hisBisim1 hxy with ⟨hbarb1, hredL1, _⟩
    rcases hisBisim2 hyz with ⟨hbarb2, _, hredR2⟩
    constructor
    · intro a
      exact (hbarb1 a).trans (hbarb2 a)
    constructor
    · intro x' hx'
      rcases hredL1 x' hx' with ⟨y', hyy', hx'y'⟩
      rcases bisim_matches_star hisBisim2 hyy' hyz with ⟨z', hzz', hy'z'⟩
      exact ⟨z', hzz', ⟨y', hx'y', hy'z'⟩⟩
    · intro z' hz'
      rcases hredR2 z' hz' with ⟨y', hyy', hy'z'⟩
      rcases bisim_matches_star_right hisBisim1 hyy' hxy with ⟨x', hxx', hx'y'⟩
      exact ⟨x', hxx', ⟨y', hx'y', hy'z'⟩⟩
  exact ⟨S, hS, ⟨Q, hPQ, hQR⟩⟩

end Barb

end Azoth
