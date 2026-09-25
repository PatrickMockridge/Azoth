/-
The process layer's claim, as `docs/src/calculus/process.md` states it: a unit operation
whose only free names are its declared channels, whose body is a total function of its
inlets and which performs no other interaction, is barbed-congruent to the pure evaluation
of that function.

**The design pass settled the statement, and this is what it settled on.** The page's own
warning is that the machinery has no value semantics: `Barb` is `Process → Name → Prop`
(`Azoth/Barb.lean`), so a barb records a channel and not a magnitude, and there is no place
for "a function of the inlets" to land as *values*. What the calculus can state is therefore
the extensional half - **the declaration determines the observable behaviour** - and that is
the half a port declaration is for: it is "the whole of what an implementation has to agree
with", in the page's words.

So a unit operation is encoded as the process that *offers* exactly its declared ports, and
the theorem relates two such processes rather than a process and a function value:

* **`Declares P C`** is "only its declared channels" - the barbs of `P` are exactly `C`.
* **`Inert P`** is "total, and no other interaction" - `P` has no step to take, so it cannot
  fail on a value its type admits and cannot act outside its ports. Both halves of the third
  hypothesis land here, which is why one predicate carries two of the three.
* **`unit_op_is_extensional`** is then the claim: two inert processes declaring the same
  channels are barbed-bisimilar.

**The non-vacuity witness is part of the result, not a footnote.** `barbed_bisim_refl` proves
`BarbedBisim P P` for every `P`, so a statement whose two sides were definitionally equal
would be true for free and would say nothing. `swapped_declarations_are_congruent_not_equal`
exhibits two declarations of the *same* channels that are congruent by the theorem and are
not equal as terms, which is what rules that out.

**What is not proved here** is the pi/rho reading of the same claim. `Rho.md` and
`docs/src/calculus/rho.md` state that a reflected loop has the same fixed point, and
`Azoth.Process.recycle_is_a_fixed_point` - the other claim `process.md` states - is
**characterised and not specified**, which is a change from when this note was first written and
the reason is worth recording.

**The first reason given was the wrong one.** It said the claim "needs a notion of a loop's
transfer function being a contraction, which this development has no metric for". A metric is
not the binding gap. `BarbedBisim` is `IsBisim`, and `IsBisim` is barb-set equality together
with the two reduction-closure clauses - so a statement whose conclusion is `BarbedBisim` is a
statement about **which channels a process communicates on and how it steps**, and the five
constructors of `Process` carry no value at all. No two processes in this development differ in
what they compute, so a loop and a process emitting its fixed point barb the same channel
whatever the loop converges to: the claim's conclusion holds for **every** solution and not only
the unique one, and uniqueness is the half `process.md` says matters.

A metric would let the contraction hypothesis be *stated*. It would not give the conclusion
anything to fail against, which is what separates a claim this repository can guard from one it
can only describe - the division `docs/src/calculus/index.md` draws between Specified and
Characterised. So the claim is characterised, and per that page's own rule it is **not dressed
as a theorem**: a lemma restating that congruence sees only barbs and steps would be the
vacuity the rule exists to prevent.
-/

import Azoth.Barb

namespace Azoth

namespace Process

open Rho
open Azoth.Barb

/-- Which way a declared port points. -/
inductive Polarity where
  /-- An inlet. -/
  | input
  /-- An outlet. -/
  | output
  deriving DecidableEq, Repr

/-- One declared port: a channel and its polarity.

A polarity is an annotation and not an intrinsic - the pi-calculus has no direction - so it
is said once here, at the port, which is where `docs/src/calculus/process.md` says it is. -/
structure Port where
  /-- The channel the port names. -/
  channel : Name
  /-- Which way it points. -/
  polarity : Polarity

/-- The action a port declares: a receive on an inlet, a send on an outlet.

An output sends the channel's own name, because the calculus records *that* a barb happened
at a channel and not what crossed it. That is the design pass's first consequence: a port
declaration encodes to a process, and the process carries a channel set rather than a value
domain. -/
def Port.action (p : Port) : Process :=
  match p.polarity with
  | .input => .in p.channel .nil
  | .output => .out p.channel p.channel .nil

/-- **The encoding of a port declaration into `Rho.Process`**: the parallel composition of
its ports' actions. -/
def declaration : List Port → Process
  | [] => .nil
  | p :: ps => .par p.action (declaration ps)

/-- The channels a declaration names. -/
def declared : List Port → Name → Prop
  | [], _ => False
  | p :: ps, a => a = p.channel ∨ declared ps a

/-- A receive barbs on its channel and on nothing else. -/
theorem barb_in {a b : Name} {P : Process} : Barb (Process.in a P) b ↔ b = a := by
  constructor
  · intro h
    cases h
    rfl
  · intro h
    rw [h]
    exact Barb.in

/-- A send barbs on the channel it sends on and on nothing else. -/
theorem barb_out {a b c : Name} {P : Process} : Barb (Process.out a b P) c ↔ c = a := by
  constructor
  · intro h
    cases h
    rfl
  · intro h
    rw [h]
    exact Barb.out

/-- **A declaration's barbs are exactly the channels it declares.** This is what makes the
encoding an encoding: it is the first hypothesis of `unit_op_is_extensional`, discharged for
any declaration rather than per example. -/
theorem barb_declaration (ports : List Port) (a : Name) :
    Barb (declaration ports) a ↔ declared ports a := by
  induction ports with
  | nil =>
      constructor
      · intro h
        exact absurd h (barb_nil a)
      · intro h
        exact h.elim
  | cons p ps ih =>
      obtain ⟨channel, polarity⟩ := p
      cases polarity
      · show Barb (Process.par (Process.in channel Process.nil) (declaration ps)) a ↔
          (a = channel ∨ declared ps a)
        rw [barb_par, barb_in, ih]
      · show Barb (Process.par (Process.out channel channel Process.nil) (declaration ps)) a ↔
          (a = channel ∨ declared ps a)
        rw [barb_par, barb_out, ih]

/-- **Inert**: a process that takes no step.

This is where "total" and "no other interaction" both land. A total function of its inlets
cannot fail, so it has no reduction to take; and a process whose only actions are its ports'
barbs cannot act outside them, so it has no reduction either. -/
def Inert (P : Process) : Prop := ∀ Q, ¬ Reduces P Q

/-- **Declares**: a process offers exactly the channels `C` and no others. -/
def Declares (P : Process) (C : Name → Prop) : Prop := ∀ a, Barb P a ↔ C a

/-- A receive is inert: `Reduces` has no constructor for one. -/
theorem inert_in (a : Name) (P : Process) : Inert (Process.in a P) := by
  intro Q h
  cases h

/-- A send is inert on its own: a communication needs a partner in parallel. -/
theorem inert_out (a b : Name) (P : Process) : Inert (Process.out a b P) := by
  intro Q h
  cases h

/-- The inactive process is inert. -/
theorem inert_nil : Inert Process.nil := by
  intro Q h
  cases h

/-- **A dropped *free* name is inert, and a dropped quote is not.** `Reduces` has exactly one
constructor whose subject is a `drop` - `drop_quote P : Reduces (drop (quote P)) P` - so a
dropped name takes no step unless it is a quoted process, which is why the theorem is stated
for `.free n` and not for an arbitrary `Name`. -/
theorem inert_drop_free (n : Nat) : Inert (Process.drop (.free n)) := by
  intro Q h
  cases h

/-- **Two receives in parallel are inert, and this is the case with content.** A parallel
composition is *not* inert in general: `par (out a b P) (in a Q)` communicates, which is what a
unit operation's own inlets and outlets do to each other. Inertness therefore holds for a
declaration only while no channel carries a send and a receive together, and that is a property
of the declaration rather than of parallel composition.

This is the design pass's second consequence: the encoding's ports must be on channels that do
not talk to each other, or the declaration is not a declaration of an inert interface. -/
theorem inert_two_inputs (a b : Name) :
    Inert (Process.par (Process.in a Process.nil)
      (Process.par (Process.in b Process.nil) Process.nil)) := by
  intro R h
  cases h with
  | par_left h' => cases h'
  | par_right h' =>
      cases h' with
      | par_left h'' => cases h''
      | par_right h'' => cases h''


/-- **A unit operation is extensional.**

Two inert processes that declare the same channels are barbed-congruent: the declaration is
the whole of what an implementation has to agree with, which is `docs/src/calculus/process.md`'s
adequacy claim.

The witness relation is "the same barbs, and both inert". Carrying inertness in the relation
is what makes the two reduction clauses vacuous *for every pair the relation relates* and not
only for the pair the theorem is stated about - a relation that only knew about `P` and `Q`
would not be a bisimulation. -/
theorem unit_op_is_extensional {P Q : Process} {C : Name → Prop}
    (hP : Declares P C) (hQ : Declares Q C) (iP : Inert P) (iQ : Inert Q) :
    BarbedBisim P Q := by
  refine ⟨fun A B => (∀ a, Barb A a ↔ Barb B a) ∧ Inert A ∧ Inert B, ?_, ?_⟩
  · intro A B hAB
    refine ⟨hAB.1, ?_, ?_⟩
    · intro A' hred
      exact absurd hred (hAB.2.1 A')
    · intro B' hred
      exact absurd hred (hAB.2.2 B')
  · exact ⟨fun a => (hP a).trans (hQ a).symm, iP, iQ⟩

/-- Two inlets on distinct channels. -/
def twoInlets : List Port := [⟨.free 0, .input⟩, ⟨.free 1, .input⟩]

/-- The same two ports, the other way round. -/
def twoInletsSwapped : List Port := [⟨.free 1, .input⟩, ⟨.free 0, .input⟩]

/-- The channel set both declarations name. -/
def twoChannels (a : Name) : Prop := a = .free 0 ∨ a = .free 1

/-- Both declarations name the same channels. -/
theorem twoInlets_declare_twoChannels :
    (∀ a, Barb (declaration twoInlets) a ↔ twoChannels a) ∧
    (∀ a, Barb (declaration twoInletsSwapped) a ↔ twoChannels a) := by
  constructor <;> intro a <;> rw [barb_declaration]
  · simp only [twoInlets, declared, twoChannels, or_false]
  · simp only [twoInletsSwapped, declared, twoChannels]
    constructor
    · intro h
      rcases h with h | h
      · exact Or.inr h
      · rcases h with h | h
        · exact Or.inl h
        · exact h.elim
    · intro h
      rcases h with h | h
      · exact Or.inr (Or.inl h)
      · exact Or.inl h

/-- Both declarations are inert. -/
theorem twoInlets_inert :
    Inert (declaration twoInlets) ∧ Inert (declaration twoInletsSwapped) := by
  constructor <;>
    simp only [declaration, twoInlets, twoInletsSwapped, Port.action] <;>
    exact inert_two_inputs _ _

/-- **The two declarations are different processes.**

`Process.par` is a constructor, so two parallel compositions are equal only when their halves
are, and the halves here are receives on different channels. -/
theorem twoInlets_ne_swapped : declaration twoInlets ≠ declaration twoInletsSwapped := by
  simp only [declaration, twoInlets, twoInletsSwapped, Port.action]
  intro h
  injection h with hhead _
  injection hhead with hname _
  injection hname with hn
  exact Nat.noConfusion hn

/-- **The theorem is not reflexivity, and this is the witness.**

The two declarations are congruent by `unit_op_is_extensional` and they are *not equal*, so
the congruence did not come from `barbed_bisim_refl` - which is the vacuity the design pass
had to rule out, and which a statement relating a process to itself would have fallen into. -/
theorem swapped_declarations_are_congruent_not_equal :
    BarbedBisim (declaration twoInlets) (declaration twoInletsSwapped) ∧
    declaration twoInlets ≠ declaration twoInletsSwapped :=
  ⟨unit_op_is_extensional (twoInlets_declare_twoChannels.1) (twoInlets_declare_twoChannels.2)
      twoInlets_inert.1 twoInlets_inert.2,
    twoInlets_ne_swapped⟩

end Process

end Azoth
