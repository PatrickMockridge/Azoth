/-
The connection judgment: what a graph's edges are allowed to join.

`docs/src/calculus/session.md` states the middleware's layer, and its first claim is about the
graph a document's connections make. A port declares a **channel type** - a finite record of
fields, each at a dimension and a shape - and two ports may be joined when the records agree:
the same field names, each at the same dimension. `crates/azoth-process/src/check.rs` runs that
rule on every connection of a flowsheet, and its verdict is `Mismatch`, which names the field
that differs.

# The judgment

`check` below is that rule as an algorithm: equal length, then every field of the first record
found in the second. The theorem is **soundness** - accepted implies supplied - and it is not a
projection of its own definition: it says something about a `Bool` that a sentence about `∈`
does not, which is the difference between a decision procedure and a restatement of what one
would decide.

# Why the record is a list of pairs and not a lookup

The Rust carrier is a `BTreeMap`, so a field's name is a key and a record cannot carry one twice.
A list can, and the first draft of this module wrote the judgment as `find? name record = some
field` - first match wins - and then could not prove that a found field is a member: a duplicate
name later in the record is shadowed, and the lemma is false. So the judgment asks whether the
*whole pair* occurs, which is the same question on any record whose names are distinct, and the
length test is what makes the two records' name sets agree. The shadowing case is not modelled
because the Rust cannot express it, and saying so is better than a lemma that is false for
inputs the tree cannot produce.

# Two witnesses, because a soundness proof can be vacuous

`check_refuses_a_missing_field` shows the judgment can refuse, so `check_sound` is not a
statement about a function that always accepts. And `compatibility_is_not_equality` shows two
records in a different order that pass and are **not** equal - the guard
`Azoth.Process.swapped_declarations_are_congruent_not_equal` is the same idea one layer up: a
proof whose two sides were definitionally equal would be true for free.

# What is not here

The rule cannot fire on the shipped palette, and that is measured rather than assumed: all 77
of its port field-sets declare the same five dimensions, so no two of its ports mismatch. This
module is about the rule a *caller's* palette reaches, and the witness above is a record no
palette in this tree declares.
-/

import Mathlib.Data.List.Basic

namespace Azoth

namespace Graph

/-- One field of a channel: the dimension it carries, and whether it is a vector.

The dimension is a name rather than an exponent vector here, and deliberately: what
`check.rs` compares is the *exponent tuple* the vocabulary resolves that name to, so a model
that carried the exponents would be modelling the vocabulary's table rather than the rule. -/
structure Field where
  dimension : String
  vector : Bool
deriving DecidableEq, Repr

/-- A channel type: the fields a port declares.

**A list of pairs whose names are distinct**, which is what the Rust's `BTreeMap` is. Nothing
in the type says so; the length test in `check` is what makes two records' name sets agree, and
the shadowing a duplicate would cause is not modelled. -/
abbrev Chan := List (String × Field)

/-- **Why two channel types do not fit together.** The Lean reading of `check::Mismatch`,
which is what a canvas marks. -/
inductive Mismatch where
  /-- The records are different sizes, so there is nothing to line up. -/
  | differentFields
  /-- One record declares a field the other does not. -/
  | missingField (field : String)
  /-- Both declare the field, at different dimensions or shapes. -/
  | differentField (field : String)
deriving DecidableEq, Repr

/-- **The judgment.** Equal length, then every field of `a` present in `b`.

This is `check::compatible`'s algorithm, including its length test: two records that agree field
by field but differ in size would be a document where one port declares a field twice. -/
def check (a b : Chan) : Bool :=
  (a.length == b.length) && a.all fun p => b.any fun q => q == p

/-- **Soundness: accepted implies supplied.** Every field the first record declares, the second
declares too, at the same dimension and shape.

This is the claim `docs/src/calculus/session.md` makes about a graph: an edge the checker
accepts carries a stream whose every field the consumer's own declaration names. -/
theorem check_sound (a b : Chan) (h : check a b = true) : ∀ p ∈ a, p ∈ b := by
  intro p hp
  rw [check, Bool.and_eq_true] at h
  obtain ⟨-, all⟩ := h
  rw [List.all_eq_true] at all
  have found := all p hp
  rw [List.any_eq_true] at found
  obtain ⟨q, hq, heq⟩ := found
  rw [beq_iff_eq] at heq
  rwa [heq] at hq

/-- **A record agrees with itself**, which is what makes the judgment total on the documents a
palette can produce. -/
theorem check_refl (a : Chan) : check a a = true := by
  rw [check, Bool.and_eq_true]
  refine ⟨beq_iff_eq.mpr rfl, ?_⟩
  rw [List.all_eq_true]
  intro p hp
  rw [List.any_eq_true]
  exact ⟨p, hp, beq_iff_eq.mpr rfl⟩

/-- **The judgment refuses.** A witness that `check_sound` is not a sentence about a function
that always answers `true`: the second record carries none of the first's fields. -/
theorem check_refuses_a_missing_field :
    check [("P", Field.mk "pressure" false)] [("T", Field.mk "thermodynamic_temperature" false)] = false := by
  decide

/-- **Compatibility is not equality.** Two records that declare the same fields in a different
order pass the judgment and are not the same record - so it is the *record* the rule is about
and not the way it was written down, and `check_sound` is not `rfl` in disguise. -/
theorem compatibility_is_not_equality :
    check [("P", Field.mk "pressure" false), ("T", Field.mk "thermodynamic_temperature" false)]
        [("T", Field.mk "thermodynamic_temperature" false), ("P", Field.mk "pressure" false)] = true ∧
      [("P", Field.mk "pressure" false), ("T", Field.mk "thermodynamic_temperature" false)] ≠
        [("T", Field.mk "thermodynamic_temperature" false), ("P", Field.mk "pressure" false)] := by
  constructor
  · decide
  · simp

end Graph

end Azoth
