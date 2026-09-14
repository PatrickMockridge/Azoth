/-
The dimension group, as `docs/src/calculus/dimensions.md` states it.

A dimension is an element of the free `ℚ`-module on the seven SI base dimensions,
which is what `lean-units` carries it as: a finitely supported function from a
symbol into `ℚ`. That page is the specification and this module is its proof; where
the two disagree, the page is right and this is a bug.

The carrier is `ℚ` rather than `ℤ` for the reason the page gives: a fractional
power is a dimension, and `uom`'s type-level arithmetic is over integers, so the
vocabulary's exponents - which both implementations have to be able to express -
are the integer sublattice of what this module is about.
-/

import LeanUnits.Systems.SI
import LeanUnits.Framework.Dimensions.Lemmas

namespace Azoth

/-- The seven base dimensions, in the order every exponent vector in this
repository is written in.

`uom::si::ISQ`'s order, and the vocabulary table's `slots`. The order is not
cosmetic: every exponent vector in `Vocabulary.lean` and every generated Rust tuple
is written this way, so a different order permutes all seven coordinates at once
and compares equal to any other permuted vector.
-/
def slots : List String := ["L", "M", "T", "I", "Θ", "N", "J"]

namespace Dim

open Units (Dimension)

/-- A dimension from an exponent vector, one exponent per base dimension listed.

`base` is an argument rather than always `slots` so that the definition is a plain
map-and-sum and the lemmas below are `List` lemmas rather than inductions:
`ofExponents` is the only caller that passes `slots`, and it is the only one
anybody should use.

An exponent list shorter than `base` is read against `base`'s leading entries and a
longer one is truncated, because `List.zip` is. Neither is a state the vocabulary
can be in - the generator writes one exponent per slot - so the width is held below
rather than by the type.
-/
def ofExponentsOn (base : List String) (e : List ℚ) : Dimension :=
  ((base.zip e).map (fun p => p.2 • Dimension.ofString p.1)).sum

/-- A dimension from an exponent vector in `slots` order.

The map from the representation the vocabulary and the two implementations use - a
vector of exponents - into the group this module is about.
-/
def ofExponents (e : List ℚ) : Dimension := ofExponentsOn slots e

/-- The exponent of each base dimension in a dimension, in `slots` order.

The inverse of `ofExponents` on vectors of one entry per slot, which is
`exponents_ofExponents` below. -/
def exponents (d : Dimension) : List ℚ := slots.map (fun s => d._impl s)

/-! ## The group laws

`lean-units` gives `Dimension` an `AddCommGroup` and a `Module ℚ` instance, so
associativity, commutativity, the identity and inverses are the instances' and are
not reproved here. What these say is that `ofExponents` *respects* them, which is
the content of "multiplication of dimensions is addition of exponent vectors". -/

/-- The empty exponent vector names the dimensionless dimension. -/
theorem ofExponentsOn_nil (base : List String) : ofExponentsOn base [] = 0 := by
  simp [ofExponentsOn]

/-- A vector of no exponents, in `slots` order. -/
theorem ofExponents_nil : ofExponents [] = 0 := by
  simp [ofExponents, ofExponentsOn]

/-- A vector with one entry names one base dimension, scaled by that entry.

The unit case of the homomorphism, and the one a generator emits a vector for:
every canonical unit's exponents decompose into these. -/
theorem ofExponentsOn_singleton (b : String) (q : ℚ) :
    ofExponentsOn [b] [q] = q • Dimension.ofString b := by
  simp [ofExponentsOn]

/-! ## Reading a dimension back

`exponents` reads a dimension's exponents out, so the two directions are the
representation the vocabulary table, the generated Rust and the generated Python
all exchange and the group element it means. That the two are inverse is what makes
the representation faithful, and it is the claim in this module a change could
break - silently, because every one of the twenty-four unit theorems in
`Vocabulary.lean` would still agree with a table that meant nothing.

The four lemmas below say what the projection `_impl` does to the operations, and
each is `rfl`: `Dimension`'s `AddCommGroup` and `Module` instances are the
`DFinsupp` ones carried across the equivalence `_impl`. They are stated because
`simp` and `rw` will not use a definitional equality on their own, so without them
the projection never reaches the arithmetic. -/

/-- The projecton of a sum is the sum of the projections. -/
theorem add_impl (d₁ d₂ : Dimension) (s : String) :
    (d₁ + d₂)._impl s = d₁._impl s + d₂._impl s := rfl

/-- The projection of a scaling is the scaling of the projection. -/
theorem smul_impl (q : ℚ) (d : Dimension) (s : String) :
    (q • d)._impl s = q • d._impl s := rfl

/-- The projection of the dimensionless dimension is zero, at every slot. -/
theorem zero_impl (s : String) : (Dimension._impl (0 : Dimension)) s = 0 := rfl

/-- A base dimension carries one in its own slot and nothing in any other.

The `if` is on the slot names and both are literals at every use, so the branch
resolves when the goal is evaluated rather than needing a side condition. -/
theorem ofString_impl (s t : String) :
    (Dimension.ofString s)._impl t = if s = t then 1 else 0 := by
  by_cases h : s = t
  · subst h
    simp [Dimension.ofString]
  · simp [Dimension.ofString, h]

/-- **The exponent vector determines the dimension, and back.**

`ofExponents` and `exponents` are mutually inverse on vectors of one entry per
slot: a vector survives the trip into the group and out again unchanged.

Seven arguments rather than a list, because seven is what `slots` has and what
every vector in this repository is. A general statement over a list of arbitrary
length would need `base.Nodup` as a hypothesis - with a repeated slot the map is
not injective and the claim is false - and no vector this repository produces is
of any other length, so the hypothesis would be carried for nobody. -/
theorem exponents_ofExponents (a b c d e f g : ℚ) :
    exponents (ofExponents [a, b, c, d, e, f, g]) = [a, b, c, d, e, f, g] := by
  simp only [exponents, ofExponents, ofExponentsOn, slots,
    List.zip_cons_cons, List.zip_nil_right, List.map_cons, List.map_nil,
    List.sum_cons, List.sum_nil, add_impl, smul_impl, zero_impl, ofString_impl]
  simp

end Dim

end Azoth
