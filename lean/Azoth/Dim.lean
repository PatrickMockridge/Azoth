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

The inverse of `ofExponents` on vectors of the right width. That the two are
inverse is the claim `docs/src/calculus/dimensions.md` makes when it says the
exponents determine the dimension and no information is lost; it is stated as
`exponents_ofExponents` in `Azoth.Units`, where the vocabulary is in scope.
-/
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

end Dim

end Azoth
