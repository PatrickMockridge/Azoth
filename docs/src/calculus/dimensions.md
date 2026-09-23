# Dimensions

A dimension is an exponent for each of the seven SI base dimensions. That is the
whole of the definition, and everything else on this page is a consequence of it.

## The group

Fix the basis, in the order the vocabulary table writes it:

```
slots = [L, M, T, I, Th, N, J]      length, mass, time, current,
                                    thermodynamic temperature,
                                    amount of substance, luminous intensity
```

A **dimension** assigns an exponent to each of the seven. Write `D` for the set of
them. Multiplication of dimensions is addition of exponent vectors, so:

- the zero vector is the **dimensionless** dimension;
- `D` is closed under addition, which is associative and commutative;
- every element has an inverse, the negation of its exponents.

So `D` is an additive abelian group, and it is *free* on `slots`: every element is
a unique combination of the seven with no relation among them. Freedom is what
makes the exponents decidable — two dimensions are equal exactly when their seven
exponents are, and no cleverness about which of them are "really" base is required
or permitted.

### Rational exponents, and where the integers are

The exponents are **rationals in the development and integers in the vocabulary**,
and the difference is not a detail.

`lean-units` carries a dimension as a finitely supported function into `ℚ`, which
is what a fractional power needs, so the Lean development is over the free
`ℚ`-module on `slots` rather than over a free abelian group. That is strictly more
general, and it means `D` is a vector space: `m**0.5` is a dimension, and any
theorem stated over `ℚ` holds for the integer exponents too.

The **vocabulary's** exponents are integers, because `uom`'s type-level arithmetic
is over integers — a `uom` quantity's exponents are type parameters like `P1` and
`N2`, with no fractional member — and a unit whose dimension had a fractional
exponent would be one no `uom` type could name. `specs/schema/vocabulary.schema.json`
enforces `integer` for that reason, and the table is therefore a submodule of the
development's carrier rather than all of it.

So the honest statement is: the development is over `ℚ`; the vocabulary is the
`ℤ`-sublattice that both implementations can express; and the tables
[here](./vocabulary.md) and in `specs/` are integers throughout. A future unit
with a genuinely fractional dimension would be expressible in the calculus and
unrepresentable in Rust, which is the kind of gap the vocabulary table exists to
make visible rather than to discover.

**The order matters and is not a convention.** Every exponent vector in this
repository is written in `slots` order, and the order is `uom::si::ISQ`'s. A
different order permutes all seven coordinates at once, and a permuted vector
compares equal to any other permuted vector — so an ordering mistake is invisible
from either side alone and visible only against the table. That is why the order
is recorded in the table rather than assumed, and why the Rust and Python sides
compare their slot lists.

## Quantities, and the dimension map

A **quantity** is a real number together with a dimension:

```
Q  ≅  ℝ × D
```

Multiplication of quantities multiplies the magnitudes and adds the dimensions,
and taking the dimension is a group homomorphism:

```
dim : (Q, ·) → (D, +)        dim(a, δ) = δ
```

`dim` is surjective — every dimension is some quantity's — and its kernel is the
dimensionless quantities. So `Q` is an extension of `D` by the dimensionless
numbers, and the entire content of "dimensional analysis" is that `dim` is
well-defined and cannot be inverted.

**Claim (dimension is the universal invariant).** Let `φ : Q → G` be a
homomorphism into any group that is invariant under multiplication and powers —
that is, `φ` respects the algebra. Then `φ` factors uniquely through `dim`: there
is exactly one `φ̄ : D → G` with `φ = φ̄ ∘ dim`.

This is the precise sense in which *the dimension is the only thing the algebra
preserves*. It is the formal counterpart of the observation that a conversion
between units is the identity on dimensions, and it is what makes a unit choice
free: units move the magnitude and leave `dim` alone.

*Status: **characterised**. `Azoth.Dim.dimension_is_the_universal_invariant` is not
proved, and deliberately: the group exists and this statement is true of it, but
nothing this repository can do would violate it, so a proof would add a theorem to
a gate whose value comes from every entry in it being one a change could break.
What is proved about the group instead is the part a change *could* break — that
the exponent vector determines the dimension and back, below.*

## The exponents determine the dimension, and back

The representation everything in this repository actually exchanges is the
**exponent vector**, and the group element is what it means. So the claim that
matters for the code is that the two are the same data:

```
Dim.exponents (Dim.ofExponents e)  =  e        for e of `slots.length` entries
```

**Claim (the representation is faithful).** `ofExponents` from exponent vectors in
`slots` order into the dimension group, and `exponents` back, are mutually inverse
on vectors of one entry per slot.

*Status: **proved**, as `Azoth.Dim.exponents_ofExponents`.*

This is the one statement in this section a change could break, and breaking it
would be silent. The vocabulary table, the generated Rust, the generated Python and
the thirty-seven unit theorems all carry vectors; if a vector could name two
different dimensions, or two vectors the same one, every one of those would agree
with a table that meant nothing. The three theorems beside it —
`ofExponentsOn_nil`, `ofExponents_nil` and `ofExponentsOn_singleton` — are the
cases that proof is built from, and each says what `ofExponents` does to a vector
of one element or none.

## Named kinds, and why the refinement is not injective

A dimension is not a name. `J/mol` is the dimension of both an enthalpy and a
Gibbs energy; `J/(mol*K)` is the dimension of both a molar entropy and a molar
heat capacity. The second is not an accident of notation — `azoth-core`'s own
source says a `MolarHeatCapacity` is the carrier for a molar entropy, and gives
the reason: the two are dimensionally identical and a second quantity type would
be a second name for one dimension.

So `Kind → Dim` **cannot be injective**, and a calculus that claims it is has
claimed something false. A **kind** is instead a pair:

```
Kind  ≅  D × Tag
```

with `Tag` a closed set whose default member is `none`. The tag is what separates
two kinds that share a dimension — exactly the mechanism `uom`'s own `quantity!`
macro provides for `AngleKind`, `SurfaceTensionKind` and
`KinematicViscosityKind`, and adopted here for the same reason.

**Claim (the refinement is faithful).** The map `refine : Kind → D × Tag` is
injective. Its first projection is well-defined and surjective onto the named
dimensions, and two distinct kinds are never conflated by the dimension map —
`Entropy` and `MolarHeatCapacity` share a dimension and differ in kind, and
`Enthalpy` and `GibbsEnergy` likewise, or they would not be distinct kinds.

*Status: **specified**. There is no `Kind` in this repository — not in Lean, not in
Rust, not in Python — so nothing could falsify the statement either way, and
`Azoth.Kind.refinement_is_faithful` is not written. What the page fixes is the
*shape* the layer has to have: a pair of a dimension and a tag, not a bare
dimension, because the injectivity a bare dimension would claim is false. The
tranche that builds the kinds proves it, and this paragraph is what it will be
checked against.*

The kinds are the layer at which a unit operation's arithmetic is *typed*, and the
dimensions are the layer at which a conversion is *checked*. A candidate value
enters as a quantity with a dimension; the kind layer is what says that adding a
heat capacity to an entropy is not a permissible operation even though the group
permits the addition. That is a deliberate cost: `Cp = T·dS/dT` is a true identity
that the kind layer rejects without an explicit bridge lemma. The trade is made in
favour of catching `Cp`-for-`S`, which is a real class of error and a silent one.

## What the layer does not do

- **It does not make Rust's type system stronger.** `generic_const_exprs` is
  unstable, so an exponent tuple cannot be an arithmetic const parameter, and
  `uom` closes its quantity list at macro expansion — a quantity added outside
  `uom::si` does not interoperate with the ones inside it. The value of this layer
  is not that more is caught at compile time; it is that **every dimension becomes
  a machine-checked table entry**, so the places the type system cannot reach are a
  listed, checked set rather than an accident. `uom: null` in the table is that
  list.
- **It does not re-derive any physics.** The dimension of a quantity is not a
  physical claim; it is a claim about which group element the unit names, and it
  is checked against the units libraries rather than against a standard.
