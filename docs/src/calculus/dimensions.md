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

A **dimension** is a function from `slots` to the integers. Write `D` for the set
of them. Multiplication of dimensions is addition of exponent vectors, so:

- the zero vector is the **dimensionless** dimension;
- `D` is closed under addition, which is associative and commutative;
- every element has an inverse, the negation of its exponents.

So `D` is an additive abelian group, and it is *free* on `slots`: every element is
a unique integer combination of the seven, with no relation among them. Freedom is
what makes the exponents decidable — two dimensions are equal exactly when their
seven exponents are, and no cleverness about which of them are "really" base is
required or permitted.

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

*Status: not yet proved. Will be `Azoth.Dim.dimension_is_the_universal_invariant`.*

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

*Status: not yet proved. Will be `Azoth.Kind.refinement_is_faithful`, with the
non-injectivity of `Kind → Dim` as `Azoth.Kind.kind_to_dim_is_not_injective`, so
that the true statement and the false one are both in the development and only one
of them is proved.*

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
