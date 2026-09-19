# Raw and normalised variables

An equation of state is written in mole *fractions* and a flash is written in mole
*numbers*. The same mixture, the same residual Helmholtz energy, and two functions of
two different arguments:

```
f(x)               the molar residual A/(RT), at mole fractions x
Phi(n) = (sum n) f(n / sum n)      the extensive A/(RT), at mole numbers n
```

Every composition derivative in the library is one of these or the other, and they are
not the same derivative. `Phi`'s is what a flash wants, because a flash varies mole
numbers at a fixed volume; `f`'s is what an equation of state gives, because it is
written in fractions.

## The conversion

**Claim (the extensive derivative).** For `f` differentiable at the normalised
composition and `sum x ≠ 0`:

```
d Phi/d n_j  =  f(x)  +  d f/d x_j  −  sum_k x_k (d f/d x_k)
```

Three terms and not one. The extensive derivative carries the molar *value*, plus the
molar partial at the normalised point, minus the composition-weighted sum of the molar
partials.

Written as linear maps, with `M` the derivative of `f` and `x` the normalised point:

```
fderiv Phi = M + (f(x) − M x) (x) sum
```

where `sum` is the functional `v ↦ sum_k v_k`. The middle `x` in the second term is
the tensor, not a multiplication.

**Why the last term is there.** `Phi(c n) = c Phi(n)`: the extensive function is
degree-one homogeneous, so Euler's theorem ties its derivatives together and the
normalised partials cannot be read off independently. A kernel that omits the sum is
computing the derivative of a function that is *not* degree-one homogeneous — and it
will agree with any reference that omits it too, which is what makes the error silent.
It has been made twice in this project, once in a `ln phi` and once in a root
sensitivity, and in both cases the check that was supposed to catch it compared two
places that had both dropped the term.

## What is proved, and what that is worth

*Status: **proved**. `Azoth/Normalisation.lean` formalises `extensive` — the molar
function recovered as an extensive one — and proves `Azoth.Normalisation.hasFDerivAt_
extensive`, the conversion as a statement about `fderiv`; its basis form
`Azoth.Normalisation.fderiv_extensive_single`, which is how it is written in code;
and `Azoth.Normalisation.clm_apply_eq_sum_single`, the expansion that makes `M x` the
composition-weighted sum of the molar partials it is claimed to be. Nothing here is
*specified*.*

**A proof prevents no code defect, and this one is not an exception.** Every theorem
below rests on `propext`, `Classical.choice` and `Quot.sound` and nothing else, so it
has no gap — but a kernel can still compute the wrong derivative, and eight times in
this tranche the question "is the algebra wrong or the code wrong?" was answered "the
code". What the theorem removes is the question. It fixes the *statement*, so the code
is wrong rather than the mathematics, and a reader comparing an implementation against
the claim is comparing against something that has been checked rather than against
another implementation.
