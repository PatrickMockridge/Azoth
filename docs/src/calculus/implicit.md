# The sensitivity of a solution

A flash solves

```
F(X(p), p) = 0
```

for the state `X` at the parameters `p` — a temperature, a pressure, a feed. Then it wants
the *derivative* of that solution: how the state moves when a parameter does. That
derivative is not free, and it is not read off the residual.

## The first derivative

**Claim (the sensitivity).** For `X` differentiable where it solves `F`, with the
residual's derivative splitting into a state part `J` and a parameter part `F_p`, and `J`
invertible:

```
X'(p) = −J⁻¹ ∘ F_p
```

Read as an equation in continuous linear maps. `J` is `∂F/∂X` at the solution and `F_p` is
`∂F/∂p` at fixed state; the conclusion says the sensitivity is **minus** the inverse of the
first composed with the second.

**Why it is worth a theorem.** Both halves are easy to get wrong and neither is visible in
the answer. A dropped sign turns a converging Newton step into a diverging one, and it does
so on the states where the answer moves fastest — which is exactly where a caller looks. A
`J⁻¹` replaced by a division by some scalar is not a smaller error, it is a different
formula, and a solver built on it converges on the problems where the Jacobian happens to
be diagonal. Eight times in this port the question "is the algebra wrong or the code wrong?"
was answered "the code", and this is the other formula beside the raw-versus-normalised
conversion that a kernel writes by hand.

**What the theorem does not do** is prove any code correct. It fixes the statement, so that
an implementation can be read against something that has been checked rather than against
another implementation that may share the mistake.

## The second derivative

The plan of record pairs this with `X''(p) = −J⁻¹(F_pp + (∂J/∂p ∘ X') ∘ X')`, the
second-order sensitivity a Hessian-based stability test wants.

*Status: **specified**, not proved. `Azoth/Implicit.lean` proves the first-order claim —
`Azoth.Implicit.hasFDerivAt_of_solution` — and nothing here is stated about the second. The
gap is real and it is an API gap rather than a mathematical one: the second-order statement
needs the derivative of a map into the space of linear maps, which is
`iteratedFDeriv`'s territory and a different set of lemmas from the ones this file uses. A
proof of a weakened version would have no gap and still not be the claim, so the claim is
left stated rather than approximated.*

## What the first-order proof rests on

*Nothing unusual.* `Azoth/Implicit.lean` uses `fderiv` and the chain rule and no inverse
function theorem: the existence of `X` is a hypothesis (`hF`), not a conclusion, because
what a solver needs is the derivative of a solution it already has. The only hypothesis
beyond differentiability is that the state part of the residual's derivative **is** an
equivalence — which is what makes `J⁻¹` a map rather than a formal inverse, and is the
hypothesis a flash checks by refusing when its Jacobian is singular.
