/-
Powers with a fractional exponent: what the dimension does, and what the value means.

`docs/src/calculus/numerics.md` states the rule; this module is its proof. Two things are
at stake and they are separate:

**The dimension.** A power scales every exponent by the exponent, so `x^q` has dimension
`q • dim x` - a statement over `ℚ`, which is why `Azoth.Dim`'s carrier is the free
`ℚ`-module rather than a free abelian group. `exponents_smul` is the law, and
`not_integer_smul` is the consequence that matters: a `q` that is not an integer puts the
dimension outside the vocabulary's integer lattice, so no `uom` type can name it.

**The value.** `x^q` for `x < 0` has no single meaning, and the three conventions in play
each return a different number:

| who | `(-8)^(1/3)` | rule |
|---|---|---|
| C, and Rust's `f64::powf` | `NaN` | `exp (q log x)` with a negative argument |
| Lean, `Real.rpow` | `1` | `exp (log x · q) · cos (q π)`, the real part of the principal complex power |
| the real odd root | `-2` | `exp (log |x| · q)` with the sign carried separately |

The theorems below establish the second column, and `the_conventions_disagree` says the
second and third can never agree on the negative axis: one has the sign of `cos (q π)` and
the other the sign of the base. That is not a rounding difference between two
implementations of one function - it is two different functions, and a kernel that does
not say which it means is not implementing either.

**Where the three agree is exactly the integral exponents**, which is `rpow_natCast`
below and the reason `powi` is safe where `powf` is not. `Azoth.Pow.rpow_natCast` states
it as this development's anchor rather than leaving a reader to find Mathlib's.

# How many roots, which is not two cases

**`x^(p/q)` is a `q`-valued expression**, and the number of *real* values is not settled by
the parity of `q` alone. In lowest terms: one for odd `q` whatever the sign; and for even
`q`, **two** when `x > 0` and none when `x < 0`, since even `q` forces odd `p` and the
radicand then carries the sign of `x`. Over the complexes the count is always `q`, for any
`q` - which is the form of the statement that matters, because `q` is a denominator and not
a flag.

So a single-valued form is a **branch selection**, not a formula. For odd `q` there is
nothing to select - `(-8)^(1/3) = -2` and `(-8)^(2/3) = +4` are each the only real value.
For even `q` and `x > 0` there are two, and `sign(x)^p |x|^(p/q)` returns the non-negative
one: `4^(1/2)` is `±2`, and a kernel returning `+2` has made a convention decision where the
mathematics has not. `even_power_has_two_roots` is the fact underneath: an even power is not
injective, so its inverse is not a function.
-/

import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Azoth.Dim

namespace Azoth

open Units (Dimension)

namespace Dim

/-- **A fractional power scales every exponent by the exponent.**

A dimension is an exponent per base, and raising a quantity to the `q`-th power multiplies
each of them by `q`. That is the whole dimensional content of `x^q`, and it is why the
carrier is a `ℚ`-module: `m**0.5` is a dimension, and its exponents are the halves of the
original's. -/
theorem exponents_smul (q : ℚ) (d : Dimension) :
    exponents (q • d) = (exponents d).map (fun e => q * e) := by
  simp [exponents, smul_impl]

/-- **A dimension the vocabulary cannot name.**

For a rational `q` that is no integer, `q` times a base dimension is not `n` times it for
any integer `n`. So a unit whose dimension carries a fractional exponent has no
representative in the vocabulary's integer lattice, and `specs/schema/vocabulary.schema.json`
is right to require integers: `uom`'s exponent parameters have no fractional member.

The consequence this is here for is `A_phi`, the Debye-Huckel parameter, whose dimension is
`(mol/kg)^-1/2` - a half. It is carried as a bare `f64` for that reason and not by
oversight, and the check that would have said so is this theorem rather than a reviewer.

Freeness is what makes it true: the projection onto one slot reads the exponent off, and two
dimensions equal in the group are equal exponent by exponent. -/
theorem not_integer_smul (q : ℚ) (hq : ∀ n : ℤ, q ≠ n) :
    ∀ n : ℤ, q • Dimension.ofString "L" ≠ (n : ℚ) • Dimension.ofString "L" := by
  intro n h
  have h3 : exponents (q • Dimension.ofString "L")
      = exponents ((n : ℚ) • Dimension.ofString "L") := by rw [h]
  rw [exponents_smul, exponents_smul] at h3
  have h4 : (exponents (Dimension.ofString "L")).map (fun e => q * e)
      = (exponents (Dimension.ofString "L")).map (fun e => (n : ℚ) * e) := h3
  have h5 : (exponents (Dimension.ofString "L")).head? = some 1 := by
    simp [exponents, Azoth.slots, Dimension.ofString]
  have h6 := congrArg (fun l : List ℚ => l.head?) h4
  simp only [List.head?_map, h5, Option.map_some] at h6
  have : q * 1 = (n : ℚ) * 1 := Option.some.inj h6
  exact hq n (by simpa using this)

end Dim

namespace Pow

/-- **The three conventions agree on an integral exponent, and only there.**

`Real.rpow x (n : ℝ) = x ^ n`, so a whole-number exponent is the one case where Lean's
principal-complex definition, C's `pow` and the elementary power are the same function.
It is why `powi` is the safe form and `powf` the one that needs a stated domain. -/
theorem rpow_natCast (x : ℝ) (n : ℕ) : x ^ (n : ℝ) = x ^ n := Real.rpow_natCast x n

/-- **Mathlib's convention for a negative base is the real part of the complex power.**

`x ^ y = exp (log x * y) * cos (y * pi)` for `x < 0`. `exp` is strictly positive, so the
sign of the whole is the sign of `cos (y * pi)` - which for `y = 1/3` is `+1/2`, so
`(-8) ^ (1/3)` comes out **positive** where the real cube root is negative.

This is not a defect in Mathlib: it is the principal branch, chosen and documented. It is
a statement about what a kernel inherits when it uses `rpow` without saying which
convention it means. -/
theorem rpow_third_of_neg_pos {x : ℝ} (hx : x < 0) : 0 < x ^ ((1 : ℝ) / 3) := by
  rw [Real.rpow_def_of_neg hx]
  apply mul_pos (Real.exp_pos _)
  rw [show (1 : ℝ) / 3 * Real.pi = Real.pi / 3 by ring, Real.cos_pi_div_three]
  norm_num

/-- The same for `2/3`, where the cosine is `-1/2` and the principal part is therefore
negative - the opposite sign to the real odd root, which is positive there. -/
theorem rpow_two_thirds_of_neg_neg {x : ℝ} (hx : x < 0) : x ^ ((2 : ℝ) / 3) < 0 := by
  rw [Real.rpow_def_of_neg hx]
  apply mul_neg_of_pos_of_neg (Real.exp_pos _)
  rw [show (2 : ℝ) / 3 * Real.pi = 2 * (Real.pi / 3) by ring, Real.cos_two_mul,
      Real.cos_pi_div_three]
  norm_num

/-- **An even power is not injective, so its inverse has two branches.**

`(-y)^n = y^n` for even `n`. The `n`-th root of a positive number is therefore two real
numbers, not one, and any single-valued root - `Real.sqrt`, or a `sign(x)^p |x|^(p/q)`
formula - is selecting one of them. That selection is a convention the kernel owes its
reader, and it is what `4^(1/2) = +2` is doing rather than a consequence of the arithmetic.

This is also why `q` odd is the *easy* case: with an odd denominator there is one real root
whatever the sign of the radicand, so there is nothing for a formula to choose.-/
theorem even_power_has_two_roots (n : ℕ) (hn : Even n) (y : ℝ) : (-y) ^ n = y ^ n :=
  hn.neg_pow y

/-- **The two conventions can never agree on the negative axis, for a third-power
exponent.**

The principal branch has the sign of `cos (pi/3) = 1/2`, so it is positive; the real odd
root has the sign of the base, so it is negative. They differ by sign and not by rounding,
and for every negative base rather than at some awkward point. -/
theorem the_conventions_disagree {x : ℝ} (hx : x < 0) :
    x ^ ((1 : ℝ) / 3) ≠ -((-x) ^ ((1 : ℝ) / 3)) := by
  have hpos : 0 < x ^ ((1 : ℝ) / 3) := rpow_third_of_neg_pos hx
  have hroot : 0 < (-x) ^ ((1 : ℝ) / 3) := Real.rpow_pos_of_pos (by linarith : (0 : ℝ) < -x) _
  exact (ne_of_lt (lt_trans (neg_lt_zero.mpr hroot) hpos)).symm

end Pow

end Azoth
