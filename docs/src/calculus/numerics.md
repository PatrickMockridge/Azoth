# Numerical safety

A calculation in this library is a real-valued function on a domain, computed in `f64`.
The two are not the same object, and this page is the rule for keeping them together.

**A partial function is not a function.** `sqrt`, `ln`, `powf` with a non-integral
exponent and `/` are each undefined somewhere, and `f64` gives them a value there anyway —
`NaN`, `±inf`, or a finite number from the wrong branch. A kernel that evaluates one
without having established its domain returns a number, and a number is what a caller
checks for, so the failure is silent. That is the whole reason for the rule.

## The rule

Every use of a partial function states **which of the two it is**:

1. **Total by construction.** The argument's domain is established by what it is, not by a
   test — a density floored at `700 kg/m3` before its square root, a mole fraction known to
   be in `[0, 1]`. Stated in a comment naming the reason.
2. **Guarded, with a stated policy.** An explicit test, and the out-of-domain branch does
   one of three *named* things: **refuse** (an error the caller sees), **clamp** (a value
   with a stated bound), or **propagate** (`NaN`, deliberately). Which one is written down,
   because the three are different and all three appear in this library.

An unguarded evaluation is a defect even when the current callers cannot reach the
undefined region. "It cannot happen" is a claim about today's callers, and the region is
reachable the moment one of them moves.

## Fractional exponents

A non-integral exponent is where the rule bites hardest, because the *value* is not merely
undefined at the edge — for a negative base it changes identity.

`x^q` for `x < 0` has three conventions in use, and they return three different numbers:

| who | `(-8)^(1/3)` | rule |
|---|---|---|
| C, and Rust's `f64::powf` | `NaN` | `exp (q log x)`, undefined for `x < 0` |
| Lean, `Real.rpow` | `1` | `exp (log x · q) · cos (q π)` — the real part of the principal complex power |
| the real odd root | `-2` | `exp (log |x| · q)`, sign carried separately |

**The second and third never agree on the negative axis.** `Azoth.Pow.the_conventions_disagree`
proves it: for `q = 1/3` the principal branch has the sign of `cos (π/3)`, so it is positive
where the real root is negative. That is a difference of sign, for every negative base, and
not a rounding difference between two implementations of one function.

So the rule has two parts:

- **A rational exponent is written as a rational.** `1.0 / 3.0`, not `0.3333`. A decimal is
  a different exponent: measured, `x^0.3333` against `x^(1/3)` differs by `4.6e-4` at
  `x = 10^6` and the error grows with `x`. Where a source writes the decimal — NeqSim's
  `ComponentGEWilson` carries `Math.pow(x, 0.3333)` from Morgan–Kobayashi — the decimal is
  **transcribed, and the transcription is a finding to record**, because the model is then
  not the correlation whose exponents are the rationals.
- **An integral exponent uses `powi`.** `x.powi(2)` is total and exact; `x.powf(2.0)` goes
  through `exp (2 log x)` and is `NaN` for `x < 0`. `Azoth.Pow.rpow_natCast` is the theorem:
  the integral exponents are exactly where the conventions coincide.

### How many roots

**`x^(p/q)` is a `q`-valued expression, and any single number returned for it is a branch
selection.** In lowest terms, the number of *real* values is:

| `x^p` | `q` odd | `q` even |
|---|---|---|
| `> 0` | 1 | **2** |
| `< 0` | 1 | 0 |

`q` even forces `p` odd, so `x^p` carries the sign of `x`; the table's second column is
therefore "`x > 0`" and its third is "`x < 0`". Over the complexes the count is always `q`,
for any `q` — which is the form of the statement that matters: the number of roots is a
property of the exponent's denominator, not something a two-case rule settles.

**So `sign(x)^p · |x|^(p/q)` is not a formula, it is a choice.** For `q` odd it picks the
one real root and there is nothing to choose — `(-8)^(1/3) = -2` and `(-8)^(2/3) = +4`, both
unique. For `q` even and `x > 0` it picks the non-negative one and **silently discards the
other**: `4^(1/2)` is `±2`, and a kernel that returns `+2` without saying so has made a
convention decision where the mathematics has not.

`Azoth.Pow.even_power_has_two_roots` is the fact underneath: `(-y)^n = y^n` for even `n`,
so the `n`-th power is not injective and its inverse is not a function. A kernel computing a
root of an even power states which branch it wants — `sqrt` returns the non-negative one,
and any use of it is that choice.

## Where a fractional dimension goes

`x^q` has dimension `q • dim x`, and that is `Azoth.Dim.exponents_smul`: the power scales
every exponent by `q`.

**A unit with a fractional dimension is not in the vocabulary.**
`Azoth.Dim.not_integer_smul` proves it for any rational that is no integer: `q • L` is not
`n • L` for any `n : ℤ`, so no `uom` type can name it and `specs/schema/vocabulary.schema.json`
is right to require integers. The live case is `A_phi`, the Debye–Hückel parameter, whose
dimension is `(mol/kg)^-1/2` — a half. It is a bare `f64` for that reason, and the theorem
is what says so rather than a reviewer.

## What this page does not do

It does not make the kernels total. A guard is a claim about a region, and the claim can be
wrong — `powf`'s guard at `x > 1e-12` in NeqSim's Pitzer `g` function sits seven decades
below where the expression `1 - (1 + x)e^-x` stops returning a meaningful value, because
the numerator is `x^2/2` computed as a cancellation of two numbers near one. Measured,
`g(1e-8) = 2.22` and `g(1e-9) = 0.0` against a true value of `1`.

**A guard is only as good as the arithmetic it protects**, which is the reason the policy
above is stated per call site rather than as a wrapper: the wrapper cannot know.
