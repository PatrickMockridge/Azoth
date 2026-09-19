/-
The raw-versus-normalised conversion, as `docs/src/calculus/normalisation.md` states it.

An equation of state is written in mole fractions and a flash is written in mole
numbers. The two are the same function of two different arguments - the molar
residual Helmholtz energy `f(x)`, and the extensive `n * f(n / sum n)` - and their
derivatives do not line up: the extensive one carries the molar value and subtracts
the composition-weighted sum of the molar derivatives.

That subtraction is the whole content of this file, and it is here because getting it
wrong is silent. It was written out by hand twice in one session of this port, once in
a kernel's `ln phi` and once in a root sensitivity, and in both cases the code agreed
with the reference it was compared against, because the same omission was in both.
A theorem cannot agree with a reference that has the same bug.
-/

import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.FDeriv.Add
import Mathlib.Analysis.Calculus.FDeriv.Comp
import Mathlib.Analysis.Calculus.FDeriv.Const
import Mathlib.Analysis.Calculus.FDeriv.Mul
import Mathlib.Analysis.Calculus.FDeriv.Basic

open scoped BigOperators

namespace Azoth
namespace Normalisation

variable {n : ℕ}

/-- The sum functional, `m ↦ sum_i m_i`, as a continuous linear map.

The derivative of the total mole number, and the `1` in the `x (x) 1` that the
conversion's second term carries. -/
noncomputable def sumCLM : ((Fin n → ℝ) →L[ℝ] ℝ) := ∑ i, ContinuousLinearMap.proj i

theorem sumCLM_apply (v : Fin n → ℝ) : sumCLM v = ∑ i, v i := by
  simp [sumCLM, ContinuousLinearMap.sum_apply, ContinuousLinearMap.proj_apply]

/-- The molar function recovered as an extensive one: `Phi(m) = (sum m) * f(m / sum m)`.

`Phi` is the residual Helmholtz energy of a mixture as a function of its mole numbers,
built from the molar function of its mole fractions. It is degree-one homogeneous, and
its derivative is what a flash wants while `f`'s is what an equation of state gives. -/
noncomputable def extensive (f : (Fin n → ℝ) → ℝ) : (Fin n → ℝ) → ℝ :=
  fun m => (∑ i, m i) * f (fun i => m i / (∑ j, m j))

/-- **The conversion.** For a molar function differentiable at the molar composition,
the extensive function built from it has derivative
`M + (f(x) - M x) (x) sum`, where `M` is `f`'s derivative and `x` is the normalised
composition.

Read on a component it says `d Phi/d n_j = f(x) + d f/d x_j - sum_k x_k d f/d x_k`:
the extensive derivative carries the molar *value* and removes the composition-weighted
sum of the molar derivatives. Dropping either half is the error this file exists for,
and it is silent - a kernel that omits the sum agrees with a reference that omits it
too, and both agree with nothing. -/
theorem hasFDerivAt_extensive (f : (Fin n → ℝ) → ℝ) (M : (Fin n → ℝ) →L[ℝ] ℝ)
    (x : Fin n → ℝ) (hx : ∑ i, x i ≠ 0)
    (hf : HasFDerivAt f M (fun i => x i / (∑ j, x j))) :
    HasFDerivAt (extensive f)
      (M + (f (fun i => x i / (∑ j, x j)) - M (fun i => x i / (∑ j, x j))) • sumCLM) x := by
  -- `m ↦ sum_i m_i` is the sum of the coordinate projections.
  have hsum : HasFDerivAt (fun m : Fin n → ℝ => ∑ i, m i) (sumCLM (n := n)) x := by
    have h : HasFDerivAt (sumCLM (n := n) : (Fin n → ℝ) → ℝ) (sumCLM (n := n)) x :=
      ContinuousLinearMap.hasFDerivAt (sumCLM (n := n))
    have hfun : (⇑(sumCLM (n := n)) : (Fin n → ℝ) → ℝ) = fun m => ∑ i, m i := by
      funext m
      simp [sumCLM, ContinuousLinearMap.sum_apply, ContinuousLinearMap.proj_apply]
    simpa only [hfun] using h
  -- `m ↦ (sum m)⁻¹`, and then the normalised composition `m ↦ (sum m)⁻¹ • m`.
  have hinv : HasFDerivAt (fun m : Fin n → ℝ => (∑ i, m i)⁻¹)
      ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (-((∑ i, x i) ^ 2)⁻¹)).comp
            (sumCLM (n := n))) x :=
    HasFDerivAt.comp (f := fun m : Fin n → ℝ => ∑ i, m i)
      (f' := sumCLM (n := n)) (x := x)
      (hasFDerivAt_inv (x := ∑ i, x i) hx) hsum
  have hnorm : HasFDerivAt (fun m : Fin n → ℝ => (∑ i, m i)⁻¹ • m)
      (((∑ i, x i)⁻¹ : ℝ) • (1 : (Fin n → ℝ) →L[ℝ] (Fin n → ℝ))
        + ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (-((∑ i, x i) ^ 2)⁻¹)).comp
            sumCLM).smulRight x) x :=
    hinv.smul (hasFDerivAt_id x)
  -- The derivative is read at the normalised composition, which is `(sum x)⁻¹ • x`.
  have hpoint : ∀ y : Fin n → ℝ, (fun i => y i / (∑ j, y j)) = (∑ i, y i)⁻¹ • y := by
    intro y
    funext i
    simp [div_eq_mul_inv, smul_eq_mul, mul_comm]
  have hf' : HasFDerivAt f M ((∑ i, x i)⁻¹ • x) := by
    simpa only [hpoint x] using hf
  have hcomp : HasFDerivAt (fun m : Fin n → ℝ => f ((∑ i, m i)⁻¹ • m))
      (M.comp (((∑ i, x i)⁻¹ : ℝ) • (1 : (Fin n → ℝ) →L[ℝ] (Fin n → ℝ))
        + ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (-((∑ i, x i) ^ 2)⁻¹)).comp
            sumCLM).smulRight x)) x :=
    HasFDerivAt.comp (f := fun m : Fin n → ℝ => (∑ i, m i)⁻¹ • m)
      (f' := ((∑ i, x i)⁻¹ : ℝ) • (1 : (Fin n → ℝ) →L[ℝ] (Fin n → ℝ))
        + ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (-((∑ i, x i) ^ 2)⁻¹)).comp
            (sumCLM (n := n))).smulRight x) (x := x) hf' hnorm
  have hprod : HasFDerivAt (fun m : Fin n → ℝ => (∑ i, m i) * f ((∑ i, m i)⁻¹ • m))
      (((∑ i, x i) : ℝ) • (M.comp (((∑ i, x i)⁻¹ : ℝ) • (1 : (Fin n → ℝ) →L[ℝ] (Fin n → ℝ))
        + ((ContinuousLinearMap.smulRight (1 : ℝ →L[ℝ] ℝ) (-((∑ i, x i) ^ 2)⁻¹)).comp
            sumCLM).smulRight x)) + (f ((∑ i, x i)⁻¹ • x)) • sumCLM) x := by
    simpa only [Pi.mul_apply] using hsum.mul hcomp
  -- The extensive function is that product, and the derivative simplifies to the claim.
  have hsame : (fun m : Fin n → ℝ => (∑ i, m i) * f ((∑ i, m i)⁻¹ • m)) = extensive f := by
    funext m
    simp only [extensive]
    rw [hpoint m]
  rw [hsame] at hprod
  convert hprod using 1
  ext v
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply,
    ContinuousLinearMap.comp_apply, ContinuousLinearMap.smulRight_apply,
    ContinuousLinearMap.one_apply, sumCLM_apply, map_add, map_smul]
  -- One form for the normalised composition, so `f` and `M` are applied to the same
  -- term on both sides before any arithmetic starts.
  rw [hpoint x]
  -- `simp` expands `M` over the sum and the scalar on the right and leaves the one
  -- occurrence on the left alone; that occurrence is the whole difference between the
  -- sides, and `x` is what the right side has.
  rw [map_smul]
  -- Everything is a real now, so `•` is multiplication and what is left is arithmetic
  -- whose only fact beyond `ring` is `(sum x) * (sum x)⁻¹ = 1`.
  simp only [smul_eq_mul]
  field_simp
  ring

/-- A linear functional is its values on the coordinate basis: `M y = sum_k y_k M e_k`.

`M x` in the conversion is the composition-weighted sum of the molar partials, and that
is this expansion at `y = x`, so the corollary below can be read as the textbook one. -/
theorem clm_apply_eq_sum_single (M : (Fin n → ℝ) →L[ℝ] ℝ) (y : Fin n → ℝ) :
    M y = ∑ k, y k * M (Pi.single k (1 : ℝ) : Fin n → ℝ) := by
  have hbasis : y = ∑ k, y k • (Pi.single k (1 : ℝ) : Fin n → ℝ) := by
    funext i
    simp [Finset.sum_apply, Pi.single_apply]
  conv_lhs => rw [hbasis]
  simp [map_sum, map_smul, smul_eq_mul]

/-- **The conversion on one component**, which is how it is written in code:
`d Phi/d n_j = f(x) + d f/d x_j - sum_k x_k d f/d x_k`.

The middle term is the molar partial at the normalised composition and the last is the
composition-weighted sum of all of them, so the extensive derivative is the molar one
plus the value minus its own weighted sum - and dropping the last term is the silent
error this file exists for. -/
theorem fderiv_extensive_single (f : (Fin n → ℝ) → ℝ) (M : (Fin n → ℝ) →L[ℝ] ℝ)
    (x : Fin n → ℝ) (hx : ∑ i, x i ≠ 0) (j : Fin n)
    (hf : HasFDerivAt f M (fun i => x i / (∑ j, x j))) :
    fderiv ℝ (extensive f) x (Pi.single j (1 : ℝ) : Fin n → ℝ)
      = M (Pi.single j (1 : ℝ) : Fin n → ℝ)
        + (f (fun i => x i / (∑ j, x j)) - M (fun i => x i / (∑ j, x j))) := by
  rw [(hasFDerivAt_extensive f M x hx hf).fderiv]
  rw [clm_apply_eq_sum_single]
  simp [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, sumCLM_apply,
    Pi.single_apply, Finset.sum_ite_eq', Finset.mem_univ]

end Normalisation
end Azoth
