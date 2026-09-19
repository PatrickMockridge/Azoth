/-
The sensitivity of a solution, as `docs/src/calculus/implicit.md` states it.

A flash solves `F (X p) p = 0` for the state `X` at the parameters `p`, and then wants
the *derivative* of that solution - how the state moves when a parameter does. The
derivative is not free: it is `-J⁻¹ F_p`, the inverse of the Jacobian composed with the
parameter derivative, and both halves are easy to get wrong. A sign is easy to drop and
an inverse is easy to replace with a division by something that is not a number.

`Azoth.Normalisation` fixes the raw-versus-normalised conversion, which was made twice in
this project; this fixes the sensitivity, which is the other formula a kernel writes by
hand and cannot check against a reference that made the same mistake.
-/

import Mathlib.Analysis.Calculus.FDeriv.Add
import Mathlib.Analysis.Calculus.FDeriv.Comp
import Mathlib.Analysis.Calculus.FDeriv.Const
import Mathlib.Analysis.Calculus.FDeriv.Bilinear
import Mathlib.Analysis.Calculus.FDeriv.Prod

open scoped BigOperators

namespace Azoth
namespace Implicit

variable {E P : Type*}
variable [NormedAddCommGroup E] [NormedSpace ℝ E]
variable [NormedAddCommGroup P] [NormedSpace ℝ P]

/-- The two projections out of a product, as continuous linear maps. -/
noncomputable def fstCLM : (E × P) →L[ℝ] E := ContinuousLinearMap.fst ℝ E P

noncomputable def sndCLM : (E × P) →L[ℝ] P := ContinuousLinearMap.snd ℝ E P

/-- **A solution's sensitivity is minus the inverse Jacobian times the parameter
derivative.**

`F` is the residual, `X` the solution, `p` the parameters. The hypothesis `hF` says `X`
solves it everywhere; `hpair` says the residual's derivative at `(X p, p)`, read as a
function of the *pair*, splits into a state part and a parameter part - which is what
`J` and `Fp` are; and `hJ` says the state part is invertible, which is the only reason
`J⁻¹` appears at all.

The conclusion is the whole content: **the sign is negative and the inverse is a real
inverse**. -/
theorem hasFDerivAt_of_solution
    (F : E → P → E) (X : P → E) (p : P)
    (J : E ≃L[ℝ] E) (Fp : P →L[ℝ] E) (X' : P →L[ℝ] E)
    (hF : ∀ q, F (X q) q = 0)
    (hX : HasFDerivAt X X' p)
    (hpair : HasFDerivAt (fun z : E × P => F z.1 z.2)
      ((J : E →L[ℝ] E).comp (fstCLM : (E × P) →L[ℝ] E)
        + Fp.comp (sndCLM : (E × P) →L[ℝ] P)) (X p, p)) :
    X' = -((J.symm : E →L[ℝ] E).comp Fp) := by
  -- The composite `q ↦ F (X q) q` is the constant zero, so its derivative is zero.
  have hzero : fderiv ℝ ((fun z : E × P => F z.1 z.2) ∘ fun q => (X q, q)) p = 0 := by
    rw [show ((fun z : E × P => F z.1 z.2) ∘ fun q => (X q, q)) = fun _ => (0 : E) by
      funext q
      exact hF q,
      fderiv_const_apply]
  have hpair' : HasFDerivAt (fun q => (X q, q)) (X'.prod (ContinuousLinearMap.id ℝ P)) p :=
    hX.prodMk (hasFDerivAt_id p)
  have hcomp := HasFDerivAt.comp (x := p) (f := fun q => (X q, q)) hpair hpair'
  have hderiv : (((J : E →L[ℝ] E).comp fstCLM + Fp.comp sndCLM).comp
      (X'.prod (ContinuousLinearMap.id ℝ P))) = 0 := by
    rw [← hcomp.fderiv, hzero]
  -- Read that identity on the two components: `J X' + Fp = 0`.
  have hsplit : (J : E →L[ℝ] E).comp X' = -Fp := by
    ext v
    have h := congrArg (fun L : P →L[ℝ] E => L v) hderiv
    simp only [ContinuousLinearMap.comp_apply, ContinuousLinearMap.prod_apply,
      ContinuousLinearMap.add_apply, ContinuousLinearMap.id_apply, fstCLM,
      sndCLM] at h
    exact eq_neg_of_add_eq_zero_left h
  -- Compose with the inverse on the left, and the inverse cancels.
  have hfin : X' = -((J.symm : E →L[ℝ] E).comp Fp) := by
    ext v
    have h := congrArg (fun L : P →L[ℝ] E => (J.symm : E →L[ℝ] E) (L v)) hsplit
    simpa [ContinuousLinearMap.comp_apply, map_neg, ContinuousLinearEquiv.symm_apply_apply]
      using h
  exact hfin

end Implicit
end Azoth
