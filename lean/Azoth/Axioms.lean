/-
Every theorem this development claims, so that a proof which quietly acquired a
`sorry` - here or in a vendored dependency - fails the build rather than being
trusted.

`tools/check_lean_axioms.py` runs this file through `lake env lean` and refuses any
axiom set outside `propext`, `Classical.choice` and `Quot.sound`. It reports the
*transitive* axiom set of each proof term, so a `sorry` anywhere in the dependency
chain shows up as `sorryAx` - including one inside `lean/vendor/lean-units/`, which
a scan over `Azoth/` cannot see. It also catches `admit`, a hand-applied `sorryAx`
and an `axiom` declaration, none of which a text search finds.

That is a completeness gate and not a correctness one: it proves a proof has no
gap, not that it is the theorem a reader expects. The guard for the second is
`python/tests/test_lean_claims.py`, which holds the theorem names below to the ones
`docs/src/calculus/` states.

Add a line here in the same commit as the theorem it names, and add the theorem to
that page too. A theorem this file does not name is a theorem nothing gates.
-/

import Azoth.Dim
import Azoth.Vocabulary
import Azoth.Rho
import Azoth.Barb
import Azoth.Capability
import Azoth.Implicit
import Azoth.Normalisation
import Azoth.Pow

#print axioms Azoth.Dim.ofExponentsOn_nil
#print axioms Azoth.Dim.ofExponents_nil
#print axioms Azoth.Dim.ofExponentsOn_singleton
#print axioms Azoth.Dim.exponents_ofExponents
#print axioms Azoth.Dim.exponents_smul
#print axioms Azoth.Dim.not_integer_smul
#print axioms Azoth.Pow.rpow_natCast
#print axioms Azoth.Pow.rpow_third_of_neg_pos
#print axioms Azoth.Pow.rpow_two_thirds_of_neg_neg
#print axioms Azoth.Pow.the_conventions_disagree
#print axioms Azoth.Rho.quote_injective
#print axioms Azoth.Rho.drop_injective
#print axioms Azoth.Rho.round_trip
#print axioms Azoth.Rho.out_injective
#print axioms Azoth.Rho.in_injective
#print axioms Azoth.Rho.subst0_var_zero
#print axioms Azoth.Rho.subst0_var_succ
#print axioms Azoth.Rho.subst0_under_binder
#print axioms Azoth.Rho.comm_reduction
#print axioms Azoth.Rho.round_trip_congr
#print axioms Azoth.Rho.round_trip_par
#print axioms Azoth.Barb.reduces_star_refl
#print axioms Azoth.Barb.reduces_star_trans
#print axioms Azoth.Barb.barbed_bisim_refl
#print axioms Azoth.Barb.barbed_bisim_symm
#print axioms Azoth.Barb.barbed_bisim_trans
#print axioms Azoth.Barb.barb_par
#print axioms Azoth.Barb.barb_nil
#print axioms Azoth.Barb.barb_drop

-- The keycard as a capability, `docs/src/calculus/capability.md`'s three claims and the
-- witness that keeps the second from being vacuous. Naming the witness theorems here as
-- well as the claims is deliberate: a gate whose only gated theorems are the ones it
-- satisfies is a gate that cannot fail.
#print axioms Azoth.Capability.run_deterministic
#print axioms Azoth.Capability.run_exists_unique
#print axioms Azoth.Capability.run_derives_in_grant
#print axioms Azoth.Capability.witness_datum_is_outside
#print axioms Azoth.Capability.the_gate_refuses
#print axioms Azoth.Capability.the_gate_can_succeed

-- The raw-versus-normalised conversion, `docs/src/calculus/normalisation.md`'s claim. The
-- basis expansion is named here as well as the conversion: it is what makes the `M x` in
-- the claim the composition-weighted sum of molar partials a reader takes it for.
#print axioms Azoth.Normalisation.hasFDerivAt_extensive
#print axioms Azoth.Normalisation.fderiv_extensive_single
#print axioms Azoth.Normalisation.clm_apply_eq_sum_single
#print axioms Azoth.Normalisation.sumCLM_apply

-- The sensitivity of a solution, `docs/src/calculus/implicit.md`'s first claim. The second
-- is stated on the page and not proved here; a page that said otherwise would be the
-- failure this gate exists to catch.
#print axioms Azoth.Implicit.hasFDerivAt_of_solution
