/-
Every theorem this development claims that is not a unit's dimension or a guard. Those two
groups are gated by the generated `Azoth/Gate.lean` and `Azoth/GuardGate.lean`, one line per
theorem, because a hand-maintained list of names goes stale the first time one is added. This
file is hand-written and carries the rest, so that a proof which quietly acquired a `sorry` -
here or in a vendored dependency - fails the build rather than being trusted.

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
that page too; the two generated gates add their own. A theorem no gate file names is
a theorem nothing gates.
-/

import Azoth.Dim
import Azoth.Vocabulary
import Azoth.Rho
import Azoth.Barb
import Azoth.Process
import Azoth.Capability
import Azoth.Implicit
import Azoth.Normalisation
import Azoth.Pow
import Azoth.View
import Azoth.Graph
import Azoth.Session

#print axioms Azoth.Dim.ofExponentsOn_nil
#print axioms Azoth.Dim.ofExponents_nil
#print axioms Azoth.Dim.ofExponentsOn_singleton
#print axioms Azoth.Dim.exponents_ofExponents
#print axioms Azoth.Dim.exponents_smul
#print axioms Azoth.Dim.not_integer_smul
#print axioms Azoth.Pow.rpow_natCast
#print axioms Azoth.Pow.even_power_has_two_roots
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

-- The process layer's adequacy claim, `docs/src/calculus/process.md`'s second claim. Its
-- statement is settled by the design pass recorded in `Azoth/Process.lean`: with a barb that
-- records a channel and not a magnitude, the claim's provable half is that the *declaration*
-- determines the observable behaviour, and the witness below is what keeps that from being
-- reflexivity in disguise.
#print axioms Azoth.Process.barb_declaration
#print axioms Azoth.Process.unit_op_is_extensional
#print axioms Azoth.Process.swapped_declarations_are_congruent_not_equal

-- The rendering, `docs/src/calculus/view.md`'s claims. Three groups: the editor's ticket
-- discipline (`keep` and its fold), the gesture that is the one piece of local state, and the
-- display unit's conversion. `keep_keeps_the_newer` and `one_point_does_not_determine_a_shift`
-- are the witnesses beside them - the first says `keep` is not a projection, the second says a
-- single point cannot tell a scale from a shifted one, which is why the units test asks twice.
#print axioms Azoth.View.keep_ticket
#print axioms Azoth.View.keep_mem
#print axioms Azoth.View.keep_keeps_newest
#print axioms Azoth.View.foldr_keep_ticket
#print axioms Azoth.View.foldr_keep_mem
#print axioms Azoth.View.keep_keeps_the_newer
#print axioms Azoth.View.drawn_after_release
#print axioms Azoth.View.a_drag_moves_one_node
#print axioms Azoth.View.a_drag_in_flight_is_not_the_document
#print axioms Azoth.View.Shift.fromSi_toSi
#print axioms Azoth.View.Shift.toSi_fromSi
#print axioms Azoth.View.celsius_is_not_kelvin
#print axioms Azoth.View.one_point_does_not_determine_a_shift
#print axioms Azoth.View.two_points_determine_a_shift

-- The connection judgment, `docs/src/calculus/session.md`'s first claim. `check_sound` is the
-- claim; `check_refuses_a_missing_field` and `compatibility_is_not_equality` are the witnesses
-- beside it - the first says the judgment can refuse, so soundness is not a sentence about a
-- function that always accepts, and the second says two records that differ in order pass and
-- are not equal, so it is not reflexivity in disguise.
#print axioms Azoth.Graph.check_refl
#print axioms Azoth.Graph.check_sound
#print axioms Azoth.Graph.check_refuses_a_missing_field
#print axioms Azoth.Graph.compatibility_is_not_equality

-- The session, `docs/src/calculus/session.md`'s other three. The door table is data and the
-- lemma reads it; the two witnesses are that its rows differ in their status and that an edit
-- really does break the invariant `dirty` marks.
#print axioms Azoth.Session.every_door_uses_the_shared_transition
#print axioms Azoth.Session.the_doors_do_not_agree_on_a_status
#print axioms Azoth.Session.a_run_restores_it
#print axioms Azoth.Session.an_edit_leaves_the_values_behind
#print axioms Azoth.Session.a_failed_run_keeps_none
#print axioms Azoth.Session.reading_leaves_the_state_alone
