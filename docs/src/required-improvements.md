# Required improvements

What is wrong with this library, with the evidence and the cause. Feedstock for the next
remediation plan.

Every entry has the same shape, so the page can be read as a set rather than as prose:

```
### <title>            severity: blocking | major | minor
Status:     open | fixed at <commit>
Symptom     what is wrong, in one or two sentences
Evidence    the exact reproduction and the measured numbers
Root cause  as determined by the root-cause pass and tested by an adversarial reviewer
Remediation one sentence, enough to size it
```

**How a root cause gets here.** One pass determines a cause per entry from evidence in the
tree; a second reader is then tasked with *refuting* each one. A cause that survives is
recorded as determined. A refuted cause is re-derived. A cause neither can settle is
recorded as unresolved **and names the experiment that would settle it**. Where the two
readings genuinely differ, both are recorded and the entry says so.

An entry leaves its severity section when it is fixed, and the fix carries a test. It moves to
"Fixed" rather than being deleted: a defect that was real and is closed is a different fact
from one that was never recorded.

See [Test plan](./test-plan.md) for what is tested and how a run is performed.

---

## Blocking

### `eos.ph_flash` returns a wrong temperature over most of its range

**Severity: blocking. Status: open.**

**Symptom.** `eos.ph_flash` inverts an enthalpy by bisecting the temperature. It declares
convergence on the width of the *temperature* bracket rather than on the enthalpy residual,
so where the function it inverts is discontinuous the bracket collapses onto the jump and
the call **returns normally** with a temperature that does not have the requested enthalpy.
The discontinuity itself comes from `eos.pt_flash`, which returns a different compressibility
root either side of a threshold while the feed is single-phase on both sides.

**Evidence.** Methane/n-butane, `z = [0.6, 0.4]`, 20 bar. Sweeping `H` from 800 to 1900 J/mol
in steps of 10 — **105 of 111 points return the same temperature, 386.6944 K**, with residuals
from 3.1e-3 to 1.275:

```
H = 1277.0   ->  T = 386.6944 K   residual = 4.6109e-01
H = 1800.0   ->  T = 386.6944 K   residual = 3.656e-02
H = 1870.0   ->  T = 386.8115 K   residual = 2.272e-10     (correct)
```

0.461 against the declared `algorithm.tolerance` of 1e-8 is **7.66 orders** outside, and
nothing raises. `iterations = 19` of a 200 budget shows the loop stopped on the bracket test,
not on the cap, so the non-convergence error at `flash_property.rs:268` is never reached.

**Root cause — two layers, and the second is upstream of the first.**

*Layer 1, where the wrong number is produced:* `crates/azoth-eos/src/flash_property.rs:256`
stops on `(hi - lo) <= algorithm.tolerance * mid`. The comment at `:252-255` justifies it —
"the bracket is narrowed on the *temperature*, not on the property … the temperature interval
is well behaved" — and that reasoning is sound **for a continuous `H(T)`**, where the bracket
width bounds the residual through `|H(mid) − target| ≈ |H′|·(hi−lo)/2`. The residual is then
computed at `:286` and returned, but never tested against the tolerance.

*Layer 2, where the discontinuity comes from:* `eos.pt_flash` at 20 bar switches phase
behaviour at a threshold — the temperature where Wilson's K for n-butane crosses 1, so
`rachford_rice_bounds` returns `None`:

| T (K) | phase | beta | `z_vapour` |
|---|---|---|---|
| 385.0 | all_vapour | **26.57** | 0.6726 |
| 386.5 | all_vapour | **92.96** | 0.6636 |
| 387.0 | all_vapour | None | **0.9189** |

Two separate defects are visible in that table. `z_vapour` flips 0.6636 → 0.9189 between two
temperatures where the feed is single-phase **on both sides**, so `H(T)` as
`flash_property.rs:79-84` builds it is discontinuous — and the jump it spans, ~814 to ~1876
J/mol, is the whole range that returns 386.6944. And **`beta` is reported as 92.96 on an
`all_vapour` phase**: a vapour fraction outside `[0, 1]`, `Some` at 386.5 K and `None` at
387.0 K, with no error and no warning. Beta is read at `flash_property.rs:78` as
`flash.beta.unwrap_or(0.0)`.

**Status of the cause.** Recorded as determined for the behaviour — the 105-of-111 measurement
is reproduced. The *mechanism* is recorded as two layers rather than one because a first
reading blamed only the stopping criterion, and that reading does not survive: fixing the
stopping criterion converts 105 silent wrong answers into 105 `SolverNotConverged` errors and
still has no root to find. The upstream layer has to be fixed for the function to be
invertible at all.

**Why 13 passing tests do not see it.** `python/tests/models/test_ph_flash.py` passes: both of
its spec cases round-trip through the same function, so a wrong `H` agrees with itself. The
defect was found through `process.heater`, not through the flash's own tests.

**Remediation.** Three parts, in this order: `pt_flash` must not hand back a root that jumps
between single-phase states; `beta` must not be reported outside `[0, 1]`; and `solve_temperature`
must test the residual against `algorithm.tolerance` before declaring convergence, so that a
case like this fails loudly instead of quietly. The first is the real fix and needs a decision
on what the flash should return across the Wilson threshold — which is a modelling choice, not
a patch.

---

## Major

### `crates/azoth-python` has no tests

**Severity: major. Status: open.**

**Symptom.** The PyO3 binding layer — every bridge function, every result conversion, the
error mapping — is exercised only by Python calling into it. There is no `tests/` directory
and no `#[cfg(test)]` anywhere in the crate.

**Evidence.** `grep -rn '#\[test\]' crates/azoth-python/src/` returns nothing;
`crates/azoth-python/tests/` does not exist.

**Root cause.** The binding is tested from the only side that can call it, so a defect that
is symmetric across the boundary — a conversion applied in both directions, a field mapped
to the wrong name on both sides — passes. Python's introspection tests
(`test_result_shapes_agree_across_languages`) compare Python's fields to Rust's *reported*
fields, which is exactly the check a self-consistent binding passes.

**Remediation.** Rust-side tests that build a result and assert the Python-visible shape
without going through Python.

### `azoth-cli`'s `main.rs` and `report.rs` are untested

**Severity: major. Status: open.**

**Symptom.** `crates/azoth-cli/tests/pipe.rs` exercises the `pipe` module. `main.rs`
(argument parsing) and `report.rs` (output formatting) are covered by nothing.

**Evidence.** Twelve tests in `tests/pipe.rs`, none naming either file.

**Root cause.** The CLI was built with the computational path first, and the compute path is
what has tests. Argument parsing and report formatting are the parts a user meets and the
parts no numerical test can see.

**Remediation.** Tests over `main`'s argument handling and `report`'s formatting, plus one
end-to-end run of the built binary.

### A keycard's `fittings` and `fluids` replace the shipped file instead of extending it

**Severity: major. Status: open.**

**Symptom.** `tools/gen_user_data.py` writes `data/fittings/crane_k_factors.csv` from the
card's `fittings` section and nothing else, so a card carrying one fitting compiles to a
one-row file and the seven shipped rows are gone. The same holds for `fluids`. But the
format's own text says the opposite in two places: `specs/schema/keycard.schema.json`
describes both sections as "overriding or **extending**", and `docs/src/keycard.md`'s
section table says "Overrides" beside a paragraph explaining that everything is by name.

**Evidence.** `tools/gen_user_data.py:337-356` — `plan()` calls `render_fittings(document
["fittings"])`, which renders the card's rows and no others. The file is then written in
place, so the shipped rows are not merged with the card's, they are replaced by them.

**Root cause.** The section was specified as an overlay and implemented as a whole-file
substitution. All six other sections are overrides applied *per name* at the point of
lookup, so nothing in the design anticipated a section whose unit is the file rather than
the row.

**Remediation.** Merge by `id` for fittings and by name for fluids, so a one-row card
yields seven rows and the six shipped placeholders keep their `estimated_dummy` status.
The banner logic already handles a mixed file, and `Fitting::is_estimated()` is per row, so
the placeholder warning survives the merge. A test asserts that a card with one fitting
compiles to the shipped count.

### 42% of the shipped `kij` table is zeros

**Severity: major. Status: open.**

**Symptom.** `data/components/kij.csv` ships 516 rows and **217 of them carry exactly
`0.0`**. A zero is the ideal-mixture default — the value the mixing rule uses when a pair
is absent — so those rows state "no interaction parameter" in the same form as a fitted
one. Nothing in the file distinguishes a pair someone fitted to zero from a pair nobody
fitted at all.

**Evidence.** Measured over NeqSim's `INTER.csv` restricted to the 173 vendored
components: 516 pairs, 299 with a non-zero `KIJPR`. `tools/gen_databank.py:171-173`
skips a row whose `KIJPR` is empty and writes one whose text is `"0"`, so all 217 come
through. The magnitudes that are real span −0.276 to 1.0.

**Root cause.** Zero is a legitimate value — `eos.components.kij_for` documents that
overriding a fitted pair *back* to zero is a deliberate act — so the generator is right
not to drop it. The defect is one level up: the file has no way to say whether a value
was fitted, so it cannot distinguish the two cases it is storing identically.

**Remediation.** Drop the zero rows. It is a no-op behaviourally: `kij_for` records only
non-zero pairs, and `mixture()` already defaults an absent pair to zero, so the compiled
answers are identical and the row count starts meaning something. A test asserting no
shipped row is zero is what would keep it true. Changing the file is a separate change
set from the vendoring audit, which is why this is recorded rather than fixed here.

### The citations across the registry are unconfirmed

**Severity: major. Status: open.**

**Symptom.** Nearly every spec records the same thing: the equation is not in doubt, and the
*citation* is. Nobody has opened the source and confirmed the equation number, the page, or
the form the relation is printed in. `source.equation` is left unstated rather than guessed
at, which is the right call and leaves the claim unverified.

**Evidence.** The report appears in eighteen specs: `eos.pr_kappa`, `eos.pr_alpha_ab`,
`eos.pr_departure`, `eos.prsv_kappa`, `eos.vdw1f_mix_binary`, `eos.rachford_rice_binary`,
`eos.critical_point`, `eos.stability_test`, `hydraulics.darcy_weisbach`,
`hydraulics.reynolds_number`, `hydraulics.friction_factor_haaland`,
`hydraulics.friction_factor_swamee_jain`, `hydraulics.orifice_flow`, `hydraulics.pump_power`,
`hydraulics.choked_flow_area`, `hydraulics.crane_k_factors`, `thermal.conduction_plane_wall`
and `process.heater`.

**Root cause.** The sources are paywalled standards and papers and this project has no
subscription to them, so the honest record was a gap rather than an equation number recalled
from memory.

**Remediation.** Not a software task. It needs one reader with access to the standards. It is
the largest body of unverified claims in the registry, and **consolidating it is what makes it
countable**: it was eighteen paragraphs each saying "this is unconfirmed", which reads as
diligence and is one fact stated eighteen times.

### Nothing checks that a separator's feed is stable before it is split

**Severity: major. Status: open.**

**Symptom.** `process.separator` splits a feed by a flash's answer without asking whether the
feed is stable as a single phase. An unstable feed — the case where a split is *required*
rather than merely permitted — is split by whatever the flash converged to.

**Evidence.** `specs/models/process/separator.yaml`, the assumptions list. The test that
answers the question already exists and is registered: `eos.stability_test`, Michelsen's
tangent-plane test.

**Root cause.** The port carried NeqSim's `Separator.run` (lines 674-782), which calls a flash
and not a stability test, so the gap was inherited rather than introduced.

**Remediation.** Either the separator runs `eos.stability_test` first, or its contract states
that an unstable feed is the caller's to detect and is not diagnosed here. Both are
defensible; leaving it unsaid is not.

### `process.pump` accepts a vapour at its inlet and returns a plausible pressure rise

**Severity: major. Status: open.**

**Symptom.** Nothing checks that the inlet is a liquid. A pump run on a vapour returns a
pressure rise and a power that look entirely reasonable.

**Evidence.** `specs/models/process/pump.yaml`; `python/src/azoth/process/reference/pump.py`.

**Root cause.** The port follows NeqSim's `Pump`, whose liquid assumption is implicit in the
equipment type rather than checked.

**Remediation.** A range check on the inlet phase, or an explicit statement that the phase is
the caller's responsibility.

### `process.pump` raises `SolverNotConvergedError` for ordinary states

**Severity: major. Status: open.**

**Symptom.** A pump run at ordinary conditions fails. 20 to 60 bar at 300 K raises; so does
40 to 50 bar on a leaner feed; at 60 bar the failing temperature is 410 K.

**Evidence.** Measured on the shipped pump. One temperature on the isentropic search's scan
fails to converge, and one failing scan point is fatal to the whole model.

**Root cause.** `eos.ps_flash`'s scan skips a temperature with no admissible root but does not
skip one where the flash fails to settle. Those are different conditions and only the first
was considered.

**Remediation.** A decision, not a patch. Skipping loses the sign change that brackets the
answer; aborting makes the model unusable at ordinary states. Decide, then make the scan's
rule match the decision.

### A separator's duty branch measures enthalpy from the feed's own state

**Severity: major. Status: open.**

**Symptom.** Where a duty is supplied, the model measures the enthalpy change from the feed's
own state. NeqSim reads `getEnthalpy()`, which is history-dependent: the value depends on the
path the stream was taken through, not only on where it is.

**Evidence.** `specs/models/process/separator.yaml`.

**Root cause.** A port of a state-function call from a library whose streams carry their
history. The two agree whenever the fluid is single-phase throughout, which is the common
case, so the divergence has no symptom in ordinary use.

**Remediation.** State the difference and the condition under which the two agree. A
documentation fix unless a case can be built where they diverge.

---

## Minor

### `UnverifiedCalculationError` is raised by nothing

**Severity: minor. Status: open.**

**Symptom.** The error is defined in the Rust enum, mapped across the binding and
re-exported from Python. No code path raises it, in either language.

**Evidence.** `grep -rn UnverifiedCalculationError crates/ python/src/` finds only its
definition, its `to_pyerr` mapping, the name table entry, and the re-exports
(`crates/azoth-python/src/errors.rs:36,99`, `python/src/azoth/core/errors.py:115`).

**Root cause.** The `source_needed` status that would have raised it was removed by a policy
decision (`azoth-source-status-policy`), which left the error with no producer and no
consumer.

**Remediation.** Either a code path raises it or it is deleted. Dead surface a caller can
catch and never will is worse than no surface.

### `InvalidInputError` is missing from the same-class-object identity test

**Severity: minor. Status: open.**

**Symptom.** `test_errors_are_the_same_class_object` checks three error classes. The binding
registers **eight**, so five are unchecked — and the real gap is that nothing asserts the
*set* of registered names at all.

**Evidence.** `python/tests/test_cross_impl.py:277` is a literal
`("OutOfRangeError", "UnknownFittingError", "SolverNotConvergedError")`.
`crates/azoth-python/src/errors.rs:88-105` `register()` binds `AzothError`,
`InvalidInputError`, `OutOfRangeError`, `PropertyUnavailableError`, `SolverNotConvergedError`,
`UnitMismatchError`, `UnknownFittingError` and `UnverifiedCalculationError`, each by
`errors.getattr(name)`, so identity holds by construction for anything it binds.

**Root cause.** The list was written when three classes crossed the boundary. The same
"a list that was complete when it was written, and a second thing that grew beside it"
failure `test_registry_contract.py` exists to catch.

**Remediation.** Derive the list from `register()` rather than naming three of them — the
useful assertion is that the names bound into the extension module are exactly the names
exported from `azoth.core.errors`, which covers all eight and cannot go stale.

### `gen_models.py` and `gen_databank.py` are checked by no test

**Severity: minor. Status: open.**

**Symptom.** A generator that fails for a spec shape not currently in the tree is caught by
nothing. `gen_databank.py` has no test **and no gate at all**.

**Evidence.** `gen_docs` and `gen_registry` are subprocessed with `--check` by
`test_registry_contract.py`; `gen_user_data` is loaded as a module by `test_gen_user_data.py`;
`spec_lint` is run by `test_spec_lint.py`. `gen_models`, `gen_stub` and `gen_databank` are
loaded or called by no test.

Gating is a separate axis and does not close the gap: the `docs-drift` job regenerates
`gen_docs`, `gen_registry`, `gen_models` and `gen_stub` and fails on any diff
(`.github/workflows/ci.yml:170-174`), so those three are at least *gated*. **`gen_databank`
appears in no workflow** (`grep -rn gen_databank .github/` is empty) — and cannot, because
`tools/gen_databank.py:214` requires a NeqSim checkout as an argument that CI has not got.

**Correction.** An earlier version of this entry claimed `test_registration_completeness.py`
invokes `gen_stub` (it names it only inside a failure message, `:258`) and that
`test_gen_user_data.py` calls `gen_registry` (it loads `gen_user_data`). Both were wrong, and
the entry had been used to "correct" the test plan's gap table, which was right to begin with.

**Root cause.** A `--check` proves the current tree regenerates identically. It says nothing
about a spec shape that is not in the tree, and no gate can cover a generator whose input CI
does not have. The uniform-conversion risk in `gen_databank` is partly covered another way —
by its own `check_round_trip` (`tools/gen_databank.py:228-241`) and by
`python/tests/eos/test_components.py`, which holds the generated table to published constants.

**Remediation.** A test per generator over a synthetic spec exercising a shape the tree does
not contain, for the three that need no external input. `gen_databank` needs a decision:
a recorded fixture standing in for the NeqSim tree, or documented as uncheckable.

### `eos.expander` and `process.pump` have one spec case each

**Severity: minor. Status: open.**

**Symptom.** Both declare a single worked case. Every other id has two or more.

**Evidence.** `grep -c '^  - id:'` is 1 for `specs/models/process/expander.yaml` and
`specs/models/process/pump.yaml`.

**Root cause.** Both were written to close the Pareto argument's coverage of the three
machines, and one case each was enough to demonstrate the shared procedure.

**Remediation.** A second case each, at a state that exercises a different branch — an
expander with `efficiency = 1`, or a pump whose inlet is at a different phase.

### `pt_flash.yaml` explains the residual with a claim that is false by ten orders

**Severity: minor. Status: open.**

**Symptom.** The spec explains why `residual` is not pinnable with a claim about magnitude
that is wrong by about ten orders of magnitude — and the explanation it gives is not the
reason, which matters because the real reason is the one that generalises.

**Evidence.** `specs/models/eos/pt_flash.yaml:148-155` (the "# Correction 4" block; the entry's
title named `ph_flash.yaml`, which contains no instance of the word) says two implementations
"differ in its last few *ulps*".

Measured on `pt_flash`/330 K: values `4.658240e-11` against `4.658231e-11`, absolute difference
`9.32e-17`, relative `2.0e-06`. The ulp of `4.658e-11` is `6.46e-27`, so **the difference is
1.4e10 ulps of the residual — ten orders, not "a few"**. It is ~0.4 ulp *of the order-1
quantity inside* the residual, which is the reason no relative tolerance can be met on the
field and which the text does not say.

**Root cause.** Written by analogy with a genuinely ulp-scale difference, without measuring.

**What the block gets right, and why it matters.** The same block already concludes that
"`residual` is not a value a case can pin". So the spec recorded the finding — in a `notes`
block that renders to the model's documentation page — and the cross-implementation harness
pinned it anyway, because `assert_results_equal` compares every field of the dataclass rather
than every output the spec declares. **A finding recorded where nothing reads it is the same
failure as one recorded in a test docstring**, and this one survived long enough to fail
fourteen tests.

**Remediation.** Correct the arithmetic and state the real reason. The harness half is fixed
(see "Fixed in this change set"); the text half is not.

### A separator routes the flash's `trivial` verdict to the gas outlet

**Severity: minor. Status: open.**

**Symptom.** A feed whose flash reports `trivial` is sent to the gas outlet. That is a
convention, not a checked result.

**Evidence.** `specs/models/process/separator.yaml`.

**Root cause.** `trivial` means the two trial phases converged to the same composition. It
carries no information about which outlet the stream belongs in.

**Remediation.** Route by the feed's own state, or say in the contract that `trivial` goes to
the gas outlet by convention.

### A separator models no entrainment carry-over

**Severity: minor. Status: open.**

**Symptom.** The separator is ideal: each phase goes wholly to its own outlet, with no liquid
carried into the gas or the reverse.

**Evidence.** `specs/models/process/separator.yaml`.

**Root cause.** Entrainment needs a correlation and a geometry. Neither is in the Pareto set,
and NeqSim's own entrainment models are out of scope by S8.

**Remediation.** A sentence in the contract. Nothing to implement.

### A unit operation inherits the flash's confidence and adds none

**Severity: minor. Status: open.**

**Symptom.** A separator's split is the flash's split, and nothing the separator does
independently verifies the flash. The heater's duty is the same: it is the flash's answer,
scaled.

**Evidence.** `specs/models/process/separator.yaml`, `specs/models/process/heater.yaml`.

**Root cause.** A unit operation is a flash call plus arithmetic — the Pareto argument for the
scope of this port — so the unit cannot be more trustworthy than the flash inside it.

**Remediation.** Nothing to fix. It is a statement about what a passing `process.separator`
test does and does not tell a reader, and it belongs in the contract where a reader meets it.

### The molar-versus-total enthalpy translation from NeqSim's joule sum

**Severity: minor. Status: open.**

**Symptom.** NeqSim sums *total* enthalpy in joules across a mixer's inlets. These models work
in molar enthalpy and multiply by the flow. The translation is exact only when the streams
share a molar mass.

**Evidence.** `specs/models/process/mixer.yaml`, `specs/models/process/heater.yaml`.

**Root cause.** NeqSim's streams expose a total-enthalpy accessor; azoth's results carry molar
properties, because a molar property is what the equations of state produce.

**Remediation.** State the translation and construct the case where it is not exact —
differing molar masses across a mixer's inlets is the one to try. If they diverge, this is a
wrong answer rather than a documentation gap and moves up a severity.

### A cubic equation of state is not trustworthy at a high-ratio compressor's outlet

**Severity: minor. Status: open.**

**Symptom.** `process.compressor`'s largest assumption is that Peng-Robinson describes the
outlet of a high-pressure-ratio machine. Nothing tests it.

**Evidence.** `specs/models/process/compressor.yaml`.

**Root cause.** A property of the equation of state rather than of the model. Testing it needs
a real machine's data, which this library does not ship and could not check.

**Remediation.** None in software. It belongs in the model's contract as a stated limit.

---

## Fixed

Entries that have left the sections above. Each names the commit that closed it and the test
that would have caught it. Four were closed at `ba98b99`, in the change set that first reached
a green baseline; `f8b5a70` closed the one that makes that baseline true rather than merely
claimed; and `01763f5` closed the prose linter's wiring.

**A note on evidence, for the first four.** They were uncommitted work when they were fixed,
so `ba98b99^` is not their pre-fix state and the failures cannot be re-derived from git. The
evidence for each is the captured run: the exact assertion messages and exit codes from the
session that made the change, quoted below. Where a fix is checkable from the tree — the
`critical_point` names, the mutation test — that is stated and is reproducible. The fifth is
not in that position: it is reproducible from wheel artifacts at any time, and the
reproduction is in its own entry below.

### `tools/prose_lint.py` was wired into nothing

**Severity: major while it stood. Status: fixed at `01763f5`.**

**Symptom.** P5 states the prose rules and names `tools/prose_lint.py` as their mechanism. The
tool appeared in no workflow and no hook, so it had never run in CI — and it could not pass on
the tree anyway, so wiring it in was not the one-line fix it looked like.

**Evidence.** Before the fix, `python tools/prose_lint.py` exited 1 with nine hits, **seven of
them the tool matching its own `HISTORY_PHRASES` list**: `sources()` included `("tools",
"*.py")` and nothing exempted the file doing the measuring. Its docstring claimed each phrase
"occurs zero times in this tree outside the change set", which was false for the same reason.
`grep -rn prose_lint .github/` returned nothing.

**Root cause.** Two, and the second is why the entry was not a formality. The tool was never
wired in — that is the symptom's cause. It was also never *run*, which is why nobody noticed it
reported itself: **a check that is never executed is a check that cannot be debugged.**

**Remediation, and what it cost.** Excluding the tool's own phrase list from its scan, then
adding `python tools/prose_lint.py` to the `spec-validate` job beside `spec_lint.py`. It
passes on this tree and fails the build on a planted violation.

**The part that is not fixed, and why it stays.** A phrase list cannot judge prose, which the
tool's own docstring concedes. It is a backstop, not the standard — and the phrases that
matter most (`deliberately`, `the whole point`) are exactly the ones it cannot catch, because
they are legitimate in some sentences and narration in others. Using the corpus as a licence
to exempt them is the scope-shrinking failure in tool form. **P2's second reader is the real
mechanism, and this tool does not substitute for it.**

### `assert_results_equal` compared a solver diagnostic on an impossible tolerance

**Severity: blocking while it stood (14 failing cases). Status: fixed at `ba98b99`.**

**Symptom.** Eleven cross-implementation cases failed on a comparison no two implementations
could satisfy, because the harness applied the case's *relative* answer tolerance to
`residual`.

**Evidence.** `python/tests/_helpers.py:210`. A residual is `|S − 1|` with `S ≈ 1`, so the
subtraction pins its absolute error at the rounding of `S` (~1e-16) however small the
residual becomes. Attainable relative accuracy is `2e-16/residual` and diverges as the
residual falls. No relative tolerance can be met, at any value.

**Root cause.** Every result field was compared the same way, and `residual` is not a result
— it is the solver's evidence that it converged. The bound that applies is the solver's own,
which the specs already declare (`pr_z_factor.yaml:284`, "a residual within 1e-12") and which
four model test files already assert per implementation. The cross-language comparison was
adding nothing.

**Fix and its test.** Diagnostic fields are compared on the model's declared
`algorithm.tolerance`, with NaN equal to NaN. A mutation test confirms it still fails on a
divergence beyond that bound, that a physical field perturbed still fails, and that one side
failing to converge while the other succeeds is reported.

### `_model_kwargs` fed a pure component's scalars to the mixture path

**Severity: major (2 failing cases). Status: fixed at `ba98b99`.**

**Symptom.** `eos.pure_saturation` cases raised `TypeError: 'float' object is not iterable`.

**Evidence.** `python/tests/test_cross_impl.py:119` iterated `inputs["Tc"]`, which
`pure_saturation` declares as a scalar.

**Root cause.** The branch tested `"Tc" in declared`, which is true of a mixture and of a
pure component alike. `kij` is what distinguishes them: it is the only declared input that
cannot be passed through, and a pure component has no pair to give.

**Fix and its test.** The branch keys on `kij`; the two `pure_saturation` cases pass.

### A matrix result fell through to exact equality

**Severity: major (1 failing case). Status: fixed at `ba98b99`.**

**Symptom.** `eos.stability_test::a_two_phase_feed_is_unstable` reported a 1.4e-15 rounding
difference in `w` as a divergence in the physics.

**Evidence.** `python/tests/_helpers.py:252-263` tested one nesting level
(`all(isinstance(item, (int, float)))`), so a matrix — a tuple of tuples — failed the vector
branch and reached `assert a == b`.

**Root cause.** The vector branch was written for a vector. A matrix output was added later
and the branch grew no depth, and the fall-through comparison was exact.

**Fix and its test.** A recursive numeric-nested comparison; the `stability_test` case passes.

### `eos.critical_point` declared outputs that name no field

**Severity: major. Status: fixed at `ba98b99`.**

**Symptom.** The spec declared `Tc`, `Pc`, `Vc` and `Z_c` — which name no field on
`CriticalPointResult` (`tc`, `pc`, `vc`, `z_c`) — and omitted `iterations` and `residual`
entirely. Two fields were compared without being specified anywhere.

**Evidence.** `specs/models/eos/critical_point.yaml` outputs against
`python/src/azoth/_core.pyi:98-105`. Nothing caught it:
`test_result_shapes_agree_across_languages` compares Python's fields to **Rust's** fields, and
`test_model_contract.py` checked only that the schema *requires* an `outputs` key.

**Root cause.** Two names are load-bearing for different consumers, and this model conflated
them. `valid_range[].quantity` is resolved by `apply_checks` (`crates/azoth-core/src/range.rs:290`)
and is emitted from the spec's `valid_range` block (`tools/gen_models.py:196`) — **not** from
`outputs`. The closure at `critical_point.rs:415-421` matches lowercase `"tc"`/`"pc"`/`"z_c"`,
so against a `valid_range` declaring `Tc`/`Pc`/`Z_c` every one resolved to `None`.

That did not remove the checks silently — `apply_checks` pushes a `RangeCheckSkipped` warning
for each (`range.rs:299-302`). So the bounds on the critical temperature, pressure and
compressibility were **skipped and reported as skipped, and no test reads that warning**.
Separately, the output keys are a contract with the contract test, the case `expected` keys
and the generated docs page, and they matched no field at all. Sixteen other models had all
of these names equal; this one had none of them equal.

**Fix and its test.** The output keys, the `valid_range` quantities and the case `expected`
keys are all the field names. `test_model_contract.py::test_every_result_field_is_declared_in_the_spec`
now holds every result field to its spec's `outputs`, in both directions — the check that
found this. **What it does not yet cover** is the `valid_range`-to-result-name agreement that
was the other half, which is why the skipped-check warning went unnoticed; that is an open
remediation item.

### An installed azoth cannot be imported

**Severity: blocking. Status: fixed at `f8b5a70`.**

**Symptom.** `pip install azoth` succeeded, and the first `import azoth` raised
`ModuleNotFoundError`. The package imports `keycard` at module scope, `keycard` imports
`yaml` at module scope, and `pyyaml` was declared only as a development extra.

**Evidence.** Reproduced on a wheel built at `65c2c67` and installed into a fresh venv:

```
$ maturin build --release -o /tmp/dist && python -m venv /tmp/wv
$ /tmp/wv/bin/pip install /tmp/dist/*.whl
$ cd /tmp && /tmp/wv/bin/python -c "import azoth"
  File ".../azoth/eos/components.py", line 49, in <module>
    from azoth import keycard
  File ".../azoth/keycard.py", line 50, in <module>
    import yaml
ModuleNotFoundError: No module named 'yaml'
```

`python/src/azoth/__init__.py:78` imports `keycard` in the module-scope list;
`python/src/azoth/keycard.py:50` is a top-level `import yaml`; `pyproject.toml:43` listed
`pyyaml>=6` under `[project.optional-dependencies]` `dev`, while `[project] dependencies`
held only `pint`.

**This made CI red at `ba98b99`.** The gate that catches it exists and is not weak:
`tools/check_wheel_data.py:43` does `import azoth` and `from azoth.hydraulics import …` from
outside the repository, and the `wheel-data` job builds a wheel, installs it into
`/tmp/wheel-venv` with **no dev extras**, and runs that script from `/tmp`. The run was red
there, so the baseline commit was not a baseline — see the note under "Reconciled at
`ba98b99`".

**Root cause.** The defect was fixed in `aa4a705` ("Declare pyyaml, because the wheel could
not be imported") and reverted by `37e6a34`, whose entire message is the stock two-line
revert text — it does not say why. The revert also carried `ruff format` line-wrapping into
`python/tests/test_registration_completeness.py`, verified as pure `re.compile` reflowing, so
a formatter was allowed to move a correctness fix. The cause of the *revert* is not recorded
anywhere and is not recoverable from the commit.

**Fix and its test.** `pyyaml` moved into `[project] dependencies`. The test is the gate that
was already red and is not new: a wheel built from the fix, installed into a clean venv with
no extras, with `check_wheel_data.py` run from `/tmp`, which reports OK and resolves both
CSVs out of `site-packages`. The pre-fix wheel was also rebuilt from `65c2c67`'s manifest and
the failure reproduced from that artifact, so the reproduction is a property of the wheel
rather than of a session's environment.

Two things travelled with it, both instances of the same defect — a list written down twice.
`uv.lock` records the resolved graph, so it was stale the moment `pyyaml` changed sections;
CI's `uv sync` would have rewritten it in place instead of failing, which is how a stale lock
survives. And `CONTRIBUTING.md` hand-enumerated the dependency list and disagreed with the
manifest in two places; it now runs `uv sync --extra dev --no-install-project`, the command
the CI jobs run, so the list has one home.

**A second route into the same import.** The traceback reaches `keycard` through
`azoth/eos/components.py:49`, not only through `__init__.py:78`. Both were checked; either
alone makes `import azoth` need yaml, and the first traceback taken only showed the first.

**`jsonschema` was checked and stays a dev extra.** Nothing under `python/src/azoth` imports
it; only `tools/check_user_data.py` and the tests do.

### Three `models` fields were accepted by the schema and stored by nothing

**Severity: major. Status: fixed at `fc37cce`.**

**Symptom.** A keycard's `models` section admitted `critical_rule`, `volume_translation` and
`root_selection`. A card carrying any of them passed the schema, passed the loader, and the
value was then dropped — `Model` has no field for it and `from_model` reads only
`components`. A key someone can write and never see again is a promise the format does not
keep, and it is the same class as the silently-ignored `source_ref` fields removed earlier.

**Evidence.** `specs/schema/keycard.schema.json:159-168` (before the fix) admitted all
three, each with a one-member enum. `python/src/azoth/keycard.py:448-462` validated
`kind`/`shape`/`alpha`/`mixing_rule` and read nothing else; `_models` had no unknown-key
check, so an unrecognised key passed through `body` untouched.

**Root cause.** The model vocabulary was written when the schema listed five cubics and six
alpha functions, and the extra fields described dimensions of that larger set. They survived
the narrowing that removed the other members, because a field with one legal value still
looks meaningful. It is not: it can only ever restate the single behaviour, and
`root_selection: from_phase` is what `mixture.rs` already does unconditionally.

**Fix and its test.** All three are deleted rather than implemented — the same reasoning the
schema already applies to the removed Mathias-Copeman branches. The loader now refuses a
model key it does not know, and
`test_keycard_loader.py::test_the_model_keys_are_the_schema_s_properties` binds the loader's
accepted set to the schema's `properties` plus `required`, both directions, so the next
field added to one and not the other fails the build.

### `provenance.json` did not hash the component databank

**Severity: major. Status: fixed at `9a29f54`.**

**Symptom.** The release provenance record hashed the fittings table and the two fluid
tables, and omitted `data/components/components.csv` and `data/components/kij.csv` — the
two largest shipped data files, and the two every `eos` calculation reads. A record that
does not cover them does not describe the release: change one critical constant and every
answer the library gives changes, with nothing in the record to compare against.

**Evidence.** `tools/provenance.py:150-154` listed three paths. `provenance.json`'s `data`
array carried three entries where five files ship.

**Root cause.** A hand-written list of files, in a project whose pages and registries exist
precisely because hand-written lists go stale. Nothing compared the list to what is
actually in `data/`.

**Fix and its test.** The tuple is walked from the data directory, so a shipped file is
hashed because it ships rather than because somebody remembered it. Five files where there
were three, and `--verify` round-trips against its own record.

### `gen_databank.py` decoded its input with `errors="replace"`

**Severity: minor. Status: fixed at `9a29f54`.**

**Symptom.** A byte in NeqSim's `COMP.csv` or `INTER.csv` that is not valid UTF-8 became
U+FFFD, the row parsed, and the corruption shipped inside a component's name or formula.

**Evidence.** `tools/gen_databank.py:137` opened with `errors="replace"`. The unit
conversions are guarded by `ROUND_TRIP` against published methane constants; the encoding
was guarded by nothing, so a mangled name would have been caught only if it happened to hit
one of four checked fields.

**Root cause.** `errors="replace"` is the forgiving default when a CSV is only being
inspected. Here the rows are transcription input to a shipped artefact, where the
forgiving behaviour is the harmful one.

**Fix and its test.** The file is decoded strictly and the failure names the byte offset.
The vendored copies under `databank/sources/neqsim/` are byte-identical to upstream, so the
check is reproducible.

---

## Source prose that is a report rather than a comment

**Actioned.** This section was a backlog: defect and limitation reports living in source and
specs, which had to move here *before* the prose was removed from the specs. Each is now an
entry in the sections above, at its own severity, and the prose it came from is gone.

The backlog existed because [P5](./process.md) said a comment either explains the code in
front of it, reports something for this page, or does not belong in the file — and the third
kind had no home. It has one now.

Two of the items were not findings but *pointers*: sentences inside specs that told a reader to
go and look at a test file to find one. `specs/calcs/eos/pr_z_factor.yaml` and
`specs/models/eos/stability_test.yaml` each carried one. **A finding recorded where nothing
reads it is not recorded**, and pointing at a test file is the same failure with an extra step.
Both sentences are gone; the findings they pointed at were already here.

One reading was contested and is resolved by deletion rather than by argument: whether
`separator.rs`'s prose had been moved to the *right spec field*. The question presupposed that
a spec is a home for prose. It is not — see [Spec files](./spec-files.md).

---

## Untested

Gaps that no test fills. Reconciled against the tree at `ba98b99`; three entries in the test
plan's table are now closed and are marked so rather than deleted, because a gap that was
covered by a test that did not pass is a different fact from a gap that was closed.

| Gap | Status |
|---|---|
| `crates/azoth-process` had no Rust tests at all | **closed** — 8 tests and a doctest in `crates/azoth-process/tests/` |
| `test_cross_impl.py` walked `CALCS` only, leaving the 17 models to per-model files | **closed** — `test_python_and_rust_agree_on_models` covers every model case |
| `test_mixture_layer.py`'s Helmholtz test claimed cross-language agreement it never made | **closed** — replaced by `test_the_helmholtz_layer_matches_its_recorded_values`, whose name says what it does |
| The PyO3 binding has no Rust tests | open |
| `azoth-cli`'s `main.rs` and `report.rs` | open |
| `gen_models.py`, `gen_databank.py` checked by no test | open |
| `InvalidInputError` absent from the identity test | open |
| `UnverifiedCalculationError` raised by nothing | open |
| `eos.expander` and `process.pump` have one spec case each | open |
| `check_wheel_data.py`, `check_links.py`, `provenance.py` have no test | open |
| `_dispatch.py`, `_rust_bridge.py`, `_data.py` have no dedicated test | open |

---

## Reconciled at `ba98b99`

**The baseline commit is not green in CI, and this page said it was.** Locally every gate
passed. One gate could not be run from inside the repository — `check_wheel_data.py` refuses
to report anything when `azoth` resolves to the source tree, which it did — and that refusal
was recorded and then treated as satisfied. It is not: the wheel job installs a wheel into a
venv with no dev extras and imports the package, and **that job is red**, for the `pyyaml`
defect now closed at `f8b5a70`. The test was run afterwards and reproduces `ModuleNotFoundError`
from a clean venv. A gate that declines to prove something has proved nothing.

**`main` was knowingly red between `65c2c67` and `f8b5a70`.** The fix was held while the
remediation plan was written, by the owner's decision, so that the plan existed before any of
it was executed. The interim is recorded here rather than left to be inferred from the
history: for three commits this page described a red baseline and the tree did not disagree
with it.

**Every root cause on this page was attacked by a second reader, and five were refuted.**

| Entry | First reading | After the adversarial pass |
|---|---|---|
| `eos.ph_flash` | the stopping criterion; a narrow band at a root-structure change | stopping criterion **and** an upstream `pt_flash` defect; **105 of 111 points**, one returned temperature. The first reading alone would have converted 105 silent wrong answers into 105 errors and fixed nothing |
| An installed azoth cannot be imported | the gate that would catch it does not exist | the gate exists and is red; the revert's message does **not** state the consequence, and its reason is unrecoverable |
| `prose_lint.py` is wired into nothing | wire it into `ci.yml` | it **fails on its own phrase list** (7 of its 9 hits are itself) and cannot pass on this tree. "Wire it in" was red on arrival |
| `eos.critical_point` output keys | the output key is the load-bearing name | `valid_range[].quantity` is what `apply_checks` resolves, and against the old keys every bound was **skipped with a `RangeCheckSkipped` warning that no test reads** |
| `pt_flash.yaml`'s residual claim | wrong by thirteen orders | wrong by **ten**; and the file named in the title was wrong. The block had already concluded residual is not pinnable — the harness pinned it anyway |
| `gen_models`/`gen_databank` untested | `gen_stub` is invoked by a test | it is not; the entry's own "correction" to the test plan was wrong and the test plan was right |

Two arithmetic errors of mine were also caught — "four and a half orders" for 7.66, and
"thirteen orders" for 10.2. Both are corrected above rather than quietly.

**What this changes about the process.** The root-cause pass produced five causes that were
wrong or half-wrong, on a page whose whole purpose is being trusted. The adversarial reader
found every one of them by going to the cited line and measuring. That is the mechanism P2
already names and that no gate performs — which is the same finding as the `prose_lint` entry,
arrived at independently.
