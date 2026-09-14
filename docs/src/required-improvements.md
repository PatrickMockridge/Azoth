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

### `tools/prose_lint.py` is wired into nothing

**Severity: major. Status: open.**

**Symptom.** [The process standard](./process.md) P5 states prose rules and names
`tools/prose_lint.py` as their mechanism. The tool appears in no workflow and no hook, so it
has never run in CI — and **it cannot pass on this tree anyway**, so wiring it in is not the
one-line fix it looks like.

**Evidence.**

```
$ python tools/prose_lint.py ; echo $?
tools/prose_lint.py:48: history phrase 'formerly'   ...   (7 such lines)
crates/azoth-process/src/mixer.rs:24: history phrase 'had to change'
specs/models/process/mixer.yaml:207: history phrase 'had to change'
prose_lint: 9 line(s) narrate the code's history.
1
```

**Seven of the nine hits are the tool matching its own `HISTORY_PHRASES` list** at
`tools/prose_lint.py:47-53`. `sources()` includes `("tools", "*.py")` and `is_generated()`
does not exempt it, so the linter reports itself, exits non-zero, and would fail any build it
were added to. Its docstring's claim that each phrase "occurs zero times in this tree outside
the change set" is false — the measurement it describes did not count the file doing the
measuring.

`grep -rn prose_lint .github/` returns nothing; CI runs `spec_lint.py`, `cargo fmt` and
`ruff check`.

The docstring is self-defeating in a second way: it exempts `deliberately`, `the whole point`
and `load-bearing` because those occur elsewhere in the tree (`docs/src/spec.md:5,19`). Using
the corpus as a licence is the scope-shrinking failure in tool form.

**Root cause.** Two, and the second is why the entry is not a formality. The tool was never
wired in — that is the symptom's cause. But it was also never *run*, which is why nobody
noticed it fails on itself: a check that is never executed is a check that cannot be
debugged. The deeper cause is that a phrase list cannot judge prose, which the tool's own
docstring concedes.

**Remediation.** Three parts, and the first is not optional. Exclude `tools/` from its own
scan (or the phrase list from the scan), because until that is done the tool reports itself
and cannot be enabled. Then wire it into `ci.yml` beside `spec_lint.py` and accept it as a
narrow backstop. Then accept that the phrases that matter most — `deliberately`, `the whole
point`, `load-bearing` — are not catchable by a list, and that P2's second reader is the real
mechanism. Deleting the tool and saying so in P5 is the honest alternative to all three.

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

---

## Fixed

Entries that have left the sections above. Each names the commit that closed it and the test
that would have caught it. The first four were closed at `ba98b99`, in the change set that
first reached a green baseline; the fifth at `f8b5a70`, and it is what makes that baseline
true rather than merely claimed.

**A note on evidence, for the first four.** They were uncommitted work when they were fixed,
so `ba98b99^` is not their pre-fix state and the failures cannot be re-derived from git. The
evidence for each is the captured run: the exact assertion messages and exit codes from the
session that made the change, quoted below. Where a fix is checkable from the tree — the
`critical_point` names, the mutation test — that is stated and is reproducible. The fifth is
not in that position: it is reproducible from wheel artifacts at any time, and the
reproduction is in its own entry below.

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

---

## Source prose that is a report rather than a comment

[P5](./process.md) says a comment either explains the code in front of it, reports something
for this page, or does not belong in the file. The entries below are the third kind found in
the change set: defect and limitation reports currently living in source and specs, which
must move here before the prose is removed. The de-spam is scoped in the remediation plan.

Each is **open** and each is a MOVE from a named line range.

| Finding | Where it lives now |
|---|---|
| Separator: nothing checks that the feed is stable before splitting it | `specs/models/process/separator.yaml` assumptions |
| Separator: the `trivial` phase is routed to the gas outlet — a convention, not a checked result | `specs/models/process/separator.yaml` |
| Separator: no entrainment carry-over is modelled | `specs/models/process/separator.yaml` |
| Separator: the split is the flash's split; the flash's correctness is not independently verified | `specs/models/process/separator.yaml` |
| Separator: the duty branch measures enthalpy from the feed's own state, where NeqSim reads a history-dependent `getEnthalpy()` | `specs/models/process/separator.yaml` |
| Pump: nothing checks that the inlet is a liquid — a pump run on a vapour returns a plausible pressure rise | `specs/models/process/pump.yaml`, `python/src/azoth/process/reference/pump.py` |
| Mixer: the molar-versus-total enthalpy translation from NeqSim's joule sum | `specs/models/process/mixer.yaml` |
| Heater: the same molar-versus-total enthalpy translation | `specs/models/process/heater.yaml` |
| Heater: the spec's citations are unconfirmed — no equation number is claimed | `specs/models/process/heater.yaml`, `specs/models/eos/critical_point.yaml` references |
| Compressor: the largest assumption, that a cubic EOS is trustworthy at a high-ratio outlet, is unchecked | `specs/models/process/compressor.yaml` |
| `eos.ph_flash`'s bracket scan aborted the whole search on one inadmissible scan temperature | `specs/models/process/throttling_valve.yaml` — **fixed**, needs its entry and its test |
| `eos.ph_flash`'s `ph_flash` root-structure discontinuity, found through `process.heater` | `python/tests/models/test_heater.py` — **the blocking entry above**, filed where it cannot be found |

A pointer into a test file that carries no finding is the same class, and two exist:
`specs/calcs/eos/pr_z_factor.yaml:12-16` ends a sentence with "See
`python/tests/models/test_mixer.py` carries the state that exposed it", and
`specs/models/eos/stability_test.yaml:299` does the like. Both are open.

One reading is **contested** and recorded rather than resolved:
`specs/models/process/separator.yaml:14-20` ("This is the first unit operation ported from
NeqSim…") is arguably itself a MOVE under P5's second specific. `separator.rs` was cleaned by
moving prose *into* the spec, and whether it moved to the right spec field — `source:`,
`assumptions:` or `notes:` — is a question for the review.

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
