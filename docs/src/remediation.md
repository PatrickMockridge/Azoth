# Remediation plan

The next plan, written from [Required improvements](./required-improvements.md). That page is
the feedstock; this one orders it.

Every item names **what changes**, **the test that closes it**, and **the exit criterion**.
An item is done when its test passes at a named commit and the defect log entry moves to
fixed — not when the code changed.

Read [Process](./process.md) first. Its rules bind here: P3 (nothing red is committed), P4 (a
fix is verified by a test that would have caught it), P5 (prose in the right place).

---

## The order, and why

Four axes, applied in this order:

1. **Restore the baseline.** A red `main` makes every other item unverifiable, because a
   change set cannot be measured against a tree that does not build.
2. **Correctness before coverage.** A test written over a wrong answer records the wrong
   answer. The `eos.ph_flash` defect is the worked example: its two spec cases round-trip
   through the same wrong function, so adding a third case would have added a third
   agreement with a wrong number.
3. **Coverage before prose.** Deleting a comment is only safe once something would notice if
   the code it described changed meaning.
4. **Prose last**, because the de-spam is the largest change set and the one with no
   behaviour to verify. It needs the tests from (3) underneath it.

---

## R1. Restore the baseline

### R1.1 `pyyaml` is a runtime dependency

**Change.** Move `pyyaml` from `[project.optional-dependencies] dev` to `[project]
dependencies` in `pyproject.toml`. Check whether the same is true of `jsonschema` — the
keycard loader and the spec tools read it, but only the tools may, in which case it stays.

**Test that closes it.** The wheel job, which is already written and already red: build,
install into a venv with no extras, `import azoth`, run `check_wheel_data.py` from outside the
repository.

**Exit.** `wheel-data` green, and `ba98b99`'s claim to be a baseline becomes true.

**Note.** This is first because it is the only item that blocks a user from using the library
at all, and because the process that let it through — a gate that declined to prove something
being read as a pass — is `R4.1`'s subject.

---

## R2. Correctness

### R2.1 `eos.ph_flash` cannot invert most of its own range

**Change.** Two layers, and both are required.

*Upstream, in `eos.pt_flash`.* It returns `z_vapour` jumping 0.6636 → 0.9189 across the
temperature where Wilson's K for n-butane crosses 1, while the feed is single-phase on both
sides. It also reports `beta = 92.96` on an `all_vapour` phase. The flash has to hand back a
root that varies continuously across that threshold, and `beta` has to satisfy `0 ≤ beta ≤ 1`
or be `None` — never a value outside its own range.

*Downstream, in `solve_temperature`.* `crates/azoth-eos/src/flash_property.rs:256` declares
convergence on the width of the temperature bracket. The residual is computed at `:286` and
discarded. It must be tested against `algorithm.tolerance` before convergence is declared.

**Test that closes it.** The 111-point sweep becomes a property case: for every enthalpy in a
range spanning the defect, `ph_flash` either returns a temperature whose enthalpy is the one
asked for, within tolerance, or raises `SolverNotConverged`. Plus a case asserting every
`beta` a flash returns is within `[0, 1]` or absent.

**Exit.** No point in the sweep returns a wrong answer silently. The `beta` assertion passes
on every case in every flash spec.

**Decision required.** What `pt_flash` should return across the threshold is a modelling
choice — the two roots are a real degeneracy of the cubic at that composition. **This item
cannot be sized until that is answered**, and it is the one item here that may need its own
spec discussion rather than a patch.

### R2.2 A skipped range check is invisible

**Change.** `eos.critical_point`'s bounds were skipped for as long as its `valid_range`
declared `Tc`/`Pc`/`Z_c` while the resolver matched `tc`/`pc`/`z_c`. `apply_checks` reported
each as `RangeCheckSkipped` and **no test reads that warning**. Add the check that
`valid_range[].quantity` resolves to a name the implementation reports, alongside the
`outputs`-to-field contract test added in `ba98b99`.

**Test that closes it.** A contract test asserting that every `valid_range` quantity of every
model is resolvable — the same shape as `test_every_result_field_is_declared_in_the_spec`, for
the other name.

**Exit.** No spec has a range check that cannot run.

### R2.3 `pt_flash.yaml`'s arithmetic is wrong by ten orders

**Change.** `specs/models/eos/pt_flash.yaml:148-155` says two implementations "differ in its
last few *ulps*". The difference is 1.4e10 ulps of the residual. Replace with the measured
reason: the difference is ~0.4 ulp of the **order-1 quantity inside** the residual, which is
why no relative tolerance can be met on that field.

**Test that closes it.** None — prose. Verified by review against the measurement recorded in
the defect log.

**Exit.** The statement matches the numbers.

---

## R3. Coverage

Ordered by what a test in each place would catch. Each item is one change set.

| | Item | Test that closes it |
|---|---|---|
| R3.1 | `crates/azoth-python` has no Rust tests — the binding is exercised only from Python | Rust tests that construct results and assert the Python-visible shape without going through Python |
| R3.2 | `azoth-cli`'s `main.rs` and `report.rs` untested | Argument-parsing and report-formatting tests, plus one end-to-end run of the built binary |
| R3.3 | Five of eight registered error names are unchecked, and nothing asserts the *set* | Derive the expected names from `register()` and assert equality with `azoth.core.errors`' exports |
| R3.4 | `gen_models`, `gen_stub`, `gen_databank` exercised by no test | A synthetic spec per generator, exercising a shape the tree does not contain. `gen_databank` additionally needs a decision: a recorded fixture, or documented as uncheckable |
| R3.5 | `process.expander` and `process.pump` have one spec case each | A second case each, at a state exercising a different branch |
| R3.6 | `check_wheel_data`, `check_links`, `provenance`, `_dispatch`, `_rust_bridge`, `_data` have no dedicated test | One test each, or a recorded disposition of why not |

**Exit.** Every row of the defect log's "Untested" table is closed or dispositioned.

---

## R4. Mechanisms

### R4.1 The prose standard has no working mechanism

**Change**, in this order, because the first is a prerequisite:

1. **Exclude `tools/` from its own scan.** `prose_lint.py` currently reports itself: 7 of its
   9 hits are its own `HISTORY_PHRASES` list. It cannot pass on any tree and cannot be
   enabled as it stands. Its docstring's claim that each phrase "occurs zero times" was
   measured without counting the file doing the measuring.
2. **Wire it into `ci.yml`** beside `spec_lint.py`.
3. **Accept what it cannot do**, and say so in P5: it is a phrase list, and the phrases that
   matter most — `deliberately`, `the whole point`, `load-bearing` — are house style, so
   catching them means catching the corpus. **P2's second reader is the real mechanism for
   those**, and no gate performs it today.

**Test that closes it.** The tool exits 0 on a clean tree and 1 on a planted violation; CI
runs it.

**Exit.** P5 names a mechanism that runs, or P5 says the rule is enforced by review alone.
Shipping a tool that reads as enforcement while enforcing nothing is the defect being fixed.

### R4.2 Nothing performs P2's second reader

**Change.** P2 requires verification by something other than the author, at a named commit.
This cycle showed why: a second reader refuted five of six root causes, including one whose
remediation would have fixed nothing. Make the practice a named step rather than an intention
— a review artefact per change set, recording who read it and what they challenged.

**Exit.** A change set can be shown to have been read by someone other than its author.

---

## R5. De-spam

Scoped **whole tree, sequenced per crate**, and sized by audit rather than by comment density —
`pr_kappa.rs` is the densest file in the tree at 70% and has zero P5 violations, so density is
not the instrument.

**The standard is [P5](./process.md#p5-the-coding-standard-for-prose)**, applied per crate:
**MOVE** to the defect log or a spec's `source:`/`assumptions:` block, **DELETE**, or **KEEP**.

**Per crate, one change set:** audit → execute → verify no behaviour change (tests green,
generators `--check` clean) → commit. Order: `azoth-process` (the audit already exists) →
`azoth-eos` → `azoth-hydraulics` → `azoth-core` → `azoth-python` → `azoth-cli` →
`azoth-thermal` → specs → `python/src`.

**Two rules that are not negotiable.**

- **A deletion must not remove the only copy.** Prose deleted from a `.rs` file is added to a
  spec or the defect log **in the same commit** (P3 atomicity). A block of three lines or more
  deleted with no matching addition under `specs/` or `docs/` is a bug.
- **Do not over-delete.** KEEP is *why this expression and not another* —
  `mixer.rs:38-43`, `separator.rs:39-45`. The audit's eight named tics are the targets, not
  every comment.

**Entry.** The ~120 MOVE lines are already itemised in the defect log and must be filed before
the prose that carries them is removed.

---

## Definition of done

1. Every blocking and major entry in the defect log is `fixed`, with a named test.
2. Every minor entry is fixed or dispositioned with a reason.
3. `main` is green in CI including `wheel-data`.
4. The de-spam is complete per crate, each commit showing no behaviour change.
5. P5 names a mechanism that runs, or states honestly that it does not.

**Not in scope here:** new calculations and new unit operations. This plan returns the tree to
a state where adding them is safe, and nothing more.
