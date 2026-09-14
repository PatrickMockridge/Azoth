# Required improvements

What is wrong with this library, with the evidence and the cause. Feedstock for the next
remediation plan.

Every entry has the same shape, so the page can be read as a set rather than as prose:

```
### <title>            severity: blocking | major | minor
Filed by:   the test run | the inventory | the expert panel
Symptom     what is wrong, in one or two sentences
Evidence    the exact reproduction and the measured numbers
Root cause  as determined by the panel; contested readings recorded as contested
Remediation one sentence, enough to size it
```

Root causes are determined by a panel of five specialists — numerics, the cross-language
contract, the specification system, the public API, and test quality — each reading the
same failure set. Where the panel agrees, that is the cause. Where it splits, **both
readings are recorded and the entry says so**; a manufactured consensus is worth less than
a recorded disagreement. Where it does not know, the entry says that and names the
experiment that would settle it.

An entry leaves this page when it is fixed, and the fix carries a test.

See [Test plan](./test-plan.md) for what is tested and how a run is performed.

---

## Blocking

### `eos.ph_flash` returns a wrong answer across a root-structure change

**Severity: blocking.** Filed by: the test run.

**Symptom.** `eos.ph_flash` inverts an enthalpy by bisecting `H(T)`, which is well posed
only while `H` is monotonic. At a change in the cubic's root structure it is not, and the
bisection converges to the discontinuity and returns a temperature whose enthalpy is not
the one asked for — **reported as converged**.

**Evidence.** Methane/n-butane, `z = [0.6, 0.4]`, 20 bar:

| T (K) | 386.0 | 386.5 | 386.7 | 387.0 |
|---|---|---|---|---|
| `z_liquid` | 0.084495 | 0.918637 | 0.918637 | 0.918855 |
| `z_vapour` | 0.666606 | 0.918637 | 0.918637 | 0.918855 |
| `H` (J/mol) | 811.2 | 1813.6 | 1866.0 | 1876.7 |

The largest admissible root moves by 0.25 between 386.0 K and 386.5 K and the enthalpy
built from it jumps about 1000 J/mol. Asking for `H = 1277.0 J/mol` returns **386.694 K
with `H(T) = 1865.8`** — 589 J/mol out.

Found through the unit operation that hits it: `process.heater` at +80 kW on 10 mol/s
closes the energy balance to 7 per cent instead of to the flash's own 1e-8.

```python
r = azoth.eos.ph_flash(fluid, ideal_gas, P=Q(20.0e5, "Pa"), H=Q(1276.996, "J/mol"),
                       z=[0.6, 0.4])
r.T          # 386.694 K
r.residual   # large: the enthalpy there is 1865.8 J/mol
```

**Root cause.** _Awaiting the panel._

**Remediation.** _Awaiting the panel._

---

## Major

### `crates/azoth-process` has no Rust tests

**Severity: major.** Filed by: the inventory.

**Symptom.** The eight unit operations' Rust implementations are exercised only
indirectly, through Python calling the extension. No `crates/azoth-process/tests/`
directory exists and there is no `#[cfg(test)]` anywhere in the crate.

**Evidence.** `cargo test -p azoth-process` reports zero tests. Every other crate with
domain code has integration tests; this is the newest code in the tree and the only crate
in this position.

**Root cause.** _Awaiting the panel._

**Remediation.** _Awaiting the panel._

### The 17 models are outside the central cross-language check

**Severity: major.** Filed by: the inventory.

**Symptom.** `python/tests/test_cross_impl.py` parametrises over `CALCS` only. Model
agreement rests entirely on per-model test files, so a model added without one is covered
by nothing — and nothing reports that.

**Evidence.** `test_python_and_rust_agree` iterates `CALCS` (21 ids). The 17 ids in
`_models_gen.MODELS` appear in no central agreement test.

**Root cause.** _Awaiting the panel._

**Remediation.** _Awaiting the panel._

### A test's name claims a cross-language guarantee its body does not make

**Severity: major.** Filed by: the inventory.

**Symptom.** `python/tests/eos/test_mixture_layer.py::test_both_implementations_agree_on_the_helmholtz_layer`
never calls the Rust backend. It compares Python results against hard-coded literals
recorded from the Rust tests.

**Evidence.** The test body pins neither backend — no `use_backend`, no comparison to a
Rust result.

**Why major rather than minor.** An absent test is a known gap. This one reads as covered,
so it suppresses the question. It is the failure mode `test_registry_contract.py` was
written about: a list that was complete when it was written, and a second thing that grew
beside it.

**Root cause.** _Awaiting the panel._

**Remediation.** _Awaiting the panel._

---

## Minor

### `UnverifiedCalculationError` is raised by nothing

**Severity: minor.** Filed by: the inventory.

**Symptom.** The error exists in the Rust `AzothError` enum, is mapped across the binding
and re-exported from Python, and no code path raises it in either language.

**Evidence.** Grep finds only its definition and the `to_pyerr` mapping.

**Root cause.** _Awaiting the panel._

**Remediation.** Either a code path raises it or it is deleted. Dead surface a caller can
catch and never will is worse than no surface.

### Two unit operations have one spec case each

**Severity: minor.** Filed by: the inventory.

**Symptom.** `eos.expander` and `process.pump` declare a single worked case. Every other
id has two or more.

**Evidence.** Spec case counts by id; both are 1.

**Root cause.** _Awaiting the panel._

**Remediation.** _Awaiting the panel._

---

## Untested

Gaps that no test fills. Each is a place a defect can sit unnoticed, listed with what a
test there *would* catch. Whether any is worth filling is the remediation plan's decision.

| Gap | What a test there would catch |
|---|---|
| `azoth-cli`'s `main.rs` and `report.rs` | An argument-parsing or report-formatting fault — what a user actually types |
| The PyO3 binding layer has no Rust tests | A binding that agrees with Python only because Python introspects it |
| `gen_models.py`, `gen_stub.py`, `gen_databank.py` have no test | A generator that fails for a spec shape not currently in the tree |
| `InvalidInputError` is absent from the same-class-object identity test | A second error class that is not the same object across the boundary |
| No end-to-end test of the compiled CLI binary | The binary failing to run where the library succeeds |
| `check_wheel_data.py`, `check_links.py`, `provenance.py` have no test | A CI-only script that breaks and is noticed only in CI |
| `_dispatch.py`, `_rust_bridge.py`, `_data.py` have no dedicated test | A dispatch fault covered only incidentally by the contract tests |

---

## Contested

Root causes on which the panel split. Recorded rather than resolved.

_Empty until the panel reports._
