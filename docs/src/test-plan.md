# Test plan

What is tested, how, and — the part worth writing down — **what is not**.

This page is the entry point to a cycle: test everything, record what is found and why on
[Required improvements](./required-improvements.md), and remediate from that page. It is
written down rather than held in someone's head because a coverage gap that nobody has
written down is a gap that gets rediscovered as a bug.

## What is under test

| | |
|---|---|
| Registered calculations and models | **38** — 21 calcs, 17 models |
| Crates | 11 — `azoth-core`, four domains, the process layer, the binding, the CLI, test support |
| The Python package | `azoth` — the public API, the dispatch layer, the keycard, the batch API |
| The generators | `tools/*.py` — six that emit code or docs from the specs |
| The CLI | `azoth pipe` |

Each registered id has **two implementations that must agree** — Rust and pure Python —
plus a spec that both are held to. That is the central claim of this library and it is
what most of the levels below exist to check.

## The levels

Ten, each catching something the others cannot. They are listed in the order a failure
travels outward.

**1. Spec cases.** Every id declares worked cases with expected values. These pin what the
model *does*, and they are the weakest check here — a recorded number only says the model
still does what it did when the number was written. They catch a change, not a mistake.

**2. Cross-language agreement.** The Rust and Python implementations are run on the same
input in the same process and compared field by field. This is the only thing that makes
"two implementations, always" true rather than aspirational. It catches a divergence that
every single-language test would pass.

**3. Contract and shape.** The spec, the generated registry and the code must agree about
names: a declared output that is not a result field, an input that is not a parameter, an
id announced in one place and missing from another. Cheap, and they catch the drift that
every numerical test is blind to — the numbers agree perfectly and only the attribute
name differs.

**4. Property and conservation.** Identities that hold at *every* answer rather than at
one recorded one: moles conserved across a separator, energy conserved across a mixer, an
isentropic state genuinely isentropic. These are what make a model trustworthy away from
its worked examples.

**5. Error paths.** Each typed error is raised, in both languages, and is the *same class
object* on both sides so `except OutOfRangeError` works whichever backend answered.

**6. Range and boundary.** A value outside the range a correlation was validated in is
returned with a warning rather than refused — so the warning must actually be emitted, and
a check that could not run must report `RANGE_CHECK_SKIPPED` rather than pass silently.

**7. Data.** The component databank, the fitting tables, the shipped CSVs — parsed
identically by both languages, field for field.

**8. Tooling.** The generators are held to their own output: regenerate and diff. A spec
edited without regenerating is caught here rather than in a user's broken import.

**9. Integration.** The CLI, the batch API, the keycard loader — the pieces a user meets
before the thermodynamics.

**10. End-to-end.** A worked case through the whole stack.

## What is covered today

| | |
|---|---|
| Python | **558** test functions across 56 files |
| Rust | **302** — 240 in 30 integration files, 62 in 12 `#[cfg(test)]` modules |
| Spec cases | Every one of the 38 ids has them; 3 to 22 per id |
| CI | 10 jobs, all gates, none with `continue-on-error` |

Levels 1, 3, 4, 5, 6 and 8 are covered well. Level 2 is covered well **for calculations
and unevenly for models** — see the gaps. Levels 7, 9 and 10 are partial.

## The gaps

Ordered by what a test in each place would catch.

| Gap | What a test there would catch |
|---|---|
| **`crates/azoth-process` has no Rust tests at all** — no `tests/` directory, no `#[cfg(test)]`. Eight unit operations whose Rust implementations are exercised only indirectly, through Python | A Rust-side defect in the whole process layer. It is the only crate in the tree in this position, and it is the newest code |
| **`test_cross_impl.py` walks `CALCS` only** | A model whose two implementations diverge. The 17 models rest entirely on per-model files, so a model added without one is covered by nothing |
| **`test_mixture_layer.py::test_both_implementations_agree_on_the_helmholtz_layer` never calls the Rust backend** — it compares Python against hard-coded literals | A Helmholtz-layer divergence. The test's *name* claims a guarantee its body does not make, which is worse than the gap alone: it reads as covered |
| **`azoth-cli`'s `main.rs` and `report.rs` untested** — only `pipe::compute` is | An argument-parsing or report-formatting fault. Neither is thermodynamics, and both are what a user actually types |
| **The PyO3 binding has no Rust tests** | A binding that agrees with Python only because Python is what introspects it |
| **`UnverifiedCalculationError` is raised by nothing**, in either language | Dead surface: an error a caller can catch and never will |
| **`InvalidInputError` is missing from the same-class-object identity test** | A second error class that is not the same object across the boundary |
| **`gen_models.py`, `gen_stub.py`, `gen_databank.py` are checked by no test** — only by CI regenerating and diffing | A generator that fails only for a spec shape not currently in the tree |
| **`eos.expander` and `process.pump` have one spec case each** | A wrong answer at any state other than the one recorded |

## Entry and exit criteria

**Entry.** The tree builds, the extension is built, and the generators are current.

**Exit for one run.** Every command below has been run once, with its complete output
captured, and **nothing changed while it ran**.

```bash
cargo test --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings && cargo fmt --all --check
env -u PYTHONPATH AZOTH_REQUIRE_RUST=1 .venv/bin/python -m pytest -q
.venv/bin/ruff check && .venv/bin/ruff format --check && .venv/bin/mypy
.venv/bin/python tools/spec_lint.py
for g in gen_registry gen_models gen_docs gen_stub; do .venv/bin/python tools/$g.py --check; done
.venv/bin/python tools/check_links.py
```

Plus the two CI jobs pytest does not cover: `check_wheel_data.py` and `provenance.py`.

**Exit for the page.** [Required improvements](./required-improvements.md) is complete when
**every failure in the run's output and every gap in the table above has an entry**. That
is the criterion — not "the important ones".
