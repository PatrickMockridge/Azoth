# Contributing

Thank you for considering it. This is a calculation library, so the bar is
deliberately higher than "the tests pass": a wrong number that looks reasonable
is the failure this project is organised against, and the process is shaped
around making that hard to ship.

**It is becoming a calculation library *with a process layer*** — unit operations,
flowsheets that compose them, reports and an agent surface. That is a deliberate
crossing of a boundary the project used to hold, and it does not relax anything
below. Two things follow from it, and both are requirements rather than notes:

- **The two-implementation rule covers the new layer too.** A unit operation is a
  registered model with a spec, two implementations and a case both run; a flowsheet
  is a document both implementations execute and compare. There is no third lane
  where code is allowed to have one implementation.
- **A flowsheet has no single retraceable worked example**, and pretending otherwise
  would be the exact failure this file is organised against. What replaces it is
  stated in `docs/src/roadmap.md`: a hand-computable case for the composition, and
  conservation checks that hold at every answer.

## Getting set up

You need Rust (the version in `rust-toolchain.toml`), Python 3.12, and for the
docs, mdBook with two preprocessors.

```bash
git clone https://github.com/PatrickMockridge/Azoth
cd Azoth

uv venv --python 3.12
uv sync --extra dev --no-install-project   # the list is `pyproject.toml`'s, not this file's
maturin develop                 # builds the Rust extension into the venv
```

`maturin develop` needs `VIRTUAL_ENV` set and `CONDA_PREFIX` unset. If you have
conda active it will refuse with "Both VIRTUAL_ENV and CONDA_PREFIX are set" -
that is maturin being careful, not a broken environment.

For the docs:

```bash
cargo install mdbook mdbook-katex mdbook-mermaid
mdbook-mermaid install docs     # vendors the JS and registers additional-js
```

Versions matter here. `mdbook-katex` 0.9 does not work with mdBook 0.5; the
preprocessor protocol changed, and the failure reads as "Unable to parse the
input". `mdbook-linkcheck` 0.7 requires mdBook 0.4 and is unmaintained, which is
why link checking is `tools/check_links.py` instead.

**If your shell sets `PYTHONPATH`,** clear it for the test run:
`env -u PYTHONPATH .venv/bin/python -m pytest`. A `PYTHONPATH` pointing at other
projects makes pytest autoload plugins that fail to import, and the error looks
nothing like its cause.

## Where work lands

**Directly on `main`.** Commit, push, and let CI tell you.

There was a rule that every change went on its own branch behind its own pull
request, and it is **suspended, not repealed**. It is worth writing down why,
because the reasoning has not changed even though the practice has: a branch per
milestone is what makes `git bisect` able to separate "this equation is wrong"
from "this tooling is wrong", and it is what keeps a reviewer from having to read
a diff of a whole subsystem at once.

What changed is the baseline. This library is part-way through porting NeqSim, and
while that is true a milestone is not a unit anybody can review on its own — the
pieces only make sense together, and the branch machinery was costing more than it
bought. When the port is done and the surface stops moving, the rule comes back.

Two things survive regardless, because they were never about branches:

- **A commit is still the unit of reasoning.** One concern per commit, and a
  message that says *why*. `git bisect` and `git log` read the commit sequence, not
  the pull request.
- **The generated documentation is the review artefact.** Reading `docs/src/eos/`
  *is* reviewing the physics, because the docs are rendered from the specs the
  implementations read. That was true when a PR carried it and it is true now.

## The checks

All of these run in CI. Run them before pushing.

```bash
cargo test && cargo clippy --all-targets -- -D warnings && cargo fmt --all --check

env -u PYTHONPATH AZOTH_REQUIRE_RUST=1 .venv/bin/python -m pytest
.venv/bin/ruff check && .venv/bin/ruff format --check && .venv/bin/mypy

.venv/bin/python tools/spec_lint.py
.venv/bin/python tools/gen_registry.py --check
.venv/bin/python tools/gen_models.py --check
.venv/bin/python tools/gen_docs.py --check
.venv/bin/python tools/gen_stub.py --check
.venv/bin/python tools/check_links.py
mdbook build docs
```

`AZOTH_REQUIRE_RUST=1` makes a missing extension a failure rather than a skip.
Without it a broken build makes every test pass on the Python path while the
cross-language agreement is verified by nothing.

If you have supplied your own licensed data, one more:

```bash
.venv/bin/python tools/check_user_data.py keycard.yaml
.venv/bin/python tools/gen_user_data.py   keycard.yaml --check
```

Both need the file to exist, so neither is in CI — `keycard.yaml` is
gitignored and never committed. The tooling itself is covered by
`python/tests/test_gen_user_data.py`, which runs against the template.

## Adding data, rather than code

Most of what a user wants to change is not a calculation. **A keycard is a YAML file
that overrides or extends what the library ships** — components, interaction
parameters, coefficients, fittings, fluids and named model variants — with no code at
all, in either language. If what you are adding is a value rather than an equation, it
belongs there, and it costs a file rather than a contribution. See
[The keycard](docs/src/keycard.md).

## Adding a calculation

**Five new files, and nine edits to existing ones** — plus eight more if the calc is
the first in a new namespace.

That number is smaller than it was, and it is still falling. It was fourteen. Five of
those are gone, and the shape of how they went is worth naming:

- the bridge's `id → function` table and the `RESULT_TYPES` table are both **derived**
  — from the id, and from each implementation's return annotation — so there is no list
  to keep in step with another list;
- the README's table and this book's index are **generated into marked blocks**, so the
  two places a *reader* meets a new calculation no longer need editing;
- `_core.pyi` is **generated** by `tools/gen_stub.py` from the result dataclasses and
  the specs, so the type stub is no longer a second place to declare a function by hand.

What remains is mostly the Rust boilerplate: a name has to be attached to a function
before a Python caller can reach it, and that attachment is still typed by hand.

**This number is measured, not guessed**, and the last measurement was by adding a
calculation and counting. It has been wrong before — the section used to claim "the
rest follows", which was never true, and later gave a count that predated the batch
API. Don't trust a smaller number without re-measuring.

The paths below use `hydraulics` as the worked example. Substitute your own
namespace everywhere, and see the next section if it does not exist yet.

### The five new files

**1. The spec** — `specs/calcs/<namespace>/<name>.yaml`. This is the real work.
It must carry an equation, a `latex` form for the docs, a source, inputs and
outputs with units, a valid range with a reason for every bound, the assumptions
that are *not* checked, a worked example, and the tests.

**2. Python** — `python/src/azoth/hydraulics/reference/<name>.py`.

**3. Rust** — `crates/azoth-hydraulics/src/<name>.rs`.

**4. Tests, Python** — `python/tests/hydraulics/test_<name>.py`.

**5. Tests, Rust** — `crates/azoth-hydraulics/tests/<name>.rs`.

Both test files are driven by the spec's `tests` list, so they follow from the
spec's contents rather than being written against the implementation.

### The nine edits

Listed because "the rest follows" was a claim nobody had checked, and because a
forgotten one fails in a different way in each case:

| File | What to add |
|---|---|
| `crates/azoth-<ns>/src/lib.rs` | `pub mod`, the re-export, the crate docstring's list |
| `crates/azoth-<ns>/src/results.rs` | the result struct, its `CalcResult` impl, the two test tables |
| `crates/azoth-python/src/<ns>.rs` | the `#[pyfunction]` wrapper |
| `crates/azoth-python/src/results.rs` | the `Py*Result` transport class, `result_fields`, `calc_ids` |
| `crates/azoth-python/src/lib.rs` | `add_class`, `add_function` |
| `crates/azoth-python/src/batch.rs` | the batch arm — the Rust half of the batch API |
| `python/src/azoth/<ns>/__init__.py` | the dispatch wrapper, `__all__`, the id constant, the docstring list |
| `python/src/azoth/<ns>/reference/__init__.py` | the import and `__all__` |
| `python/src/azoth/batch/<ns>.py` | the batch wrapper — the Python half of the batch API |

Five of the files that used to appear here no longer do. `_rust_bridge.py` and
`core/result.py` are unchanged by a new calculation, because the tables that lived in
them are derived from the id and from the return annotation. `_core.pyi` is generated
by `tools/gen_stub.py`. `README.md` and `docs/src/index.md` *are* changed by one — but
by `tools/gen_docs.py`, into a marked block, and a hand-edit there fails the drift
check.

**The batch API is two edits, in two languages, and it is not optional.** There is
a Rust arm and a Python module because the batch path loops over the *scalar*
implementation on each side rather than introducing a third one, so both halves
have to learn the new calc. `python/tests/test_batch.py` enforces this: its
`test_the_excluded_set_is_exactly_crane_k_factors` asserts that the only calc
without a batch form is `crane_k_factors`, so a new calc that skips either edit
fails there rather than raising `NotImplementedError` at call time.

Then:

```bash
.venv/bin/python tools/gen_registry.py && .venv/bin/python tools/gen_docs.py
```

which produces the registries, the calc's doc page and the book's contents.

### If it is the first calc in a new namespace

Eight more, all of which a namespace needs once rather than per calc:

| File | What to add |
|---|---|
| `crates/azoth-<ns>/Cargo.toml` | the crate manifest (new file) |
| `Cargo.toml` | the crate in `members` and in `[workspace.dependencies]` |
| `crates/azoth-python/Cargo.toml` | the dependency |
| `crates/azoth-python/src/lib.rs` | `mod <ns>;`, alongside the per-calc `add_class`/`add_function` |
| `python/src/azoth/__init__.py` | `from azoth import ... <ns> ...`, and `__all__` |
| `python/src/azoth/batch/__init__.py` | the import and `__all__` |
| `tools/gen_docs.py` | the namespace's display title in `NAMESPACES` |
| `tools/provenance.py` | the namespace's shared files in `NAMESPACE_SUPPORT` |

The two generator dicts are the ones that fail loudly: `gen_docs` exits naming the
namespace it has no title for, and `provenance` would silently hash nothing for it.

**`azoth/__init__.py` is the one that fails quietly**, and it has now been missed
once — `thermal` was unreachable as `azoth.thermal` from a bare `import azoth` for
as long as that namespace existed, because the submodule was never imported. Add
the import, not just the `__all__` entry.

### How you find out you forgot one

Run the suite. Each omission has its own failure:

| Omission | What fails |
|---|---|
| the result dataclass, the stub, `calc_ids()`, the bridge function | `test_registration_completeness.py`, which names what is missing |
| the function signatures, the declared outputs, the result fields | `test_registry_contract.py` |
| either half of the batch API | `test_batch.py::test_the_excluded_set_is_exactly_crane_k_factors` |
| the generated block in `README.md` or `docs/src/index.md` | `test_the_announced_list_is_generated_rather_than_hand_edited`, by running the generator |
| the `NAMESPACES` entry | `gen_docs.py` exits, naming the namespace |

The two front-page files are also where the prose about what is *not* implemented
lives, and no test can read that for you: if you add the first relief-valve
calculation, go and fix the paragraphs around the generated block by hand.

Read an existing calc end to end first — `reynolds_number` is the simplest, and
`eos.pr_kappa` is the most recent and has the fewest moving parts.

## The rules that are not negotiable

Seven, and each is either enforced by a test or is a constraint from outside the
project. Anything that is merely good advice is under
[What reviewers will look at](#what-reviewers-will-look-at) instead, so that this
list stays short enough to actually read.

**Never reproduce copyrighted tables or text.** This is why the fitting coefficients
are placeholders rather than sourced: Crane TP-410 is copyrighted, and transcribing
its tables is not something this project does. Cite the equation; work your own
arithmetic. A *fact* is not copyrightable, so a single coefficient recorded with a
citation is fine; a *table* is not.

[Copyright and licensed data](docs/src/copyright.md) lists every place this rule has
changed what the library does, and the three patterns that resolve all of them: a
single coefficient becomes a function argument (`Cd`, `eta`, `f_t` — no file, no
licensing question), a *set* of coefficients becomes keycard data, and *code under a
permissive licence* may be ported with attribution. Read it before designing a
calculation around a value you cannot ship.

**A port is accepted on this library's tests, never on its provenance.** That a
well-known library implements something is evidence that it *can* be implemented. It
is not evidence that it is right. The concrete case: NeqSim's `CriticalPointFlash`
implements Heidemann & Khalil correctly by inspection and validates the result
**nowhere** — no pure-component check, no mixture-locus check, a silent `break` on
`NaN`. Porting it carries the algorithm and none of the confidence, so every port
needs a test that would fail if the port were wrong, and that test has to be ours.

Attribution for anything ported lives in [`NOTICE`](NOTICE), once, rather than in a
per-calc field — a block repeated in twenty specs is a block nobody reads. The spec
still cites two things: its `references` name the paper the method comes from, and its
`notes` say what was **changed** and why, because a port is never a transcription and
the differences are the part a reader cannot recover from either source.

**A worked example must be retraceable by hand.** Write the substitution out. An
example nobody can follow is a number somebody typed.

Be aware of what is **not** checked. Nothing verifies the arithmetic in a
`derivation`, because deciding which `=` in a paragraph is a claim turns out to need
more discrimination than prose carries: `1 US gallon = 3.785411784e-3 m**3` is a unit
conversion, `rho * g = 998 * 9.80665 = 9787.036699999999 W/(m**3/s)` is a claim with a
unit on the end, and a checker that cannot tell those apart flagged thirty-six of them
on this registry alone. It was withdrawn for being too noisy, because a check that
noisy teaches people to skim its output. So the arithmetic is a reviewer's job —
**recompute your derivation rather than restating what you expect it to be.**

**Units-safe signatures, no bare floats for physical quantities.** Genuinely
dimensionless quantities (`f`, `Re`, `epsilon/D`) are plain floats; everything
else is a `pint` quantity in Python and a `uom` quantity in Rust.

**A temperature difference is not an absolute temperature, and the spec has to
say which.** They share a dimension, so nothing in the units can separate them:
`pint` converts an absolute `Q(30, "degC")` to 303.15 K, and a calc that meant "a
30 kelvin difference" returns a plausible answer ten times too large. Mark such an
input `interval: true` in the spec, and convert it with `input_to_si(spec, name,
value)` rather than `to_si(...)` so the flag is read rather than restated - a flag
the spec holds and the call site ignores is a check that exists on paper only.
Rust needs no equivalent: `TemperatureInterval` and `ThermodynamicTemperature` are
different types and the absolute one does not fit. `spec_lint` errors on the flag
used away from `K`, and warns when a kelvin-dimensioned input named like a
difference lacks it; `test_every_interval_input_in_the_registry_is_honoured_on_both_backends`
walks the registry and checks that each marked input actually refuses an offset
unit on both backends, which is what catches a spec and an implementation that
disagree.

**Errors are typed; range violations are warnings.** They are different things
and the library keeps them different. A value outside its validated range is
still a value - say so, do not refuse it.

**Justify new dependencies.** The tree is deliberately small. Say in the pull request
what the dependency buys and why something already present will not do.

**Do not hand-edit generated files.** `spec_gen.rs`, `model_gen.rs`,
`_registry_gen.py`, `_models_gen.py` and every page under `docs/src/eos/`,
`docs/src/hydraulics/` and `docs/src/thermal/` are output. Each carries a banner
saying so. Change the spec and regenerate.

## What reviewers will look at

- **Is the source real?** Not "does it look plausible" - could you point someone
  at it. An equation number that is wrong is worse than one marked TODO.
- **Is the valid range honest?** Bounds that cannot be evaluated are worse than
  no bounds, because they read as validation. `spec_lint` checks that every
  range check names a quantity the calc can actually compute.
- **Are the assumptions complete?** This is where the real gaps live. The Darcy-
  Weisbach spec records that it assumes Mach < 0.3 and explicitly says it cannot
  check it.
- **Would the tests fail if the calc were wrong?** A test that compares an
  implementation against itself proves nothing. Watch for tolerance choices that
  would pass for an answer that is off by more than the source claims.

## Commits

Signed. The repository is configured for SSH signing; `git config commit.gpgsign`
should be true and `git verify-commit HEAD` should report a good signature.

Small commits with a message that explains *why* rather than restating the diff.
If you made a judgement call, put the reasoning in the message - the next person
needs to know whether a bound was chosen deliberately or copied.

## Licence

Contributions are under AGPL-3.0-or-later for code and CC-BY-4.0 for
documentation and data. See [LICENSE](LICENSE) and
[LICENSE-CC-BY-4.0](LICENSE-CC-BY-4.0).
