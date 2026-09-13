# Contributing

Thank you for considering it. This is a calculation library, so the bar is
deliberately higher than "the tests pass": a wrong number that looks reasonable
is the failure this project is organised against, and the process is shaped
around making that hard to ship.

## Getting set up

You need Rust (the version in `rust-toolchain.toml`), Python 3.12, and for the
docs, mdBook with two preprocessors.

```bash
git clone https://github.com/PatrickMockridge/Azoth
cd Azoth

uv venv --python 3.12
uv pip install maturin pytest ruff mypy pyyaml jsonschema "pint>=0.24"
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

## The checks

All of these run in CI. Run them before opening a pull request.

```bash
cargo test && cargo clippy --all-targets -- -D warnings && cargo fmt --all --check

env -u PYTHONPATH AZOTH_REQUIRE_RUST=1 .venv/bin/python -m pytest
.venv/bin/ruff check && .venv/bin/ruff format --check && .venv/bin/mypy

.venv/bin/python tools/spec_lint.py
.venv/bin/python tools/gen_registry.py --check
.venv/bin/python tools/gen_docs.py --check
.venv/bin/python tools/check_links.py
mdbook build docs
```

`AZOTH_REQUIRE_RUST=1` makes a missing extension a failure rather than a skip.
Without it a broken build makes every test pass on the Python path while the
cross-language agreement is verified by nothing.

## Adding a calculation

Four files, and the rest follows.

**1. The spec** — `specs/calcs/<namespace>/<name>.yaml`. This is the real work.
It must carry an equation, a `latex` form for the docs, a source, inputs and
outputs with units, a valid range with a reason for every bound, the assumptions
that are *not* checked, a worked example, and the tests.

**2. Python** — `python/src/azoth/hydraulics/reference/<name>.py`.

**3. Rust** — `crates/azoth-hydraulics/src/<name>.rs`, plus a result struct in
`results.rs`.

**4. Tests** — `python/tests/hydraulics/test_<name>.py` and
`crates/azoth-hydraulics/tests/<name>.rs`, both driven by the spec's `tests`
list.

Then:

```bash
.venv/bin/python tools/gen_registry.py && .venv/bin/python tools/gen_docs.py
```

The docs, the range checks both implementations enforce, and the test cases both
implementations run are all generated from that one YAML file.

Read an existing calc end to end first — `reynolds_number` is the simplest.

## The rules that are not negotiable

**No calculation without a source and a worked example.** If you cannot find a
source, set `verification.status: source_needed` in the spec and skip the worked
example test with a recorded reason. A gap that is visible is fine; a citation
nobody checked is not.

**Never reproduce copyrighted tables or text.** This is why the fitting
coefficients are labelled placeholder rather than sourced: Crane TP-410 is
copyrighted, and transcribing its tables is not something this project does. Cite
the equation, work your own arithmetic. *Facts* are not copyrightable, so a
single coefficient recorded with a citation is fine; a table is not.

**Every number in a worked example must be derivable and shown.** Write the
substitution out. A worked example nobody can retrace is a number somebody typed.

**Units-safe signatures, no bare floats for physical quantities.** Genuinely
dimensionless quantities (`f`, `Re`, `epsilon/D`) are plain floats; everything
else is a `pint` quantity in Python and a `uom` quantity in Rust.

**Errors are typed; range violations are warnings.** They are different things
and the library keeps them different. A value outside its validated range is
still a value - say so, do not refuse it.

**Justify new dependencies.** The tree is deliberately small. Say in the pull
request what the dependency buys and why something already present will not do.

**Do not hand-edit generated files.** `spec_gen.rs`, `_registry_gen.py` and
everything under `docs/src/hydraulics/` are output. Each carries a banner saying
so. Change the spec and regenerate.

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
