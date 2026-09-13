## What this changes

<!-- One paragraph. If it adds or changes a calculation, say which. -->

## If this adds a calculation

Every calculation must ship with all of these, and CI enforces most of them:

- [ ] A spec under `specs/calcs/` with an equation, a source, a valid range, assumptions, a worked example, and tests
- [ ] The LaTeX form of the equation (`latex`), not just the evaluable form
- [ ] One Python function and one Rust function
- [ ] The worked example's expected values derived and shown, not copied from a source's answers
- [ ] `python tools/gen_registry.py && python tools/gen_docs.py` run and committed

## If a source is uncertain

Say so rather than guessing. Set `verification.status` in the spec and either
find a source or skip the test with a recorded reason:

- [ ] `verification.status` reflects what is actually known
- [ ] Any test that cannot run is `skipped` with a `skip_reason`

Do not reproduce tables or text from a standard, textbook, or paper. Cite the
equation; write your own worked example.

## Dependencies

<!-- The project keeps its dependency tree small. Any new dependency needs a
     sentence here saying why it earns its place, and why something already
     present would not do. -->

None, or:

- `crate-or-package` — why

## Checklist

- [ ] `cargo test && cargo clippy --all-targets -- -D warnings && cargo fmt --all --check`
- [ ] `pytest` (with the extension built: `maturin develop && AZOTH_REQUIRE_RUST=1 pytest`)
- [ ] `ruff check && ruff format --check && mypy`
- [ ] `python tools/spec_lint.py`
- [ ] `mdbook build docs` and `python tools/check_links.py`
- [ ] Commits are signed

## Anything a reviewer should look at first

<!-- The part you are least sure about. -->
