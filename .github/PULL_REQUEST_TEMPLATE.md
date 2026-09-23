## What this changes

<!-- One paragraph. If it adds or changes a calculation, say which. -->

## If this adds a calculation

The contract in full is [the specification](docs/src/architecture/specification.md#how-a-calculation-is-ported).
The short version:

- [ ] A spec under `specs/calcs/` — equation, LaTeX form, source, inputs and outputs
      with units, valid range, assumptions, a worked example, and tests
- [ ] One Python file and one Rust file
- [ ] The worked example's expected values derived and shown, not copied from a
      source's answers
- [ ] The generators run and committed: `python tools/gen_registry.py &&
      python tools/gen_docs.py`

**Nothing is registered.** There is no dispatch table, no `__all__` and no
registration call to update, because a calculation's id *is* its address: the module
path and the function name follow from it by convention. If you found yourself
editing a list to say "this calculation exists", something is wrong — say so in the
PR rather than doing it.

## If a source is uncertain

Say so in words, in the spec's `notes`, rather than guessing.

There is deliberately no status field to set. This library does not record whether a
person checked an attribution, because no tool can verify that a person did — so the
field would be a form to fill in rather than a fact, and a form teaches people to
fill it in. Write down what is confirmed and what is not, and a reader can judge it.

- [ ] The spec's `notes` say what has and has not been confirmed
- [ ] Any test that cannot run is `skipped`, with a `skip_reason`

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
