#!/usr/bin/env python3
"""Hold both implementations to the refusals their spec declares, and the spec to theirs.

# Why this exists

A spec says which values its inputs carry and which the port does not. That fact was
prose: an `assumptions` sentence, a `notes` paragraph, and an error message assembled at
each refusal site. Prose in two languages drifts, and it did - measured,
`billet_schultes_1999` was refused by the Python half, carried by the Rust half, and
declared ported by the spec, and no test in the tree could see any of it, because the
cross-implementation test compares *cases* and no case states that value.

Three spellings of one fact is the defect. The fix is one spelling: the spec's
`[[unported]]` array, read by `crates/azoth-process/src/unported.rs` and
`python/src/azoth/process/reference/_unported.py`, and named by one string literal at
each refusal site.

# What it checks

**Every key the declaration carries is refused by both implementations, and every key
either implementation refuses is declared.** The comparison is on the *key* - a string
literal at the call site - so nothing here reads a sentence, and a refusal that went
around the helper is invisible to this rule rather than silently agreeing with it.

**A row names something the spec has.** `parameter` is a declared input; a `value` is a
member of that input's `values` where the input enumerates them; a `capture` is a file
that exists. `spec_lint.py` holds the same three against the JSON Schema, and that is
deliberate duplication rather than a second source of truth: this tool has to run without
`jsonschema` so the test suite can point it at a synthetic tree.

# The comparison is per spec, by file stem

A refusal site is attributed to a spec by the stem of the file it sits in -
`crates/azoth-process/src/kernels/rate_based_packed_column.rs` and
`python/src/azoth/process/reference/rate_based_packed_column.py` are both
`rate_based_packed_column`, and `specs/models/**/rate_based_packed_column.toml` is the
spec they answer to. A stem that names more than one model spec is an error, because the
attribution would be a guess.

Usage:
    python tools/check_unported.py            # check; 0 ok, 1 a refusal and its spec disagree
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

MODEL_DIR = ROOT / "specs" / "models"

#: A call site, in either language: the helper's name and one string literal. The literal
#: is the row's key, so a site that names a row differently is a site this reads
#: differently - which is the point rather than a limitation. The comma is optional because
#: `rustfmt` breaks a long argument onto its own line and leaves a trailing one, which is
#: how two rows went unread on this checker's first run over the converted tree.
CALL = re.compile(r"\b(?:unported::refuse|_unported\.refuse)\(\s*\"([^\"]+)\"\s*,?\s*\)")

#: Where a call site is written, in each language.
RUST_SITES = ROOT / "crates"
PYTHON_SITES = ROOT / "python" / "src"


def unported_key(row: dict[str, object]) -> str:
    """The canonical key a row is declared under, spelled the way both helpers spell it."""
    parameter = str(row["parameter"])
    if "value" in row:
        return f"{parameter}={row['value']}"
    conditions = sorted(
        f"{condition['parameter']}={condition['value']}"  # type: ignore[index]
        for condition in row["when"]  # type: ignore[union-attr]
    )
    return f"{parameter}@{'&'.join(conditions)}"


def model_specs(model_dir: Path = MODEL_DIR) -> list[tuple[Path, dict[str, object]]]:
    """Every model spec under a tree, parsed, with its path."""
    specs: list[tuple[Path, dict[str, object]]] = []
    for path in sorted(model_dir.rglob("*.toml")):
        specs.append((path, tomllib.loads(path.read_text(encoding="utf-8"))))
    return specs


def declared(model_dir: Path = MODEL_DIR) -> tuple[dict[str, set[str]], list[str]]:
    """Every declared row as `model id -> its keys`, and every well-formedness failure.

    **Keyed by model rather than by key**, because two ids legitimately refuse the same
    thing: `process.absorption_column` and `process.distillation_column` both refuse the
    eight `ColumnSolverFactory` strategies neither carries, and a key-only map called that
    a duplicate on this checker's first run over the converted tree.

    `model_dir` is a parameter so the test suite can point this at a synthetic tree - the
    reason `spec_lint.py` takes `--spec-dir`, and the only way to prove the check fails on
    a disagreement rather than on the tree it happens to be reading.
    """
    failures: list[str] = []
    rows: dict[str, set[str]] = {}
    for path, spec in model_specs(model_dir):
        where = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        inputs = spec.get("inputs", {})
        assert isinstance(inputs, dict)
        mine = rows.setdefault(str(spec["id"]), set())
        for row in spec.get("unported", []) or []:
            assert isinstance(row, dict)
            key = unported_key(row)
            if key in mine:
                failures.append(f"{where}: declares {key!r} twice")
            mine.add(key)

            parameter = str(row["parameter"])
            declared_input = inputs.get(parameter)
            if not isinstance(declared_input, dict):
                failures.append(f"{where}: {key!r} names {parameter!r}, which is not an input")
            else:
                values = declared_input.get("values")
                if "value" in row and isinstance(values, list) and row["value"] not in values:
                    failures.append(
                        f"{where}: {key!r} refuses a value {parameter!r} does not carry - its "
                        f"`values` are {values}"
                    )
            capture = row.get("capture")
            if capture is not None and not (ROOT / str(capture)).exists():
                failures.append(f"{where}: {key!r} names the capture {capture!r}, which is absent")
    return rows, failures


def sites(root: Path, suffix: str) -> dict[str, set[str]]:
    """Every refusal key written under a tree, by the file stem it is written in."""
    found: dict[str, set[str]] = {}
    for path in sorted(root.rglob(f"*{suffix}")):
        if "/target/" in str(path) or "/vendor/" in str(path):
            continue
        keys = set(CALL.findall(path.read_text(encoding="utf-8")))
        if keys:
            found.setdefault(path.stem, set()).update(keys)
    return found


def check(
    model_dir: Path = MODEL_DIR,
    rust_root: Path = RUST_SITES,
    python_root: Path = PYTHON_SITES,
) -> list[str]:
    """Every disagreement between the declaration and the two implementations."""
    declared_rows, failures = declared(model_dir)

    stems: dict[str, list[str]] = {}
    for path, spec in model_specs(model_dir):
        stems.setdefault(path.stem, []).append(str(spec["id"]))
    ambiguous = {stem: ids for stem, ids in stems.items() if len(set(ids)) > 1}
    for stem, ids in sorted(ambiguous.items()):
        failures.append(
            f"specs/models: the stem {stem!r} names {sorted(set(ids))}, so a refusal site "
            f"in a file of that name cannot be attributed to one of them"
        )

    for language, root, suffix in (
        ("rust", rust_root, ".rs"),
        ("python", python_root, ".py"),
    ):
        for stem, keys in sorted(sites(root, suffix).items()):
            model_ids = stems.get(stem)
            if model_ids is None:
                failures.append(
                    f"{language}: {stem} refuses {sorted(keys)}, and no model spec of that "
                    f"name declares anything - a refusal with nothing to hold it"
                )
                continue
            model = model_ids[0]
            undeclared = sorted(keys - declared_rows.get(model, set()))
            if undeclared:
                failures.append(
                    f"{language}: {stem} refuses {undeclared}, which {model} does not declare"
                )

    # **And the other direction**, which is what catches a declaration the code has not
    # kept: every row of a model both implementations must refuse, and one of them does
    # not. Compared per model, because a key belongs to the spec that declares it and two
    # specs can declare the same one.
    for language, root, suffix in (
        ("rust", rust_root, ".rs"),
        ("python", python_root, ".py"),
    ):
        written = {key for keys in sites(root, suffix).values() for key in keys}
        for model, keys in sorted(declared_rows.items()):
            missing = sorted(keys - written)
            if missing:
                failures.append(
                    f"{language}: {model} declares {missing} and no {language} site refuses "
                    f"them, so the row is a claim this implementation has not kept"
                )
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)

    failures = check()
    if failures:
        for failure in failures:
            print(f"  ERROR  {failure}", file=sys.stderr)
        print(f"\ncheck_unported: FAILED with {len(failures)} problem(s)", file=sys.stderr)
        return 1

    rows, _ = declared()
    total = sum(len(keys) for keys in rows.values())
    print(
        f"check_unported: OK ({total} declared row(s) over "
        f"{len([m for m, k in rows.items() if k])} model(s), refused in both languages)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
