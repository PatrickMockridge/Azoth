#!/usr/bin/env python3
"""Fail on a JSON document that carries the same key twice.

# Why this exists

`json.loads` keeps the last of two identical keys and discards the first without
saying so, and `jsonschema`'s `check_schema` does not look for them either. So a
schema can carry a duplicated keyword and pass every gate in this repository.

It did. `specs/schema/keycard.schema.json` had `"additionalProperties": false`
twice at the document root and twice inside `models`. Both pairs happened to
agree, which is the only reason nothing broke - the value that survived was the
intended one, by luck. A pair that disagreed would have resolved silently to the
second, and a keyword whose meaning depends on which of two identical lines a
maintainer reads last is a rule nobody can rely on.

# What it checks, and what it does not

Every `.json` file under `specs/schema/` and `validation/`. Both are documents a
person writes and edits in place, which is the shape that produces a duplicate.
Generated JSON is not scanned, because no person is editing it.

This is a backstop, not a linter. It has one rule and reports one thing. It does
not validate the document - `spec_lint.py` and `check_user_data.py` do that
against the schemas themselves.

Usage:
    python tools/check_json_keys.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: The directories whose JSON a person writes. See the docstring for why
#: generated documents are excluded.
SCANNED = ("specs/schema", "validation")


def duplicate_keys(path: Path) -> list[str]:
    """Every key `path` defines more than once, in the order they were met.

    All of them rather than the first, because the file this exists for carried
    two - a reporter that stopped at one would have to be run twice to say so.
    `json` calls the hook once per object after reading it whole, so a key
    repeated in two *different* objects is one entry here; that is rare enough
    that the line list below is still the useful part of the message.
    """
    found: list[str] = []

    def hook(items: list[tuple[str, object]]) -> dict[str, object]:
        seen: dict[str, object] = {}
        for key, value in items:
            if key in seen and key not in found:
                found.append(key)
            seen[key] = value
        return seen

    json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=hook)
    return found


def lines_of(path: Path, key: str) -> list[int]:
    """Where a key appears, as `"key":` at the start of a line.

    Anchored to the line rather than searched anywhere in it, because a
    description mentioning a keyword in a sentence is not a definition of it.
    Every document scanned here is pretty-printed one key per line, so this is
    exact for them; a document that is not would report fewer lines than it has,
    which makes the message less useful and not wrong.
    """
    pattern = re.compile(rf'^\s*"{re.escape(key)}"\s*:')
    return [
        number
        for number, text in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if pattern.match(text)
    ]


def describe(path: Path) -> list[str]:
    """One message per duplicated key in `path`; empty if it has none."""
    where = path.relative_to(ROOT)
    messages = []
    for key in duplicate_keys(path):
        found = lines_of(path, key)
        at = f"lines {found[0]} and {found[1]}" if len(found) >= 2 else "twice"
        messages.append(
            f"{where}: {at} both define {key!r}. A parser keeps the last and discards "
            f"the first without saying so, so which value a reader sees depends on "
            f"which line they read last."
        )
    return messages


def main() -> int:
    roots = [ROOT / name for name in SCANNED]
    missing = [name for name in SCANNED if not (ROOT / name).is_dir()]
    if missing:
        sys.exit(f"check_json_keys: no such directory: {', '.join(missing)}")

    paths = sorted(path for root in roots for path in root.rglob("*.json"))
    if not paths:
        # An empty scan passes for the same reason an empty glob does: the check
        # did not run, and a check that did not run has not passed.
        sys.exit(f"check_json_keys: no JSON found under {', '.join(SCANNED)}")

    errors = [message for path in paths for message in describe(path)]

    if errors:
        for error in errors:
            print(f"  ERROR  {error}", file=sys.stderr)
        print(
            f"\ncheck_json_keys: FAILED with {len(errors)} duplicated key(s)",
            file=sys.stderr,
        )
        return 1

    print(f"check_json_keys: OK ({len(paths)} document(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
