#!/usr/bin/env python3
"""Did the tests that ran execute the lines that changed?

    python tools/check_coverage.py --check --base origin/main --json coverage.json

**Why this exists.** `--verify-selection` proves every emitted pytest target collects a test.
That is not the same claim as the tests *running the changed line*: a target can collect forty
tests and touch none of the code the diff changed. Scoping a push to what it reaches is only
sound if "reaches" means executed, and the only way to know is to measure it.

**A threshold is deliberately not offered.** A percentage is a number to game - raise it by
testing the easy modules, and the file nobody touched keeps its zero. What is asserted here is
one fact per changed file: **at least one line of it ran.** A changed file with no executed
line is a file the gate never looked at, whatever the total says.

Exit codes follow the tree's convention: 1 is a finding, 2 is could-not-decide.

The decision is a pure function of the coverage document and the changed paths, so both are
tested from fixtures; only `main` reads git and the filesystem.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: What the instrumented run covers. The Python half only: `cargo-llvm-cov` builds the whole
#: workspace with instrumentation, which is minutes per push for a signal the Rust side already
#: has from its own per-crate tests - every crate the diff reaches is emitted by name, and
#: `--verify-selection` proves those targets collect.
PREFIX = "python/src/azoth/"


def executed_files(document: dict) -> set[str]:
    """Every file the coverage run reports as having run at least one line."""
    files = document.get("files")
    if not isinstance(files, dict):
        raise ValueError(f"no `files` mapping in the coverage document: {sorted(document)}")
    return {
        # Coverage records paths as the interpreter resolved them, which is absolute.
        os.path.relpath(str(name), ROOT)
        for name, data in files.items()
        if isinstance(data, dict) and data.get("executed_lines")
    }


def files_with_no_line_run(changed: list[str], executed: set[str]) -> list[str]:
    """The changed files no line of which ran. Empty is the only good answer."""
    return sorted(path for path in changed if path not in executed)


def changed_python(base: str) -> list[str]:
    """The `python/src/azoth` files this diff changed and that still exist."""
    result = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=d", f"{base}...HEAD", "--", PREFIX],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return sorted(line for line in result.stdout.split("\n") if line.strip())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="the only mode there is")
    parser.add_argument("--base", default="origin/main", help="the revision the diff is against")
    parser.add_argument("--json", default="coverage.json", help="the coverage document")
    parser.add_argument(
        "--stdin", action="store_true", help="read the changed list from stdin, for a test"
    )
    args = parser.parse_args()

    path = ROOT / args.json
    if not path.exists():
        print(
            f"check_coverage: could not decide - no {args.json}. The instrumented run writes it; "
            "if pytest-cov did not run, this check has measured nothing and says so rather than "
            "passing."
        )
        return 2
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        executed = executed_files(document)
        changed = (
            [line.strip() for line in sys.stdin if line.strip()]
            if args.stdin
            else changed_python(args.base)
        )
    except (ValueError, subprocess.CalledProcessError, OSError) as exc:
        print(f"check_coverage: could not decide - {exc}")
        return 2

    missing = files_with_no_line_run(changed, executed)
    if missing:
        print(f"check_coverage: {len(missing)} changed file(s) no line of which ran:")
        for path_ in missing:
            print(f"    {path_}")
        print(
            "A file the diff changed and the gate never executed is a file nothing checked. "
            "Either the scoped set does not reach it - which is a map to fix in "
            "tools/test_scope.py - or it is generated, and a generated file's own generator is "
            "what covers it."
        )
        return 1
    print(f"check_coverage: every one of {len(changed)} changed file(s) ran at least one line")
    return 0


if __name__ == "__main__":
    sys.exit(main())
