#!/usr/bin/env python3
"""Check `databank/manifest.toml` against the files it describes.

# Why this exists

The manifest is where a decision about vendored data is written down: this column
was carried across, that one was left because no model reads it, the other one
because a model does not exist yet. A record like that is worth having only if it
cannot rot, and rot here is silent in both directions - a column added to the
vendoring script and not the manifest, or removed from the slice and left in the
manifest. Neither breaks anything. Both make the manifest a description of a tree
that no longer exists.

# What it checks

The rules are named in `tools/manifest.py`, which does the work; this is the
command. They are, in one line each: every vendored file exists and is declared;
the manifest's `used` columns and the vendored header agree in both directions;
the row counts agree; a `not-yet` column naming a spec id names one that exists.

# What it cannot check

Whether the vendored slice is *current* against NeqSim. That needs the upstream
checkout, which CI does not have, so it runs only under
`python tools/gen_databank.py --check <checkout>`. The manifest carries the commit
it was last checked against, and nothing here can tell you a newer NeqSim exists.

Usage:
    python tools/check_manifest.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import manifest as manifest_module


def main() -> int:
    if not manifest_module.MANIFEST.is_file():
        sys.exit(f"check_manifest: {manifest_module.MANIFEST} does not exist")

    found, problems = manifest_module.read()
    problems.extend(manifest_module.validate(found))
    # **The two citation checks: the commit a provenance names, and the version it must not.**
    # A version identifies neither of NeqSim's trees, so a port citing one cites a revision
    # nobody can check it against - which is what this caught the first time it ran.
    problems.extend(manifest_module.citation_problems(found))
    problems.extend(manifest_module.version_problems(found))

    if problems:
        for problem in problems:
            print(f"  ERROR  {problem}", file=sys.stderr)
        print(f"\ncheck_manifest: FAILED with {len(problems)} problem(s)", file=sys.stderr)
        return 1

    columns = sum(len(f.columns) for f in found.files())
    files = found.files()
    carried = [c for f in files for c in f.columns if c.disposition in manifest_module.CARRIED]
    used = sum(1 for c in carried if c.disposition == "used")
    unstated = sum(1 for c in carried if c.unit == manifest_module.NO_STATED_UNIT)
    print(
        f"check_manifest: OK ({len(files)} vendored file(s), {columns} column(s), "
        f"{len(carried)} carried of which {used} read, "
        f"{len(found.not_vendored)} not-vendored entr(ies))"
    )
    # The two numbers that answer the question this file exists for: is the data here,
    # and does anything use it. A carried column with no stated unit is neither - it is
    # present and unmeaning, so it is counted rather than left to look finished.
    print(f"  {len(carried) - used:>3}  carried, nothing reads it yet")
    print(f"  {unstated:>3}  carried with no unit NeqSim states ({manifest_module.NO_STATED_UNIT})")

    # Printed rather than left to `grep`: a reason is a quoted flow mapping, so
    # `grep 'reason: not-yet'` matches nothing. See `manifest.reasons`.
    for prefix, count in manifest_module.reasons(found).items():
        print(f"  {count:>3}  {prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
