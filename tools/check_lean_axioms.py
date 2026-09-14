#!/usr/bin/env python3
"""Run Lean's `#print axioms` over every theorem this project claims, and refuse a gap.

Reads `lean/Azoth/Axioms.lean`, which is a hand-written file whose entire content
is one `#print axioms <theorem>` per claim, runs it through `lake env lean`, and
fails unless every reported axiom set is a subset of the three Lean permits:

    propext, Classical.choice, Quot.sound

**Why this rather than a search for `sorry`.** A `sorry` elaborates to `sorryAx`,
and `#print axioms` reports the *transitive* axiom set of a proof term - so a
theorem that depends on a `sorry` anywhere in its dependency chain reports
`sorryAx`, including one inside `lean/vendor/lean-units/`, which a scan over
`lean/Azoth/` cannot see. It also catches what a text search cannot: `admit`,
`sorryAx` applied by hand, and an `axiom` declaration, which appears by name. A
search over the sources would miss all four.

**And what it does not do.** The gate proves a proof is *complete*; it does not
prove the theorem is the one wanted. A lemma with a weakened hypothesis is a proof
with no gaps and still not the claim a reader expects. The guard for that is not
here: each claim in `docs/src/calculus/` names the theorem it states, and
`python/tests/test_lean_claims.py` asserts that every theorem the book names is
named in `Axioms.lean` - so deleting one fails a check rather than passing quietly.

Usage:
    python tools/check_lean_axioms.py            # check
    python tools/check_lean_axioms.py --quiet    # only report problems

Exit status is non-zero if a theorem is missing from the gate or rests on
something outside the three.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEAN_DIR = ROOT / "lean"
AXIOMS = LEAN_DIR / "Azoth" / "Axioms.lean"

#: The axioms a proof may rest on. `propext` and `Quot.sound` are the two Lean
#: includes that are not definitional, and `Classical.choice` is what makes
#: classical reasoning available; all three are axioms of the standard library
#: rather than of this development. `sorryAx` and `Lean.ofReduceBool` are the two
#: that mean something is unfinished or unsound, and are named here to be refused
#: with a message that says which rather than merely "not in the list".
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}

#: `'Azoth.Dim.foo' depends on axioms: [propext, Quot.sound]`
_DEPENDS = re.compile(r"^'?(?P<name>[A-Za-z0-9_.«»]+)'? depends on axioms: \[(?P<axioms>.*)\]$")

#: Every `#print axioms` line in the gate file, which is what it exists for.
_PRINTED = re.compile(r"^\s*#print axioms\s+(?P<name>[A-Za-z0-9_.]+)\s*$", re.MULTILINE)


def fail(message: str) -> None:
    print(f"check_lean_axioms: {message}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="only report problems")
    args = parser.parse_args()

    if not AXIOMS.exists():
        fail(f"{AXIOMS.relative_to(ROOT)} does not exist, so nothing is gated")
        return 1

    source = AXIOMS.read_text(encoding="utf-8")
    claimed = _PRINTED.findall(source)
    if not claimed:
        # A file that prints nothing would otherwise pass every check below by
        # having nothing to check - the vacuous-gate failure this whole tool is
        # written against.
        fail(f"{AXIOMS.relative_to(ROOT)} contains no `#print axioms` line, so the gate is empty")
        return 1

    proc = subprocess.run(
        ["lake", "env", "lean", "Azoth/Axioms.lean"],
        cwd=LEAN_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        fail(f"`lake env lean Azoth/Axioms.lean` failed:\n{proc.stdout}{proc.stderr}")
        return 1

    reported: dict[str, list[str]] = {}
    for line in proc.stdout.splitlines():
        match = _DEPENDS.match(line.strip())
        if match:
            reported[match.group("name")] = [
                a.strip() for a in match.group("axioms").split(",") if a.strip()
            ]

    problems: list[str] = []

    # Both directions. A theorem in the file that Lean did not report means the
    # `#print` line did not take effect - and a gate that silently stopped gating
    # one theorem is worse than no gate, because the file still looks complete.
    for name in claimed:
        if name not in reported:
            problems.append(
                f"{name} is named in {AXIOMS.name} but Lean reported nothing for it, so "
                f"it is not actually gated"
            )
    for name in reported:
        if name not in claimed:
            problems.append(
                f"Lean reported {name}, which {AXIOMS.name} does not name - the parse "
                f"and the file disagree"
            )

    for name in sorted(set(claimed) & set(reported)):
        extra = sorted(set(reported[name]) - ALLOWED)
        if extra:
            problems.append(
                f"{name} depends on {extra}, which is outside {sorted(ALLOWED)}. A "
                f"sorryAx here means the proof has a gap; any other name means it "
                f"rests on something this project did not agree to."
            )

    if problems:
        for problem in problems:
            fail(problem)
        return 1

    if not args.quiet:
        print(
            f"check_lean_axioms: {len(claimed)} theorem(s) gated, every one resting only "
            f"on {sorted(ALLOWED)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
