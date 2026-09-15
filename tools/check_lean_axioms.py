#!/usr/bin/env python3
"""Run Lean's `#print axioms` over every theorem this project claims, and refuse a gap.

Reads the gate files under `lean/Azoth/` - hand-written prose and a generated list,
each carrying one `#print axioms <theorem>` per claim - runs each through
`lake env lean`, and fails unless every reported axiom set is a subset of the three
Lean permits:

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
with no gaps and still not the claim a reader expects, and no check closes that -
it is what review is for. `python/tests/test_lean_claims.py` covers the two things
that *are* mechanical: that the gate is not empty, and that every name in it is a
declaration some file actually makes.

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

#: The files carrying `#print axioms` lines, and the reason there are two.
#:
#: `Axioms.lean` is hand-written and holds the general theorems, beside the prose
#: explaining what the gate is for. `Gate.lean` is generated from the vocabulary
#: table and holds one line per canonical unit - generated because a
#: hand-maintained list of twenty-four names goes stale the first time a unit is
#: added, and it goes stale *silently*: the theorem is proved and nothing gates it.
GATES = (
    LEAN_DIR / "Azoth" / "Axioms.lean",
    LEAN_DIR / "Azoth" / "Gate.lean",
)

#: The axioms a proof may rest on. `propext` and `Quot.sound` are the two Lean
#: includes that are not definitional, and `Classical.choice` is what makes
#: classical reasoning available; all three are axioms of the standard library
#: rather than of this development. `sorryAx` and `Lean.ofReduceBool` are the two
#: that mean something is unfinished or unsound, and are named here to be refused
#: with a message that says which rather than merely "not in the list".
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}

#: `'Azoth.Dim.foo' depends on axioms: [propext, Quot.sound]`
_DEPENDS = re.compile(r"^'?(?P<name>[A-Za-z0-9_.«»]+)'? depends on axioms: \[(?P<axioms>.*)\]$")

#: `'Azoth.Rho.quote_injective' does not depend on any axioms` - the empty case,
#: which `#print axioms` spells differently from a zero-entry list.
_EMPTY = re.compile(r"^'?(?P<name>[A-Za-z0-9_.«»]+)'? does not depend on any axioms$")

#: Every `#print axioms` line in the gate file, which is what it exists for.
_PRINTED = re.compile(r"^\s*#print axioms\s+(?P<name>[A-Za-z0-9_.]+)\s*$", re.MULTILINE)


def fail(message: str) -> None:
    print(f"check_lean_axioms: {message}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="only report problems")
    args = parser.parse_args()

    claimed: list[str] = []
    reported: dict[str, list[str]] = {}
    for gate in GATES:
        if not gate.exists():
            fail(f"{gate.relative_to(ROOT)} does not exist, so what it gates is not gated")
            return 1

        claimed.extend(_PRINTED.findall(gate.read_text(encoding="utf-8")))

        proc = subprocess.run(
            ["lake", "env", "lean", str(gate.relative_to(LEAN_DIR))],
            cwd=LEAN_DIR,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            fail(
                f"`lake env lean {gate.relative_to(LEAN_DIR)}` failed:\n{proc.stdout}{proc.stderr}"
            )
            return 1
        for line in proc.stdout.splitlines():
            if match := _DEPENDS.match(line.strip()):
                reported[match.group("name")] = [
                    a.strip() for a in match.group("axioms").split(",") if a.strip()
                ]
            elif match := _EMPTY.match(line.strip()):
                reported[match.group("name")] = []

    if not claimed:
        # Files that print nothing would otherwise pass every check below by
        # having nothing to check - the vacuous-gate failure this whole tool is
        # written against.
        fail("no gate file names a theorem, so the gate is empty")
        return 1

    problems: list[str] = []

    # Both directions. A theorem in the file that Lean did not report means the
    # `#print` line did not take effect - and a gate that silently stopped gating
    # one theorem is worse than no gate, because the file still looks complete.
    for name in claimed:
        if name not in reported:
            problems.append(
                f"{name} is named in a gate file but Lean reported nothing for it, so "
                f"it is not actually gated"
            )
    for name in reported:
        if name not in claimed:
            problems.append(
                f"Lean reported {name}, which no gate file names - the parse and the files disagree"
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
