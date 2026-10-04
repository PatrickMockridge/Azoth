#!/usr/bin/env python3
"""Does every edge-case class have an owner, and is the owner there?

`docs/src/calculus/numerics.md` requires every partial function to be total, clamped or
deliberately `NaN`; `lean/guards.toml` declares which of those the code states, who holds the
statement, and what each class of `valid_range` row is held by. This is what makes the declaration
true rather than remembered:

* every `guarded:` marker in the kernels has an entry, and every entry has a marker - two
  directions, because a one-way check passes when a site is deleted;
* every shape a `valid_range` row in `specs/` takes has a class entry, so a new shape is a
  decision rather than an omission;
* the generated `Azoth/Guards.lean` carries a theorem per `owner = "lean"` entry and the generated
  `Azoth/GuardGate.lean` a `#print axioms` line per theorem, so a theorem the gate does not name
  is caught here rather than being a theorem nothing gates;
* every owner is one of the five the vocabulary allows, and every path an entry names resolves.

Stdlib only and no Lean: it answers "is the declaration complete", which is what the build needs
to know before Lean is asked whether the claims are true. That is `tools/check_lean_axioms.py`.
"""

from __future__ import annotations

import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import guards

ROOT = pathlib.Path(__file__).resolve().parent.parent
THEOREMS = ROOT / "lean" / "Azoth" / "Guards.lean"
GATE = ROOT / "lean" / "Azoth" / "GuardGate.lean"


def main() -> int:
    try:
        declared, classes = guards.manifest()
        marked = guards.markers()
        shapes = guards.class_keys()
    except guards.Refusal as refusal:
        print(f"check_guards: {refusal}", file=sys.stderr)
        return 1

    problems: list[str] = []

    # --- the markers and the manifest ------------------------------------------------------
    declared_keys = {guard.key for guard in declared}
    for key in sorted(set(marked) - declared_keys):
        where, line = marked[key]
        problems.append(f"{where}:{line}: `guarded: {key}` has no entry in {guards.MANIFEST.name}")
    for key in sorted(declared_keys - set(marked)):
        problems.append(f"{key}: declared in {guards.MANIFEST.name} and no call site carries it")

    # --- the classes ------------------------------------------------------------------------
    owned = {entry["key"] for entry in classes}
    for shape in sorted(shapes - owned):
        problems.append(
            f"the specs declare a `{shape}` bound, which no class in {guards.MANIFEST.name} owns"
        )

    # --- the generated files ----------------------------------------------------------------
    lean = [guard for guard in declared if guard.owner == "lean"]
    if not THEOREMS.exists() or not GATE.exists():
        problems.append(
            "lean/Azoth/Guards.lean or GuardGate.lean does not exist - run "
            "tools/gen_guard_theorems.py"
        )
    else:
        body = THEOREMS.read_text(encoding="utf-8")
        gate = GATE.read_text(encoding="utf-8")
        named = set(re.findall(r"^#print axioms Azoth\.Guards\.(\w+)", gate, re.M))
        for guard in lean:
            if f"theorem {guard.name} " not in body:
                problems.append(f"{guard.key}: no theorem `{guard.name}` in Guards.lean")
            if guard.name not in named:
                problems.append(f"{guard.key}: `{guard.name}` is not gated in GuardGate.lean")
        for extra in sorted(named - {guard.name for guard in lean}):
            problems.append(
                f'GuardGate.lean gates `{extra}`, which no `owner = "lean"` entry declares'
            )

    if problems:
        for problem in problems:
            print(f"check_guards: {problem}", file=sys.stderr)
        return 1
    print(
        f"check_guards: {len(declared)} guard(s) marked and declared, "
        f"{len(shapes)} bound shape(s) owned, {len(lean)} theorem(s) gated"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
