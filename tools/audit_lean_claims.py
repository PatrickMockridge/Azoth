#!/usr/bin/env python3
"""Every claim the calculus pages make about the Lean tree, checked against the tree.

`python/tests/test_lean_claims.py` states its own limit in its docstring: "the pages carry the
statement in prose and each names its theorem, so the two are held together by a reader rather
than by a machine". This closes the part of that which is mechanical, and says which part is not.

Three rules, and the third is the one that was missing:

* **A module a page names exists.** The layer table in `docs/src/calculus/index.md` gives every
  layer a Lean module; a row naming a module that is not there is a row pointing at nothing.
* **A reference in a claim the page calls proved resolves.** Every `Azoth.Module.name` in the two
  normative trees must be a declaration under `lean/Azoth/`, or a module that exists - *unless* the
  paragraph carrying it says **specified** or **characterised**. Those two statuses mean the layer
  does not exist, so a name in one is a name the tree does not have yet, deliberately. Reading the
  status is the whole of the classification: a rule that ignored it would report every
  specification as a missing theorem, and a rule that ignored the names would check nothing.
* **A name a page calls proved is gated.** `Azoth.Capability.run_deterministic` named in a
  `*Status: **proved**, by ...*` block must have a `#print axioms` line. `check_lean_axioms.py`
  proves a *gated* theorem is complete; nothing until now said that the theorems the pages call
  proved are the gated ones.

**What this still does not do, and cannot**: check that a theorem is the claim the page states. A
proof of a weakened hypothesis has no gap and is still not the claim a reader expects. This tool
finds a name that is missing and a claim that is ungated; whether the statement matches the prose
is what review is for, and saying so here is better than a check that appears to do it.

    python tools/audit_lean_claims.py            # the audit, as a report
    python tools/audit_lean_claims.py --check    # fail on a claim the tree does not hold
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEAN = ROOT / "lean" / "Azoth"
GATES = (LEAN / "Axioms.lean", LEAN / "Gate.lean", LEAN / "GuardGate.lean")

#: The two trees that are normative for the types and the surface. The same pair
#: `check_doc_claims.swept_pages` sweeps, and for the same reason.
TREES = ("calculus", "architecture")

#: A theorem reference: `Azoth.Dim.exponents_smul`, `Azoth.View.Shift.toSiFromSi`.
_REFERENCE = re.compile(r"`(Azoth(?:\.[A-Za-z_][\w']*)+)`")

#: A module path: `Azoth/Dim.lean`.
_MODULE = re.compile(r"`(Azoth/[A-Za-z_]+\.lean)`")

#: The statuses that mean the layer does not exist, so a name under one is a name the tree
#: does not have yet. Everything else - including no status at all - is required to resolve.
#:
#: The word is matched inside a bold span rather than as a bare `**specified**`, because the
#: pages write it as a sentence: `index.md` says "**The balance is specified.**" and a rule
#: looking for the two asterisks alone reads that block as claiming a proof.
_UNPROVED = re.compile(r"\*\*[^*]*\b(?:specified|characterised)\b[^*]*\*\*", re.IGNORECASE)

#: What a page may name. `theorem` and `lemma` are what a claim is *proved* by; `def`, `abbrev`,
#: `structure` and `inductive` are what it is *about* - `capability.md` names `Answer` and
#: `Computation`, both structures, and the pages are allowed to talk about the objects a claim
#: ranges over as well as the theorem that discharges it.
_DECLARED = re.compile(
    r"^\s*(?:theorem|lemma|def|abbrev|structure|inductive)\s+(?P<name>[A-Za-z0-9_'À-ÿ.]+)",
    re.MULTILINE,
)
_PRINTED = re.compile(r"^\s*#print axioms\s+(?P<name>[A-Za-z0-9_.]+)\s*$", re.MULTILINE)


def declared() -> set[str]:
    """Every theorem or definition under `lean/Azoth/`, fully qualified.

    A namespace stack rather than a variable: `Dim.lean` nests `Azoth` and `Dim`, and a single
    variable would attribute `ofExponents_nil` to `Azoth`.
    """
    names: set[str] = set()
    for path in sorted(LEAN.glob("*.lean")):
        open_namespaces: list[str] = []
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("namespace "):
                open_namespaces.append(stripped.removeprefix("namespace ").strip())
                continue
            if stripped.startswith("end ") and open_namespaces:
                open_namespaces.pop()
                continue
            match = _DECLARED.match(line)
            if match:
                names.add(".".join([*open_namespaces, match.group("name")]))
    return names


def gated() -> set[str]:
    """Every theorem the gate files ask Lean to report the axioms of."""
    names: set[str] = set()
    for gate in GATES:
        names.update(_PRINTED.findall(gate.read_text(encoding="utf-8")))
    return names


def pages() -> list[Path]:
    out: list[Path] = []
    for tree in TREES:
        out += sorted((ROOT / "docs" / "src" / tree).glob("*.md"))
    return out


def paragraphs(text: str) -> list[str]:
    """The page as blank-line-separated blocks, which is the unit a status applies to."""
    return re.split(r"\n\s*\n", text)


def audit() -> tuple[list[str], int]:
    """`(failures, checked)`. Every reference in the two normative trees, classified."""
    names = declared()
    gated_names = gated()
    modules = {path.stem for path in LEAN.glob("*.lean")}
    failures: list[str] = []
    checked = 0

    for page in pages():
        where = page.relative_to(ROOT).as_posix()
        text = page.read_text(encoding="utf-8")

        for module in _MODULE.findall(text):
            checked += 1
            if not (ROOT / "lean" / module).exists():
                failures.append(f"{where}: names the module `{module}`, which does not exist")

        for block in paragraphs(text):
            unproved = bool(_UNPROVED.search(block))
            proved = "**proved**" in block
            references = _REFERENCE.findall(block)
            checked += len(references)
            for reference in references:
                leaf = reference.removeprefix("Azoth.")
                # A module reference (`Azoth.Dim`) is the module; a dotted one is a declaration.
                if reference in names or leaf in modules:
                    continue
                if unproved:
                    continue
                failures.append(
                    f"{where}: `{reference}` is named and no file under lean/Azoth/ declares "
                    f"it, in a block that claims neither **specified** nor **characterised**"
                )
            # **A block that claims a proof has to name the proof.** Not "every name in it is
            # gated": a proved claim legitimately names the *definitions* it is about beside the
            # theorem that proves it - `view.md`'s drag claim names `setPosition`, which is a
            # `def`. What must not happen is a block calling a claim proved and naming no gated
            # theorem at all, which is a proof pointed at nothing.
            #
            # **What counts as a claim** is a status block or a paragraph that names something.
            # `index.md`'s status *vocabulary* - `| **proved** | in lean/Azoth/, and the axiom
            # gate covers it |` - is a paragraph carrying the word and making no claim, and a
            # rule that counted it would report the definitions table as an unheld proof.
            claim = block.lstrip().startswith("*Status:") or bool(references)
            if proved and claim and not any(reference in gated_names for reference in references):
                failures.append(
                    f"{where}: a block claims **proved** and names no gated theorem - the "
                    f"`#print axioms` lines in {', '.join(g.name for g in GATES)} cover none of "
                    f"{references}"
                )
    return failures, checked


def main() -> None:
    failures, checked = audit()
    if not checked:
        raise SystemExit(
            "audit_lean_claims: no claim was read - this is looking in the wrong place"
        )
    if "--check" in sys.argv:
        if failures:
            raise SystemExit(
                "audit_lean_claims: the pages and the tree disagree\n  " + "\n  ".join(failures)
            )
        print(f"audit_lean_claims: {checked} reference(s) checked, every one held by the tree")
        return
    print(f"{checked} reference(s) across {len(pages())} page(s)")
    for failure in failures:
        print(f"  {failure}")
    print(f"{len(failures)} finding(s)")


if __name__ == "__main__":
    main()
