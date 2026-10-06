"""The Lean tree and the pages that describe it are held to each other, and the auditor is real.

`tools/audit_lean_claims.py` closes the gap `test_lean_claims.py` states in its own docstring -
that the correspondence between a theorem and the prose claim it answers is "held together by a
reader rather than by a machine". The audit is the machine half; this file is what says it is
still running and still finding things.

Two assertions, and the second is the one this project always needs: an audit that read nothing
would pass every check in it by having nothing to check.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_the_audit_finds_the_pages_and_the_tree_agreeing() -> None:
    """The gate `spec-validate` runs, run here too so a push is not the first thing to see it."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "audit_lean_claims.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the pages and the Lean tree disagree\n{result.stdout}{result.stderr}"
    )


def test_the_audit_reads_something() -> None:
    """**A vacuous audit passes by checking nothing.**

    Every rule in the tool is satisfied by an empty page set and an empty declaration index, so
    the counts are asserted here rather than trusted. The declaration side matters as much as the
    reference side: a namespace parser that stopped working would make every reference
    unresolvable, which fails loudly - but one that returned *everything* would make every
    reference resolve, which is the quiet direction.
    """
    audit = _tools_module("audit_lean_claims")

    failures, checked = audit.audit()
    assert checked > 50, (
        f"{checked} reference(s) read - the tool believes it checked the pages, and a walk that "
        f"found nothing would pass every rule it has"
    )
    assert not failures, failures

    assert audit.declared(), "no declaration was found under lean/Azoth/ - the parser is broken"
    assert audit.gated(), "no `#print axioms` line was read - the gate files are not being read"
    assert "Azoth.Capability.run_deterministic" in audit.declared()
    assert "Azoth.Capability.run_deterministic" in audit.gated()
