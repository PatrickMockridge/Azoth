"""`tools/pin_impact.py`: the trees it walks, and a fix it has to see.

The tool's first question is which NeqSim classes changed, and that is a function of the tree
list. **The list excluded `process` and `fluidmechanics` on the belief that azoth wrote its own
unit operations and hydraulics**; since P11/P12 it ports the palette and the executor, and since
P1 the packing geometry's hydraulics, so a process fix reported zero changed classes. The
Kirkbride fix (#4079, PR #4108) is the regression that shows it: one class changed, it is cited
by this tree, and before the list moved the tool saw none of it.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKOUT = Path(os.environ.get("NEQSIM_CHECKOUT", "/tmp/neqsim-check/neqsim"))
KIRKBRIDE_MERGE = "1f7581f914892878e69957fb37a49a54965b77d4"


def pin_impact() -> Any:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("pin_impact")
    finally:
        sys.path.pop(0)


def test_the_ported_process_and_packing_trees_are_in_scope() -> None:
    """A fix under `process` or `fluidmechanics` is a fix this tree can have to re-measure."""
    trees = pin_impact().TREES
    assert "src/main/java/neqsim/process" in trees
    assert "src/main/java/neqsim/fluidmechanics" in trees


@pytest.mark.skipif(
    not (CHECKOUT / ".git").exists(),
    reason="needs the NeqSim checkout; CI has none, so this runs where the oracle does",
)
def test_the_kirkbride_fix_is_seen() -> None:
    """PR #4108 changed `ShortcutDistillationColumn`, which this tree cites; the tool names it."""
    done = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools" / "pin_impact.py"),
            "--checkout",
            str(CHECKOUT),
            "--from",
            f"{KIRKBRIDE_MERGE}^",
            "--to",
            KIRKBRIDE_MERGE,
            "--all",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode == 0, done.stderr
    assert "ShortcutDistillationColumn" in done.stdout
