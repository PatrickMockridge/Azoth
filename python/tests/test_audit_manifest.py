"""The manifest's unread rows are read by the tree, and the auditor is real.

`tools/audit_manifest.py` answers the one question about `databank/manifest.toml` that
`tools/manifest.py` cannot: a column recorded `not-ported` or `unread-upstream` claims *nothing
reads this*, and the tree is what decides. Four rows have been corrected by hand after going false -
`PARACHOR`, `LJEPS` and the three `LIQUIDCONDUCTIVITY` columns - which is why the claim has a
reader now.

Three assertions, and the last two are the ones this repository always needs: an auditor that
mapped no alias would pass every rule by having nothing to check, and one that mapped every alias
to the wrong field would too.
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


def test_the_gate_spec_validate_runs() -> None:
    """The check CI runs, run here too so a push is not the first thing to see it."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "audit_manifest.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the manifest says nothing reads these\n{result.stdout}{result.stderr}"
    )


def test_the_audit_reads_a_mapping_and_a_manifest() -> None:
    """**A walk that found nothing passes every rule it has.**

    Both directions: the alias-to-field mapping has to be non-empty, and the manifest has to
    contain rows that claim to be unread. Either one empty and the audit reports success while
    having checked nothing - which is the failure `check_json_keys.py` calls a check that did not
    run.
    """
    audit = _tools_module("audit_manifest")

    mapping = audit.carried_as()
    assert len(mapping) > 50, (
        f"{len(mapping)} alias(es) mapped from the readers - the walk is looking in the wrong "
        f"place, or the readers have been restructured"
    )
    # The two the tool exists for: an alias landing in a field whose name is *not* the alias, one
    # from each of the two shapes the walk has to handle. `ljeps` is the single-line form and
    # `liquidconductivity1` is inside the array form, where the field opens a block.
    assert mapping.get("ljeps") == "lennard_jones_energy"
    assert mapping.get("liquidconductivity1") == "liquid_conductivity"

    rows, read = audit.verdicts()
    assert read, "no column claims to be unread, so the manifest's rows are not being read"
    assert len(rows) == read


def test_the_rows_that_went_false_are_no_longer_claims() -> None:
    """**The correction is that the row stops making the claim**, and this is what says so.

    `LJEPS`, the three `LIQUIDCONDUCTIVITY` columns and `PARACHOR` each said `not-ported` while
    the tree read them. Their reasons now name the reader, which takes them out of this audit's
    scope - it reads a *claim* of unreadness, and they no longer make one. So the ratchet is that
    they are absent from the rows, which fails if the correction is reverted.
    """
    audit = _tools_module("audit_manifest")
    rows, _ = audit.verdicts()
    claimed = {row["where"] for row in rows}

    for where in (
        "neqsim/COMP.csv.LJEPS",
        "neqsim/COMP.csv.LIQUIDCONDUCTIVITY1",
        "neqsim/COMP.csv.LIQUIDCONDUCTIVITY2",
        "neqsim/COMP.csv.LIQUIDCONDUCTIVITY3",
        "neqsim/COMP.csv.PARACHOR",
    ):
        assert where not in claimed, (
            f"{where} claims to be unread again, and the tree reads it - the audit's `--check` "
            f"fails on this, and so does the reason it was corrected for"
        )

    # And the tool is not passing by reporting everything as read: a column that genuinely is
    # unread is still there, parsed into no field at all.
    unread = {row["where"]: row for row in rows}
    assert unread["neqsim/COMP.csv.criticalViscosity"]["verdict"] == "unparsed", (
        "criticalViscosity is the friction-theory viscosity method, which is not ported and is "
        "not parsed into the entry - a verdict of `READ` here means the search is too loose"
    )
