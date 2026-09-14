"""prose_lint, and the file it used to report: itself.

P5 names this tool as the mechanism for the prose standard, and the tool was wired into
nothing - so nothing ran it, and when it was run it failed. `SEARCH` included `tools/`,
and `HISTORY_PHRASES` is a literal list of the phrases being searched for, so the file
matched every entry in its own data block: seven of its nine hits were itself. Since it
exits non-zero on any hit, it could not have been enabled as it stood.

These tests pin both halves: the phrase list is not scanned as prose, and a violation
planted in a tree is still caught.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
PROSE_LINT = REPO_ROOT / "tools" / "prose_lint.py"

#: A line of history prose to plant, assembled rather than written out. This file is
#: itself under `python/tests`, which the tool reads, so a literal phrase here would be a
#: violation in the very tree `test_this_tree_is_clean` asserts is clean - which is how
#: this was found.
PLANTED_PHRASE = "was " + "renamed"
PLANTED_LINE = f"// the accumulator {PLANTED_PHRASE} blend\n"


def lint(tree: Path) -> subprocess.CompletedProcess[str]:
    """Run prose_lint over a tree, returning the completed process."""
    return subprocess.run(
        [sys.executable, str(PROSE_LINT), "--path", str(tree), "--quiet"],
        capture_output=True,
        text=True,
        check=False,
    )


def test_this_tree_is_clean() -> None:
    """The gate itself, run over the repository it guards."""
    result = lint(REPO_ROOT)
    assert result.returncode == 0, f"{result.stdout}{result.stderr}"


def test_a_planted_violation_is_caught(tmp_path: Path) -> None:
    """A tree with history prose in it fails, and says where.

    Without this the test above would pass just as happily on a tool that had been
    deleted, or on one whose phrase list had been emptied.
    """
    planted = tmp_path / "crates" / "planted.rs"
    planted.parent.mkdir(parents=True)
    planted.write_text(PLANTED_LINE, encoding="utf-8")

    result = lint(tmp_path)

    assert result.returncode == 1, "a planted history phrase was not caught"
    assert PLANTED_PHRASE in result.stdout, result.stdout
    assert "planted.rs" in result.stdout, result.stdout


def test_the_phrase_list_is_not_scanned_as_prose() -> None:
    """The tool's own data block is not in the set of files it reads.

    The measurement the docstring describes - that each phrase occurs nowhere in the
    tree - was made without counting the file the list is written in.
    """
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        tool: ModuleType = importlib.import_module("prose_lint")
    finally:
        sys.path.pop(0)

    assert tool.SELF not in tool.sources(tool.ROOT), (
        "prose_lint is scanning its own HISTORY_PHRASES literal, so it matches every "
        "phrase it looks for and fails on any tree"
    )
