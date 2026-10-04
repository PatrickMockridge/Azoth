"""The executor's id-to-kernel table is generated, and this is the ratchet that says so.

`tools/gen_dispatch.py` emits `crates/azoth-process/src/executor/dispatch_gen.rs` from the
palette's own ids and the adapter functions `executor/dispatch.rs` declares. The table was a third
written copy of a join of two things already written, so a palette entry added without its line was
a flowsheet that failed at run time with "no kernel".

`docs-drift` regenerates and diffs it; `test_the_committed_table_is_current` runs the same `--check`
in the suite, and `test_every_palette_entry_is_dispatched_or_refused` drives the derivation, because
the count that matters is the palette's.
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


def test_the_committed_table_is_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_dispatch.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the dispatch table is out of date - run `python tools/gen_dispatch.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_table_is_the_palette_less_what_it_refuses() -> None:
    """Every palette entry is either dispatched or named in `UNRUNNABLE`, and not both.

    Driven rather than read from the generated file: the claim is that the *derivation* accounts
    for the palette, so a palette entry whose leaf stops matching an adapter shows up here rather
    than only as a byte diff.
    """
    generator = _tools_module("gen_dispatch")
    ids = generator.palette_ids()
    assert ids, "the palette parsed to no ids"

    dispatched = {
        line.split('"')[1] for line in generator.emit().splitlines() if line.startswith('    ("')
    }
    assert dispatched <= set(ids), "the table names something the palette does not declare"

    refused = {
        line.split('"')[1]
        for line in (REPO_ROOT / "crates" / "azoth-process" / "src" / "executor" / "dispatch.rs")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.startswith('    "unit_ops.')
    }
    assert not dispatched & refused, f"both dispatched and refused: {sorted(dispatched & refused)}"
    assert dispatched | refused == set(ids), (
        "the palette has entries the executor neither runs nor refuses: "
        f"{sorted(set(ids) - dispatched - refused)}"
    )
