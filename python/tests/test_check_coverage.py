"""The one fact diff-coverage asserts, from fixtures rather than from an instrumented run.

`tools/check_coverage.py` answers a question `--verify-selection` cannot: the emitted targets
may all collect tests and still touch none of the code the diff changed. What it asserts is not
a percentage but one fact per changed file - at least one line of it ran - and that decision is
pure, so it is pinned here.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

TOOLS = Path(__file__).resolve().parents[2] / "tools"


def tool() -> ModuleType:
    """`tools/check_coverage.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_coverage")
    finally:
        sys.path.pop(0)


def document(**files: list[int] | None) -> dict[str, object]:
    """A coverage document with one entry per named file; `None` means it never ran."""
    return {
        "files": {
            name: {"executed_lines": lines or [], "summary": {"covered_lines": len(lines or [])}}
            for name, lines in files.items()
        }
    }


def test_a_file_that_ran_no_line_is_the_whole_finding() -> None:
    executed = tool().executed_files(document(**{"python/src/azoth/eos/flash.py": [1, 2, 3]}))
    assert executed == {"python/src/azoth/eos/flash.py"}


def test_a_file_with_an_empty_executed_list_counts_as_unrun() -> None:
    """The shape that matters: present in the document, and none of it executed."""
    assert tool().executed_files(document(**{"python/src/azoth/eos/flash.py": []})) == set()


def test_a_changed_file_that_did_not_run_is_reported() -> None:
    changed = ["python/src/azoth/eos/flash.py", "python/src/azoth/eos/other.py"]
    executed = {"python/src/azoth/eos/flash.py"}
    assert tool().files_with_no_line_run(changed, executed) == ["python/src/azoth/eos/other.py"]


def test_nothing_missing_is_the_good_answer() -> None:
    changed = ["python/src/azoth/eos/flash.py"]
    assert tool().files_with_no_line_run(changed, {"python/src/azoth/eos/flash.py"}) == []


def test_a_document_with_no_files_is_not_a_pass() -> None:
    """A malformed document must raise rather than make every file look covered."""
    with pytest.raises(ValueError):
        tool().executed_files({"totals": {}})


def test_a_run_that_did_not_happen_is_two(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """**An absent document is not a covered diff.** pytest-cov not having run has measured
    nothing, and returning 0 would make the whole check a no-op that looks green."""
    module = tool()
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["check_coverage.py", "--check", "--json", "missing.json"])
    assert module.main() == 2


def test_a_written_document_with_a_missing_file_is_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The finding itself, through `main`: one changed file, none of it executed."""
    module = tool()
    (tmp_path / "coverage.json").write_text(
        json.dumps(document(**{"python/src/azoth/eos/flash.py": [1]})), encoding="utf-8"
    )
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["check_coverage.py", "--check", "--json", "coverage.json", "--stdin"],
    )
    monkeypatch.setattr(
        sys, "stdin", type("S", (), {"__iter__": staticmethod(lambda: iter(["a.py\n"]))})()
    )
    assert module.main() == 1


def test_a_type_stub_is_not_demanded() -> None:
    """**A `.pyi` has no executable line, so "a line of it ran" is unsatisfiable.**

    Measured on the first push that regenerated `_core.pyi`: the diff-coverage step failed on the
    stub and on no code at all. The filter is a pure function here so the rule has a test.
    """
    assert tool().executed_only(
        ["python/src/azoth/eos/flash.py", "python/src/azoth/_core.pyi"]
    ) == ["python/src/azoth/eos/flash.py"]
