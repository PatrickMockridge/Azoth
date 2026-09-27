"""The front end's fixture recapture: its rules, with the CLI faked.

`tools/gen_ui_fixtures.py` shells out to `cargo run`, so a test that ran it for real would
build the workspace in the pytest job to re-measure what `--check` measures in `rust-test` and
what `ui/test/currency.test.ts` measures again in `ui`. So what is tested here is the part that
is *this tool's*: which exit codes it accepts, that it refuses a capture that wrote nothing,
and that `--check` compares rather than writes.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def tool() -> ModuleType:
    """`tools/gen_ui_fixtures.py`, imported by name. The tools are not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_ui_fixtures")
    finally:
        sys.path.pop(0)


def canned(monkeypatch: pytest.MonkeyPatch, code: int, out: str) -> None:
    """The CLI, as a process that answered with `code` and `out`."""

    def run(*_args: Any, **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(args=[], returncode=code, stdout=out, stderr="")

    monkeypatch.setattr(subprocess, "run", run)


def test_a_refused_document_is_a_capture_and_not_a_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """**`2` is `broken.json`'s whole content**, so a tool that refused it could not write the
    one fixture whose point is that a document can be broken."""
    module: Any = tool()
    assert 2 in module.ACCEPTED
    canned(monkeypatch, 2, '{"ok": false}\n')
    assert module.capture(["edit", "--json"]) == '{"ok": false}\n'


def test_another_exit_code_is_a_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    module: Any = tool()
    canned(monkeypatch, 101, "")
    with pytest.raises(SystemExit, match="neither a document nor a refused one"):
        module.capture(["forms"])


def test_a_capture_that_wrote_nothing_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    """A capture that wrote nothing would make every fixture empty and every mirror agree."""
    module: Any = tool()
    canned(monkeypatch, 0, "  \n")
    with pytest.raises(SystemExit, match="wrote nothing"):
        module.capture(["forms"])


def test_check_compares_against_the_committed_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module: Any = tool()
    monkeypatch.setattr(module, "FIXTURES", tmp_path)
    monkeypatch.setattr(module, "CASES", (("one.json", ["forms"]),))
    (tmp_path / "one.json").write_text("the same\n", encoding="utf-8")

    canned(monkeypatch, 0, "the same\n")
    monkeypatch.setattr(sys, "argv", ["gen_ui_fixtures.py", "--check"])
    module.main()

    canned(monkeypatch, 0, "something else\n")
    with pytest.raises(SystemExit) as raised:
        module.main()
    assert raised.value.code == 1
    # **And `--check` wrote nothing**, which is the property that lets it run in a job that
    # must not leave the tree dirty for the step after it.
    assert (tmp_path / "one.json").read_text(encoding="utf-8") == "the same\n"


def test_a_missing_fixture_is_stale_rather_than_an_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The first capture is the tool's to write, so a fixture that is not there yet is a
    difference to report and not a crash before the reason is printed."""
    module: Any = tool()
    monkeypatch.setattr(module, "FIXTURES", tmp_path)
    monkeypatch.setattr(module, "CASES", (("never-written.json", ["forms"]),))
    canned(monkeypatch, 0, "anything\n")
    monkeypatch.setattr(sys, "argv", ["gen_ui_fixtures.py", "--check"])
    with pytest.raises(SystemExit) as raised:
        module.main()
    assert raised.value.code == 1
