"""The checks that only ever run in CI.

`tools/check_links.py` and `tools/provenance.py` are both wired into workflows and
neither has a test, which is the shape of a check that breaks silently: it is not run
by the suite, so the first person to see it fail is whoever pushes next.

What is tested here is the part that can be *wrong* rather than the part that is
obviously right - link resolution, which has three cases that are easy to conflate
(a page that exists, a page that does not, and a link that is not a page at all) -
and the shape of the provenance record, which is consumed by a release rather than by
a reader.
"""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


def link_checker() -> ModuleType:
    """`tools/check_links.py`, imported by name.

    The tools are not a package and are outside mypy's `files`, so the import is
    dynamic - the same arrangement `test_prose_lint.py` uses.
    """
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_links")
    finally:
        sys.path.pop(0)


def check(tmp_path: Path, text: str) -> tuple[int, list[str]]:
    """Run the link checker over one page, rooted at a temporary directory.

    The root is redirected because the error message is relative to it; a page under
    `/tmp` against the real root would raise from `relative_to` rather than report.
    """
    # `Any` at the point the tool is reached, rather than a widened signature:
    # `ModuleType` declares no `ROOT` or `check_file`, and the tools are outside
    # mypy's `files`.
    tool: Any = link_checker()
    tool.ROOT = tmp_path
    page = tmp_path / "page.md"
    page.write_text(text, encoding="utf-8")
    errors: list[str] = []
    checked: int = tool.check_file(page, errors)
    return checked, errors


def test_a_link_to_a_page_that_exists_is_checked_and_passes(tmp_path: Path) -> None:
    (tmp_path / "target.md").write_text("# Target\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [the target](./target.md)\n")

    assert checked == 1
    assert errors == []


def test_a_link_to_a_page_that_does_not_exist_is_reported(tmp_path: Path) -> None:
    checked, errors = check(tmp_path, "see [nothing](./nowhere.md)\n")

    assert checked == 1, "a broken link is still a link the checker examined"
    assert len(errors) == 1, errors
    assert "nowhere.md" in errors[0], errors[0]
    assert "page.md:1" in errors[0], f"the line number is part of the report: {errors[0]}"


def test_a_link_into_a_directory_resolves_to_its_index(tmp_path: Path) -> None:
    """`./guide/` is a page in a book, and mdBook renders its `index.md`."""
    (tmp_path / "guide").mkdir()
    (tmp_path / "guide" / "index.md").write_text("# Guide\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [the guide](./guide/)\n")

    assert checked == 1
    assert errors == []


def test_external_links_and_bare_anchors_are_not_resolved(tmp_path: Path) -> None:
    """Two different cases, and the count distinguishes them.

    An external link is skipped outright - fetching it would make the docs build
    depend on someone else's server. A bare anchor is *examined* and then needs no
    file, because it points at a heading on the page it is written on. Neither can
    produce an error, and only the second is counted as a link looked at.
    """
    checked, errors = check(
        tmp_path,
        "see [the docs](https://example.invalid/x) and [below](#a-heading)\n",
    )

    assert checked == 1, "the same-page anchor is examined; the external link is skipped"
    assert errors == []


def test_a_link_with_an_anchor_resolves_on_its_file(tmp_path: Path) -> None:
    (tmp_path / "target.md").write_text("# Target\n", encoding="utf-8")

    checked, errors = check(tmp_path, "see [a heading](./target.md#a-heading)\n")

    assert checked == 1
    assert errors == []


def test_the_link_checker_passes_on_this_tree() -> None:
    """The gate itself, run where it runs in CI."""
    result = subprocess.run(
        [sys.executable, str(TOOLS / "check_links.py")],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, f"{result.stdout}{result.stderr}"


def test_provenance_records_a_verifiable_shape(tmp_path: Path) -> None:
    """The record a release is made from, written and read back."""
    out = tmp_path / "provenance.json"

    result = subprocess.run(
        [sys.executable, str(TOOLS / "provenance.py"), "--out", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"{result.stdout}{result.stderr}"

    record = json.loads(out.read_text(encoding="utf-8"))
    for key in ("schema_version", "git", "calcs", "models", "data", "licences"):
        assert key in record, f"the record has no {key!r}: {sorted(record)}"
    assert record["git"]["commit"], "a record with no commit names no revision"
    assert record["calcs"], "the record lists no calculations"

    # And it must verify against the tree it was written from, or the check the
    # release job runs would fail on the record that job just produced.
    verified = subprocess.run(
        [sys.executable, str(TOOLS / "provenance.py"), "--verify", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert verified.returncode == 0, f"{verified.stdout}{verified.stderr}"
