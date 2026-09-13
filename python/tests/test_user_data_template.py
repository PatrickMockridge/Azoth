"""The user data template must be valid, and must stay a template.

``azoth-data.example.yaml`` is the file a user copies to supply values from a
standard they licensed, and ``tools/check_user_data.py`` is what tells them
whether they filled it in correctly. Both are documentation as much as they are
code, which is what makes them worth testing.

The failures these tests exist for are not crashes:

* A template that does not pass its own checker teaches the wrong shape. Someone
  copies it, sees six errors they did not cause, and concludes the tool is broken
  rather than their file.
* A template row that claims ``verified`` is a citation nobody has checked, sitting
  in the one file whose stated purpose is to make that impossible. The schema and
  the checker both allow ``verified``; only a test can hold the *template* to
  placeholders, because for a real user's file ``verified`` is exactly right.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKER = REPO_ROOT / "tools" / "check_user_data.py"
TEMPLATE = REPO_ROOT / "azoth-data.example.yaml"


def run_checker(path: Path) -> subprocess.CompletedProcess[str]:
    """Run the checker over a file, returning the completed process."""
    return subprocess.run(
        [sys.executable, str(CHECKER), str(path)],
        capture_output=True,
        text=True,
        check=False,
    )


def template_document() -> dict[str, Any]:
    """The template, parsed."""
    document: dict[str, Any] = yaml.safe_load(TEMPLATE.read_text(encoding="utf-8"))
    return document


def test_the_template_exists() -> None:
    """A missing template is the failure mode this file is named for."""
    assert TEMPLATE.is_file(), (
        "azoth-data.example.yaml is missing. It is what docs/src/copyright.md tells "
        "a user to copy; without it that page describes a file that is not there."
    )


def test_the_template_passes_the_checker() -> None:
    """It validates as it stands.

    Enforced rather than assumed, because the template is the worked example of the
    format: if it is wrong, every file derived from it starts wrong.
    """
    result = run_checker(TEMPLATE)
    assert result.returncode == 0, (
        f"the template does not pass its own checker, so it teaches the wrong shape:\n"
        f"{result.stdout}{result.stderr}"
    )


def test_every_row_in_the_template_is_a_placeholder() -> None:
    """No row may look sourced, because no row is.

    A template row cannot honestly cite anything - the reader has not supplied the
    source yet. The marker is now the citation itself: a row whose citation does not
    say it is a placeholder is a row that reads as real, and the worst case is one
    real enough to copy without reading.

    This used to be a `verify_status` field on every row. That field is gone: it
    asked the reader to make a claim about their own diligence and could not check
    the answer, which is a form rather than a provenance record. The citation says
    the same thing and is the thing a reader actually reads.
    """
    document = template_document()
    citations = [row["citation"] for row in document["fittings"]]
    for rows in document["fluids"].values():
        citations.extend(row["citation"] for row in rows)

    unmarked = [c for c in citations if "DUMMY" not in c.upper()]
    assert not unmarked, (
        f"the template carries {len(unmarked)} row(s) that do not say they are "
        f"placeholders: {unmarked}. Every template row must say so in its citation; "
        f"show the shape of a finished row in a comment instead."
    )


def test_the_checker_is_not_vacuous(tmp_path: Path) -> None:
    """A broken file must fail, or the passing test above proves nothing.

    Sabotages the template the way a user would break it by accident - a typo in the
    status column - and asserts the checker says so.

    The sabotage used to be promoting a row to `verified` without a `source_ref`. That
    requirement is gone: a row now names its source in its citation, in prose, and
    there is nothing for a tool to check about it. What is still checkable is that the
    status is one the pipeline understands.
    """
    document = template_document()
    document["fittings"][0]["verify_status"] = "probably_fine"

    broken = tmp_path / "broken.yaml"
    broken.write_text(yaml.safe_dump(document), encoding="utf-8")

    result = run_checker(broken)
    assert result.returncode != 0, "an unrecognised verify_status must fail"
    assert "verify_status" in result.stderr


def test_the_checker_rejects_an_unknown_schema_version(tmp_path: Path) -> None:
    """The version gate refuses rather than guessing.

    A future format may mean something different by the same key names, so reading
    it on the assumption that it does not is how a value gets silently reinterpreted.
    """
    document = template_document()
    document["schema_version"] = 99

    future = tmp_path / "future.yaml"
    future.write_text(yaml.safe_dump(document), encoding="utf-8")

    result = run_checker(future)
    assert result.returncode != 0
    assert "schema_version" in result.stderr


@pytest.mark.parametrize("section", ["fittings", "fluids"])
def test_the_template_documents_every_section_it_carries(section: str) -> None:
    """Each section exists, so a user knows where their data goes."""
    assert section in template_document(), (
        f"the template has no '{section}' section, so a user has nowhere to put that kind of data"
    )
