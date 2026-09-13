"""Run every validation case under ``validation/``.

These are the *external* checks. The tests elsewhere in the suite are generated
from the specs and check that the code does what the spec says; these check
whether the spec was right, because they come from outside the registry and can
disagree with it. See ``validation/README.md``.

A mismatch fails the build. That is the whole point: a validation case nobody
runs is a claim in a JSON file.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import _helpers as h
from azoth import hydraulics
from azoth._registry_gen import BY_ID

REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"

REQUIRED_KEYS = {"id", "calc", "source", "inputs", "expected", "tolerance"}

#: Statuses that oblige the case to explain itself. A case whose arithmetic is
#: ours must say so, or a passing test reads as an external validation it is not.
MUST_BE_EXPLAINED = {"unverified", "source_needed"}


def _case_paths() -> list[Path]:
    return sorted(p for p in VALIDATION_DIR.rglob("*.json"))


def _load(path: Path) -> dict[str, Any]:
    # Annotated on its own line: `json.loads` returns Any, and returning it
    # directly is an implicit Any at the boundary of every test here.
    case: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return case


CASES = _case_paths()


def test_there_is_at_least_one_case() -> None:
    """Guard against a refactor that silently stops running any of them."""
    assert CASES, f"no validation cases found under {VALIDATION_DIR}"


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_case_is_well_formed(path: Path) -> None:
    """A case must carry what the runner needs, and name a real calc.

    Checked separately from running it, so a malformed case reports as a
    malformed case rather than as an arithmetic failure.
    """
    case = _load(path)
    where = path.relative_to(REPO_ROOT)

    missing = REQUIRED_KEYS - set(case)
    assert not missing, f"{where} is missing {sorted(missing)}"
    assert case["calc"] in BY_ID, (
        f"{where} names calc {case['calc']!r}, which is not in the registry"
    )
    assert case["tolerance"] > 0, f"{where} has a non-positive tolerance"

    source = case["source"]
    assert source.get("verification") in {"verified", "unverified", "source_needed"}, (
        f"{where} has an unrecognised source.verification"
    )


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_unconfirmed_sources_are_explained(path: Path) -> None:
    """A case that is not backed by a confirmed source has to say why.

    Without this, `source.verification: source_needed` is a field nobody reads,
    and a green test tick reads as "validated against the standard" when it
    means "the arithmetic we wrote matches the arithmetic we wrote".
    """
    case = _load(path)
    where = path.relative_to(REPO_ROOT)
    status = case["source"]["verification"]

    if status in MUST_BE_EXPLAINED:
        notes = str(case["source"].get("notes", "")).strip()
        assert notes, (
            f"{where}: source.verification is {status!r} but there are no notes "
            f"explaining what is unconfirmed. A reader has to be able to tell this "
            f"from a case validated against a source."
        )


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_case_matches_implementation(path: Path) -> None:
    """Run the case and compare, output by output.

    Compares each named output separately so a failure says which one disagreed,
    by how much, and against what - rather than reporting that two dictionaries
    differ. The source's verification status is deliberately *not* a skip
    condition: the arithmetic is worth checking either way, and a skip would
    report a case that ran as one that did not.
    """
    case = _load(path)
    calc = BY_ID[case["calc"]]
    where = path.relative_to(REPO_ROOT)

    _, _, function_name = case["calc"].rpartition(".")
    result = getattr(hydraulics, function_name)(**h.kwargs_for(calc, case["inputs"]))

    for name, want in case["expected"].items():
        assert hasattr(result, name), (
            f"{where}: expected output {name!r} is not a field of {type(result).__name__}"
        )
        got = getattr(result, name)
        if hasattr(got, "magnitude"):
            got = got.magnitude
        h.assert_close(float(got), float(want), float(case["tolerance"]), f"{where} ({name})")
