"""spec_lint's rule about verification status and runnable worked examples.

The defect these tests exist for: ``source_needed`` was a documented,
schema-legal, contributor-facing status that no spec could actually hold.

``check_tests`` demanded an active ``worked_example`` test unconditionally, and
the rule immediately below it forbade an active ``worked_example`` for exactly
the calcs that are permitted to omit one. The two rules were contradictory, so
every ``source_needed`` spec failed lint - and because no spec in the registry
used the status, nothing ever noticed. CONTRIBUTING.md told contributors to
reach for a state the tooling rejected.

A test that only asserted "a source_needed spec lints clean" would pass just as
happily if the worked-example requirement had been deleted outright. So these
pin both directions: the status is reachable, **and** the check it relaxes still
bites for every other status, and still bites for the contradiction it was
guarding against.

The specs are built by mutating a copy of a real spec that already passes, so a
failure here is about the rule under test rather than about some unrelated field
being malformed.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC_LINT = REPO_ROOT / "tools" / "spec_lint.py"
SOURCE_SPEC = REPO_ROOT / "specs" / "calcs" / "hydraulics" / "reynolds_number.yaml"

Mutator = Callable[[dict[str, Any]], None]


def lint(spec_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run spec_lint over a spec tree, returning the completed process."""
    return subprocess.run(
        [sys.executable, str(SPEC_LINT), "--spec-dir", str(spec_dir), "--quiet"],
        capture_output=True,
        text=True,
        check=False,
    )


def spec_tree(tmp_path: Path, mutate: Mutator | None = None) -> Path:
    """Lay out one spec where the linter expects to find it, optionally mutated.

    The filename has to match the id's final segment or `check_identity` objects,
    so the copy keeps the original name under its own namespace directory.
    """
    spec: dict[str, Any] = yaml.safe_load(SOURCE_SPEC.read_text(encoding="utf-8"))
    if mutate is not None:
        mutate(spec)
    directory = tmp_path / "hydraulics"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / SOURCE_SPEC.name).write_text(
        yaml.safe_dump(spec, sort_keys=False), encoding="utf-8"
    )
    return tmp_path


def skip_worked_example(spec: dict[str, Any], reason: str) -> None:
    """Mark the worked-example test skipped, as a source-less calc must."""
    for test in spec["tests"]:
        if test["type"] == "worked_example":
            test["status"] = "skipped"
            test["skip_reason"] = reason


def as_source_needed(spec: dict[str, Any]) -> None:
    """No adequate source found, so the worked example is skipped."""
    spec["verification"]["status"] = "source_needed"
    skip_worked_example(spec, "source_needed: no adequate source located for this calc")


def as_unverified_without_an_example(spec: dict[str, Any]) -> None:
    """The status that does *not* excuse a missing runnable example."""
    spec["verification"]["status"] = "unverified"
    skip_worked_example(spec, "source_needed: this excuse does not apply here")


def as_source_needed_with_an_active_example(spec: dict[str, Any]) -> None:
    """The contradiction: no source found, yet claiming a passing example."""
    spec["verification"]["status"] = "source_needed"


def as_source_needed_with_nothing_active(spec: dict[str, Any]) -> None:
    """Every test skipped, so the spec is exercised by nothing at all."""
    spec["verification"]["status"] = "source_needed"
    for test in spec["tests"]:
        test["status"] = "skipped"
        test["skip_reason"] = "source_needed: no adequate source located for this calc"


def test_source_needed_with_a_skipped_example_lints_clean(tmp_path: Path) -> None:
    """The status is reachable: this is the regression test for the defect.

    Under the old rules this spec failed, because a skipped worked example left
    `active_worked_example` false and nothing consulted the status before
    demanding one. `source_needed` was therefore write-only documentation.
    """
    result = lint(spec_tree(tmp_path, as_source_needed))
    assert result.returncode == 0, (
        f"a source_needed spec with a skipped worked example must lint clean, "
        f"but spec_lint failed:\n{result.stdout}{result.stderr}"
    )


def test_relaxing_source_needed_does_not_relax_anything_else(tmp_path: Path) -> None:
    """The exemption is narrow: `unverified` still requires a runnable example.

    Without this, deleting the worked-example requirement entirely would leave
    the test above passing.
    """
    result = lint(spec_tree(tmp_path, as_unverified_without_an_example))
    assert result.returncode != 0, (
        "an unverified spec with no active worked example must still fail; the "
        f"source_needed exemption has leaked:\n{result.stdout}{result.stderr}"
    )
    assert "no active worked_example test" in result.stderr


def test_source_needed_with_an_active_example_still_fails(tmp_path: Path) -> None:
    """The other half of the contradiction is still guarded.

    `source_needed` says no adequate source was found. An active worked example
    asserts the opposite, and the pair must not both be true.
    """
    result = lint(spec_tree(tmp_path, as_source_needed_with_an_active_example))
    assert result.returncode != 0, (
        "a source_needed spec claiming an active worked example must fail:\n"
        f"{result.stdout}{result.stderr}"
    )
    assert "must not claim a passing example" in result.stderr


def test_source_needed_with_no_active_test_is_reported(tmp_path: Path) -> None:
    """A calc nothing exercises should be loud, though it is not an error.

    Skipping the worked example is legal under `source_needed`, and
    `tests: minItems: 1` is satisfied by the skipped test alone - so a spec can
    reach the registry with nothing running at all. That is a warning rather than
    a failure, because there are honest reasons to record a calc in that state,
    but it must not pass silently.
    """
    result = lint(spec_tree(tmp_path, as_source_needed_with_nothing_active))
    assert result.returncode == 0, (
        f"an untested source_needed spec is a warning, not an error:\n"
        f"{result.stdout}{result.stderr}"
    )
    assert "no test of any kind is active" in result.stderr


def test_the_harness_lints_the_unmutated_spec_clean(tmp_path: Path) -> None:
    """Sanity check on the harness itself: the unmutated spec still passes.

    Guards against these tests passing for the wrong reason. If copying a spec
    into a temporary tree failed lint on its own - a mangled round-trip, say -
    then a non-zero exit above would prove nothing about the mutation.
    """
    result = lint(spec_tree(tmp_path))
    assert result.returncode == 0, (
        f"an unmutated copy of a passing spec must lint clean, so a failure above "
        f"is attributable to the mutation:\n{result.stdout}{result.stderr}"
    )
