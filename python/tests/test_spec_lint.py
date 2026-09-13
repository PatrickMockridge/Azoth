"""spec_lint's worked-example rule, which is the one requirement it kept.

The library used to carry a ``verification.status`` — ``verified``, ``unverified`` or
``source_needed`` — and a rule coupling it to whether a calc's worked example could run.
The status is gone: provenance is the engineer's business, not a YAML linter's, and a
field nobody can enforce honestly trains people to fill it in rather than to know the
answer.

What survives is the rule underneath it, which was never about provenance at all:
**a calculation ships a runnable worked example.** It costs a dozen lines, it is the
cheapest check in the registry, and it is what pins a number to something a reader can
retrace. Skipping it is allowed, but the test has to say why.

These tests pin both directions. A test that only asserted "a skipped example lints
clean" would pass just as happily if the requirement had been deleted outright, so the
last case asserts the requirement still bites when nothing is declared at all.

The specs are built by mutating a copy of a real spec that already passes, so a failure
here is about the rule under test rather than about some unrelated field being
malformed.
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
    """Mark the worked-example test skipped, as a calc with nothing to check must."""
    for test in spec["tests"]:
        if test["type"] == "worked_example":
            test["status"] = "skipped"
            test["skip_reason"] = reason


def as_skipped_with_a_reason(spec: dict[str, Any]) -> None:
    skip_worked_example(spec, "the arithmetic is not reproducible on a machine")


def as_skipped_without_a_reason(spec: dict[str, Any]) -> None:
    for test in spec["tests"]:
        if test["type"] == "worked_example":
            test["status"] = "skipped"
            test.pop("skip_reason", None)


def as_having_no_worked_example(spec: dict[str, Any]) -> None:
    spec["tests"] = [t for t in spec["tests"] if t["type"] != "worked_example"]


def test_an_active_worked_example_lints_clean(tmp_path: Path) -> None:
    """The baseline. Every other case here is a mutation of a spec that passes."""
    result = lint(spec_tree(tmp_path))
    assert result.returncode == 0, (
        f"the unmutated spec must lint clean:\n{result.stdout}{result.stderr}"
    )


def test_skipping_the_example_with_a_reason_is_allowed(tmp_path: Path) -> None:
    """A calc with nothing checkable may say so and move on.

    This is the case the old ``source_needed`` status existed to permit and then
    failed to, because two rules contradicted each other and no spec ever used the
    status, so nothing noticed. The permission is real now: the reason is the whole
    requirement.
    """
    result = lint(spec_tree(tmp_path, as_skipped_with_a_reason))
    assert result.returncode == 0, (
        f"a skipped worked example with a recorded reason must lint clean:\n"
        f"{result.stdout}{result.stderr}"
    )


def test_skipping_the_example_without_a_reason_fails(tmp_path: Path) -> None:
    """Skipping is legal; a silent skip is not. The schema requires the reason."""
    result = lint(spec_tree(tmp_path, as_skipped_without_a_reason))
    assert result.returncode != 0, (
        f"a skipped worked example with no skip_reason must fail:\n{result.stdout}{result.stderr}"
    )


def test_a_calc_with_no_worked_example_at_all_fails(tmp_path: Path) -> None:
    """The requirement still bites, which is what makes the permission mean anything.

    Without this case the whole file would pass if the check had been deleted rather
    than relaxed - the failure mode the old test was written to avoid and which is
    just as easy to hit now.
    """
    result = lint(spec_tree(tmp_path, as_having_no_worked_example))
    assert result.returncode != 0, (
        f"a calc with no worked_example test of any kind must fail:\n{result.stdout}{result.stderr}"
    )
    assert "worked_example" in result.stdout + result.stderr, (
        "the failure must name the missing worked example rather than being a generic schema error"
    )
