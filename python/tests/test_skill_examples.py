"""The eighteen example scripts the skills ship, run as programs.

Every skill page shows a worked example, and eighteen of them ship a runnable script beside the
prose. **Nothing ran them.** `pytest skills` collects the skills' own `test_*.py` files and not
these, so a page could document an API that no longer exists and the only thing that would
notice is a reader.

Each runs in its own interpreter, because that is how a reader runs it - what this catches is a
script that does not work, not a module that does not import.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS = REPO_ROOT / "skills"
EXAMPLES = sorted(SKILLS.rglob("examples/*.py"))


def example_id(path: Path) -> str:
    """`eos/cubic/pr_cubic` - the skill's own path, so a failure names the page it is on."""
    parts = path.relative_to(SKILLS).parts
    return "/".join((*parts[:-2], parts[-1].removesuffix(".py")))


def test_there_are_examples_left_to_run() -> None:
    """A glob that stopped matching would make every test below it pass by not existing."""
    assert len(EXAMPLES) >= 18, [str(path.relative_to(REPO_ROOT)) for path in EXAMPLES]


@pytest.mark.parametrize("path", EXAMPLES, ids=example_id)
def test_the_example_runs(path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python" / "src")},
    )
    assert result.returncode == 0, f"{path.name} failed:\n{result.stdout}\n{result.stderr}"
    assert result.stdout.strip(), f"{path.name} printed nothing, so it did nothing"
