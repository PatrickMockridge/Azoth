"""The databank README's tally is what the checker prints, not a copy of it.

`databank/README.md` shows `tools/check_manifest.py`'s output as a sample, and the
numbers in it are the page's central claim: how much of NeqSim's data is vendored,
and how much of it anything reads. Nothing checked them, and they had gone stale -
the page said "1459 carried of which 14 read" with 1445 unread where the tool
printed 1463, 981 and 482. A count is cheap to check and expensive to notice, which
is the argument `test_registry_contract.py` already makes for the specification's
id count and the front page's model count.

**The tool is run, not reimplemented.** The page is a sample *of its output*, so the
only honest check is to produce the output and compare the numbers in order. A test
that recomputed the tally here would agree with a README that had drifted from the
tool as long as the two drifted together, which is the failure mode this is for.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
README = REPO_ROOT / "databank" / "README.md"

#: The fenced block holding the sample.
_BLOCK = re.compile(r"```\n(?P<body>check_manifest: OK.*?)```", re.DOTALL)

#: Every integer in the sample, in the order it appears.
_NUMBERS = re.compile(r"\d+")


def sample_numbers() -> list[int]:
    """The numbers the page shows, in order."""
    match = _BLOCK.search(README.read_text(encoding="utf-8"))
    assert match, (
        f"{README.relative_to(REPO_ROOT)} no longer carries a run of `check_manifest`'s "
        f"output in a fenced block. Either the page changed and this regex did not, or "
        f"the sample was removed - and a claim nothing checks is the one that goes stale."
    )
    return [int(n) for n in _NUMBERS.findall(match.group("body"))]


def test_the_readme_shows_what_the_checker_prints() -> None:
    result = subprocess.run(
        [sys.executable, "tools/check_manifest.py"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    printed = [int(n) for n in _NUMBERS.findall(result.stdout)]
    assert printed, "check_manifest printed no numbers, so this test cannot compare"

    assert sample_numbers() == printed, (
        f"{README.relative_to(REPO_ROOT)} shows {sample_numbers()} and check_manifest "
        f"prints {printed}. Update the sample from a fresh run."
    )
