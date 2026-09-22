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

**The prose above the sample is a third place the same counts are written**, and the
sample's check does not cover it: the sample holds `1497` and `37`, so a fixer who changes
the manifest updates the block and leaves the sentence. The sentence said "35 files" while
the manifest declared 37 and enumerated the columns of 31, and nothing read it - which is
the failure this module's first paragraph describes, one paragraph further down the page.
The last test here reads the sentence.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
README = REPO_ROOT / "databank" / "README.md"

#: The fenced block holding the sample.
_BLOCK = re.compile(r"```\n(?P<body>check_manifest: OK.*?)```", re.DOTALL)

#: Every integer in the sample, in the order it appears.
_NUMBERS = re.compile(r"\d+")

#: The page's prose claim, above the sample: how many files have their columns enumerated
#: and how many columns that is in total.
_PROSE_COUNTS = re.compile(
    r"(?P<files>\d+) files, (?P<columns>[\d,]+) columns"
)


def manifest_tool() -> ModuleType:
    """`tools/manifest.py`, imported by name because `tools/` is not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("manifest")
    finally:
        sys.path.pop(0)


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


def test_the_prose_counts_are_the_manifests_counts() -> None:
    """The sentence above the sample, held to the manifest it describes.

    It counted files while the manifest grew, and the sample below it did not cover the
    number: the block holds the *total* file count and never the enumerated subset, so
    the sentence is the only place the second number appears and nothing read it.
    """
    match = _PROSE_COUNTS.search(README.read_text(encoding="utf-8"))
    assert match, (
        f"{README.relative_to(REPO_ROOT)} no longer states its file and column counts in "
        f"the form this test reads. Either the sentence changed and this regex did not, "
        f"or the claim was removed - and a claim nothing checks is the one that goes stale."
    )

    files = manifest_tool().read()[0].files()
    enumerated = [f for f in files if f.columns_declared]
    columns = sum(len(f.columns) for f in enumerated)

    assert int(match.group("files")) == len(enumerated), (
        f"{README.relative_to(REPO_ROOT)} says {match.group('files')} files have their "
        f"columns enumerated; the manifest enumerates {len(enumerated)} of its "
        f"{len(files)}."
    )
    assert int(match.group("columns").replace(",", "")) == columns, (
        f"{README.relative_to(REPO_ROOT)} says {match.group('columns')} columns; the "
        f"manifest's enumerated files declare {columns}."
    )
