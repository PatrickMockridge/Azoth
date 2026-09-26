"""A column the component parser indexes is `used`, and for three tranches it was not.

`used` and `vendored` differ by one claim - something reads it - and that claim is the one
a manifest cannot check for itself. The parser can: `parse_components` holds a closed
allow-list of the names it looks up in every record, so being in it *is* the claim.

**Four columns were `vendored` while the parser indexed all four**, and the tally printed
them as carried-and-unread. `LJDIAMETER`, which both Furst phases read for the Born radius
and the shielding parameter; and `SCHWARTZENTRUBER1`-`3`, the fitted parameters water's
alpha moves 2.6% on, which three registered specs read - one of them naming those columns
as its source. The three had been that way since the alpha terms were ported and every gate
was green, because nothing compared the two sides of the claim.

**One file, one direction.** The allow-list belongs to `parse_components`, whose compiled
output is one table, and the reverse direction is false: `cas`,
`liquid_density_kg_per_m3` and `viscosity_correction_factor` are `used` and read by
something other than this parser. The tests below hold both edges - the sabotage has to
fire, and a column that no reader indexes has to be left alone.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]

RUST = Path("crates") / "azoth-eos" / "src" / "databank.rs"

#: The column the sabotage tests point the rule at. `PARACHOR` is genuinely unread by the
#: parser - the two parachor ids take it as an *input* rather than indexing it - which is why
#: it is still `vendored`, so a scratch allow-list naming it is the only thing that can make
#: it a violation.
#:
#: It was `LJEPS` until `eos.phase_transport` started indexing it, and then `PARACHOR` until the
#: rate-based packed column's segment model did: the sabotage has to name a column whose marking
#: is genuinely wrong, so a parser that starts reading one retires its own example.
UNREAD = "neqsim/COMP.csv.TRIPLEPOINTDENSITY"

#: What the scratch parser has to hold for the rule to find it. The include is the anchor
#: that pairs the allow-list with the table it reads; the two paths resolve inside whichever
#: root the rule is given, which is what lets these tests use a scratch tree.
INCLUDE = 'const COMPONENTS_CSV: &str = include_str!("../../../data/components/components.csv");'


def manifest_tool() -> ModuleType:
    """`tools/manifest.py`, imported by name because `tools/` is not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("manifest")
    finally:
        sys.path.pop(0)


def _manifest() -> object:
    found, problems = manifest_tool().read()
    assert not problems, problems
    return found


def scratch(root: Path, body: str) -> None:
    """A tree holding only the parser the rule reads."""
    source = root / RUST
    source.parent.mkdir(parents=True, exist_ok=True)
    (root / "data" / "components").mkdir(parents=True, exist_ok=True)
    source.write_text(body, encoding="utf-8")


def parser(*names: str) -> str:
    """A `parse_components` whose allow-list is exactly `names`."""
    listed = "".join(f'    "{name}",\n' for name in names)
    return f"{INCLUDE}\nfn parse_components() {{\n  for name in [\n{listed}  ] {{}}\n}}\n"


def test_the_repository_marks_every_indexed_column_as_used() -> None:
    """The live check: nothing the parser indexes is still calling itself unread."""
    problems = manifest_tool().indexed_column_problems(_manifest(), REPO_ROOT)
    assert not problems, "\n".join(problems)


def test_an_indexed_column_marked_vendored_is_reported(tmp_path: Path) -> None:
    """Sabotage, against the real manifest and a parser that reads one unread column.

    The root is a scratch tree rather than the repository, so this is a test of the rule
    and not of the tree it happens to be sitting in.
    """
    # `triplepointdensity`, and it has to be a column the manifest still calls `vendored`:
    # **`parachor` was this test's example until the rate-based packed column's segment model
    # started reading it**, at which point the sabotage stopped firing - which is the rule
    # working, and a stale test rather than a stale manifest.
    scratch(tmp_path, parser("triplepointdensity"))

    problems = manifest_tool().indexed_column_problems(_manifest(), tmp_path)
    assert len(problems) == 1, problems
    assert UNREAD in problems[0]
    assert "vendored" in problems[0]

    # And the same tree with a name no column claims is left alone, which is what stops the
    # rule being one that fires on everything.
    scratch(tmp_path, parser("nosuchcolumn"))
    assert manifest_tool().indexed_column_problems(_manifest(), tmp_path) == []


def test_a_parser_the_rule_cannot_read_is_reported_rather_than_skipped(tmp_path: Path) -> None:
    """A rule that cannot find its subject must say so, because silence looks like agreement.

    Three ways to lose it: the file is gone, the change that reads the table is renamed, or
    the allow-list is written some other way. Each has to be a message, not an empty list.
    """
    missing = manifest_tool().indexed_column_problems(_manifest(), tmp_path)
    assert len(missing) == 1 and str(RUST) in missing[0], missing

    scratch(tmp_path, "fn parse_components() {}\n")
    unnamed = manifest_tool().indexed_column_problems(_manifest(), tmp_path)
    assert len(unnamed) == 1 and "COMPONENTS_CSV" in unnamed[0], unnamed

    scratch(tmp_path, f"{INCLUDE}\nfn parse_components() {{}}\n")
    unlisted = manifest_tool().indexed_column_problems(_manifest(), tmp_path)
    assert len(unlisted) == 1 and "parse_components" in unlisted[0], unlisted


def test_the_rule_reads_its_root_and_not_the_working_directory(tmp_path: Path) -> None:
    """A scratch root missing the parser is one problem, never a pass.

    The first draft resolved the include against the process's working directory, so a
    caller passing any other root got a path outside its own tree and the rule reported a
    resolution failure instead of the missing file - right answer, wrong reason, and the
    wrong reason is what a later change would build on.
    """
    problems = manifest_tool().indexed_column_problems(_manifest(), tmp_path)
    assert len(problems) == 1
    assert str(RUST) in problems[0]
    assert "outside the tree" not in problems[0]
