"""Every document naming a NeqSim revision names the one vendored under `sources/`.

A spec's `references` says which upstream revision its model was ported from, and the
manifest says which revision is in `databank/sources/`. Those are the same revision, so a
document naming a different one cites a file this project does not hold.

**They had come apart, and nothing could see it.** The refresh that moved the manifest to
`805cf0f` left seventeen specs and `NOTICE` citing `dedba873`. Neither hash was wrong on
its own terms, so no gate failed: the difference was between two documents rather than
between a document and a file. Checking afterwards, the two revisions differ in one of the
twenty-eight cited sources - so the stale citation named a revision whose source says
something else.

What is tested here is mostly the *scanner*, because which commits in a document are
NeqSim's is the part with a judgement in it: `NOTICE` cites two upstreams, and the
`lean-units` commit is not this manifest's business.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

VENDORED = "f0c7436c6923766b1e22957b7075f650600457a7"
STALE = "dedba8735d030c6411e09b6fd7e69f6c4a136114"


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


def test_a_document_with_two_upstreams_yields_only_neqsims_commit() -> None:
    """`NOTICE` gives each upstream a section, and only one of them is this manifest's."""
    text = f"""
azoth
Copyright

--------
lean-units
--------
  Project:   lean-units - https://github.com/ecyrbe/lean-units
  Version:   commit 5322fc77dd4e3ec8d31af2c757a47abbeb94d500

--------
NeqSim
--------
  Project:   NeqSim - https://github.com/equinor/neqsim
  Version:   master (commit {VENDORED})
"""
    assert manifest_tool().neqsim_citations(text) == [VENDORED]


def test_a_one_letter_line_is_not_a_heading() -> None:
    """`databank/README.md` draws its stages with a lone `v` on a line of its own.

    Read as a heading, that `v` would close the page's one passage and silently stop the
    commit two lines further down being checked at all - a check that quietly stops
    covering a document is worse than one that never covered it.
    """
    text = f"""
databank/sources/    upstream files, whole, at a named revision    EXISTS
        |  compile             tools/gen_databank.py
        v
data/components/     the files both languages read                 EXISTS

`manifest.toml` records the NeqSim commit it was last checked against
(`{VENDORED}`, master), and nothing here can tell you a newer NeqSim exists.
"""
    assert manifest_tool().neqsim_citations(text) == [VENDORED]


def test_the_repository_names_the_vendored_revision_everywhere() -> None:
    """The live check, over the real specs, `NOTICE` and the databank page."""
    problems = manifest_tool().citation_problems(_manifest(), REPO_ROOT)
    assert not problems, "\n".join(problems)


def test_a_spec_citing_another_revision_is_reported(tmp_path: Path) -> None:
    """Sabotage: the check fires on the one thing it exists for.

    The root is a scratch tree rather than the repository, so this is a test of the rule
    and not of the tree it happens to be sitting in.
    """
    spec = tmp_path / "specs" / "models" / "eos"
    spec.mkdir(parents=True)
    (spec / "water_phase.toml").write_text(
        f'references = ["NeqSim, `Iapws_if97.java`, Apache-2.0, commit {STALE}."]\n',
        encoding="utf-8",
    )

    problems = manifest_tool().citation_problems(_manifest(), tmp_path)
    assert len(problems) == 1, problems
    assert STALE in problems[0]
    assert str(Path("specs/models/eos/water_phase.toml")) in problems[0]

    # And the same spec against the vendored revision is left alone, which is what stops
    # the check being one that fires on everything.
    (spec / "water_phase.toml").write_text(
        f'references = ["NeqSim, `Iapws_if97.java`, Apache-2.0, commit {VENDORED}."]\n',
        encoding="utf-8",
    )
    assert manifest_tool().citation_problems(_manifest(), tmp_path) == []


@pytest.mark.parametrize("document", ["NOTICE", "databank/README.md"])
def test_the_two_pages_that_state_the_revision_are_scanned(document: str) -> None:
    """A page whose passage the scanner never enters would pass this check vacuously."""
    text = (REPO_ROOT / document).read_text(encoding="utf-8")
    assert manifest_tool().neqsim_citations(text) == [VENDORED]
