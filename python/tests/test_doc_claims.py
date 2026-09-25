"""The doc-claim gate: the declaration's own rules, and the ways it can rot.

`tools/check_doc_claims.py` runs in CI over a tree that is currently clean, which is the
shape of a check that can stop working without anyone noticing. So the rules are tested
against a synthetic tree rather than only against the repository.

**The rule that matters most is that the declaration holds no measured value.** It is what
stops a failing run being silenced by editing `docs/claims.toml` instead of the page, and
it is the one rule whose violation would look like a fix. So it is tested, and so is every
way a declaration can go stale: a template that stops matching, a probe that measures
nothing, an escape that exempts nothing, a skip whose sentence was rewritten.

The last test runs the real tool over the real tree. A check that is only tested against
fakes is a check whose own subject is unverified.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def gate() -> ModuleType:
    """`tools/check_doc_claims.py`, imported by name. The tools are not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("check_doc_claims")
    finally:
        sys.path.pop(0)


def rooted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *files: str) -> Any:
    """The tool pointed at a synthetic root holding `files`, each `path\\ntext`."""
    tool: Any = gate()
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    for entry in files:
        name, _, body = entry.partition("\n")
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body + "\n", encoding="utf-8")
    return tool


def specs(calculations: int, models: int) -> list[str]:
    """The spec trees a count probe reads, and a page to point a claim at."""
    return [f"specs/calcs/eos/c{i}.toml\nid = 'eos.c{i}'" for i in range(calculations)] + [
        f"specs/models/eos/m{i}.toml\nid = 'eos.m{i}'" for i in range(models)
    ]


# --- the declaration's own rules ---------------------------------------------


def test_a_claim_that_agrees_with_the_tree_passes(tmp_path, monkeypatch):
    tool = rooted(
        tmp_path,
        monkeypatch,
        *specs(3, 2),
        "docs/page.md\nthe library has 3 calculations and 2 models.",
    )
    claim = tool.Claim(
        {
            "page": "docs/page.md",
            "text": "the library has {calcs} calculations and {models} models.",
            "measure": {"calcs": "specs.calcs", "models": "specs.models"},
        },
        0,
    )
    assert claim.check() == []


def test_a_claim_the_tree_contradicts_fails_naming_both(tmp_path, monkeypatch):
    tool = rooted(
        tmp_path,
        monkeypatch,
        *specs(3, 2),
        "docs/page.md\nthe library has 5 calculations and 2 models.",
    )
    claim = tool.Claim(
        {
            "page": "docs/page.md",
            "text": "the library has {calcs} calculations and {models} models.",
            "measure": {"calcs": "specs.calcs", "models": "specs.models"},
        },
        0,
    )
    failures = claim.check()
    assert len(failures) == 1
    assert "calcs = 5" in failures[0] and "the tree is 3" in failures[0]
    # The message a reader acts on names the command that reproduces the measurement.
    assert "--probe specs.calcs" in failures[0]


def test_a_declaration_that_carries_a_value_is_refused(tmp_path, monkeypatch):
    """The rule that stops the check being tuned until it passes."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *specs(3, 2),
        "docs/page.md\nthe library has 3 calculations.",
        """docs/claims.toml
[[claim]]
page = "docs/page.md"
text = "the library has {3} calculations."
measure = { 3 = "specs.calcs" }
""",
    )
    with pytest.raises(tool.ProbeError, match="carries a value"):
        tool.load_claims()


def test_a_template_that_stops_matching_is_a_broken_declaration(tmp_path, monkeypatch):
    tool = rooted(
        tmp_path,
        monkeypatch,
        *specs(3, 2),
        "docs/page.md\nthe library now says three calculations.",
    )
    claim = tool.Claim(
        {
            "page": "docs/page.md",
            "text": "the library has {calcs} calculations",
            "measure": {"calcs": "specs.calcs"},
        },
        0,
    )
    failures = claim.check()
    assert len(failures) == 1
    assert "matches its page 0 time(s)" in failures[0]


def test_a_probe_that_measures_nothing_is_a_failure_not_a_pass(tmp_path, monkeypatch):
    """A check that did not run is not a check that has passed."""
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nthe library has 0 calculations.")
    claim = tool.Claim(
        {
            "page": "docs/page.md",
            "text": "the library has {calcs} calculations",
            "measure": {"calcs": "specs.calcs"},
        },
        0,
    )
    with pytest.raises(tool.ProbeError):
        tool.specs_calcs()
    # And an empty measurement is not silently a match for a page that says 0.
    failures = claim.check()
    assert len(failures) == 1
    assert "could not measure" in failures[0]


def test_an_unknown_probe_is_named_with_the_ones_that_exist(tmp_path, monkeypatch):
    tool = rooted(
        tmp_path,
        monkeypatch,
        *specs(3, 2),
        "docs/page.md\nthe library has 3 calculations.",
    )
    claim = tool.Claim(
        {
            "page": "docs/page.md",
            "text": "the library has {calcs} calculations",
            "measure": {"calcs": "specs.nonsense"},
        },
        0,
    )
    failures = claim.check()
    assert "not a known probe" in failures[0]


# --- reading a captured hole --------------------------------------------------


def test_a_hole_reads_digits_separators_and_number_words():
    tool: Any = gate()
    assert tool.as_int("184") == 184
    assert tool.as_int("1,299,006") == 1_299_006
    assert tool.as_int("**44**") == 44
    assert tool.as_int("Twenty-seven") == 27
    assert tool.as_int("six") == 6
    with pytest.raises(ValueError):
        tool.as_int("about forty")


def test_a_hole_captures_one_token_so_it_cannot_reach_back_into_prose():
    """A hole allowed spaces would let the leftmost match start in the sentence before."""
    tool: Any = gate()
    pattern = tool.template_to_regex("{declared} unit operations declared on typed channels")
    page = "- **Unit operations** — the palette: 29 unit operations declared on typed channels"
    assert pattern.search(page).group("declared") == "29"


# --- the path sweep -----------------------------------------------------------


def test_the_sweep_names_an_inline_code_path_that_does_not_exist(tmp_path, monkeypatch):
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nsee `crates/gone/src/lib.rs` for it")
    checked, failures, escapes = tool.sweep_paths([tmp_path / "docs" / "page.md"])
    assert checked == 1 and escapes == []
    assert len(failures) == 1
    assert "crates/gone/src/lib.rs" in failures[0]


def test_an_escape_exempts_only_the_path_that_is_absent(tmp_path, monkeypatch):
    """The usual case: one line names a path because it is gone and another that is there."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        "docs/page.md",
        "crates/real/src/lib.rs",
    )
    page = tmp_path / "docs" / "page.md"
    page.write_text(
        "`crates/gone/` does not exist; `crates/real/src/lib.rs` does "
        "<!-- doc-claims-ok: the sentence says the first is absent -->\n",
        encoding="utf-8",
    )
    _, failures, escapes = tool.sweep_paths([page])
    assert failures == []
    assert escapes == ["docs/page.md:1 'crates/gone/'"]


def test_an_escape_that_exempts_nothing_is_stale(tmp_path, monkeypatch):
    tool = rooted(tmp_path, monkeypatch, "docs/page.md", "crates/real/src/lib.rs")
    page = tmp_path / "docs" / "page.md"
    page.write_text(
        "`crates/real/src/lib.rs` <!-- doc-claims-ok: leftover -->\n", encoding="utf-8"
    )
    _, failures, _ = tool.sweep_paths([page])
    assert len(failures) == 1
    assert "exempts nothing" in failures[0]


def test_an_escape_without_a_reason_is_refused(tmp_path, monkeypatch):
    tool = rooted(tmp_path, monkeypatch, "docs/page.md")
    page = tmp_path / "docs" / "page.md"
    page.write_text("`crates/gone/` <!-- doc-claims-ok: -->\n", encoding="utf-8")
    _, failures, _ = tool.sweep_paths([page])
    assert len(failures) == 1
    assert "no reason" in failures[0]


def test_the_sweep_leaves_upstream_names_and_globs_alone(tmp_path, monkeypatch):
    """The prefix rule is what keeps NeqSim's vocabulary out without an allowlist."""
    tool = rooted(tmp_path, monkeypatch, "docs/page.md")
    page = tmp_path / "docs" / "page.md"
    page.write_text(
        "`PhaseGEUniquac` and `thermo/util/leachman/` and `src/main` are NeqSim's;\n"
        "`specs/calcs/<namespace>/<id>.toml` is a template.\n",
        encoding="utf-8",
    )
    checked, failures, _ = tool.sweep_paths([page])
    assert (checked, failures) == (0, [])


# --- skips --------------------------------------------------------------------


def test_a_skip_goes_stale_when_its_sentence_leaves_the_page(tmp_path, monkeypatch):
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nthe page was reworded entirely.")
    skips = [
        {
            "page": "docs/page.md",
            "text": "NeqSim master is 1,299,006 lines",
            "reason": "no checkout in CI",
        }
    ]
    failures = tool.check_skips(skips, [])
    assert len(failures) == 1
    assert "stale" in failures[0]


# --- and the tree this all runs against ---------------------------------------


def test_the_shipped_tree_passes(monkeypatch):
    """The real declaration against the real tree - the check's own subject."""
    tool: Any = gate()
    monkeypatch.setattr(tool, "ROOT", REPO_ROOT)
    monkeypatch.setattr(sys, "path", sys.path + [str(REPO_ROOT / "tools")])
    assert tool.main([]) == 0
