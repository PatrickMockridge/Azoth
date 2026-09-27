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


def test_a_claim_that_agrees_with_the_tree_passes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_claim_the_tree_contradicts_fails_naming_both(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_declaration_that_carries_a_value_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_template_that_stops_matching_is_a_broken_declaration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_probe_that_measures_nothing_is_a_failure_not_a_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_an_unknown_probe_is_named_with_the_ones_that_exist(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_hole_reads_digits_separators_and_number_words() -> None:
    tool: Any = gate()
    assert tool.as_int("184") == 184
    assert tool.as_int("1,299,006") == 1_299_006
    assert tool.as_int("**44**") == 44
    assert tool.as_int("Twenty-seven") == 27
    assert tool.as_int("six") == 6
    with pytest.raises(ValueError):
        tool.as_int("about forty")


def test_a_hole_captures_one_token_so_it_cannot_reach_back_into_prose() -> None:
    """A hole allowed spaces would let the leftmost match start in the sentence before."""
    tool: Any = gate()
    pattern = tool.template_to_regex("{declared} unit operations declared on typed channels")
    page = "- **Unit operations** — the palette: 29 unit operations declared on typed channels"
    assert pattern.search(page).group("declared") == "29"


# --- the path sweep -----------------------------------------------------------


def test_the_sweep_names_an_inline_code_path_that_does_not_exist(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nsee `crates/gone/src/lib.rs` for it")
    checked, failures, escapes = tool.sweep_paths([tmp_path / "docs" / "page.md"])
    assert checked == 1 and escapes == []
    assert len(failures) == 1
    assert "crates/gone/src/lib.rs" in failures[0]


def test_an_escape_exempts_only_the_path_that_is_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_an_escape_that_exempts_nothing_is_stale(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = rooted(tmp_path, monkeypatch, "docs/page.md", "crates/real/src/lib.rs")
    page = tmp_path / "docs" / "page.md"
    page.write_text("`crates/real/src/lib.rs` <!-- doc-claims-ok: leftover -->\n", encoding="utf-8")
    _, failures, _ = tool.sweep_paths([page])
    assert len(failures) == 1
    assert "exempts nothing" in failures[0]


def test_an_escape_without_a_reason_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = rooted(tmp_path, monkeypatch, "docs/page.md")
    page = tmp_path / "docs" / "page.md"
    page.write_text("`crates/gone/` <!-- doc-claims-ok: -->\n", encoding="utf-8")
    _, failures, _ = tool.sweep_paths([page])
    assert len(failures) == 1
    assert "no reason" in failures[0]


def test_the_sweep_leaves_upstream_names_and_globs_alone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
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


def test_a_skip_goes_stale_when_its_sentence_leaves_the_page(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nthe page was reworded entirely.")
    skips = [
        {
            "page": "docs/page.md",
            "text": "NeqSim master is 1,299,006 lines",
            "reason": "no checkout in CI",
        }
    ]
    failures = tool.check_skips(skips)
    assert len(failures) == 1
    assert "stale" in failures[0]


# --- the enforcement sweep ----------------------------------------------------
#
# Two trees: the nine layers in `docs/src/calculus/`, and the surface a front-end binds in
# `docs/src/architecture/`. Each page states once where its claim is enforced. The rule worth
# testing is the one that makes the sweep a gate rather than a formality: **a page may claim
# `nothing` only where no status claims a proof.** A proved claim whose implementation half
# enforces nothing is the finding the sweep exists for.


def calculus(*pages: str) -> list[str]:
    """Synthetic calculus pages, each `name\\ntext`."""
    return [f"docs/src/calculus/{entry}" for entry in pages]


def architecture(*pages: str) -> list[str]:
    """Synthetic architecture pages. The same rules, the other normative tree."""
    return [f"docs/src/architecture/{entry}" for entry in pages]


def test_an_architecture_page_with_no_marker_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The surface a front-end binds makes the same kind of claim, so it is swept the same way.

    Without this the sweep reads one tree and reports OK while the pages that describe the
    editor and the middleware state their enforcement nowhere.
    """
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus("process.md\n*Status: **proved**.*\n\n*Enforcement: check — `x.rs`.*"),
        *architecture("middleware.md\nIt holds no copy of the flowsheet."),
        "x.rs\n",
    )
    checked, failures, _ = tool.sweep_enforcement()
    assert checked == 1
    assert len(failures) == 1
    assert "docs/src/architecture/middleware.md" in failures[0]


def test_an_enforcement_naming_a_status_the_page_does_not_state_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A marker about a status this page never states is one written for another page's claim.

    The shape a contradiction takes: `process.md` says the balance is **specified** while three
    other pages group it with the recycle as **characterised**, and no rule could see it because
    each page's own statuses are consistent. This catches the half that is decidable - a marker
    whose status word is nowhere on the page it sits under.
    """
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "process.md\n*Status: **specified**.*\n\n"
            "*Enforcement: check — `x.rs`, and the claim stays **characterised**.*"
        ),
        "x.rs\n",
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "'characterised'" in failures[0]
    assert "this page does not make" in failures[0]


def test_a_marker_naming_a_status_the_page_does_state_passes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The rule above is about a stray word, not about bold text in a marker."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "process.md\n*Status: **specified**.*\n\n*Status: **characterised**.*\n\n"
            "*Enforcement: check — `x.rs`, and the second claim stays **characterised**.*"
        ),
        "x.rs\n",
    )
    _, failures, _ = tool.sweep_enforcement()
    assert failures == []


def test_a_proved_claim_with_nothing_enforcing_it_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The whole point: proved in Lean, enforced nowhere, and the page says so - a failure."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "dimensions.md\n*Status: **proved**.*\n\n*Enforcement: nothing — it is a proof.*"
        ),
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "**proved**" in failures[0]


def test_a_specified_claim_may_say_nothing_enforces_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """And the same page passes once the status says the layer is not there yet."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "barbs.md\n*Status: **specified**.*\n\n*Enforcement: nothing — the layer is specified.*"
        ),
    )
    checked, failures, unenforced = tool.sweep_enforcement()
    assert (checked, failures) == (1, [])
    assert len(unenforced) == 1


def test_a_page_with_no_status_may_not_claim_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Nothing then says the absence is deliberate rather than unexamined."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus("numerics.md\n*Enforcement: nothing — nobody looked.*"),
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "deliberate" in failures[0]


def test_an_enforcement_naming_a_path_that_does_not_resolve_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`construction` and `check` are claims about the tree, so they name where."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "process.md\n*Status: **proved**.*\n\n"
            "*Enforcement: check — `crates/nowhere/src/check.rs`.*"
        ),
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "none of them exists" in failures[0]


def test_an_enforcement_naming_no_path_at_all_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Prose where a path belongs is a claim a reader cannot follow."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "rho.md\n*Status: **proved**.*\n\n*Enforcement: construction — everywhere, obviously.*"
        ),
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "names no path" in failures[0]


def test_a_page_with_no_marker_or_two_of_them_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = rooted(tmp_path, monkeypatch, *calculus("vocabulary.md\n*Status: **proved**.*"))
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "0 `*Enforcement:*` marker(s)" in failures[0]

    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "vocabulary.md\n*Status: **proved**.*\n\n"
            "*Enforcement: check — `tools/check_numerics.py`.*\n\n"
            "*Enforcement: check — `tools/check_doc_claims.py`.*"
        ),
        "tools/check_numerics.py\n",
        "tools/check_doc_claims.py\n",
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "2 `*Enforcement:*` marker(s)" in failures[0]


def test_the_marker_goes_after_the_last_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """It summarises the page's claims, so one above a status is a marker about the wrong page."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        *calculus(
            "implicit.md\n*Enforcement: check — `tools/check_numerics.py`.*\n\n"
            "*Status: **specified**.*"
        ),
        "tools/check_numerics.py\n",
    )
    _, failures, _ = tool.sweep_enforcement()
    assert len(failures) == 1
    assert "above the last" in failures[0]


def test_a_tree_with_no_calculus_pages_is_not_a_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The sweep would agree with itself at zero if the directory moved."""
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nnothing here")
    with pytest.raises(tool.ProbeError):
        tool.sweep_enforcement()


# --- and the tree this all runs against ---------------------------------------


def test_the_shipped_tree_passes(monkeypatch: pytest.MonkeyPatch) -> None:
    """The real declaration against the real tree - the check's own subject."""
    tool: Any = gate()
    monkeypatch.setattr(tool, "ROOT", REPO_ROOT)
    monkeypatch.setattr(sys, "path", [*sys.path, str(REPO_ROOT / "tools")])
    assert tool.main([]) == 0


# --- the port's sources ------------------------------------------------------


#: A spec of the shape the source sweep reads: an id and a `[source]` block.
def a_spec(name: str, source: str) -> str:
    return f"specs/models/eos/{name}.toml\nid = 'eos.{name}'\n[source]\n{source}\n"


def test_a_spec_naming_its_neqsim_class_passes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A `.java` path is the shape the probe can hold a port to."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec("thing", 'standard = "NeqSim process/equipment/Thing.java"'),
    )
    checked, failures, markers = tool.sweep_sources()
    assert (checked, failures, markers) == (1, [], [])


def test_a_spec_naming_a_paper_or_a_standard_passes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """And so are the three other kinds of place a port's arithmetic comes from."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec("paper", 'standard = "Chung, T.-H.; Ajlan, M. (1988)"'),
        a_spec("standard", 'standard = "Crane TP-410"'),
        a_spec("textbook", 'standard = "Standard thermodynamics, as in Smith, Van Ness & Abbott"'),
    )
    checked, failures, markers = tool.sweep_sources()
    assert checked == 3, checked
    assert failures == [], failures
    assert markers == []


def test_a_bare_neqsim_is_not_a_source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """**The failure this sweep exists for.** A library is not a class, and nothing about
    "NeqSim" says what was ported or what to compare against."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec("vague", 'standard = "NeqSim, developed at NTNU and maintained by Equinor"'),
    )
    checked, failures, _ = tool.sweep_sources()
    assert checked == 1
    assert len(failures) == 1, failures
    assert "names a class, method, paper, standard or author" in failures[0]


def test_a_spec_with_no_source_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A spec nobody can check is refused rather than passed over."""
    # One spec with a source, so the probe has something to measure - a tree where *nothing*
    # names a source is a broken probe rather than a finding.
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec("named", 'standard = "NeqSim process/Thing.java"'),
        "specs/calcs/eos/bare.toml\nid = 'eos.bare'\n",
    )
    _, failures, _ = tool.sweep_sources()
    assert len(failures) == 1, failures
    assert "has no `[source]`" in failures[0]


def test_the_marker_excuses_a_source_that_cannot_be_named(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A definition has no paper, and the marker says so rather than the sweep assuming it."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec(
            "defined",
            'standard = "The definition of molar mass"\n'
            "# source-ok: molar mass is a definition, and there is no paper to cite for one.",
        ),
    )
    checked, failures, markers = tool.sweep_sources()
    assert (checked, failures) == (1, [])
    assert len(markers) == 1 and "molar mass" in markers[0]


def test_an_empty_source_block_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`edition` may carry the citation, but `standard` may not be missing."""
    tool = rooted(
        tmp_path,
        monkeypatch,
        a_spec("empty", 'edition = "`Thing.java` at line 4"'),
    )
    _, failures, _ = tool.sweep_sources()
    assert len(failures) == 1, failures
    assert "is empty" in failures[0]


def test_no_spec_at_all_is_a_broken_probe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A tree that moved would otherwise report zero findings and pass."""
    tool = rooted(tmp_path, monkeypatch, "docs/page.md\nnothing here")
    with pytest.raises(tool.ProbeError):
        tool.sweep_sources()
