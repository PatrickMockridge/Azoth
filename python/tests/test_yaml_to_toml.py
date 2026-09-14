"""The YAML-to-TOML converter, held to the two parsers it stands between.

`tools/yaml_to_toml.py` rewrites every data file in this repository, and the thing that
makes a rewrite of that size trustworthy is not reading the converter - it is running
both parsers over the same document and comparing the trees. So that is what this file
does, over **every YAML file in the repository**: the original through PyYAML, the
conversion through `tomllib`, and the two trees asserted equal. One implementation of
"read this document", two parsers, real data - the shape this repository uses wherever a
contract has two implementations.

What the corpus cannot reach is pinned by name below. A spec file's prose exercises the
common spellings and nothing else, so the awkward values - a triple quote inside a
paragraph, a value whose second line is indented, a path ending in a backslash - are
written out here with the rule each one holds in place. Each was found rather than
imagined: a randomised sweep over the writer's own alphabet turned up a lone `"` reaching
a one-quote-wide delimiter, and a value whose *first* line is indented losing that
indentation to a continuation backslash. Both are cases below, and both are also shown to
be load-bearing by disabling the guard that stops them and watching the same value come
back wrong.

One file cannot be converted while this is written: `specs/vocabulary/vocabulary.yaml`
holds six `uom: null`s and **TOML has no null**. That is a policy decision for the stage
that converts the specs, so the converter refuses the file and this test asserts the
refusal - a file skipped without the refusal being checked is a file whose failure mode
nobody has seen.
"""

from __future__ import annotations

import datetime
import importlib
import os
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Directories the corpus walk does not enter. `.github` is the interesting one: it holds
#: YAML, but it is GitHub's format carrying GitHub's data - workflows and issue templates
#: - rather than this repository's, and converting it would be an edit to how CI is
#: declared for no gain. The rest is build output and the vendored Lean toolchain, which
#: `git ls-files '*.yaml'` does not see either.
_SKIP_DIRS = frozenset(
    {
        ".git",
        ".github",
        ".mypy_cache",
        ".pytest_cache",
        ".venv",
        "__pycache__",
        "build",
        "dist",
        "lean",
        "node_modules",
        "target",
        "venv",
    }
)

#: Files the converter cannot write, with the reason it must give. The refusal is
#: asserted, not assumed: see `test_a_table_that_cannot_be_written_is_refused`.
REFUSED: dict[str, str] = {
    "specs/vocabulary/vocabulary.yaml": "TOML has no null",
}

#: Values that exercise one spelling decision of the writer each, named for the decision
#: rather than for the value so that a failure says which rule broke.
VALUES: dict[str, str] = {
    "prose": "The friction factor is solved by iteration.\n",
    "two_paragraphs": "First paragraph.\n\nSecond paragraph.\n",
    "trailing_whitespace": "a line ending in a space \nand the line after it\n",
    "latex": "\\kappa = 0.37464 + 1.54226\\omega - 0.26992\\omega^{2}",
    "non_ascii": "Theorie analytique - 30 degrees C, theta",
    "markdown": "**bold**, `code`, [a link](https://example.invalid)\n\nsecond\n",
    "double_space_inside": "a run  of two spaces, kept as two",
    "windows_path": "C:\\srv\\spec\\pr_kappa.yaml",
    "ends_in_a_backslash": "a Windows root, C:\\",
    "ends_in_a_quote": "the component is named 'methane'",
    "triple_quote_inside": 'the note reads """ and stops',
    "triple_single_quote_inside": "it reads ''' three quotes ''' and carries on",
    "quote_and_triple_single_quote": "the note reads \" and then ''' and stops",
    "indented_after_a_triple_quote": (
        "prose holding ''' which rules out the verbatim spelling\n    indented by four spaces\n"
    ),
    "leading_space_on_the_first_line": " leading space\nsecond line\n",
    "tab_indented": "step one\n\tstep two\n",
}

#: A document rather than a value: the shapes a table has, including the one the format
#: change is *for*. `factor` is written unquoted and read back as a float - under YAML
#: 1.1 the same literal is a **string**, which is why this repository carries a lint rule
#: that TOML makes unnecessary.
DOCUMENT: dict[str, Any] = {
    "schema_version": 2,
    "published": True,
    "factor": 1.0e12,
    "small": 1.0e-12,
    "retrieved": datetime.date(2026, 9, 13),
    "coefficients": {"hydraulics.orifice_flow": {"Cd": {"value": 0.61, "unit": "dimensionless"}}},
    "components": ["methane", "n-butane"],
    "rows": [{"a": 1, "b": "x"}, {"a": 2, "b": "y"}],
    "blocks": [{"id": "first"}, {"id": "second"}],
    "long_blocks": [
        {"id": "n-hexadecane", "reason": "Long enough to need a line of its own, " * 4}
    ],
    "empty_table": {},
    "empty_list": [],
}


def converter() -> ModuleType:
    """The tool, importable as a module.

    `tools/` is not a package and not on the path (see `pyproject.toml`), so it is
    inserted for the duration of the import rather than left there.
    """
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("yaml_to_toml")
    finally:
        sys.path.pop(0)


def relative(path: Path) -> str:
    """A path as this file names it: relative to the repository root, with slashes."""
    return path.relative_to(REPO_ROOT).as_posix()


def yaml_files() -> list[Path]:
    """Every YAML file in the repository, in path order.

    Walked rather than globbed, because `lean/.lake` holds several gigabytes of Mathlib
    and a glob would read every one of them to discard the result.
    """
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dirnames[:] = [name for name in dirnames if name not in _SKIP_DIRS]
        found.extend(Path(dirpath) / name for name in filenames if name.endswith(".yaml"))
    return sorted(found)


def revisit(document: Any) -> Any:
    """A document through the converter and back, as the comparison the test makes."""
    text = converter().convert(document)
    return tomllib.loads(text)


def test_the_corpus_is_the_repositorys_data_files() -> None:
    """The walk finds the files it is supposed to, and does not wander into the rest.

    Without this, a typo in the skip list turns the round-trip below into a comparison
    of whatever is left - which is the check-that-cannot-fail this repository refuses
    everywhere else.
    """
    found = {relative(path) for path in yaml_files()}

    assert {
        "databank/manifest.yaml",
        "specs/calcs/eos/pr_departure.yaml",
        "specs/models/eos/pt_flash.yaml",
        "specs/vocabulary/vocabulary.yaml",
    } <= found
    assert not any(name.startswith(("lean/", ".github/")) for name in found), sorted(found)


@pytest.mark.parametrize(
    "path",
    [path for path in yaml_files() if relative(path) not in REFUSED],
    ids=lambda path: relative(path),
)
def test_every_yaml_file_converts_losslessly(path: Path) -> None:
    """The document survives the conversion, value for value.

    Both parsers run over the same document and the trees are compared whole, so a
    reflowed paragraph, an escaped backslash and a table emitted after a sub-table of
    its own all fail the same way: the tree that comes back is not the tree that went
    in.
    """
    original = yaml.safe_load(path.read_text(encoding="utf-8"))
    converted = converter().convert(original)

    assert tomllib.loads(converted) == original, f"{relative(path)} does not survive"


@pytest.mark.parametrize("name", sorted(REFUSED), ids=sorted(REFUSED))
def test_a_table_that_cannot_be_written_is_refused(name: str) -> None:
    """A file TOML cannot express fails the conversion, and says why.

    The alternative - skipping it quietly - is the failure this test exists to prevent:
    a file whose absence from the comparison nobody would notice.
    """
    original = yaml.safe_load((REPO_ROOT / name).read_text(encoding="utf-8"))

    with pytest.raises(ValueError, match=REFUSED[name]):
        converter().convert(original)


def test_every_file_is_either_compared_or_refused() -> None:
    """The corpus is exactly the parametrised list plus the refusals.

    `REFUSED` is written by hand and the corpus is discovered, so this is the join
    between them: a file that the converter quietly cannot write stays a failure until
    somebody writes down why, rather than leaving the parametrised list one file short.
    """
    refused: set[str] = set()

    for path in yaml_files():
        original = yaml.safe_load(path.read_text(encoding="utf-8"))
        try:
            converter().convert(original)
        except ValueError:
            refused.add(relative(path))

    assert refused == set(REFUSED), (
        "the files the converter refuses are not the files this test declares it "
        "refuses - a refusal is a policy decision and belongs in REFUSED with its reason"
    )


@pytest.mark.parametrize("name", sorted(VALUES), ids=sorted(VALUES))
def test_the_awkward_values_survive(name: str) -> None:
    """Every string the writer has a rule for comes back as it went in.

    A value is the whole of what these pin: the writer chooses between a literal string,
    a wrapped basic string and a verbatim block, and each choice is only right if the
    parsed value is unchanged. `tomllib` is the judge, so the spelling is free to be
    ugly where that is what exactness costs.
    """
    value = VALUES[name]

    assert revisit({"k": value}) == {"k": value}


def test_a_scalar_is_spelled_the_way_a_reader_expects() -> None:
    """A short value is an ordinary quoted string unless escaping it would be a loss.

    A literal string is unusual enough that it should carry a reason, and it has one:
    this repository's values include LaTeX and citations with quoted titles, where the
    basic spelling doubles every backslash and escapes every quote. Everywhere else the
    ordinary spelling is what a reader of a TOML file is looking for - and this is forty
    files of hand-reviewed prose, not generated output.
    """
    module = converter()

    assert module.toml_string("hydraulics.reynolds_number") == '"hydraulics.reynolds_number"'
    assert module.toml_string("Reynolds number for pipe flow") == (
        '"Reynolds number for pipe flow"'
    )
    assert module.toml_string("the value is named 'methane'") == (
        "\"the value is named 'methane'\""
    )
    # The two that make the literal spelling worth having.
    assert module.toml_string(r"Re = \frac{\rho v D}{\mu}") == (r"'''Re = \frac{\rho v D}{\mu}'''")
    assert module.toml_string('Reynolds, O. (1883). "An experimental investigation"') == (
        "'''Reynolds, O. (1883). \"An experimental investigation\"'''"
    )


def test_the_document_shapes_survive() -> None:
    """A whole table converts, including the shapes that make ordering matter.

    `coefficients` is a case's shape - a name with dots in it, which is one key and
    would silently become three levels of table unquoted. `blocks` is a list of records,
    which is `[[blocks]]`. `long_blocks` is a record too long to be an inline table, so
    it is the one that has to be a block rather than a line. And `factor` is the number
    the format change is for.
    """
    converted = converter().convert(DOCUMENT)

    assert tomllib.loads(converted) == DOCUMENT
    assert "factor = 1000000000000.0" in converted, (
        "an unsigned exponent must reach TOML as a number, not as a quoted string"
    )
    assert '"hydraulics.orifice_flow"' in converted, (
        "a key containing a dot has to be quoted, or it becomes a path"
    )


@pytest.mark.parametrize(
    "name",
    ["leading_space_on_the_first_line", "indented_after_a_triple_quote"],
    ids=["leading_space", "indented_after_a_triple_quote"],
)
def test_the_layout_escaping_is_what_keeps_the_indentation(
    name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With the layout escaping disabled, the same value comes back wrong.

    The wrapped spelling relies on a trailing backslash, and TOML trims a continuation's
    leading whitespace - so indentation that belongs to the value has to be written as
    escapes. Both values here reach the wrapped spelling, and neither has anything else
    holding its indentation up, so disabling that one function is enough to lose it.
    """
    module = converter()
    value = VALUES[name]
    assert revisit({"k": value}) == {"k": value}

    monkeypatch.setattr(module, "_escaped_layout", lambda line: line)

    assert tomllib.loads(module.convert({"k": value})) != {"k": value}


def test_the_quote_escaping_is_what_stops_a_second_delimiter(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A lone `"` is a delimiter in the one-quote-wide spelling, and must be escaped.

    The multi-line spelling needs only runs of three escaped, and using that escaper
    where the delimiter is one quote wide leaves a stray quote in the middle of a
    document - which is a parse error, not a wrong value, and so is caught here by
    `tomllib` rather than by the comparison.
    """
    module = converter()
    value = "the note reads \" and stops, with ''' to force a quoted spelling"
    assert revisit({"k": value}) == {"k": value}

    monkeypatch.setattr(module, "_escaped_every_quote", module._escaped_quotes)

    with pytest.raises(tomllib.TOMLDecodeError):
        tomllib.loads(module.convert({"k": value}))
