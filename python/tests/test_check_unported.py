"""`tools/check_unported.py` fails on the disagreement it exists to find.

The tool is green against the tree, and green is the state that proves nothing on its
own: a check whose patterns have stopped matching agrees with a tree that has stopped
having rows. So every rule is exercised here against a synthetic tree that breaks it -
the same reason `spec_lint.py` takes a `--spec-dir`.
"""

from __future__ import annotations

import importlib.util
import tomllib
from pathlib import Path
from typing import NamedTuple, cast

import pytest

ROOT = Path(__file__).resolve().parents[2]

_spec = importlib.util.spec_from_file_location(
    "check_unported", ROOT / "tools" / "check_unported.py"
)
assert _spec and _spec.loader
check_unported = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_unported)

SPEC = """
id = "process.widget"
name = "Widget"
kind = "procedure"
inputs = { mode = { type = "enum", values = ["fast", "slow"] } }
outputs = {}
implementations = {}
[[unported]]
parameter = "mode"
value = "slow"
class = "Widget.Mode.SLOW"
"""

RUST_SITE = 'Err(unported::refuse("mode=slow"))'
PYTHON_SITE = '_unported.refuse("mode=slow")'


class Tree(NamedTuple):
    """The synthetic tree's three roots, and the three files a test rewrites."""

    models: Path
    rust_root: Path
    python_root: Path
    spec: Path
    rust: Path
    python: Path

    def check(self) -> list[str]:
        # The tool is loaded by path, so its functions are `Any` to mypy.
        return cast(
            "list[str]", check_unported.check(self.models, self.rust_root, self.python_root)
        )


@pytest.fixture
def tree(tmp_path: Path) -> Tree:
    """A widget spec declaring one row, and one file in each language."""
    spec = tmp_path / "specs" / "models" / "widget.toml"
    spec.parent.mkdir(parents=True)
    spec.write_text(SPEC, encoding="utf-8")
    rust = tmp_path / "crates" / "azoth-process" / "src" / "widget.rs"
    rust.parent.mkdir(parents=True)
    rust.write_text(RUST_SITE, encoding="utf-8")
    python = tmp_path / "python" / "src" / "azoth" / "process" / "reference" / "widget.py"
    python.parent.mkdir(parents=True)
    python.write_text(PYTHON_SITE, encoding="utf-8")
    return Tree(
        models=spec.parent,
        rust_root=tmp_path / "crates",
        python_root=tmp_path / "python" / "src",
        spec=spec,
        rust=rust,
        python=python,
    )


def test_the_declaration_parses_to_the_key_both_languages_spell() -> None:
    """The key is the whole interface between the spec and a call site."""
    row = tomllib.loads(SPEC)["unported"][0]
    assert check_unported.unported_key(row) == "mode=slow"


def test_agreeing_tree_is_clean(tree: Tree) -> None:
    assert tree.check() == []


def test_a_language_that_refuses_what_the_spec_carries_fails(tree: Tree) -> None:
    """**The defect this tool was built for.** Python refuses a value neither the spec
    nor Rust does, which is exactly how `billet_schultes_1999` sat in the tree."""
    tree.python.write_text(PYTHON_SITE + '\n_unported.refuse("mode=fast")', encoding="utf-8")
    failures = tree.check()
    assert any("python" in f and "mode=fast" in f for f in failures), failures


def test_a_declared_row_no_implementation_refuses_fails(tree: Tree) -> None:
    """A row is a claim the code has kept - so a row with no site is a claim it has not."""
    tree.python.write_text("pass\n", encoding="utf-8")
    failures = tree.check()
    assert any("python" in f and "mode=slow" in f for f in failures), failures


def test_a_row_naming_an_undeclared_input_fails(tree: Tree) -> None:
    tree.spec.write_text(
        SPEC + '[[unported]]\nparameter = "speed"\nvalue = "slow"\nclass = "Widget.Speed"\n',
        encoding="utf-8",
    )
    failures = tree.check()
    assert any("not an input" in f for f in failures), failures


def test_a_row_naming_an_absent_capture_fails(tree: Tree) -> None:
    tree.spec.write_text(
        SPEC.replace(
            'class = "Widget.Mode.SLOW"',
            'class = "Widget.Mode.SLOW"\ncapture = "validation/neqsim/captures/nope.tsv"',
        ),
        encoding="utf-8",
    )
    failures = tree.check()
    assert any("absent" in f for f in failures), failures


def test_a_refusal_with_no_spec_of_its_name_fails(tree: Tree) -> None:
    """A refusal in a file no spec answers to is a fact nothing holds."""
    gadget = tree.rust_root / "azoth-process" / "src" / "gadget.rs"
    gadget.write_text(RUST_SITE, encoding="utf-8")
    failures = tree.check()
    assert any("nothing to hold it" in f for f in failures), failures
