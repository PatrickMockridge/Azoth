"""`tools/check_claims.py` fails on every shape it exists to find.

The tool is green against the tree, and green proves nothing on its own: a pattern that
has stopped matching agrees with a tree that has stopped carrying claims. So each rule is
exercised here against a synthetic tree that breaks it.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import NamedTuple, cast

import pytest

ROOT = Path(__file__).resolve().parents[2]

_spec = importlib.util.spec_from_file_location("check_claims", ROOT / "tools" / "check_claims.py")
assert _spec and _spec.loader
check_claims = importlib.util.module_from_spec(_spec)
# Registered before execution: `@dataclass` resolves its module through `sys.modules`,
# and a module loaded by path alone is absent from it.
sys.modules[_spec.name] = check_claims
_spec.loader.exec_module(check_claims)

RUST_TEST = "fn the_test() {}"
PYTHON_TEST = "def test_the_widget() -> None:\n    assert True\n"


class Tree(NamedTuple):
    """The synthetic checkout's roots, and the two source files a test rewrites."""

    root: Path
    rust: Path
    python: Path

    def check(self) -> list[str]:
        # The tool is loaded by path, so its functions are `Any` to mypy. The roots are
        # passed explicitly: a default would sweep this checkout rather than the tree
        # under test, which is the one way these cases could pass without checking anything.
        return cast("list[str]", check_claims.check(self.root, self.root, self.root))


@pytest.fixture
def tree(tmp_path: Path) -> Tree:
    """A tree with a Rust and a Python source file, and a test in each language."""
    rust_dir = tmp_path / "crates" / "widget" / "src"
    rust_dir.mkdir(parents=True)
    rust = rust_dir / "widget.rs"
    rust.write_text("pub fn widget() -> u32 { 1 }\n", encoding="utf-8")
    rust_test = tmp_path / "crates" / "widget" / "tests" / "widget.rs"
    rust_test.parent.mkdir(parents=True)
    rust_test.write_text(RUST_TEST, encoding="utf-8")
    python_dir = tmp_path / "python" / "src" / "azoth"
    python_dir.mkdir(parents=True)
    python = python_dir / "widget.py"
    python.write_text("def widget() -> int:\n    return 1\n", encoding="utf-8")
    python_test = tmp_path / "python" / "tests" / "test_widget.py"
    python_test.parent.mkdir(parents=True)
    python_test.write_text(PYTHON_TEST, encoding="utf-8")
    return Tree(tmp_path, rust, python)


def test_an_unmarked_panic_is_a_finding(tree: Tree) -> None:
    tree.rust.write_text('fn widget() { panic!("never here"); }\n', encoding="utf-8")
    messages = tree.check()
    assert len(messages) == 1
    assert "carries no `claim-ok" in messages[0]


def test_a_marker_naming_an_existing_test_is_clean(tree: Tree) -> None:
    tree.rust.write_text(
        'fn widget() {\n    // claim-ok: the_test\n    panic!("never here");\n}\n',
        encoding="utf-8",
    )
    assert tree.check() == []


def test_a_marker_naming_no_test_is_a_finding(tree: Tree) -> None:
    tree.rust.write_text(
        'fn widget() {\n    // claim-ok: no_such_test\n    panic!("never here");\n}\n',
        encoding="utf-8",
    )
    messages = tree.check()
    assert len(messages) == 1
    assert "names no test" in messages[0]


def test_a_marker_covers_its_block_to_the_previous_blank_line(tree: Tree) -> None:
    tree.rust.write_text(
        'fn widget() {\n    // claim-ok: the_test\n    let why = 1;\n    panic!("{why}");\n}\n',
        encoding="utf-8",
    )
    assert tree.check() == []


def test_a_blank_line_ends_the_markers_block(tree: Tree) -> None:
    tree.rust.write_text(
        'fn widget() {\n    // claim-ok: the_test\n\n    panic!("never here");\n}\n',
        encoding="utf-8",
    )
    messages = tree.check()
    # Both directions speak: the claim has no marker in its block, and the marker above
    # the blank line has no claim.
    assert len(messages) == 2
    assert any("carries no `claim-ok" in message for message in messages)
    assert any("names no claim" in message for message in messages)


def test_the_cfg_test_module_is_not_swept(tree: Tree) -> None:
    tree.rust.write_text(
        "fn widget() {}\n\n#[cfg(test)]\nmod tests {\n    #[test]\n    fn it_works() {\n"
        "        assert!(true);\n    }\n}\n",
        encoding="utf-8",
    )
    assert tree.check() == []


def test_a_generated_file_is_not_swept(tree: Tree) -> None:
    tree.rust.write_text(
        '// GENERATED FILE - DO NOT EDIT BY HAND\nfn widget() { panic!("generated"); }\n',
        encoding="utf-8",
    )
    assert tree.check() == []


def test_a_comment_quoting_a_panic_is_not_a_claim(tree: Tree) -> None:
    tree.rust.write_text(
        '// panic!("this is prose about the rule")\nfn widget() {}\n', encoding="utf-8"
    )
    assert tree.check() == []


def test_an_expect_needs_a_guarantee_phrase(tree: Tree) -> None:
    tree.rust.write_text(
        'fn widget(held: Option<u32>) -> u32 {\n    held.expect("the value")\n}\n',
        encoding="utf-8",
    )
    assert tree.check() == []
    tree.rust.write_text(
        'fn widget(held: Option<u32>) -> u32 {\n    held.expect("checked above")\n}\n',
        encoding="utf-8",
    )
    assert len(tree.check()) == 1


def test_a_python_assert_is_a_claim_and_a_typed_refusal_is_not(tree: Tree) -> None:
    tree.python.write_text(
        "def widget(held):\n    assert held is not None\n    return held\n",
        encoding="utf-8",
    )
    messages = tree.check()
    assert len(messages) == 1
    assert "widget.py:2" in messages[0]
    tree.python.write_text(
        "def widget(held):\n    if held is None:\n"
        '        raise InvalidInputError("widget", "held is required")\n'
        "    return held\n",
        encoding="utf-8",
    )
    assert tree.check() == []


def test_a_raise_assertion_error_is_a_claim(tree: Tree) -> None:
    tree.python.write_text(
        'def widget(held):\n    if held is None:\n        raise AssertionError("never here")\n',
        encoding="utf-8",
    )
    assert len(tree.check()) == 1


def test_a_marker_with_no_claim_is_a_finding(tree: Tree) -> None:
    tree.rust.write_text("fn widget() {}\n// claim-ok: the_test\n", encoding="utf-8")
    messages = tree.check()
    assert len(messages) == 1
    assert "names no claim in its block" in messages[0]


def test_the_test_support_crate_is_not_swept(tree: Tree) -> None:
    support = tree.root / "crates" / "azoth-test-support" / "src" / "lib.rs"
    support.parent.mkdir(parents=True)
    support.write_text('pub fn verdict() { panic!("the test failed"); }\n', encoding="utf-8")
    assert tree.check() == []
