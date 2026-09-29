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
    """The synthetic tree's roots, and the files a test rewrites."""

    specs: Path
    models: Path
    rust_root: Path
    python_root: Path
    spec: Path
    rust: Path
    python: Path

    def check(self) -> list[str]:
        # The tool is loaded by path, so its functions are `Any` to mypy. `specs_root` and
        # `source_roots` are passed explicitly: a default would sweep this checkout instead of
        # the tree under test, which is the one way these cases could pass without checking
        # anything.
        return cast(
            "list[str]",
            check_unported.check(
                self.models,
                self.rust_root,
                self.python_root,
                self.specs,
                (self.rust_root, self.python_root),
            ),
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
        specs=tmp_path / "specs",
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


def test_two_models_may_refuse_the_same_value(tree: Tree) -> None:
    """**A key belongs to the spec that declares it, and two may declare one.**

    `process.absorption_column` and `process.distillation_column` both decline the eight
    `ColumnSolverFactory` strategies neither carries, so a map keyed by key alone called
    that a duplicate the first time this checker ran over the converted tree.
    """
    (tree.models / "gadget.toml").write_text(SPEC.replace("widget", "gadget"), encoding="utf-8")
    (tree.rust_root / "azoth-process" / "src" / "gadget.rs").write_text(RUST_SITE, encoding="utf-8")
    (tree.python_root / "azoth" / "process" / "reference" / "gadget.py").write_text(
        PYTHON_SITE, encoding="utf-8"
    )
    assert tree.check() == []


# --- the closure rule: a porting claim may not stand in prose ---------------------
#
# 43 strings over 30 specs carried the claim when the rule was written, in four senses - a
# declaration, a correction, an argument, and a count of a set declared elsewhere - and 38
# of them were about a class that is *absent*, which no row can hold, because a row is only
# available where the port refuses a value of a declared input. Those went to `ROADMAP.md`,
# which is the document whose charter is naming the class that would close a gap.


def test_a_porting_claim_in_a_spec_sentence_is_refused(tree: Tree) -> None:
    """The rule. Every one of the 43 was this shape."""
    tree.spec.write_text(
        SPEC.replace('name = "Widget"', 'name = "a widget, not ported"'), encoding="utf-8"
    )
    failures = tree.check()
    assert any("not ported" in f and "in prose" in f for f in failures), failures


def test_the_field_name_itself_is_not_a_claim(tree: Tree) -> None:
    """**`[[unported]]` is the declaration's own name, and the rule may not fire on it.**

    A rule that read keys rather than values would forbid the field it exists to protect.
    """
    assert "unported" in SPEC
    assert tree.check() == []


def test_a_claim_in_the_palette_tree_is_refused(tree: Tree) -> None:
    """The sweep is over every spec tree, not just the models one - `unit_ops` is where
    four of the six stale notes this tranche found had been sitting unread."""
    entry = tree.specs / "unit_ops" / "utility" / "widget.toml"
    entry.parent.mkdir(parents=True)
    entry.write_text(
        'id = "unit_ops.widget"\nname = "Widget"\ndescription = "is not carried"\n',
        encoding="utf-8",
    )
    failures = tree.check()
    # The message names the substring the pattern matched, not the whole phrase: `not
    # carried`, not `is not carried`.
    assert any("unit_ops" in f and "not carried" in f for f in failures), failures


def test_the_rule_reaches_a_nested_string(tree: Tree) -> None:
    """A string inside an array of tables is a string in the document like any other, and
    every one of the 43 findings was exactly that - an `assumptions` entry or a nested
    `source` field."""
    tree.spec.write_text(
        SPEC.replace('name = "Widget"', 'name = "Widget"\nassumptions = ["not ported"]'),
        encoding="utf-8",
    )
    failures = tree.check()
    assert any("assumptions[0]" in f and "not ported" in f for f in failures), failures


@pytest.mark.parametrize(
    "claim",
    ["None is ported", "none ported", "Neither is ported", "neither are ported"],
)
def test_the_alternation_is_the_whole_claim(tree: Tree, claim: str) -> None:
    """**The first draft of this rule read `not ported|unported|is not carried` and missed
    19 strings the tree was already carrying.**

    Sixteen `process.*` entries said "NeqSim's class also carries a `runTransient` and a
    `MechanicalDesign`. **None is ported**", and three more said "neither is ported" or
    "none ported". A phrase list that stops at one wording is one a spec dodges by
    rewording - which is worse than no list, because it reads as coverage.
    """
    tree.spec.write_text(
        SPEC.replace('name = "Widget"', f'name = "Widget"\nnote = "{claim}"'), encoding="utf-8"
    )
    failures = tree.check()
    assert any(claim.lower() in f.lower() for f in failures), failures


# --- the same rule, over code: the sentinel is confined --------------------------


def test_a_refusal_spelled_as_a_message_is_refused(tree: Tree) -> None:
    """**The rule the sentinel is named for.** A refusal goes *through* the helper, because that
    is what a declared row can hold; a sentence written beside it is held to nothing."""
    tree.rust.write_text(RUST_SITE + '\n\nfn f() { let _ = "this is not ported here"; }\n')
    failures = tree.check()
    assert any("in a message" in f and "not ported" in f for f in failures), failures


def test_a_refusal_as_a_message_is_refused_in_python_too(tree: Tree) -> None:
    tree.python.write_text(PYTHON_SITE + '\nraise ValueError("this is not carried either")\n')
    failures = tree.check()
    assert any("in a message" in f for f in failures), failures


def test_the_phrase_in_a_rust_doc_comment_is_documentation(tree: Tree) -> None:
    """**The carve-out that makes the rule usable.** `///` and `//` are where this phrase
    legitimately lives - it is how the tree explains a refusal - and neither is a string."""
    tree.rust.write_text(
        RUST_SITE + "\n/// A wetting floor the class is not ported with.\n// not ported either.\n"
    )
    assert tree.check() == []


def test_the_phrase_in_a_python_docstring_is_documentation(tree: Tree) -> None:
    """The Python half of the same carve-out. This one needs an **AST** rather than lines: a
    docstring is a string literal, and the tree has twenty-odd of them explaining refusals."""
    tree.python.write_text(
        PYTHON_SITE + '\n\ndef f() -> None:\n    """A rule that is not ported here."""\n'
    )
    assert tree.check() == []


def test_the_field_name_is_not_a_claim_in_code_either(tree: Tree) -> None:
    """**`UNPORTED_SOLVERS` and `'unported': [...]` are a name and a key.**
    The boundaries are what separate the word from the fragment of an identifier."""
    tree.python.write_text(PYTHON_SITE + '\n__all__ = ["UNPORTED_SOLVERS"]\n')
    assert tree.check() == []


def test_a_generated_file_is_not_swept(tree: Tree) -> None:
    """A `*_gen.py` is emitted from the specs, so its prose is the specs' prose and the `unported`
    key it carries is the field's own name. Checking it would ask a generator to stop emitting
    the name it must emit."""
    generated = tree.python_root / "azoth" / "_models_gen.py"
    generated.write_text('MODELS = [{"unported": [{"parameter": "mode"}]}]\n', encoding="utf-8")
    assert tree.check() == []


def test_the_rule_sees_a_message_continued_across_lines(tree: Tree) -> None:
    """**The correction this rule needed, and why its first draft was worse than no rule.**

    Every refusal in this tree is a `format!` whose text is backslash-continued onto the next
    line, so the phrase sits on a line holding no `"` at all. The per-line scan this replaced
    found **none** of the three Rust sites the rule exists to confine, while still passing a
    sabotage probe written on one line - and it missed a seventh site outright
    (`crates/azoth-eos/src/databank.rs`), where the phrase is wrapped between its own two words.
    """
    tree.rust.write_text(RUST_SITE + '\nfn f() { let _ = "this is not \\\n  ported here"; }\n')
    failures = tree.check()
    assert any("in a message" in f and "not ported" in f for f in failures), failures


def test_a_trailing_comment_after_code_is_still_a_comment(tree: Tree) -> None:
    """The case a quote count was only a proxy for: a `//` on a line that also holds a literal is
    a comment, not a message. The lexer decides that now rather than the parity of the quotes."""
    tree.rust.write_text(RUST_SITE + '\nfn f() { let x = "ok"; } // this is not ported\n')
    assert tree.check() == []


def test_the_phrase_in_a_python_comment_is_documentation(tree: Tree) -> None:
    """Python's `#` needs no carve-out of its own - an AST never sees a comment - which is the
    whole reason that side reads the tree instead of the text."""
    tree.python.write_text(PYTHON_SITE + "\n# this is not ported either\n")
    assert tree.check() == []


def test_the_two_sentinels_are_the_only_files_allowed_to_say_it() -> None:
    """**The exemption, proved against the real tree, because nothing synthetic can prove it.**

    `SENTINELS` holds repo-relative paths, so no fixture can name one: from inside a synthetic
    tree, the exemption and a sweep that never opens those files look identical. So the assertion
    is the pair - the sentinels *do* spell the sentence they are sentinel for, and the real sweep
    is clean - which is what makes the exemption load-bearing rather than vacuous.
    """
    for sentinel in check_unported.SENTINELS:
        text = (ROOT / sentinel).read_text(encoding="utf-8")
        assert "not ported" in text, f"{sentinel} no longer builds the sentence it stands for"
    assert check_unported.source_prose() == []
