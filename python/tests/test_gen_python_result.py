"""The result dataclasses are generated, and this is the ratchet that says so.

`tools/gen_python_result.py` emits `python/src/azoth/core/result_gen.py` - one frozen dataclass per
registered id, its fields in `CalcResult::FIELDS` order, and `python/src/azoth/core/result.py`
re-exports them beside the hand-written enums, mixin and record in `core/result_base.py`.
`docs-drift` regenerates and diffs it; the `--check` below runs the same gate in the suite.

`test_registry_contract.py` already holds each class's *fields* to `FIELDS`, and it does so through
the imported classes, so it now ratchets the generated file without having been changed. What is
left for this file is what that one cannot see: that the committed file is what the generator
emits, that the registry and the file account for each other, and that the prose in it is the
spec's rather than a paraphrase.
"""

from __future__ import annotations

import ast
import importlib
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULT_GEN = REPO_ROOT / "python" / "src" / "azoth" / "core" / "result_gen.py"


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _classes() -> dict[str, ast.ClassDef]:
    tree = ast.parse(RESULT_GEN.read_text(encoding="utf-8"))
    return {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}


def test_the_committed_results_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_result.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the result dataclasses are out of date - run `python tools/gen_python_result.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_file_and_the_registry_account_for_each_other() -> None:
    """One class per registered id, and no id without one.

    `test_registry_contract` compares *fields*, so a class that lost its `CALC_ID` or an id that
    never got a class is outside what it looks at.
    """
    rust_index = _tools_module("rust_index")
    registered = {calc_id for calc_id, _, _ in rust_index.implementations()}
    declared = {}
    for name, node in _classes().items():
        declared[name] = _calc_id(node)

    assert len(declared) == len(_classes()), (
        f"{len(_classes())} class(es) and {len(declared)} CALC_ID(s) - a result without one is not "
        f"reachable from its spec: {sorted(set(_classes()) - set(declared))}"
    )
    assert set(declared.values()) == registered, (
        "the generated results and the registry do not account for each other: "
        f"{sorted(registered - set(declared.values()))} have no class, "
        f"{sorted(set(declared.values()) - registered)} have no id"
    )


def test_every_field_doc_is_the_specs_own() -> None:
    """**A generated file may not paraphrase.** Each `#:` line is the spec's `description`.

    The Rust transport makes the same rule for its struct fields, and states why: a field's doc
    comment written here would be a second copy of prose that already has one home. The one field
    no spec declares is `warnings`, whose line is the module's.
    """
    generator = _tools_module("gen_python_result")
    all_specs = generator.specs()
    source = RESULT_GEN.read_text(encoding="utf-8")

    checked, wrong = 0, []
    for name, node in _classes().items():
        outputs = all_specs[_calc_id(node)].get("outputs") or {}
        for field, written, expected in _documented_fields(node, source, outputs):
            checked += 1
            if written != expected:
                wrong.append(f"{name}.{field}: wrote {written[:60]!r}")

    assert checked, "no field doc was compared - the parser found none"
    assert not wrong, (
        f"{wrong} carry a doc that is not the spec's `description`, so the prose has two homes"
    )


def _calc_id(node: ast.ClassDef) -> str:
    """One generated class's `CALC_ID`, which is also its key in the spec table."""
    for statement in node.body:
        if (
            isinstance(statement, ast.AnnAssign)
            and isinstance(statement.target, ast.Name)
            and statement.target.id == "CALC_ID"
            and isinstance(statement.value, ast.Constant)
        ):
            return str(statement.value.value)
    raise AssertionError(f"{node.name} declares no CALC_ID")


def _documented_fields(
    node: ast.ClassDef, source: str, outputs: Mapping[str, Any]
) -> list[tuple[str, str, str]]:
    """`(field, the doc written above it, the doc the spec declares)` for one class."""
    lines = source.split("\n")
    found = []
    for statement in node.body:
        if not isinstance(statement, ast.AnnAssign) or not isinstance(statement.target, ast.Name):
            continue
        field = statement.target.id
        if ast.unparse(statement.annotation) == "ClassVar[str]":
            continue
        # The `#:` block immediately above the field. `lineno` is 1-based and absolute, so the
        # field's own line is `lineno - 1` and the block is read upwards from the line before it.
        index = statement.lineno - 1
        block: list[str] = []
        while index > 0 and lines[index - 1].strip().startswith("#:"):
            block.insert(0, lines[index - 1].strip().removeprefix("#:").strip())
            index -= 1
        if field == "warnings":
            expected = "Caveats."
        else:
            text = str((outputs.get(field) or {}).get("description", ""))
            expected = "\n".join(line.strip() for line in text.splitlines())
        found.append((field, "\n".join(block), expected))
    return found
