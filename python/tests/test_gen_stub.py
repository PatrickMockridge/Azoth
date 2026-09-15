"""The type stub's renderers, on shapes the tree does not contain.

`gen_stub --check` regenerates the file and diffs it, which catches a spec edited
without regenerating. What it cannot catch is a renderer that is wrong for a shape no
current spec uses - and the stub is what a type checker reads, so a wrong rendering is
a wrong signature everywhere the checker looks, with nothing to compare it against.

The renderers are pure functions of a declaration or an annotation, so they can be
driven directly. Each case below is a shape the *emitters* have to handle rather than
one the tree happens to exercise: an optional input, a vector of strings, a union with
None, an empty result class.

`gen_stub.py` imports `azoth._dispatch` at module scope, so the package must be
importable when this runs - which it is, from `pythonpath` in `pyproject.toml`.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
STUB = REPO_ROOT / "python" / "src" / "azoth" / "_core.pyi"


def gen_stub() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_stub")
    finally:
        sys.path.pop(0)


@pytest.mark.parametrize(
    ("declaration", "expected"),
    [
        ({"unit": "K"}, "float"),
        ({"type": "quantity", "unit": "Pa"}, "float"),
        ({"type": "vector", "length": "one per component"}, "list[float]"),
        ({"type": "matrix", "shape": "N x N"}, "list[list[float]]"),
        ({"type": "fitting_list"}, "list[str]"),
    ],
)
def test_a_declared_input_renders_as_the_python_type_it_arrives_as(
    declaration: dict[str, Any], expected: str
) -> None:
    tool: Any = gen_stub()
    assert tool.parameter_type(declaration) == expected


def test_an_optional_input_gets_a_default_of_none() -> None:
    """A caller must be able to leave it out, which is what the spec declared."""
    tool: Any = gen_stub()
    rendered = tool.render_parameter("mu", {"type": "quantity", "optional": True})

    assert rendered == "mu: float | None = None"


def test_a_required_input_has_no_default() -> None:
    tool: Any = gen_stub()
    rendered = tool.render_parameter("T", {"type": "quantity", "unit": "K"})

    assert rendered == "T: float"
    assert "=" not in rendered


@pytest.mark.parametrize(
    ("annotation", "expected"),
    [
        (float, "float"),
        (int, "int"),
        (str, "str"),
        (bool, "bool"),
        (list[float], "list[float]"),
        (tuple[float, ...], "list[float]"),
        (None, "None"),
        (float | None, "float | None"),
        (list[float] | None, "list[float] | None"),
    ],
)
def test_an_annotation_renders_as_a_valid_stub_type(annotation: Any, expected: str) -> None:
    """Including the two shapes the tree's results do not have: a bare None and a
    union, both of which appear in the transport classes."""
    tool: Any = gen_stub()
    assert tool.stub_type(annotation) == expected


def test_a_long_signature_is_wrapped_rather_than_left_on_one_line() -> None:
    """A stub line over the line limit is unreadable and fails the formatter."""
    tool: Any = gen_stub()

    entry: dict[str, Any] = {
        "id": "eos.pt_flash",
        "inputs": {
            f"a_rather_long_input_name_{index}": {"type": "quantity", "unit": "K"}
            for index in range(8)
        },
    }
    rendered: str = tool.render_signature(entry)

    assert "\n" in rendered, rendered
    for line in rendered.splitlines():
        assert len(line) <= 100, f"{len(line)} characters: {line!r}"


def test_the_stub_on_disk_is_every_result_class_the_registry_declares() -> None:
    """The rendered stub and the registries name the same result types.

    Cheap, and it is the half a `--check` cannot state: the gate proves the file
    matches what the generator emits *now*, not that what it emits is complete.
    """
    tool: Any = gen_stub()
    rendered: str = tool.render()

    from azoth._dispatch import result_types
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    for entry in [*CALCS, *MODELS]:
        name = result_types()[entry["id"]].__name__
        assert f"class {name}:" in rendered, f"{entry['id']} renders no {name}"


def test_the_committed_stub_is_current() -> None:
    """The gate CI runs, run here too so it fails in the suite rather than on a push."""
    tool: Any = gen_stub()
    assert STUB.read_text(encoding="utf-8") == tool.render(), (
        "python/src/azoth/_core.pyi is out of date; run tools/gen_stub.py"
    )
