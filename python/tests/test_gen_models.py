"""The model registry generator, on a model the tree does not have.

`gen_registry` and `gen_docs` each have a test that runs `--check` and asserts the tree
is what they emit. `gen_models` had none - the drift job regenerates and diffs it, so
a spec edited without regenerating is caught in CI, but nothing ran it in the suite and
nothing drove its emitters with an input the tree does not contain.

That is the gap this closes. A generator that is only ever run over one shape is a
generator whose other branches are untested, and a model spec is the interface a
contributor writes against: the failure mode is a spec that validates and then
generates something that does not compile.

The emitters are pure functions of a list of model dicts, so a synthetic model can be
pushed through them without touching the tree.
"""

from __future__ import annotations

import ast
import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def gen_models() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_models")
    finally:
        sys.path.pop(0)


def real_models() -> list[dict[str, Any]]:
    from azoth._models_gen import MODELS

    return list(MODELS)


def synthetic(namespace: str, name: str, **overrides: Any) -> dict[str, Any]:
    """A model spec shaped like a real one, with a namespace and name of its own."""
    model: dict[str, Any] = {
        # `load_models` attaches the spec's path, and the emitter writes it into the
        # generated header - so a model built by hand has to carry one too.
        "_path": f"specs/models/{namespace}/{name}.yaml",
        "id": f"{namespace}.{name}",
        "name": "A synthetic model",
        "kind": "direct",
        "inputs": {"T": {"type": "quantity", "unit": "K"}},
        "outputs": {"T": {"unit": "K", "description": "the temperature back"}},
        "implementations": {
            "python": f"azoth.{namespace}.{name}",
            "rust": f"azoth_{namespace}::{name}",
        },
        "cases": [
            {
                "id": "the_only_case",
                "inputs": {"T": 300.0},
                "expected": {"T": 300.0},
                "tolerance": 1.0e-9,
            }
        ],
    }
    model.update(overrides)
    return model


def test_the_python_module_it_emits_parses() -> None:
    """A registry that does not parse would be an ImportError at the first caller."""
    tool: Any = gen_models()
    source: str = tool.emit_python(real_models())

    ast.parse(source)
    assert "MODELS" in source


def test_a_synthetic_model_reaches_the_python_registry() -> None:
    tool: Any = gen_models()
    models = [*real_models(), synthetic("synthetic", "one")]

    source: str = tool.emit_python(models)

    ast.parse(source)
    assert "synthetic.one" in source


def test_a_synthetic_model_reaches_the_rust_registry() -> None:
    tool: Any = gen_models()
    models = [*real_models(), synthetic("synthetic", "one")]

    source: str = tool.emit_rust(models, "synthetic")

    assert "synthetic.one" in source


def test_a_namespace_with_no_models_emits_nothing_to_register() -> None:
    """The generator is called once per namespace, so an empty one must not fail."""
    tool: Any = gen_models()

    source: str = tool.emit_rust(real_models(), "no_such_namespace")

    assert "no_such_namespace" not in source


@pytest.mark.parametrize("kind", ["direct", "procedure"])
def test_each_kind_generates_a_registry_entry(kind: str) -> None:
    """Both kinds must survive codegen; the schema allows both and they differ in
    whether there is an algorithm block at all."""
    tool: Any = gen_models()
    overrides: dict[str, Any] = {"kind": kind}
    if kind == "procedure":
        overrides["algorithm"] = {
            "scheme": "bisection",
            "convergence": "absolute",
            "tolerance": 1.0e-9,
            "max_iterations": 50,
        }
    models = [*real_models(), synthetic("synthetic", kind, **overrides)]

    source: str = tool.emit_python(models)

    ast.parse(source)
    assert f"synthetic.{kind}" in source


def test_the_committed_registries_are_current() -> None:
    """The gate the drift job runs, run here so it fails in the suite too.

    Through `--check` rather than by comparing the emit: the generator passes the Rust
    it emits through `rustfmt` before writing, so a comparison against the raw emitter
    output fails on formatting rather than on drift.
    """
    import subprocess

    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_models.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the generated model registries are out of date - run "
        f"`python tools/gen_models.py`\n{result.stdout}{result.stderr}"
    )
