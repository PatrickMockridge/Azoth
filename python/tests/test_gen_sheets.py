"""The calculation catalogue and dispatcher, held to the generator.

`tools/gen_sheets.py` reads the specs and the kernel signatures; the generated files are what the
UI and the wasm module carry. The check is the same `--check` CI runs, so a spec edited without
regenerating fails here rather than in a browser.
"""

from __future__ import annotations

import importlib
import shutil
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def tool() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_sheets")
    finally:
        sys.path.pop(0)


def test_the_generated_files_are_what_the_specs_emit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The check formats its Rust through `rustfmt`, which the Python job may not carry; the
    # Rust-test job runs the same check with the toolchain.
    if shutil.which("rustfmt") is None:
        pytest.skip("rustfmt is not on this runner; rust-test runs the generated-file check")
    module: Any = tool()
    monkeypatch.setattr(sys, "argv", ["gen_sheets.py", "--check"])
    module.main()


def test_the_catalogue_never_names_the_upstream_implementation() -> None:
    module: Any = tool()
    text = module.render_catalogue()
    assert "neqsim" not in text.lower()
    assert text.count("CalcForm {") == 98


def test_the_shipped_books_match_the_document_schema() -> None:
    import json
    import tomllib

    from jsonschema import Draft202012Validator

    schema = json.loads(
        (REPO_ROOT / "specs" / "schema" / "calcbook.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for path in sorted((REPO_ROOT / "specs" / "calcbooks").glob("*.toml")):
        document = tomllib.loads(path.read_text(encoding="utf-8"))
        errors = list(validator.iter_errors(document))
        assert errors == [], f"{path.name}: {errors}"


def test_every_available_calculation_has_a_dispatch_arm() -> None:
    module: Any = tool()
    text = module.render_dispatch()
    available = [
        spec["id"]
        for _, spec in module.documents()
        if not module.complex_inputs(spec) and not module.complex_outputs(spec)
    ]
    assert len(available) == 94
    for calc_id in available:
        assert f'"{calc_id}" => calc_' in text, calc_id
