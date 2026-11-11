"""The notebook gate's structural half.

`tools/check_notebooks.py` executes the worked notebooks in CI; the half tested here is the
one that needs no Jupyter and no kernel - the shape a committed notebook has to have for GitHub
to render it. Each refusal below is a state that has looked executed and was not: a cell with
no execution count, an error kept in the file, a notebook of prose with no figure.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

TOOLS = Path(__file__).resolve().parents[2] / "tools"


def tool() -> ModuleType:
    """`tools/check_notebooks.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_notebooks")
    finally:
        sys.path.pop(0)


def written(tmp_path: Path, document: dict[str, Any]) -> Path:
    path = tmp_path / "notebook.ipynb"
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def code_cell(
    *, count: int | None = 1, outputs: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    return {
        "cell_type": "code",
        "execution_count": count,
        "metadata": {},
        "outputs": outputs
        if outputs is not None
        else [
            {
                "output_type": "display_data",
                "data": {"image/png": "iVBORw0KGgo="},
                "metadata": {},
            }
        ],
        "source": "1 + 1",
    }


def markdown_cell() -> dict[str, Any]:
    return {"cell_type": "markdown", "metadata": {}, "source": "# A notebook"}


def document(*cells: dict[str, Any]) -> dict[str, Any]:
    return {"nbformat": 4, "nbformat_minor": 5, "metadata": {}, "cells": list(cells)}


def test_an_executed_notebook_with_a_figure_passes(tmp_path: Path) -> None:
    module = tool()
    path = written(tmp_path, document(markdown_cell(), code_cell()))
    assert module.structural_problems(json.loads(path.read_text()), path) == []


def test_a_cell_with_no_execution_count_is_refused(tmp_path: Path) -> None:
    module = tool()
    path = written(tmp_path, document(code_cell(count=None)))
    problems = module.structural_problems(json.loads(path.read_text()), path)
    assert any("execution count" in problem for problem in problems)


def test_an_error_output_is_refused(tmp_path: Path) -> None:
    module = tool()
    error = {
        "output_type": "error",
        "ename": "SolverNotConvergedError",
        "evalue": "",
        "traceback": [],
    }
    path = written(tmp_path, document(code_cell(outputs=[error])))
    problems = module.structural_problems(json.loads(path.read_text()), path)
    assert any("error output" in problem for problem in problems)


def test_a_notebook_with_no_figure_is_refused(tmp_path: Path) -> None:
    module = tool()
    stream = {"output_type": "stream", "name": "stdout", "text": "42\n"}
    path = written(tmp_path, document(code_cell(outputs=[stream])))
    problems = module.structural_problems(json.loads(path.read_text()), path)
    assert any("no figure" in problem for problem in problems)


def test_a_markdown_only_notebook_is_refused(tmp_path: Path) -> None:
    module = tool()
    path = written(tmp_path, document(markdown_cell()))
    problems = module.structural_problems(json.loads(path.read_text()), path)
    assert any("no code cell" in problem for problem in problems)


def test_an_oversize_notebook_is_refused(tmp_path: Path, monkeypatch: Any) -> None:
    module = tool()
    path = written(tmp_path, document(code_cell()))
    monkeypatch.setattr(module, "MAX_BYTES", 10)
    problems = module.structural_problems(json.loads(path.read_text()), path)
    assert any("2 MiB limit" in problem for problem in problems)


def test_the_text_outputs_are_the_code_cells_only(tmp_path: Path) -> None:
    module = tool()
    stream = {"output_type": "stream", "name": "stdout", "text": "hello\n"}
    document_with_stream = document(markdown_cell(), code_cell(outputs=[stream]))
    assert module.text_outputs(document_with_stream) == ["hello"]
