#!/usr/bin/env python3
"""Execute the worked notebooks and check what they commit.

    python tools/check_notebooks.py                 # execute every notebooks/*.ipynb
    python tools/check_notebooks.py notebooks/01-simple-calculations.ipynb
    python tools/check_notebooks.py --write         # refresh committed outputs
    python tools/check_notebooks.py --check-outputs # compare text outputs, locally

The committed notebook is a rendered example; this gate proves it still runs. Execution
happens in a temporary copy, so the repository is not touched unless `--write` is given.

Two halves. The structural half is stdlib only: it holds a committed notebook to the shape
GitHub renders - every code cell executed, no error output, at least one figure. The
execution half needs the `notebook` extra and starts a Jupyter kernel; a missing extra is a
skip locally and an error under `AZOTH_REQUIRE_NOTEBOOKS=1`, which is what CI sets.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Final

ROOT: Final = Path(__file__).resolve().parents[1]
NOTEBOOKS: Final = ROOT / "notebooks"

#: A committed notebook carries its figures. Past this it no longer renders quickly on GitHub.
MAX_BYTES: Final = 2 * 1024 * 1024

#: A worked example that needs more cells than this is a chapter, not a notebook.
MAX_CELLS: Final = 30

#: The reference is a result a notebook prints; the kernel needs a bounded wait for a cell.
TIMEOUT: Final = 300


class Problem(Exception):
    """One notebook did not meet its contract."""


def read_notebook(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Problem(f"cannot read the notebook: {error}") from None
    return document


def structural_problems(document: dict[str, Any], path: Path) -> list[str]:
    """What a rendered notebook must be, independent of Jupyter.

    Deliberately not a schema validator: it enforces the three properties a reader relies on
    - there is a code cell, every cell ran, and a figure came out - and refuses the states
    that make a notebook look executed when it is not.
    """
    problems: list[str] = []
    if document.get("nbformat") != 4:
        problems.append(f"nbformat is {document.get('nbformat')!r}, not 4")
    cells = document.get("cells")
    if not isinstance(cells, list):
        return [*problems, "`cells` is not a list"]
    code = [cell for cell in cells if cell.get("cell_type") == "code"]
    if not code:
        problems.append("no code cell")
    if len(cells) > MAX_CELLS:
        problems.append(f"{len(cells)} cells exceeds the {MAX_CELLS}-cell limit")
    figures = 0
    for index, cell in enumerate(code, start=1):
        if cell.get("execution_count") is None:
            problems.append(f"code cell {index} has no execution count (not executed)")
        for output in cell.get("outputs", []):
            if output.get("output_type") == "error":
                problems.append(f"code cell {index} carries an error output")
            data = output.get("data", {})
            if "image/png" in data:
                figures += 1
    if code and figures == 0:
        problems.append("no figure output")
    size = path.stat().st_size if path.exists() else 0
    if size > MAX_BYTES:
        problems.append(f"{size / 1024 / 1024:.1f} MiB exceeds the 2 MiB limit")
    return problems


def execute(path: Path, timeout: int) -> tuple[dict[str, Any], float]:
    """Run the notebook in a temporary copy and return the executed document and seconds."""
    try:
        import nbformat
        from nbclient import NotebookClient
    except ImportError as error:  # pragma: no cover - the CI job installs the extra
        raise Problem(f"the notebook extra is not installed: {error}") from None

    kernel_env = {
        # The kernelspec's argv says `python`; make sure it is this interpreter's venv.
        "PATH": os.pathsep.join([str(Path(sys.executable).parent), os.environ.get("PATH", "")]),
    }
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="azoth-notebook-") as work:
        os.environ.setdefault("MPLCONFIGDIR", str(Path(work) / "mpl"))
        notebook = nbformat.read(path, as_version=4)
        client = NotebookClient(
            notebook,
            timeout=timeout,
            startup_timeout=60,
            kernel_name="python3",
            resources={"metadata": {"path": work}},
            allow_errors=False,
            record_timing=True,
            env=kernel_env,
        )
        client.execute()
    return notebook, time.monotonic() - started


def text_outputs(document: dict[str, Any]) -> list[str]:
    """The `text/plain` answer of every code cell, for a local output comparison."""
    answers: list[str] = []
    for cell in document.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        parts = []
        for output in cell.get("outputs", []):
            data = output.get("data", {})
            if "text/plain" in data:
                parts.append("".join(data["text/plain"]))
            elif output.get("output_type") == "stream":
                parts.append("".join(output.get("text", "")))
        answers.append("\n".join(parts).strip())
    return answers


def report(path: Path, document: dict[str, Any], seconds: float) -> str:
    cells = [cell for cell in document.get("cells", []) if cell.get("cell_type") == "code"]
    figures = sum(
        1
        for cell in cells
        for output in cell.get("outputs", [])
        if "image/png" in output.get("data", {})
    )
    kib = path.stat().st_size / 1024
    return f"{len(cells):>3} cells {seconds:>5.1f} s {figures:>2} figure(s) {kib:>6.0f} KiB"


def run_one(path: Path, *, write: bool, check_outputs: bool, timeout: int) -> bool:
    name = path.relative_to(ROOT)
    try:
        committed = None if write else read_notebook(path)
        if committed is not None:
            problems = structural_problems(committed, path)
            if problems:
                raise Problem("; ".join(problems) + " (run --write to refresh)")
        executed, seconds = execute(path, timeout)
        problems = structural_problems(executed, path)
        if problems:
            raise Problem("; ".join(problems))
        if check_outputs and committed is not None:
            before, after = text_outputs(committed), text_outputs(executed)
            changed = [
                i for i, (a, b) in enumerate(zip(before, after, strict=True), start=1) if a != b
            ]
            if changed:
                raise Problem(f"text output differs in code cell(s) {changed}")
        if write:
            import nbformat

            nbformat.write(executed, path)
        print(f"  PASS  {name}  {report(path, executed, seconds)}")
        return True
    except Problem as error:
        print(f"  FAIL  {name}  {error}")
        return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "paths", nargs="*", type=Path, help="notebooks to run (default: notebooks/*.ipynb)"
    )
    parser.add_argument("--write", action="store_true", help="write the executed notebook back")
    parser.add_argument(
        "--check-outputs", action="store_true", help="compare text outputs with the committed file"
    )
    parser.add_argument(
        "--timeout", type=int, default=TIMEOUT, help=f"seconds per cell (default {TIMEOUT})"
    )
    arguments = parser.parse_args(argv)

    if arguments.write and os.environ.get("CI"):
        print("check_notebooks: --write is refused under CI; commit the refreshed outputs instead")
        return 2

    paths = arguments.paths or sorted(NOTEBOOKS.glob("*.ipynb"))
    paths = [path if path.is_absolute() else (Path.cwd() / path) for path in paths]
    if not paths:
        print("check_notebooks: no notebooks found")
        return 1

    try:
        import nbclient  # noqa: F401
    except ImportError:
        message = "check_notebooks: SKIP - install the notebook extra: uv sync --extra notebook"
        if os.environ.get("AZOTH_REQUIRE_NOTEBOOKS"):
            print(message.replace("SKIP", "FAIL"))
            return 2
        print(message)
        return 0

    print(f"check_notebooks: {len(paths)} notebook(s)")
    results = [
        run_one(
            path,
            write=arguments.write,
            check_outputs=arguments.check_outputs,
            timeout=arguments.timeout,
        )
        for path in paths
    ]
    if all(results):
        print(f"check_notebooks: PASS ({len(results)} of {len(results)})")
        return 0
    print(f"check_notebooks: FAILED ({results.count(False)} of {len(results)})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
