"""Running a calculation on a chosen implementation, and replaying its cases.

The public dispatch already knows both implementations and refuses to fall back silently,
so the bundle asks it rather than importing a reference module. The result is written with
the result's own ``to_dict``: that method appends the provenance block, and the block is
the traceability evidence the bundle exists to carry.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import azoth
from azoth import _dispatch
from azoth._dispatch import Backend
from azoth._registry_gen import spec

from . import cases as case_module
from . import compare, inputs
from .paths import repo_root

__all__ = [
    "BackendUnavailable",
    "build_identity",
    "environment",
    "match_case",
    "replay",
    "run",
    "solver_tolerance",
]


class BackendUnavailable(RuntimeError):
    """A backend the bundle needs is not present in this installation."""


def _git(root: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *args],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def solver_tolerance(calc_id: str) -> float | None:
    """The bound a solver diagnostic is compared on, if the spec declares one."""
    solver = spec(calc_id).get("solver")
    if not isinstance(solver, Mapping):
        return None
    tolerance = solver.get("tolerance")
    return float(tolerance) if isinstance(tolerance, (int, float)) else None


def build_identity() -> dict[str, Any]:
    """The stable source identity of this installation.

    The dirty flag is read with ``--untracked-files=no``: an output directory written under
    the checkout is not a change to the source, and including it would make a second run
    hash differently from the first. A wheel has no ``.git`` or ``provenance.json``, and
    those fields are null rather than a claim about a commit nobody can check.
    """
    root = repo_root()
    commit = _git(root, "rev-parse", "HEAD") if root is not None else None
    status = (
        _git(root, "status", "--porcelain", "--untracked-files=no") if root is not None else None
    )
    provenance_sha = None
    if root is not None and (root / "provenance.json").is_file():
        provenance_sha = hashlib.sha256((root / "provenance.json").read_bytes()).hexdigest()

    extension = _dispatch.extension()
    core_version = None
    if extension is not None:
        version = getattr(extension, "version", None)
        if callable(version):
            core_version = str(version())

    return {
        "git_commit": commit,
        "tracked_dirty": None if status is None else bool(status),
        "package_version": azoth.__version__,
        "python": f"{sys.version_info.major}.{sys.version_info.minor}",
        "rust_core_version": core_version,
        "provenance_sha256": provenance_sha,
        "source": "checkout" if root is not None else "wheel",
    }


def environment() -> dict[str, Any]:
    """The unhashed runtime record: this build, on this machine."""
    return {
        "platform": sys.platform,
        "machine": platform.machine(),
        "python": platform.python_version(),
        "extension": _dispatch.extension() is not None,
    }


def run(calc_id: str, kwargs: Mapping[str, Any], backend: Backend) -> dict[str, Any]:
    """One run on one implementation, as the result's own document."""
    if backend not in _dispatch.available():
        raise BackendUnavailable(
            f"the {backend!r} backend is not available; have {sorted(_dispatch.available())}"
        )
    with _dispatch.use_backend(backend):
        function = _dispatch.resolve(calc_id)
        result = function(**kwargs)
    return {"backend": backend, "result": result.to_dict()}


def match_case(
    calc_id: str, kwargs: Mapping[str, Any], case_set: case_module.CaseSet
) -> dict[str, Any] | None:
    """The declared case whose inputs these are, or ``None`` for an exploratory run."""
    encoded = inputs.encode(calc_id, kwargs)
    for case in case_set.cases:
        if inputs.encode(calc_id, case["inputs"]) == encoded:
            return case
    return None


def replay(
    calc_id: str,
    case_set: case_module.CaseSet,
    backends: Sequence[Backend],
    *,
    diagnostic_bound: float | None,
) -> list[dict[str, Any]]:
    """Run every declared case on both implementations and judge it against its expected values."""
    outputs = spec(calc_id)["outputs"]
    reports: list[dict[str, Any]] = []
    for case in case_set.cases:
        kwargs = inputs.decode(calc_id, inputs.encode(calc_id, case["inputs"]))
        by_backend: dict[str, Any] = {}
        passed = True
        for backend in backends:
            try:
                payload = run(calc_id, kwargs, backend)
            except Exception as error:
                by_backend[backend] = {
                    "status": "error",
                    "error": f"{type(error).__name__}: {error}",
                }
                passed = False
                continue
            expected_report = compare.against_expected(
                payload["result"],
                case["expected"],
                outputs,
                tolerance=float(case["tolerance"]),
                diagnostic_bound=diagnostic_bound,
            )
            by_backend[backend] = {
                "status": "ok",
                "result": payload["result"],
                "expected": expected_report,
            }
            if expected_report["disagreements"]:
                passed = False

        mirrored = None
        if len(backends) >= 2 and all(by_backend[b]["status"] == "ok" for b in backends):
            mirrored = compare.between_runs(
                by_backend[backends[0]]["result"],
                by_backend[backends[1]]["result"],
                tolerance=float(case["tolerance"]),
                diagnostic_bound=diagnostic_bound,
            )
            if mirrored["disagreements"]:
                passed = False

        reports.append(
            {
                "id": case["id"],
                "source": case["source"],
                "type": case["type"],
                "tolerance": case["tolerance"],
                "inputs": case["inputs"],
                "expected": case["expected"],
                "passed": passed,
                "backends": by_backend,
                "mirrored": mirrored,
            }
        )
    return reports
