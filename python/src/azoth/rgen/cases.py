"""The declared case set a calculation's claim is measured against.

Three sources, one list:

- the spec's top-level ``[worked_example]``;
- the spec's active ``[[tests]]`` entries that carry ``inputs`` and ``expected``;
- every external ``validation/**/*.json`` case whose ``calc`` names this id.

The spec's reference tests restate the external cases, so the union is deduplicated by
inputs and expected values and the external copy is preferred: it carries the source's
own verification flag. A property test is not a value case and is recorded separately;
it is not counted as a case that passed or failed.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from . import canonical
from .paths import repo_root

__all__ = ["CaseSet", "declared_cases", "strictest_tolerance"]


class CaseSet:
    """The value-carrying cases and the records that are not cases."""

    def __init__(self, cases: list[dict[str, Any]], not_value: list[dict[str, Any]]) -> None:
        self.cases = cases
        self.not_value = not_value

    def __len__(self) -> int:
        return len(self.cases)

    def as_dict(self) -> dict[str, Any]:
        return {"cases": self.cases, "not_value": self.not_value}


def _canonical_key(inputs: Mapping[str, Any], expected: Mapping[str, Any]) -> str:
    return canonical.hash_value({"inputs": inputs, "expected": expected})


def _case(
    *,
    identifier: str,
    source: str,
    case_type: str,
    inputs: Mapping[str, Any],
    expected: Mapping[str, Any],
    tolerance: float,
    verification: str | None = None,
    spec_ids: list[str] | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "id": identifier,
        "source": source,
        "type": case_type,
        "inputs": dict(inputs),
        "expected": dict(expected),
        "tolerance": tolerance,
    }
    if verification is not None:
        record["verification"] = verification
    if spec_ids:
        record["spec_ids"] = spec_ids
    return record


def _external(root: Path, calc_id: str) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for path in sorted((root / "validation").rglob("*.json")):
        try:
            document = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if document.get("calc") != calc_id:
            continue
        inputs = document.get("inputs")
        expected = document.get("expected")
        tolerance = document.get("tolerance")
        if not isinstance(inputs, dict) or not isinstance(expected, dict):
            continue
        if not isinstance(tolerance, (int, float)):
            continue
        source = document.get("source") or {}
        cases.append(
            _case(
                identifier=str(document.get("id", path.stem)),
                source="validation",
                case_type="reference",
                inputs=inputs,
                expected=expected,
                tolerance=float(tolerance),
                verification=source.get("verification"),
            )
        )
    return cases


def _from_spec(spec: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cases: list[dict[str, Any]] = []
    not_value: list[dict[str, Any]] = []

    worked = spec.get("worked_example")
    if (
        isinstance(worked, dict)
        and isinstance(worked.get("inputs"), dict)
        and isinstance(worked.get("expected"), dict)
    ):
        tolerance = worked.get("tolerance")
        if isinstance(tolerance, (int, float)):
            cases.append(
                _case(
                    identifier="worked_example",
                    source="spec",
                    case_type="worked_example",
                    inputs=worked["inputs"],
                    expected=worked["expected"],
                    tolerance=float(tolerance),
                )
            )

    for test in spec.get("tests", []):
        test_type = test.get("type")
        status = test.get("status")
        if test_type == "worked_example":
            continue
        if status != "active":
            not_value.append(
                {
                    "id": test.get("id"),
                    "type": test_type,
                    "status": "skipped",
                    "reason": test.get("skip_reason"),
                }
            )
            continue
        if test_type == "property":
            not_value.append(
                {
                    "id": test.get("id"),
                    "type": "property",
                    "property": test.get("property"),
                    "status": "not_a_value_case",
                }
            )
            continue
        inputs = test.get("inputs")
        expected = test.get("expected")
        tolerance = test.get("tolerance")
        if not (
            isinstance(inputs, dict)
            and isinstance(expected, dict)
            and isinstance(tolerance, (int, float))
        ):
            not_value.append(
                {
                    "id": test.get("id"),
                    "type": test_type,
                    "status": "no_case_values",
                }
            )
            continue
        cases.append(
            _case(
                identifier=str(test.get("id")),
                source="spec",
                case_type=str(test_type),
                inputs=inputs,
                expected=expected,
                tolerance=float(tolerance),
            )
        )
    return cases, not_value


def declared_cases(calc_id: str) -> CaseSet:
    """The union of the spec's and the external cases, deduplicated and ordered by id."""
    from azoth._registry_gen import spec

    spec_cases, not_value = _from_spec(spec(calc_id))

    root = repo_root()
    external = [] if root is None else _external(root, calc_id)

    merged: dict[str, dict[str, Any]] = {}
    for case in [*external, *spec_cases]:
        key = _canonical_key(case["inputs"], case["expected"])
        existing = merged.get(key)
        if existing is None:
            merged[key] = case
        elif existing["source"] == "validation" and case["source"] == "spec":
            existing.setdefault("spec_ids", []).append(case["id"])
        elif existing["source"] == "spec" and case["source"] == "validation":
            case.setdefault("spec_ids", []).append(existing["id"])
            merged[key] = case
        else:
            existing.setdefault("spec_ids", []).append(case["id"])

    cases = sorted(merged.values(), key=lambda case: str(case["id"]))
    return CaseSet(cases, not_value)


def strictest_tolerance(cases: CaseSet) -> float | None:
    """The tightest tolerance any declared case states, for an exploratory run."""
    if not cases.cases:
        return None
    return min(float(case["tolerance"]) for case in cases.cases)
