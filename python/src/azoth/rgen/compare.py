"""Field-by-field comparison of a result with its twin and with an expected case.

The tolerance rule is the one the cross-implementation tests already use: a relative
comparison against the expected value, with a tiny floor so a true zero is not a division
by zero. A mirrored pair has no expected value, so the same rule is applied in both
directions and the larger of the two relative errors decides.

Two exceptions come from the same tests. A solver diagnostic is compared on the bound the
solver declares, because it is near zero by construction and a relative comparison there
divides by the vanishing quantity. A path quantity is not compared at all: an iteration
count reports where the arithmetic stopped, not what it computed.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, Final

from azoth.core.units import quantity

from . import policy

__all__ = ["against_expected", "between_runs"]

#: The smallest positive normal double, so a zero expected value has a scale.
_TINY: Final[float] = float.fromhex("0x0.0000000000001p-1022")


def _base(magnitude: float, unit: str) -> float:
    return float(quantity(magnitude, unit).to_base_units().magnitude)


def _quantity(value: Any) -> tuple[float, str] | None:
    if isinstance(value, Mapping) and set(value) == {"magnitude", "unit"}:
        return float(value["magnitude"]), str(value["unit"])
    return None


def _numeric(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _leaf(path: str) -> str:
    """The field name a dotted/indexed path ends in, for the policy sets."""
    return path.rsplit(".", 1)[-1].split("[", 1)[0]


def _relative(left: float, right: float, scale: float, tolerance: float) -> bool:
    difference = abs(left - right)
    return difference / max(scale, _TINY) <= tolerance


def _diagnostic(left: float, right: float, bound: float, tolerance: float) -> bool:
    if math.isnan(left) or math.isnan(right):
        return math.isnan(left) and math.isnan(right)
    if math.isinf(left) or math.isinf(right):
        return left == right
    difference = abs(left - right)
    return difference <= bound or _relative(left, right, max(abs(left), abs(right)), tolerance)


def _warnings_key(value: Any) -> list[tuple[str, str | None]]:
    if not isinstance(value, list):
        return []
    keys: list[tuple[str, str | None]] = []
    for entry in value:
        if isinstance(entry, Mapping):
            code = entry.get("code")
            field = entry.get("field")
            keys.append((str(code), None if field is None else str(field)))
        else:
            keys.append((str(entry), None))
    return keys


def _scalar_disagreement(
    path: str,
    left: Any,
    right: Any,
    tolerance: float,
    *,
    diagnostic_bound: float | None,
    scale: float | None,
) -> dict[str, Any] | None:
    if left is None and right is None:
        return None
    if left is None or right is None:
        return {"path": path, "claimant": left, "validator": right, "reason": "absent on one side"}

    left_q = _quantity(left)
    right_q = _quantity(right)
    if left_q is not None or right_q is not None:
        if left_q is None or right_q is None:
            return {
                "path": path,
                "claimant": left,
                "validator": right,
                "reason": "quantity against a number",
            }
        left_number = _base(*left_q)
        right_number = _base(*right_q)
    elif _numeric(left) and _numeric(right):
        left_number = float(left)
        right_number = float(right)
    elif isinstance(left, (bool, str)) and isinstance(right, (bool, str)):
        return None if left == right else {"path": path, "claimant": left, "validator": right}
    else:
        return {"path": path, "claimant": left, "validator": right, "reason": "different kinds"}

    if diagnostic_bound is not None and _leaf(path) in policy.DIAGNOSTIC_FIELDS:
        ok = _diagnostic(left_number, right_number, diagnostic_bound, tolerance)
    else:
        scale_value = (
            scale if scale is not None else max(abs(left_number), abs(right_number), _TINY)
        )
        ok = _relative(left_number, right_number, scale_value, tolerance)
    return None if ok else {"path": path, "claimant": left, "validator": right}


def between_runs(
    claimant: Mapping[str, Any],
    validator: Mapping[str, Any],
    *,
    tolerance: float,
    diagnostic_bound: float | None,
) -> dict[str, Any]:
    """Compare the two implementations' result documents field by field."""
    disagreements: list[dict[str, Any]] = []
    uncompared: list[dict[str, Any]] = []
    compared = 0

    def walk(path: str, left: Any, right: Any) -> None:
        nonlocal compared
        if _leaf(path) in policy.UNCOMPARED_FIELDS:
            uncompared.append(
                {"path": path, "claimant": left is not None, "validator": right is not None}
            )
            return
        if path == "warnings":
            compared += 1
            if _warnings_key(left) != _warnings_key(right):
                disagreements.append({"path": path, "claimant": left, "validator": right})
            return
        # A quantity is two leaves that mean one value; recursing into them would compare
        # a millimetre magnitude with a metre one and call the units a disagreement.
        if _quantity(left) is not None or _quantity(right) is not None:
            compared += 1
            problem = _scalar_disagreement(
                path,
                left,
                right,
                tolerance,
                diagnostic_bound=diagnostic_bound,
                scale=None,
            )
            if problem is not None:
                disagreements.append(problem)
            return
        left_is_mapping = isinstance(left, Mapping)
        right_is_mapping = isinstance(right, Mapping)
        if left_is_mapping or right_is_mapping:
            if not (left_is_mapping and right_is_mapping):
                disagreements.append(
                    {
                        "path": path,
                        "claimant": left,
                        "validator": right,
                        "reason": "different kinds",
                    }
                )
                return
            for key in sorted(set(left) | set(right)):
                child = f"{path}.{key}" if path else str(key)
                walk(child, left.get(key), right.get(key))
            return
        if isinstance(left, list) or isinstance(right, list):
            if not (isinstance(left, list) and isinstance(right, list)):
                disagreements.append(
                    {
                        "path": path,
                        "claimant": left,
                        "validator": right,
                        "reason": "different kinds",
                    }
                )
                return
            if len(left) != len(right):
                disagreements.append(
                    {
                        "path": path,
                        "claimant": left,
                        "validator": right,
                        "reason": "different lengths",
                    }
                )
                return
            for index, (left_item, right_item) in enumerate(zip(left, right, strict=True)):
                walk(f"{path}[{index}]", left_item, right_item)
            return
        compared += 1
        problem = _scalar_disagreement(
            path,
            left,
            right,
            tolerance,
            diagnostic_bound=diagnostic_bound,
            scale=None,
        )
        if problem is not None:
            disagreements.append(problem)

    for key in sorted(set(claimant) | set(validator)):
        walk(str(key), claimant.get(key), validator.get(key))
    return {
        "checked": True,
        "tolerance": tolerance,
        "diagnostic_bound": diagnostic_bound,
        "fields_compared": compared,
        "disagreements": disagreements,
        "uncompared": uncompared,
    }


def against_expected(
    result: Mapping[str, Any],
    expected: Mapping[str, Any],
    outputs: Mapping[str, Any],
    *,
    tolerance: float,
    diagnostic_bound: float | None,
) -> dict[str, Any]:
    """Compare a result document with one case's expected values.

    The expected magnitudes are in the spec's declared output unit; the result's quantity
    is rebuilt from its own unit and converted to base units, so the two meet in SI without
    either side assuming the other's spelling.
    """
    disagreements: list[dict[str, Any]] = []
    compared = 0
    for name, wanted in expected.items():
        declaration = outputs.get(name, {})
        unit = declaration.get("unit", "dimensionless")
        actual = result.get(name)
        if actual is None:
            disagreements.append(
                {"path": name, "claimant": None, "validator": wanted, "reason": "result absent"}
            )
            continue
        actual_q = _quantity(actual)
        if actual_q is not None:
            actual_number = _base(*actual_q)
        elif _numeric(actual):
            actual_number = float(actual)
        elif isinstance(actual, (str, bool)):
            compared += 1
            if actual != wanted:
                disagreements.append({"path": name, "claimant": actual, "validator": wanted})
            continue
        else:
            disagreements.append(
                {"path": name, "claimant": actual, "validator": wanted, "reason": "uncheckable"}
            )
            continue
        if not _numeric(wanted):
            disagreements.append(
                {"path": name, "claimant": actual, "validator": wanted, "reason": "uncheckable"}
            )
            continue
        compared += 1
        wanted_number = float(quantity(float(wanted), unit).to_base_units().magnitude)
        if diagnostic_bound is not None and name in policy.DIAGNOSTIC_FIELDS:
            ok = _diagnostic(actual_number, wanted_number, diagnostic_bound, tolerance)
        else:
            ok = _relative(actual_number, wanted_number, abs(wanted_number), tolerance)
        if not ok:
            disagreements.append({"path": name, "claimant": actual, "validator": wanted})
    return {
        "checked": True,
        "tolerance": tolerance,
        "fields_compared": compared,
        "disagreements": disagreements,
    }
