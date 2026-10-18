"""A bound refuses what the arithmetic cannot carry, not only what is out of interval."""

from __future__ import annotations

import math

import pytest

from azoth.core.errors import OutOfRangeError
from azoth.core.range import RangeCheck, Severity


def test_a_non_finite_value_violates_even_an_unbounded_side() -> None:
    """An infinity passes a one-sided bound, and a NaN passes every comparison."""
    below = RangeCheck.from_spec(
        {
            "quantity": "T",
            "min": 0.0,
            "severity": "error",
            "rationale": "absolute temperature",
        }
    )
    assert below.violated(math.nan)
    assert below.violated(math.inf)
    assert below.violated(-math.inf)
    assert not below.violated(300.0)


def test_the_refusal_is_the_declared_error() -> None:
    check = RangeCheck.from_spec(
        {"quantity": "P", "min": 0.0, "severity": "error", "rationale": "absolute pressure"}
    )
    warnings: list[object] = []
    with pytest.raises(OutOfRangeError, match="P"):
        check.apply(math.inf, warnings)  # type: ignore[arg-type]
    assert not warnings
    assert check.severity is Severity.ERROR
