"""`process.expander`'s recovered-power constraint: the flags, the absent branch, the signature.

`capacity_exceeded` and `hard_limit_exceeded` are booleans with no `TestCase` slot, and the
no-rating branch is a `null` no case can state. The two thresholds are held by the cases;
this file reads the raw spec for the flags and pins the wrapper's whole input set, which
`_helpers.py` would otherwise drop an input through silently.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen
from azoth._dispatch import available, resolve, use_backend
from azoth.process import expander as run_expander

MODEL_ID = "process.expander"
SPEC = _models_gen.model(MODEL_ID)

FLAGGED = [
    case
    for case in SPEC["cases"]
    if "capacity_exceeded" in case["expected"] or "hard_limit_exceeded" in case["expected"]
]


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", FLAGGED, ids=lambda case: case["id"])
def test_the_capacity_flags_are_the_recorded_ones(case: dict[str, Any]) -> None:
    result = resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))
    expected = case["expected"]
    if "capacity_exceeded" in expected:
        assert result.capacity_exceeded is expected["capacity_exceeded"]
    if "hard_limit_exceeded" in expected:
        assert result.hard_limit_exceeded is expected["hard_limit_exceeded"]


@pytest.mark.parametrize("backend", sorted(available()))
def test_no_rating_has_no_per_constraint_output(backend: str) -> None:
    result = _run(_case("the_recovered_power_is_absent_without_a_rating"), backend)
    assert result.recovered_power_utilization is None
    assert result.bottleneck_constraint is None
    assert result.max_utilization == 0.0
    assert result.capacity_exceeded is False
    assert result.hard_limit_exceeded is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_two_thresholds_are_separate(backend: str) -> None:
    """1.2 kW is inside the hard band; 1.0 kW is outside both."""
    inside = _run(_case("a_rating_inside_the_hard_band_is_not_exceeded"), backend)
    assert inside.capacity_exceeded is False
    assert inside.hard_limit_exceeded is False

    outside = _run(_case("a_rating_below_the_recovered_power_trips_both_flags"), backend)
    assert outside.capacity_exceeded is True
    assert outside.hard_limit_exceeded is True


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_utilization_is_the_recovered_power_over_the_rating(backend: str) -> None:
    case = _case("a_rating_above_the_recovered_power_is_its_utilization")
    result = _run(case, backend)
    assert result.recovered_power_utilization == pytest.approx(
        case["expected"]["recovered_power_utilization"], rel=float(case["tolerance"])
    )
    assert result.bottleneck_constraint == "recoveredPower"


def test_the_public_wrapper_takes_every_declared_input() -> None:
    """The harness drops a case input the wrapper does not accept; this cannot."""
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_expander).parameters)
    assert accepted == declared, sorted(declared - accepted)
