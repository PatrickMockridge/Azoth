"""`process.filter`'s loaded state: the booleans, the defaults and the refusals.

The spec's cases pin the two implementations to the captures. These tests are for what a
case cannot say: `replacement_required` has no slot in a generated case row, the absent
loading inputs are the class's zero state, and the two setters that throw are refused.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import OutOfRangeError
from azoth.process import filter as run_filter

MODEL_ID = "process.filter"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    """The harness drops a case input the wrapper does not accept; this cannot."""
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_filter).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_no_loading_inputs_is_the_class_s_zero_state(backend: str) -> None:
    result = _run(_case("pressure_drop_1_bar"), backend)
    assert result.solids_loading.to("kg").magnitude == 0.0
    assert result.loading_fraction == 0.0
    assert result.breakthrough_fraction == 0.0
    assert result.bypass_fraction == 0.0
    assert result.unrestricted_drop.to("Pa").magnitude == result.clean_drop.to("Pa").magnitude
    assert result.differential_pressure_utilization == pytest.approx(0.2)
    assert result.replacement_required is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_non_positive_capacity_is_the_sentinel(backend: str) -> None:
    result = _run(_case("a_non_positive_capacity_is_no_capacity"), backend)
    assert result.loading_fraction == 0.0
    assert result.breakthrough_fraction == 0.0
    assert result.unrestricted_drop.to("Pa").magnitude == result.clean_drop.to("Pa").magnitude


@pytest.mark.parametrize("backend", sorted(available()))
def test_replacement_required_reaches_both_branches(backend: str) -> None:
    terminal = _run(_case("the_terminal_drop_calls_for_replacement"), backend)
    assert terminal.differential_pressure_utilization == pytest.approx(4.0)
    assert terminal.replacement_required is True

    loading = _run(_case("loading_past_capacity_breaks_through"), backend)
    assert loading.loading_fraction == pytest.approx(1.1)
    assert loading.differential_pressure_utilization < 1.0
    assert loading.replacement_required is True

    boundary = _run(_case("breakthrough_starts_above_its_boundary"), backend)
    assert boundary.replacement_required is False


def test_a_non_positive_bypass_pressure_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="bypass_cracking_pressure"):
        resolve(MODEL_ID)(**base, bypass_cracking_pressure=Q(0.0, "Pa"))


def test_a_non_positive_terminal_drop_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="terminal_pressure_drop"):
        resolve(MODEL_ID)(**base, terminal_pressure_drop=Q(0.0, "Pa"))


def test_a_non_finite_loading_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="solids_loading"):
        resolve(MODEL_ID)(**base, solids_loading=Q(float("nan"), "kg"))
