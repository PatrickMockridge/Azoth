"""`process.splitter`'s capacity map: the flag, the absent branches and the signature.

`test_model_cases` holds every carved number; ``capacity_exceeded`` is a boolean with no
slot in `TestCase`, so it is `SELF_ASSERTED_EXPECTATIONS`' and this file reads it from the
raw spec. The rest of what is here is what no case file can state: a `null` (the six
`Option`s of a map that does not hold the constraint) and the class's own bottleneck string
on the branches where the strict `>` leaves it absent.

**The signature check is the harness's blind spot.** `_helpers.py` drops a case input the
public wrapper does not accept, so a wrapper missing one of these four would pass every
generated case on defaults. Reading `inspect.signature` against the spec's inputs is the
one place that cannot.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
import azoth
from azoth import _models_gen
from azoth._dispatch import available, resolve, use_backend
from azoth.process import splitter as run_splitter

MODEL_ID = "process.splitter"
SPEC = _models_gen.model(MODEL_ID)

q = azoth.ureg.Quantity

FLAGGED = [case for case in SPEC["cases"] if "capacity_exceeded" in case["expected"]]

PER_CONSTRAINT = (
    "pressure_drop_value",
    "pressure_drop_design",
    "pressure_drop_utilization",
    "velocity_value",
    "velocity_design",
    "velocity_utilization",
)


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", FLAGGED, ids=lambda case: case["id"])
def test_the_capacity_flag_is_the_recorded_one(case: dict[str, Any]) -> None:
    """`isCapacityExceeded` is a case's own boolean, read from the raw spec."""
    result = resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))
    assert result.capacity_exceeded is case["expected"]["capacity_exceeded"]


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_map_with_no_setters_has_no_per_constraint_outputs(backend: str) -> None:
    """The class's fresh `Splitter` has no constraints, and an absent one is `None`."""
    result = _run(_case("the_capacity_map_is_absent_until_a_setter_creates_it"), backend)
    for name in PER_CONSTRAINT:
        assert getattr(result, name) is None, name
    assert result.bottleneck_constraint is None
    assert result.max_utilization == 0.0
    assert result.capacity_exceeded is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_header_diameter_alone_does_not_initialize_the_map(backend: str) -> None:
    """`setHeaderDiameter` stores a design field; the map is still empty."""
    result = _run(_case("the_header_diameter_alone_does_not_initialize_the_map"), backend)
    for name in PER_CONSTRAINT:
        assert getattr(result, name) is None, name
    assert result.bottleneck_constraint is None
    assert result.capacity_exceeded is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_non_positive_limit_suppresses_only_its_own_constraint(backend: str) -> None:
    drop = _run(_case("a_non_positive_drop_suppresses_only_its_own_constraint"), backend)
    assert drop.pressure_drop_value is None
    assert drop.pressure_drop_design is None
    assert drop.pressure_drop_utilization is None
    assert drop.velocity_value is not None
    assert drop.velocity_design is not None
    assert drop.bottleneck_constraint == "velocity"

    velocity = _run(_case("a_non_positive_velocity_suppresses_only_its_own_constraint"), backend)
    assert velocity.velocity_value is None
    assert velocity.velocity_design is None
    assert velocity.velocity_utilization is None
    assert velocity.pressure_drop_value is not None
    assert velocity.pressure_drop_design is not None
    assert velocity.bottleneck_constraint is None
    assert velocity.max_utilization == 0.0


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_bottleneck_is_the_class_key_and_the_cap_is_9_99(backend: str) -> None:
    tight = _run(_case("a_tight_velocity_is_exceeded_with_analysis_disabled"), backend)
    assert tight.bottleneck_constraint == "velocity"
    assert tight.capacity_exceeded is True

    defaults = _run(_case("enabling_analysis_initializes_the_class_defaults"), backend)
    assert defaults.bottleneck_constraint is None
    assert defaults.capacity_exceeded is False

    capped = _run(_case("the_utilization_is_capped_at_9_99"), backend)
    assert capped.velocity_utilization == pytest.approx(9.99)
    assert capped.max_utilization == pytest.approx(9.99)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_two_phase_velocity_reads_the_feeds_own_volume_flow(backend: str) -> None:
    case = _case("the_velocity_reads_the_two_phase_feeds_own_volume_flow")
    result = _run(case, backend)
    expected = case["expected"]
    tolerance = float(case["tolerance"])
    assert result.velocity_value.to("m/s").magnitude == pytest.approx(
        expected["velocity_value"], rel=tolerance
    )
    assert result.velocity_utilization == pytest.approx(
        expected["velocity_utilization"], rel=tolerance
    )
    assert result.bottleneck_constraint == "velocity"


def test_the_public_wrapper_takes_every_declared_input() -> None:
    """The harness drops a case input the wrapper does not accept; this cannot."""
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_splitter).parameters)
    assert accepted == declared, sorted(declared - accepted)
