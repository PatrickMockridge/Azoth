"""`process.splitter`'s mechanical design: the flag, the refusals and the sized records.

The nine outputs are `Option`s, so the flag-off branch is a `null` no case file can state;
it is asserted here. The refusals are the class's two combinations - a design input without
the flag, and a stated header with it - plus the two divisions the port guards. Everything
else the capture's cases hold.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
import azoth
from azoth import _models_gen
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError

MODEL_ID = "process.splitter"
SPEC = _models_gen.model(MODEL_ID)

q = azoth.ureg.Quantity

DESIGN_FIELDS = (
    "design_header_diameter",
    "design_branch_diameter",
    "design_header_length",
    "design_wall_thickness",
    "manifold_volume",
    "design_total_pressure_drop",
    "header_weight",
    "branch_weight",
    "total_weight",
)

UNITS = {
    "design_header_diameter": "m",
    "design_branch_diameter": "m",
    "design_header_length": "m",
    "design_wall_thickness": "m",
    "manifold_volume": "m**3",
    "design_total_pressure_drop": "Pa",
    "header_weight": "kg",
    "branch_weight": "kg",
    "total_weight": "kg",
}


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_design_outputs_are_absent_without_their_flag(backend: str) -> None:
    """The class does not run `calcDesign` unless it is asked; every output is `None`."""
    result = _run(_case("a_two_way_split_of_a_liquid"), backend)
    for name in DESIGN_FIELDS:
        assert getattr(result, name) is None, name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_default_design_is_the_captures(backend: str) -> None:
    """The 1 mol/s design at 15 m/s: both NPS sizes 1/2", the 3 mm floor and its weights."""
    case = _case("the_default_design_sizes_at_fifteen_metres_per_second")
    result = _run(case, backend)
    for name in DESIGN_FIELDS:
        value = getattr(result, name)
        assert value is not None, name
        assert value.to(UNITS[name]).magnitude == pytest.approx(
            case["expected"][name], rel=float(case["tolerance"])
        ), name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_design_and_the_capacity_map_read_one_header(backend: str) -> None:
    """Block 10: the velocity constraint answers against the NPS-rounded header."""
    case = _case("the_velocity_reads_the_designed_header")
    result = _run(case, backend)
    assert result.velocity_value is not None
    assert result.velocity_value.to("m/s").magnitude == pytest.approx(
        case["expected"]["velocity_value"], rel=float(case["tolerance"])
    )
    assert result.velocity_design.to("m/s").magnitude == pytest.approx(30.0)
    assert result.bottleneck_constraint == "velocity"


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_two_phase_design_reads_the_system_density(backend: str) -> None:
    case = _case("the_two_phase_design_reads_the_system_density")
    result = _run(case, backend)
    assert result.design_header_diameter.to("m").magnitude == pytest.approx(0.0381)
    assert result.design_total_pressure_drop.to("Pa").magnitude == pytest.approx(
        case["expected"]["design_total_pressure_drop"], rel=float(case["tolerance"])
    )


def test_a_design_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("a_two_way_split_of_a_liquid")["inputs"])
    with pytest.raises(InvalidInputError, match="design_velocity"):
        resolve(MODEL_ID)(**base, design_velocity=q(10.0, "m/s"))


def test_a_stated_header_with_the_flag_is_refused() -> None:
    case = _case("the_default_design_sizes_at_fifteen_metres_per_second")
    with pytest.raises(InvalidInputError, match="header_diameter"):
        resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]), header_diameter=q(0.02, "m"))


def test_a_non_positive_design_velocity_is_refused() -> None:
    """The spec's hard range fires first; the kernel guard is the dispatcher's backstop."""
    base = h.model_kwargs(SPEC, _case("a_two_way_split_of_a_liquid")["inputs"])
    with pytest.raises(OutOfRangeError, match="design_velocity"):
        resolve(MODEL_ID)(
            **base,
            mechanical_design=True,
            design_velocity=q(0.0, "m/s"),
        )
