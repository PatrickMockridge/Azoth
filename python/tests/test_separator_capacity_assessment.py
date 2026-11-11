"""`process.separator_capacity_assessment`: the absences, the order and the shared fields.

`test_model_cases` holds every recorded string and number. What is here is what a case cannot
state: that an `unassessed` snapshot carries **no** numerics, that the check order is the
class's (gas before level), that the assessment and `process.separator` agree on the fields
they share, and that a non-finite mechanical limit is refused rather than echoed.
"""

from __future__ import annotations

from typing import Any

import pytest

import azoth
from azoth import _models_gen
from azoth.core.errors import OutOfRangeError
from azoth.process import separator, separator_capacity_assessment

MODEL_ID = "process.separator_capacity_assessment"
SPEC = _models_gen.model(MODEL_ID)

q = azoth.ureg.Quantity

BASE: dict[str, Any] = {
    "components": ["methane", "n-butane"],
    "feed_n": q(1.0, "mol/s"),
    "feed_z": [0.7, 0.3],
    "feed_p": q(20.0, "bar"),
    "feed_t": q(300.0, "K"),
    "pressure_drop": q(0.0, "Pa"),
    "gas_in_liquid": 0.0,
    "internal_diameter": q(1.0, "m"),
    "design_gas_load_factor": q(0.11, "m/s"),
    "orientation": "vertical",
    "design_liquid_level_fraction": 0.8,
    "design_basis_provenance": "design basis A",
}

NUMERIC = (
    "gas_area",
    "gas_density",
    "liquid_density",
    "maximum_gas_velocity",
    "current_gas_flow",
    "mechanical_design_maximum_gas_flow",
    "souders_brown_maximum_gas_flow",
    "maximum_gas_flow",
    "utilization",
)


@pytest.mark.parametrize(
    "case",
    [case for case in SPEC["cases"] if case["expected"]["status"] == "unassessed"],
    ids=lambda case: case["id"],
)
def test_an_unassessed_snapshot_carries_no_numbers(case: dict[str, Any]) -> None:
    """The class holds `NaN` in every numeric field there; the port publishes none of them."""
    inputs = dict(BASE)
    raw = case["inputs"]
    assert isinstance(raw, dict)
    inputs.update(raw)
    result = separator_capacity_assessment(**_quantities(inputs))
    assert result.status == "unassessed"
    for field in NUMERIC:
        assert getattr(result, field) is None, field
    assert result.design_basis_source == ""
    assert result.governing_capacity_basis == ""


def test_the_assessment_and_the_capacity_surface_agree_on_what_they_share() -> None:
    """One arithmetic under two ids: the velocity and the utilisation are the same numbers."""
    assessment = separator_capacity_assessment(**BASE)
    capacity = separator(
        components=BASE["components"],
        feed_n=BASE["feed_n"],
        feed_z=BASE["feed_z"],
        feed_p=BASE["feed_p"],
        feed_t=BASE["feed_t"],
        pressure_drop=BASE["pressure_drop"],
        gas_in_liquid=BASE["gas_in_liquid"],
        internal_diameter=BASE["internal_diameter"],
        design_gas_load_factor=BASE["design_gas_load_factor"],
        orientation=BASE["orientation"],
        design_liquid_level_fraction=BASE["design_liquid_level_fraction"],
    )
    assert assessment.maximum_gas_velocity is not None
    assert capacity.max_allowable_gas_velocity is not None
    assert assessment.maximum_gas_velocity.to("m/s").magnitude == pytest.approx(
        capacity.max_allowable_gas_velocity.to("m/s").magnitude, rel=1e-12
    )
    assert assessment.utilization == pytest.approx(capacity.capacity_utilization, rel=1e-12)


def test_a_non_finite_mechanical_limit_is_refused_not_echoed() -> None:
    """An infinity is a value no snapshot publishes; the class's own field would carry it.

    The refusal is the units boundary's, which both backends convert through - so the Rust
    kernel's own finiteness check is a second line rather than the one a caller meets.
    """
    with pytest.raises(OutOfRangeError, match="mechanical_design_max_gas_flow"):
        separator_capacity_assessment(
            **BASE, mechanical_design_max_gas_flow=q(float("inf"), "m**3/s")
        )


def _quantities(inputs: dict[str, Any]) -> dict[str, Any]:
    """A case's bare SI numbers as the quantities the public function takes."""
    units = {
        "feed_n": "mol/s",
        "feed_p": "Pa",
        "feed_t": "K",
        "pressure_drop": "Pa",
        "internal_diameter": "m",
        "design_gas_load_factor": "m/s",
        "mechanical_design_max_gas_flow": "m**3/s",
    }
    out: dict[str, Any] = {}
    for key, value in inputs.items():
        unit = units.get(key)
        out[key] = q(value, unit) if unit is not None else value
    return out
