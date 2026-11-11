"""`process.separator`'s capacity surface: the flag, the refusals and the floor.

The five getters are cases in their own right, and `test_model_cases` holds every field but
``is_overloaded`` - a boolean has no slot in the generated `TestCase`, so it is
`SELF_ASSERTED_EXPECTATIONS`' and this file reads it from the raw spec. The rest of what is
here is what no capture can state: a vessel that is half given, an orientation that is
neither of the two, and the near-dry floor on densities no captured state reaches.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
import azoth
from azoth import _models_gen
from azoth._dispatch import resolve
from azoth.core.errors import InvalidInputError
from azoth.process import separator as run_separator
from azoth.process.reference.separator import _load_factor_liquid_density

MODEL_ID = "process.separator"
SPEC = _models_gen.model(MODEL_ID)

FLAGGED = [case for case in SPEC["cases"] if "is_overloaded" in case["expected"]]

q = azoth.ureg.Quantity


@pytest.mark.parametrize("case", FLAGGED, ids=lambda case: case["id"])
def test_the_overload_flag_is_the_recorded_one(case: dict[str, Any]) -> None:
    """`isOverloaded` is `capacity_utilization > 1.0`, exactly, and the case pins both."""
    result = resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))
    expected = case["expected"]
    assert result.is_overloaded is expected["is_overloaded"]
    assert result.capacity_utilization == pytest.approx(
        expected["capacity_utilization"], rel=float(case["tolerance"])
    )


def test_a_vessel_with_no_gas_uses_no_capacity() -> None:
    """Liquid only: utilisation zero and overload false, and no gas metric invented."""
    result = run_separator(
        ["n-butane"],
        q(1.0, "mol/s"),
        [1.0],
        q(5.0, "bar"),
        q(250.0, "K"),
        q(0.0, "Pa"),
        0.0,
        internal_diameter=q(1.0, "m"),
        design_gas_load_factor=q(0.11, "m/s"),
        orientation="vertical",
    )
    assert result.capacity_utilization == 0.0
    assert result.is_overloaded is False
    # The class's velocity getters would read phase 0, which is the liquid.
    assert result.max_allowable_gas_velocity is None
    assert result.gas_superficial_velocity is None
    assert result.gas_load_factor is None


def test_no_vessel_means_no_capacity_outputs() -> None:
    result = run_separator(
        ["methane", "n-butane"],
        q(1.0, "mol/s"),
        [0.7, 0.3],
        q(20.0, "bar"),
        q(300.0, "K"),
        q(0.0, "Pa"),
        0.0,
    )
    assert result.capacity_utilization is None
    assert result.max_allowable_gas_velocity is None
    assert result.gas_superficial_velocity is None
    assert result.gas_load_factor is None
    assert result.is_overloaded is None


def test_half_a_vessel_is_refused_both_ways() -> None:
    base = (
        ["methane", "n-butane"],
        q(1.0, "mol/s"),
        [0.7, 0.3],
        q(20.0, "bar"),
        q(300.0, "K"),
        q(0.0, "Pa"),
        0.0,
    )
    with pytest.raises(InvalidInputError, match="design_gas_load_factor"):
        run_separator(*base, internal_diameter=q(1.0, "m"))
    with pytest.raises(InvalidInputError, match="internal_diameter"):
        run_separator(*base, design_gas_load_factor=q(0.11, "m/s"))


def test_the_gas_area_inputs_need_the_vessel_they_shape() -> None:
    base = (
        ["methane", "n-butane"],
        q(1.0, "mol/s"),
        [0.7, 0.3],
        q(20.0, "bar"),
        q(300.0, "K"),
        q(0.0, "Pa"),
        0.0,
    )
    with pytest.raises(InvalidInputError, match="orientation"):
        run_separator(*base, orientation="vertical")
    with pytest.raises(InvalidInputError, match="design_liquid_level_fraction"):
        run_separator(*base, design_liquid_level_fraction=0.5)


def test_an_orientation_that_is_neither_is_refused() -> None:
    with pytest.raises(InvalidInputError, match="orientation"):
        run_separator(
            ["methane", "n-butane"],
            q(1.0, "mol/s"),
            [0.7, 0.3],
            q(20.0, "bar"),
            q(300.0, "K"),
            q(0.0, "Pa"),
            0.0,
            internal_diameter=q(1.0, "m"),
            design_gas_load_factor=q(0.11, "m/s"),
            orientation="inclined",
        )


def test_the_near_dry_floor_is_the_load_factor_s_only() -> None:
    """A constructed pair, because no captured state reaches the floor: both are dry."""
    # `rho_l - rho_g` under `10.0` substitutes `1000.0`, and exactly ten does not.
    assert _load_factor_liquid_density(1009.0, 1000.0) == 1000.0
    assert _load_factor_liquid_density(1010.0, 1000.0) == 1010.0
    # A real liquid is left alone.
    assert _load_factor_liquid_density(542.32, 20.16) == 542.32
