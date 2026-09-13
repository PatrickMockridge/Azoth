"""Fluid property providers."""

from __future__ import annotations

import pytest

from azoth import ureg
from azoth.core.errors import PropertyUnavailableError
from azoth.core.units import Q
from azoth.properties import AIR, WATER, available_fluids, provider_for
from azoth.properties.provider import STANDARD_PRESSURE_PA, PropertyProvider, TableProvider


def celsius(value: float) -> Q:
    # Annotated on its own line: pint's `Quantity` is untyped, so returning it
    # directly would be an implicit Any at the boundary of every test here.
    quantity: Q = ureg.Quantity(value + 273.15, "K")
    return quantity


def test_the_built_in_fluids_are_water_and_air() -> None:
    """Exactly two, because the brief scopes property data to the worked
    examples."""
    assert available_fluids() == ("air", "water")


def test_providers_satisfy_the_protocol() -> None:
    """The providers must be usable through the interface, not just directly."""
    for provider in (WATER, AIR):
        assert isinstance(provider, PropertyProvider)


@pytest.mark.parametrize(("provider", "celsius_value"), [(WATER, 20.0), (AIR, 20.0)])
def test_tabulated_points_are_returned_exactly(
    provider: TableProvider, celsius_value: float
) -> None:
    """A tabulated temperature must return the tabulated value, not an
    interpolation that happens to land on it."""
    point = next(p for p in provider.points if p.temperature_c == celsius_value)
    assert provider.density(celsius(celsius_value)).magnitude == point.density_kg_m3
    assert (
        provider.dynamic_viscosity(celsius(celsius_value)).magnitude == point.dynamic_viscosity_pa_s
    )


def test_interpolation_is_between_the_bracketing_points() -> None:
    """Checked against a hand-computed linear interpolation."""
    low = next(p for p in WATER.points if p.temperature_c == 20.0)
    high = next(p for p in WATER.points if p.temperature_c == 40.0)
    midpoint = 30.0
    expected = (low.density_kg_m3 + high.density_kg_m3) / 2.0

    assert WATER.density(celsius(midpoint)).magnitude == pytest.approx(expected, rel=1e-12)


def test_water_density_falls_with_temperature() -> None:
    """A physical sanity check that would catch a mis-ordered table."""
    values = [WATER.density(celsius(t)).magnitude for t in (20.0, 40.0, 60.0, 80.0)]
    assert values == sorted(values, reverse=True), f"water density is not decreasing: {values}"


def test_water_viscosity_falls_with_temperature() -> None:
    """The trend is the whole reason these are tables rather than constants:
    viscosity drops by a factor of six across 0-100 C."""
    values = [WATER.dynamic_viscosity(celsius(t)).magnitude for t in (20.0, 40.0, 60.0, 80.0)]
    assert values == sorted(values, reverse=True), f"viscosity is not decreasing: {values}"
    assert values[0] / values[-1] > 2.0, "viscosity should vary substantially over this range"


@pytest.mark.parametrize("temperature", [-10.0, 150.0])
def test_out_of_range_temperature_raises_rather_than_extrapolating(temperature: float) -> None:
    """Extrapolation would be a confident wrong number.

    Water's viscosity changes by a factor of six across the tabulated range, so a
    straight-line extension past either end is not a small error.
    """
    with pytest.raises(PropertyUnavailableError, match="does not extrapolate"):
        WATER.density(celsius(temperature))


def test_non_atmospheric_pressure_is_refused() -> None:
    """The built-in tables are at 1 atm only.

    For a gas, silently returning a density for the wrong pressure would be wrong
    by roughly the pressure ratio - and would look entirely reasonable.
    """
    with pytest.raises(PropertyUnavailableError, match="1 atm only"):
        AIR.density(celsius(20.0), ureg.Quantity(10.0, "bar"))
    # Atmospheric pressure is accepted.
    assert AIR.density(celsius(20.0), ureg.Quantity(STANDARD_PRESSURE_PA, "Pa")).magnitude > 0


def test_unknown_fluid_is_an_error_not_a_default() -> None:
    """Substituting water for an unknown fluid would produce a plausible pressure
    drop for the wrong substance."""
    with pytest.raises(PropertyUnavailableError, match="not a built-in fluid"):
        provider_for("unobtainium")


def test_fluid_names_are_normalised() -> None:
    """Case and surrounding whitespace should not decide whether a lookup works."""
    assert provider_for("  WATER ") is WATER


def test_a_bare_number_is_rejected_where_a_temperature_is_expected() -> None:
    """The units boundary applies to properties too."""
    from azoth.core.errors import UnitMismatchError

    with pytest.raises(UnitMismatchError):
        WATER.density(293.15)  # type: ignore[arg-type]
