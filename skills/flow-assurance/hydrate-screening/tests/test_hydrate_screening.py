"""The three hydrate screens `SKILL.md` documents."""

import pytest

import azoth

Q = azoth.ureg.Quantity

COMPONENTS = ["methane", "ethane", "propane", "water"]
Z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]


def test_the_pressure_is_the_temperature_curve_inverted() -> None:
    at_100_bar = azoth.eos.hydrate_formation_temperature(
        COMPONENTS, P=Q(100.0, "bar"), z=Z, eos="srk"
    )
    at_that_temperature = azoth.eos.hydrate_formation_pressure(
        COMPONENTS, T=Q(at_100_bar.temperature.to("K").magnitude, "K"), z=Z, eos="srk"
    )
    # The same point, reached from either axis, to the bisections' tolerance.
    assert at_that_temperature.pressure.to("bar").magnitude == pytest.approx(100.0, rel=1e-3)


def test_the_fraction_closes_the_material_balance() -> None:
    formed = azoth.eos.hydrate_fraction(
        COMPONENTS, T=Q(283.15, "K"), P=Q(100.0, "bar"), z=Z, eos="srk"
    )
    assert formed.beta == pytest.approx(0.1153, rel=1e-3)
    assert abs(formed.balance_error) < 1e-9
