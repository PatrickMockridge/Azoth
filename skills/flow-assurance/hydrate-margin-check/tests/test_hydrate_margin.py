"""The hydrate equilibrium temperature `SKILL.md` documents."""

import pytest

import azoth

Q = azoth.ureg.Quantity

COMPONENTS = ["methane", "ethane", "propane", "water"]
Z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]


def test_the_equilibrium_temperature_falls_with_pressure() -> None:
    high = azoth.eos.hydrate_formation_temperature(COMPONENTS, P=Q(100.0, "bar"), z=Z, eos="srk")
    low = azoth.eos.hydrate_formation_temperature(COMPONENTS, P=Q(50.0, "bar"), z=Z, eos="srk")
    assert high.temperature.to("K").magnitude == pytest.approx(293.23, rel=1e-4)
    assert low.temperature.to("K").magnitude < high.temperature.to("K").magnitude
    assert str(high.structure) == "structure_ii"


def test_a_feed_without_water_is_refused() -> None:
    with pytest.raises(azoth.InvalidInputError, match="no water"):
        azoth.eos.hydrate_formation_temperature(
            ["methane", "ethane"], P=Q(100.0, "bar"), z=[0.9, 0.1], eos="srk"
        )
