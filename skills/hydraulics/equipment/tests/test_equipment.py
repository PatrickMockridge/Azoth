"""The pump-power calc `SKILL.md` documents."""

import azoth


def test_pump_power_returns_positive() -> None:
    q = azoth.ureg.Quantity
    r = azoth.hydraulics.pump_power(q(998.0, "kg/m**3"), q(0.01, "m**3/s"), q(30.0, "m"), 0.75)
    assert r.power.magnitude > 0
