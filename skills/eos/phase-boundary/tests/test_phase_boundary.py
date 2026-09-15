"""The pure-component saturation `SKILL.md` documents."""

import azoth


def test_pure_saturation_returns_a_pressure() -> None:
    methane = azoth.eos.component("methane")
    r = azoth.eos.pure_saturation(
        Tc=methane.Tc, Pc=methane.Pc, omega=methane.omega, T=azoth.ureg.Quantity(120.0, "K")
    )
    assert r.p_sat.magnitude > 0
