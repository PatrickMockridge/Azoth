"""The plane-wall conduction `SKILL.md` documents."""

import azoth


def test_conduction_returns_heat_flow() -> None:
    r = azoth.thermal.conduction_plane_wall(
        azoth.ureg.Quantity(0.5, "W/(m*K)"),
        azoth.ureg.Quantity(2.0, "m**2"),
        azoth.ureg.Quantity(30.0, "K"),
        azoth.ureg.Quantity(0.1, "m"),
    )
    assert r.q.magnitude > 0
