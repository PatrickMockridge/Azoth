"""Steady heat flow through a plane wall."""

import azoth

r = azoth.thermal.conduction_plane_wall(
    azoth.ureg.Quantity(0.5, "W/(m*K)"),
    azoth.ureg.Quantity(2.0, "m**2"),
    azoth.ureg.Quantity(30.0, "K"),
    azoth.ureg.Quantity(0.1, "m"),
)
print("heat flow:", r.q)
