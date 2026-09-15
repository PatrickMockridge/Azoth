"""A pump's shaft power from flow, head and efficiency."""

import azoth

q = azoth.ureg.Quantity

r = azoth.hydraulics.pump_power(q(998.0, "kg/m**3"), q(0.01, "m**3/s"), q(30.0, "m"), 0.75)
print("shaft power:", r.power)
