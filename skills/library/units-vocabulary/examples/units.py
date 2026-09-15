"""The unit vocabulary: build, convert, and see a bare float refused."""

import azoth

q = azoth.ureg.Quantity
print(q(1000.0, "Pa").to("kPa"))
print(azoth.ureg.Quantity(1.5, "m/s"))

try:
    azoth.hydraulics.darcy_weisbach(0.02, 100.0, q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s"))
except azoth.UnitMismatchError as exc:
    print("refused:", exc)
