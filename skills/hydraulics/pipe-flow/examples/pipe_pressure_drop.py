"""Reynolds number, then a Darcy-Weisbach pressure drop."""

import azoth

q = azoth.ureg.Quantity

re = azoth.hydraulics.reynolds_number(
    q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
)
print("reynolds number:", round(re.re, 1), re.regime)

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
print("pressure drop:", r.dp)
print("warnings:", [str(w.code) for w in r.warnings])
