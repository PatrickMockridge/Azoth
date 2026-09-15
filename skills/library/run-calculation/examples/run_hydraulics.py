"""Run two hydraulics calculations and read their warnings.

Copyable as-is: only `azoth` is imported, and everything is computed from the
arguments rather than read from a file.
"""

import azoth

q = azoth.ureg.Quantity

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
print("pressure drop:", r.dp)
print("clean:", r.is_clean)
for warning in r.warnings:
    print("warning:", warning.code)

re = azoth.hydraulics.reynolds_number(
    q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
)
print("reynolds number:", round(re.re, 1))
print("regime:", re.regime)
