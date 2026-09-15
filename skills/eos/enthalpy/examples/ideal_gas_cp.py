"""Ideal-gas heat capacity from a constant polynomial."""

import azoth

q = azoth.ureg.Quantity
r = azoth.eos.ideal_gas_cp(
    q(30.0, "J/(mol*K)"),
    q(0.0, "J/(mol*K**2)"),
    q(0.0, "J/(mol*K**3)"),
    q(0.0, "J/(mol*K**4)"),
    q(0.0, "J/(mol*K**5)"),
    T=q(300.0, "K"),
)
print("cp:", r.cp)
