"""A pure component's saturation pressure, below its critical temperature."""

import azoth

methane = azoth.eos.component("methane")
r = azoth.eos.pure_saturation(
    Tc=methane.Tc, Pc=methane.Pc, omega=methane.omega, T=azoth.ureg.Quantity(120.0, "K")
)
print("saturation pressure:", r.p_sat)
