"""Build a mixture, look a component up, and flash it."""

import azoth

methane = azoth.eos.component("methane")
print("methane Tc:", methane.Tc)

fluid = azoth.eos.from_names(["methane", "n-butane"])
print("components:", len(fluid))

flash = azoth.eos.pt_flash(
    fluid, T=azoth.ureg.Quantity(330.0, "K"), P=azoth.ureg.Quantity(2.5e6, "Pa"), z=[0.6, 0.4]
)
print("vapour fraction:", flash.beta)
