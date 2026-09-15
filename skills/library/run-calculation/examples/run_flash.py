"""Run a two-phase flash and look a component up by name.

Copyable as-is for the flash: the mixture is built from the databank by name, so
nothing is read from a file. The keycard step reads `keycard.example.toml`, so run
it from the repository root.
"""

import azoth

q = azoth.ureg.Quantity

fluid = azoth.eos.from_names(["methane", "n-butane"])

flash = azoth.eos.pt_flash(fluid, T=q(330.0, "K"), P=q(2.5e6, "Pa"), z=[0.6, 0.4])
print("vapour fraction:", flash.beta)
print("phase:", flash.phase)

methane = azoth.eos.component("methane")
print("methane Tc:", methane.Tc)

# A keycard overrides or extends what the library ships. `load` reads a file and
# returns a value; it stores nothing.
card = azoth.keycard.load("keycard.example.toml")
overridden = azoth.eos.component("methane", card=card)
print("methane Tc (card):", overridden.Tc)
