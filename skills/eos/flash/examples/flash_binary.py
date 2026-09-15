"""Flash a binary and test its stability."""

import azoth

q = azoth.ureg.Quantity
fluid = azoth.eos.from_names(["methane", "n-butane"])

r = azoth.eos.pt_flash(fluid, T=q(330.0, "K"), P=q(2.5e6, "Pa"), z=[0.6, 0.4])
print("vapour fraction:", r.beta)
print("phase:", r.phase)

stab = azoth.eos.stability_test(fluid, T=q(330.0, "K"), P=q(2.5e6, "Pa"), z=[0.6, 0.4])
print("stability verdict:", stab.verdict)
