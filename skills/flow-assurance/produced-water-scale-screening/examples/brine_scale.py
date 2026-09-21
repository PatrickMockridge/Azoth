"""A brine's saturation ratio against halite, and how much comes out."""

import azoth

q = azoth.ureg.Quantity

BRINE = ["water", "Na+", "Cl-", "CO3--", "HCO3-"]
Z = [
    0.802568218298555,
    0.0963081861958266,
    0.0963081861958266,
    0.00321027287319422,
    0.00160513643659711,
]

pitzer = azoth.eos.pitzer_phase(BRINE, T=q(298.15, "K"), x=Z)
water = Z[0]
print(f"gamma Na+ {pitzer.gamma[1]:.6f}   gamma Cl- {pitzer.gamma[2]:.6f}")
print(f"water activity {pitzer.water_activity:.6f}")

ratio = azoth.eos.scale_saturation_ratio(
    "NaCl",
    x1=Z[1],
    x2=Z[2],
    x_water=water,
    gamma1=pitzer.gamma[1],
    gamma2=pitzer.gamma[2],
    water_activity=pitzer.water_activity,
    T=q(298.15, "K"),
    P=q(10.0, "bar"),
)
print(f"NaCl saturation ratio {ratio.saturation_ratio:.6f}")

taken = azoth.eos.salt_precipitation(BRINE, salt="NaCl", T=q(298.15, "K"), P=q(10.0, "bar"), z=Z)
print(
    f"halite taken {taken.precipitated_moles:.6f} mol/mol of feed, "
    f"extent of maximum {taken.extent_of_maximum:.4f}"
)
