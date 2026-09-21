"""The hydrate equilibrium temperature at a pressure, and the margin to it."""

import azoth

q = azoth.ureg.Quantity

components = ["methane", "ethane", "propane", "water"]
z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]

for pressure in (100.0, 50.0):
    equilibrium = azoth.eos.hydrate_formation_temperature(
        components, P=q(pressure, "bar"), z=z, eos="srk"
    )
    margin = q(285.0, "K") - equilibrium.temperature
    print(
        f"{pressure:5.0f} bar: {equilibrium.temperature.to('K').magnitude:.2f} K "
        f"({equilibrium.structure}), margin at 285 K "
        f"{margin.to('delta_degC').magnitude:.2f} K"
    )
