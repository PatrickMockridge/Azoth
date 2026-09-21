"""Screen a wet stream: where the region starts, and how much hydrate is stable inside it."""

import azoth

q = azoth.ureg.Quantity

components = ["methane", "ethane", "propane", "water"]
z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]

temperature = azoth.eos.hydrate_formation_temperature(components, P=q(100.0, "bar"), z=z, eos="srk")
pressure = azoth.eos.hydrate_formation_pressure(components, T=q(285.0, "K"), z=z, eos="srk")
print(f"at 100 bar the region starts at {temperature.temperature.to('K').magnitude:.2f} K")
print(f"at 285 K it starts at {pressure.pressure.to('bar').magnitude:.3f} bar")

formed = azoth.eos.hydrate_fraction(components, T=q(283.15, "K"), P=q(100.0, "bar"), z=z, eos="srk")
print(
    f"inside it, {formed.beta:.4f} of the feed is hydrate ({formed.structure}), "
    f"balance error {formed.balance_error:.2e}"
)
