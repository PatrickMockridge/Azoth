"""Find the wax appearance temperature by bisection, and the margin above it."""

import azoth

q = azoth.ureg.Quantity

COMPONENTS = ["methane", "n-heptane", "nc14", "nc20"]
Z = [0.7, 0.1, 0.1, 0.1]
P = q(5.0, "bar")


def wax_at(temperature: float) -> float:
    """The stable wax fraction at a temperature, zero where the two-phase answer stands."""
    result = azoth.eos.tp_multiflash_wax(COMPONENTS, T=q(temperature, "K"), P=P, z=Z, eos="srk")
    return result.wax_fraction if result.converged else 0.0


low, high = 255.0, 275.0  # wax at the bottom, none at the top
for _ in range(17):
    middle = 0.5 * (low + high)
    if wax_at(middle) > 0.0:
        low = middle
    else:
        high = middle

appearance = 0.5 * (low + high)
operating = 268.0
print(f"wax appearance temperature {appearance:.4f} K at 5 bar")
print(f"margin at {operating} K: {operating - appearance:.4f} K")
print(f"wax stable there: {wax_at(operating):.6f} of the feed")
