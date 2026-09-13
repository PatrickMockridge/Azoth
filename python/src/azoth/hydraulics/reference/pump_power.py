"""``hydraulics.pump_power`` - shaft power from flow, head and efficiency.

```text
power = rho * g * q * H / eta
```

Spec: ``specs/calcs/hydraulics/pump_power.yaml``

# ``g`` is a constant here, not an input

Standard gravity, ``9.80665 m/s**2``, is a *defined* value - the CGPM fixed it
exactly in 1901 - so it belongs to the equation in the way pi does, rather than
being a measured property assumed behind the caller's back. The spec records the
alternative of taking it as an input and why that was rejected: a caller who must
supply ``g`` can supply a wrong one, and a pump efficiency supplied from a
manufacturer's curve carries far more uncertainty than the 0.3% by which local
gravity varies.

# What ``H`` is, and what it is not

``H`` is the head the pump *delivers* - the useful head - in metres of the pumped
fluid. It is not the head the impeller generates before the pump's own internal
losses, because those losses are exactly what ``eta`` accounts for; counting them
in ``H`` as well would apply them twice and produce a power that is too high.
"""

from __future__ import annotations

from typing import Final

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PumpPowerResult
from azoth.core.units import Q, from_si, to_si
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.pump_power"

#: Standard gravity, in metres per second squared.
#:
#: Exact by definition: the CGPM fixed it in 1901, which is what makes it a
#: constant of the equation rather than an input. Named rather than written inline
#: so it is greppable and so the Rust side can be compared against it.
STANDARD_GRAVITY_M_S2: Final[float] = 9.80665


def pump_power(rho: Q, q: Q, H: Q, eta: float) -> PumpPowerResult:
    """Shaft power a pump must be supplied with.

    Args:
        rho: density of the pumped fluid.
        q: volumetric flow rate. A magnitude; a negative value is rejected.
        H: head developed, in metres of the pumped fluid. The *delivered* head.
        eta: overall pump efficiency, dimensionless, in ``(0, 1]``.

    Raises:
        OutOfRangeError: if ``rho`` is not positive, if ``q`` or ``H`` is negative,
            or if ``eta`` is outside ``(0, 1]``.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pump_power(
        ...     q(998.0, "kg/m**3"), q(0.01, "m**3/s"), q(30.0, "m"), 0.75
        ... )
        >>> round(r.power.magnitude, 3)
        3914.815
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "rho": to_si(rho, "kg/m**3", "rho"),
        "q": to_si(q, "m**3/s", "q"),
        "H": to_si(H, "m", "H"),
        # Dimensionless, so it arrives as a plain float with no unit to convert -
        # the same treatment colebrook's `re` and `relative_roughness` get.
        "eta": eta,
    }

    apply_checks(checks.on_input, values.get, warnings)

    power = values["rho"] * STANDARD_GRAVITY_M_S2 * values["q"] * values["H"] / values["eta"]

    apply_checks(checks.derived, lambda name: power if name == "power" else None, warnings)

    return PumpPowerResult(power=from_si(power, "W"), warnings=tuple(warnings))
