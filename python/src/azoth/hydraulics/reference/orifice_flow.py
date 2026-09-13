"""``hydraulics.orifice_flow`` - flow through an orifice from its pressure difference.

```text
q = Cd * (math.pi * d**2 / 4) * sqrt(2 * dP / rho)
```

Spec: ``specs/calcs/hydraulics/orifice_flow.yaml``

# The discharge coefficient is an input, deliberately

``Cd`` is supplied by the caller; this calc does not compute it. ISO 5167's
discharge-coefficient equation is a long fitted expression whose constants come from
a table of experimental results, and reproducing that table is what this project's
copyright rule forbids - "just implementing the equation" being exactly how it would
happen by accident. Taking ``Cd`` as an input keeps this calc to the relation that is
genuinely public and leaves the fitted coefficient to a caller who has the standard.

It is also the shape the rest of the registry already has: ``f`` is an input to
``darcy_weisbach``, ``f_t`` to ``crane_k_factors`` and ``eta`` to ``pump_power``. A
library that computed one of those and took the others would be inconsistent about
where its own responsibility stops.

# Which ``Cd``

This calc applies the equation exactly as written and adds no velocity-of-approach
factor. ISO 5167's coefficient already includes ``1/sqrt(1 - beta**4)``; older texts
write the two separately. The caller must supply the coefficient for *this* form -
getting the convention wrong is a silent error of a few per cent rather than a
visible one.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import OrificeFlowResult
from azoth.core.units import Q, from_si, to_si
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.orifice_flow"


def orifice_flow(d: Q, dP: Q, rho: Q, Cd: float) -> OrificeFlowResult:
    """Volumetric flow through an orifice.

    Args:
        d: orifice bore diameter. Declared in millimetres, as bores are quoted in
            piping work; any length is accepted.
        dP: pressure difference across the orifice, between the tappings ``Cd`` is
            defined against. A magnitude; a negative value is rejected.
        rho: density of the fluid.
        Cd: discharge coefficient, dimensionless. See the module docstring for which
            convention this must be in.

    Raises:
        OutOfRangeError: if ``d`` or ``rho`` is not positive, if ``dP`` is negative,
            or if ``Cd`` is outside ``(0, 1]``.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = orifice_flow(
        ...     q(50.0, "mm"), q(25000.0, "Pa"), q(998.0, "kg/m**3"), 0.62
        ... )
        >>> round(r.q.magnitude, 7)
        0.0086167
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        # Declared in mm, converted to the SI base metre here - so a 50 mm bore
        # reaches the arithmetic as 0.05.
        "d": to_si(d, "mm", "d"),
        "dP": to_si(dP, "Pa", "dP"),
        "rho": to_si(rho, "kg/m**3", "rho"),
        # Dimensionless, so it arrives as a plain float with no unit to convert.
        "Cd": Cd,
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Guarded by the checks above, so d and rho are positive and dP is not negative.
    area = math.pi * values["d"] ** 2 / 4.0
    q = values["Cd"] * area * math.sqrt(2.0 * values["dP"] / values["rho"])

    apply_checks(checks.derived, lambda name: q if name == "q" else None, warnings)

    return OrificeFlowResult(q=from_si(q, "m**3/s"), warnings=tuple(warnings))
