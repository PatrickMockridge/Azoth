"""``hydraulics.choked_flow_area`` - the throat area a choked gas flow needs.

```text
A = m_dot / (sqrt(k * rho0 * P0) * (2 / (k + 1))**((k + 1) / (2 * (k - 1))))
```

Spec: ``specs/calcs/hydraulics/choked_flow_area.toml``

# What this is, and what it deliberately is not

This is the isentropic critical-flow relation: the mass flux through a throat when
the downstream pressure is low enough that the flow reaches sonic velocity there.

It is the physical basis of relief valve sizing, and it is **not** relief valve
sizing to a standard. API 520 wraps this relation in de-rating coefficients - a
discharge coefficient, a back-pressure correction, a combination factor - whose
values are tabulated in the standard, and this library does not reproduce tables. A
caller sizing a relief valve applies those factors themselves, visibly, exactly as
the ``azoth pipe`` CLI composes fitting losses rather than a calc baking them in.
That is why the calculation is named for what it computes; see the spec.

# Why the units work out

The whole unit structure lives in ``sqrt(k rho0 P0)``, which is a mass flux:
``sqrt(Pa * kg/m**3)`` is ``kg/(m**2 s)``. Everything else is a function of ``k``
alone and is dimensionless. A formula whose units do not resolve to an area is
wrong, and this one does.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ChokedFlowAreaResult
from azoth.core.units import Q, from_si, to_si
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.choked_flow_area"


def choked_flow_area(m_dot: Q, P0: Q, rho0: Q, k: float) -> ChokedFlowAreaResult:
    """Throat area required for a choked mass flow of an ideal gas.

    Args:
        m_dot: mass flow rate the throat must pass.
        P0: stagnation (upstream total) pressure. The relation assumes the approach
            velocity is negligible, which is the usual relief-valve assumption.
        rho0: stagnation density, taken directly rather than derived from a gas
            constant and a molecular weight - form it with whichever equation of
            state is appropriate, which for a real gas is where the compressibility
            factor belongs.
        k: isentropic exponent, the ratio of specific heats. Dimensionless, and
            above 1 for every real gas.

    Raises:
        OutOfRangeError: if ``P0`` or ``rho0`` is not positive, if ``m_dot`` is
            negative, or if ``k`` is not greater than 1.

    Outside the range of real substances - ``k`` above 5/3, the monatomic limit - the
    area is still returned, carrying an ``OUT_OF_VALID_RANGE`` warning.

    The flow must be choked for this area to be the right one, and this calc cannot
    check that: the critical pressure ratio depends on ``k``, which the spec's typed
    range checks cannot express. See the spec's assumptions.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = choked_flow_area(
        ...     q(1.0, "kg/s"), q(1.0e6, "Pa"), q(10.0, "kg/m**3"), 1.4
        ... )
        >>> round(r.a.magnitude * 1e6, 1)   # in square millimetres
        461.8
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "m_dot": to_si(m_dot, "kg/s", "m_dot"),
        "P0": to_si(P0, "Pa", "P0"),
        "rho0": to_si(rho0, "kg/m**3", "rho0"),
        # Dimensionless, so it arrives as a plain float with no unit to convert.
        "k": k,
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Guarded by the checks above, so k > 1 and P0, rho0 > 0.
    #
    # The geometric factor is written as the derived expression rather than as a
    # named constant, because it depends on `k` and so cannot be one. It is the
    # `(2/(k+1))**((k+1)/(2(k-1)))` of the published form, and at k = 1.4 it is
    # exactly (5/6)**3 = 125/216 - which is what makes the spec's worked example
    # checkable by hand.
    geometric = (2.0 / (values["k"] + 1.0)) ** ((values["k"] + 1.0) / (2.0 * (values["k"] - 1.0)))
    flux = math.sqrt(values["k"] * values["rho0"] * values["P0"]) * geometric
    a = values["m_dot"] / flux

    apply_checks(checks.derived, lambda name: a if name == "a" else None, warnings)

    return ChokedFlowAreaResult(a=from_si(a, "m**2"), warnings=tuple(warnings))
