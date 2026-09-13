"""``hydraulics.reynolds_number`` - Reynolds number and flow regime.

```text
Re = rho * v * D / mu
```

Spec: ``specs/calcs/hydraulics/reynolds_number.yaml``

Every other calc in this slice consumes the Reynolds number, and two of them
consume the regime as well. The regime boundaries are the Crane/Moody convention
and are approximate: real transition depends on inlet geometry, vibration and
roughness, which is why the transitional band produces a warning rather than a
clean label.
"""

from __future__ import annotations

from chemeng._registry_gen import spec as _spec_for
from chemeng.core.range import apply_checks, checks_for
from chemeng.core.result import FlowRegime, ReynoldsNumberResult
from chemeng.core.units import Q, to_si
from chemeng.core.warnings import Warning

CALC_ID = "hydraulics.reynolds_number"


def reynolds_number(rho: Q, v: Q, D: Q, mu: Q) -> ReynoldsNumberResult:
    """Reynolds number for flow in a circular pipe.

    Returns the Reynolds number and its :class:`FlowRegime`. When the flow is
    transitional the result carries a ``TRANSITIONAL_FLOW`` warning, because the
    friction factor - and therefore any pressure drop derived from it - is
    indeterminate in that band rather than merely uncertain.

    Args:
        rho: fluid density.
        v: bulk mean velocity.
        D: internal pipe diameter.
        mu: dynamic viscosity.

    Raises:
        OutOfRangeError: if density, diameter or viscosity is not positive, or if
            velocity is negative. Those make the Reynolds number undefined or
            meaningless, as opposed to merely out of range.

    Example:
        >>> import chemeng
        >>> q = chemeng.ureg.Quantity
        >>> r = reynolds_number(
        ...     q(998.0, "kg/m**3"), q(1.5, "m/s"), q(0.1, "m"), q(1.002e-3, "Pa*s")
        ... )
        >>> round(r.re, 1)
        149401.2
        >>> str(r.regime)
        'turbulent'
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "rho": to_si(rho, "kg/m**3", "rho"),
        "v": to_si(v, "m/s", "v"),
        "D": to_si(D, "m", "D"),
        "mu": to_si(mu, "Pa*s", "mu"),
    }

    # Input checks first: they guard the arithmetic below, and their severity is
    # `error`, so a violation means we must not divide at all.
    apply_checks(checks.on_input, values.get, warnings)

    re = values["rho"] * values["v"] * values["D"] / values["mu"]

    # `re` is an output here, so the transitional band is detected at this point.
    # The spec names that check's warning code as TRANSITIONAL_FLOW, so a caller
    # gets the specific code rather than a generic out-of-range.
    apply_checks(checks.derived, lambda name: re if name == "re" else None, warnings)

    return ReynoldsNumberResult(
        re=re,
        regime=FlowRegime.from_reynolds_number(re),
        warnings=tuple(warnings),
    )
