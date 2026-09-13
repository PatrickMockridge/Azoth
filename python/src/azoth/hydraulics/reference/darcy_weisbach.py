"""``hydraulics.darcy_weisbach`` - pressure drop in a straight pipe.

```text
dP = f * (L / D) * (rho * v**2 / 2)
```

Spec: ``specs/calcs/hydraulics/darcy_weisbach.yaml``

# A note on the sign convention

``dp`` is a *drop* and is returned positive. The direction of flow is not
modelled: velocity is taken as a speed, and a negative input is rejected rather
than interpreted as reversed flow.

# Fittings are not included

This is straight pipe only. Fitting losses are a separate calc
(``crane_k_factors``) computed by a different method, and combining them is a
modelling decision the caller should make deliberately rather than have applied
invisibly. The ``azoth pipe`` CLI performs that composition and shows both
contributions separately.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DarcyWeisbachResult, FlowRegime
from azoth.core.units import Q, quantity, to_si
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.darcy_weisbach"


def darcy_weisbach(
    f: float,
    L: Q,
    D: Q,
    rho: Q,
    v: Q,
    mu: Q | None = None,
) -> DarcyWeisbachResult:
    """Pressure drop over a length of straight pipe.

    ``mu`` is optional, and the reason is worth stating. It is not needed for the
    pressure drop itself: ``f`` already carries the flow information. It is needed
    to *check* the flow regime, because the Reynolds number requires viscosity.

    When ``mu`` is supplied the result reports ``re`` and ``regime``, and a
    transitional flow produces a ``TRANSITIONAL_FLOW`` warning. When it is omitted
    the result reports neither, and carries a ``RANGE_CHECK_SKIPPED`` warning
    saying the regime went unchecked.

    That last part is the point. "Checked and fine" and "never checked" must not
    look the same to a caller, and an omitted optional input is exactly how those
    two states would otherwise become indistinguishable.

    Args:
        f: Darcy friction factor, from :func:`friction_factor_colebrook` or
            :func:`friction_factor_swamee_jain`.
        L: pipe length.
        D: internal pipe diameter.
        rho: fluid density.
        v: bulk mean velocity.
        mu: dynamic viscosity. Optional; see above.

    Raises:
        OutOfRangeError: if ``f``, ``L``, ``D``, ``rho`` or ``v`` violates a hard
            bound in the spec - a non-positive friction factor, diameter or
            density, a negative length or velocity.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = darcy_weisbach(
        ...     0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
        ... )
        >>> round(r.dp.magnitude, 1)
        22455.0
        >>> r.re is None
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "f": f,
        "L": to_si(L, "m", "L"),
        "D": to_si(D, "m", "D"),
        "rho": to_si(rho, "kg/m**3", "rho"),
        "v": to_si(v, "m/s", "v"),
    }

    apply_checks(checks.on_input, values.get, warnings)

    # Guarded by the checks above, so D is non-zero here.
    dp = values["f"] * (values["L"] / values["D"]) * (values["rho"] * values["v"] ** 2 / 2.0)

    # Reynolds number only when viscosity was supplied. Everything downstream of
    # this is Optional, which is what makes "unchecked" visible in the type rather
    # than only in the warnings.
    re: float | None = None
    regime: FlowRegime | None = None
    if mu is not None:
        mu_si = to_si(mu, "Pa*s", "mu")
        re = values["rho"] * values["v"] * values["D"] / mu_si
        regime = FlowRegime.from_reynolds_number(re)

    def resolve(name: str) -> float | None:
        if name == "re":
            return re
        if name == "dp":
            return dp
        return None

    # `re` resolves to None when mu was omitted, so the regime check reports
    # RANGE_CHECK_SKIPPED rather than quietly passing.
    apply_checks(checks.derived, resolve, warnings)

    # The spec's transitional band check carries TRANSITIONAL_FLOW, so it has
    # already warned if applicable. Adding a second warning for the same condition
    # would just train callers to ignore them.

    return DarcyWeisbachResult(
        dp=quantity(dp, "Pa"),
        f=f,
        re=re,
        regime=regime,
        warnings=tuple(warnings),
    )


def add_fitting_loss(dp_straight: float, k_total: float, rho: float, v: float) -> float:
    """Add a fitting loss to a straight-pipe pressure drop.

    ``dP_fittings = K * (rho * v**2 / 2)``, the velocity-head form of the same
    minor loss the equivalent-length method expresses as added pipe length.

    Provided as a named function rather than left to callers because the two terms
    are computed by different methods and combined in a specific way: the
    straight-pipe term uses ``f * L/D`` and the fitting term uses ``K``.
    Presenting the composition explicitly is what lets the CLI show both
    contributions instead of a single unexplained total.
    """
    return dp_straight + k_total * (rho * v**2 / 2.0)
