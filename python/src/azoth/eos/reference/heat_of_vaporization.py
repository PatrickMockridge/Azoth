"""``eos.heat_of_vaporization`` - NeqSim's pure-component heat-of-vaporisation
correlation.

```text
hov = 1e-3*c0*(1 - Tr)**(c1 + c2*Tr + c3*Tr**2)
```

Spec: ``specs/calcs/eos/heat_of_vaporization.toml``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HeatOfVaporizationResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.heat_of_vaporization"


def heat_of_vaporization(
    c0: float, c1: float, c2: float, c3: float, Tc: Q, T: Q
) -> HeatOfVaporizationResult:
    """The pure-component heat of vaporisation at a temperature, from NeqSim's correlation.

    ``c0..c3`` are the raw ``heatofvaporizationcoefs1..4`` and ``Tc`` the critical
    temperature, all caller-supplied.

    Args:
        c0: the first coefficient, NeqSim's internal scale.
        c1, c2, c3: the exponent's constant, linear and quadratic terms.
        Tc: critical temperature.
        T: absolute temperature.

    Returns:
        The heat of vaporisation, in ``J/mol``.

    Raises:
        OutOfRangeError: if ``Tc`` or ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> args = (52_100_000.0, 0.32, -0.212, 0.258)
        >>> r = heat_of_vaporization(*args, q(425.12, "K"), q(300.0, "K"))
        >>> round(r.hov.magnitude, 9)
        36147.577791
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "c0": c0,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "Tc": input_to_si(spec, "Tc", Tc),
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    tr = values["T"] / values["Tc"]
    exponent = c1 + c2 * tr + c3 * tr * tr
    hov = 1e-3 * c0 * (1.0 - tr) ** exponent

    apply_checks(checks.derived, lambda name: hov if name == "hov" else None, warnings)

    return HeatOfVaporizationResult(hov=from_si(hov, "J/mol"), warnings=tuple(warnings))
