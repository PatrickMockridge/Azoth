"""``eos.iapws_saturation_pressure`` - the IAPWS-IF97 Region 4 saturation pressure of water.

```text
p_sat = (2C / (-B + sqrt(B**2 - 4*A*C)))**4,  theta = T + n9/(T - n10)
```

Spec: ``specs/calcs/eos/iapws_saturation_pressure.toml``. Oracle:
``validation/neqsim/captures/iapws_if97_probe.tsv``.

**It is the boundary :func:`azoth.eos.reference.water_phase.water_phase` selects its region
by**, and the release states the two directions as two fits - ``p4_T`` forward and ``T4_p``
backward are not one equation inverted - so the inverse is a second id's question and not this
one's.

**A different answer from ``eos.pure_saturation`` to a question that sounds the same.** That one
is the pressure where a Peng-Robinson fluid's two roots cross, which is a property of the cubic;
this is the reference curve for water.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import IapwsSaturationPressureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.iapws_saturation_pressure"


def saturation_pressure_mpa(t: float) -> float:
    """The saturation pressure in MPa at a temperature in K, the Region 4 forward equation.

    Megapascals because that is the release's own unit and the one ``iapws_if97::
    saturation_pressure`` carries; the model converts once, at its own boundary.
    """
    theta = t - 0.23855557567849 / (t - 650.17534844798)
    a = theta * theta + 1167.0521452767 * theta - 724213.16703206
    b = -17.073846940092 * theta * theta + 12020.82470247 * theta - 3232555.0322333
    c = 14.91510861353 * theta * theta - 4823.2657361591 * theta + 405113.40542057
    return (2.0 * c / (-b + math.sqrt(b * b - 4.0 * a * c))) ** 4


def iapws_saturation_pressure(T: Q) -> IapwsSaturationPressureResult:
    """The saturation pressure of water at a temperature, by the Region 4 equation.

    Args:
        T: absolute temperature, inside the curve's own domain ``[273.15, 647.096)`` K.

    Returns:
        The pressure at which liquid water and its vapour are in equilibrium at ``T``.

    Raises:
        OutOfRangeError: if ``T`` is outside ``[273.15, 647.096)`` K - the triple point to the
            critical point, which is the curve's own domain.

    Example:
        >>> from azoth.core.units import ureg
        >>> r = iapws_saturation_pressure(ureg.Quantity(373.15, "K"))
        >>> round(r.p_sat.to("Pa").magnitude, 3)
        101417.978
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = float(input_to_si(spec, "T", T))

    apply_checks(checks.on_input, {"T": t_si}.get, warnings)

    # The equation returns the release's own unit, megapascals; the model is SI.
    p_sat = saturation_pressure_mpa(t_si) * 1.0e6

    return IapwsSaturationPressureResult(
        p_sat=from_si(p_sat, "Pa"),
        warnings=tuple(warnings),
    )
