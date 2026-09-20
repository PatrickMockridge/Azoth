"""``eos.antoine_vapor_pressure`` - NeqSim's pure-component vapour-pressure correlation.

NeqSim's ``getAntoineVaporPressure`` dispatches on a string label and evaluates one of
four correlations. The label vocabulary in the vendored data is messy - ``exp`` and
``log`` are one formula under two names, and ``loglog``/``log10`` have no branch and
fall through to Wagner - so the label is cleaned onto a closed form before the
correlation is evaluated.

Spec: ``specs/calcs/eos/antoine_vapor_pressure.toml``.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AntoineVaporPressureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.antoine_vapor_pressure"


def antoine_vapor_pressure(
    A: float, B: float, C: float, D: float, E: float, form: str, Tc: Q, Pc: Q, T: Q
) -> AntoineVaporPressureResult:
    """The pure-component vapour pressure at a temperature, from NeqSim's correlation.

    ``A``-``E`` are the raw ``ANTOINEA``-``ANTOINEE``, ``Tc``/``Pc`` the critical
    constants, all caller-supplied. ``E`` is the DIPPR-101 exponent, used by
    ``dippr101`` and ignored by the others; ``Tc`` and ``Pc`` are used only by the
    Wagner form.

    Args:
        A..E: the five Antoine coefficients, NeqSim's internal scale.
        form: one of ``"dippr101"``, ``"pow10"``, ``"pow10kpa"``, ``"exp"`` or
            ``"wagner"``.
        Tc: critical temperature.
        Pc: critical pressure.
        T: absolute temperature.

    Returns:
        The vapour pressure, in ``Pa``.

    Raises:
        OutOfRangeError: if ``T``, ``Tc`` or ``Pc`` is not positive, or if the Wagner
            form is evaluated above ``Tc``, where it is not real.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> coeffs = (5.23243, 891.0098, 332.0975, 0.0, 0.0)
        >>> tc, pc, t = q(190.56, "K"), q(4_599_000.0, "Pa"), q(300.0, "K")
        >>> r = antoine_vapor_pressure(*coeffs, "pow10", tc, pc, t)
        >>> round(r.p_sat.magnitude, 6)
        56252981.195393
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Pc": input_to_si(spec, "Pc", Pc),
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    t = values["T"]
    if form == "dippr101":
        # NeqSim returns this one in pascals already - ``exp(...) / 100000`` bar - so
        # unlike the other bar-returning forms there is no ``1e5`` factor here.
        p_sat = math.exp(A + B / t + C * math.log(t) + D * t**E)
    elif form == "pow10":
        p_sat = 1e5 * 10.0 ** (A - B / (t + C - 273.15))
    elif form == "pow10kpa":
        p_sat = 10.0 ** (A - B / (t + C))
    elif form == "exp":
        p_sat = 1e5 * math.exp(A - B / (t + C))
    else:
        # **The form is defined on `0 < T <= Tc`, and the kernel says so.** Above the
        # critical temperature `x = 1 - T/Tc` is negative and `x**1.5` is not real: Rust's
        # `powf` returns `NaN` there and Python's `**` returns a *complex*, so the two
        # kernels would disagree in kind rather than in value. The state is refused
        # instead, because there is no saturation pressure above `Tc` to report.
        x = 1.0 - t / values["Tc"]
        if x < 0.0:
            raise OutOfRangeError(
                "T",
                t,
                f"the Wagner form is defined up to the critical temperature, and this state "
                f"is above it: `1 - T/Tc` is {x}, whose 1.5 power is not a real number. "
                f"There is no saturation pressure above `Tc` to return, so the state is "
                f"refused rather than reported as `NaN`",
            )
        p_sat = math.exp((A * x + B * x**1.5 + C * x**3 + D * x**6) / (1.0 - x)) * values["Pc"]

    apply_checks(checks.derived, lambda name: p_sat if name == "p_sat" else None, warnings)

    return AntoineVaporPressureResult(p_sat=from_si(p_sat, "Pa"), warnings=tuple(warnings))
