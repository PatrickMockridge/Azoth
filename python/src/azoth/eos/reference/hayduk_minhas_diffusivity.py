"""``eos.hayduk_minhas_diffusivity`` - the liquid binary diffusivity from the
Hayduk-Minhas correlation.

Spec: ``specs/calcs/eos/hayduk_minhas_diffusivity.toml``, which records the two
solvent correlations and the clamps.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HaydukMinhasDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.hayduk_minhas_diffusivity"


def hayduk_minhas_diffusivity(form: str, VA: Q, T: Q, eta: Q) -> HaydukMinhasDiffusivityResult:
    """The binary diffusion coefficient at infinite dilution, from Hayduk-Minhas.

    ``form`` selects the solvent correlation; ``VA`` is the solute molar volume at
    the normal boiling point and ``eta`` the solvent viscosity. ``VA`` and ``eta``
    are clamped to NeqSim's `[20, 500]` cm**3/mol and `[0.1, 100]` cP before the
    correlation.

    Args:
        form: ``"paraffin"`` or ``"aqueous"``.
        VA: solute molar volume at its normal boiling point.
        T: absolute temperature.
        eta: solvent dynamic viscosity at ``T``.

    Returns:
        The binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> va = q(4.0203262233375156e-5, "m**3/mol")
        >>> eta = q(9.163064908813372e-4, "Pa*s")
        >>> r = hayduk_minhas_diffusivity("paraffin", va, q(298.15, "K"), eta)
        >>> round(r.d.magnitude, 16)
        4.3917780890442e-09
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "VA": input_to_si(spec, "VA", VA),
        "T": input_to_si(spec, "T", T),
        "eta": input_to_si(spec, "eta", eta),
    }
    apply_checks(checks.on_input, values.get, warnings)

    va_cm3 = min(500.0, max(20.0, values["VA"] * 1.0e6))
    eta_cp = min(100.0, max(0.1, values["eta"] * 1000.0))

    if form == "paraffin":
        d_cm2s = (
            13.3e-8
            * math.pow(values["T"], 1.47)
            * math.pow(eta_cp, 10.2 / va_cm3 - 0.791)
            / math.pow(va_cm3, 0.71)
        )
    else:
        vi_term = max(0.01, math.pow(va_cm3, -0.19) - 0.292)
        d_cm2s = (
            1.25e-8 * vi_term * math.pow(values["T"], 1.52) * math.pow(eta_cp, 9.58 / va_cm3 - 1.12)
        )
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return HaydukMinhasDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
