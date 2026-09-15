"""``eos.siddiqi_lucas_diffusivity`` - the liquid binary diffusivity from the
Siddiqi-Lucas correlation.

Spec: ``specs/calcs/eos/siddiqi_lucas_diffusivity.toml``, which records the two
solvent correlations and the viscosity floor.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SiddiqiLucasDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.siddiqi_lucas_diffusivity"


def siddiqi_lucas_diffusivity(
    form: str, VA: Q, VB: Q, T: Q, eta: Q
) -> SiddiqiLucasDiffusivityResult:
    """The binary diffusion coefficient at infinite dilution, from Siddiqi-Lucas.

    ``form`` selects the solvent correlation; ``VA`` and ``VB`` are the solute and
    solvent molar volumes and ``eta`` the solvent viscosity, floored at 0.01 cP.

    Args:
        form: ``"aqueous"`` or ``"organic"``.
        VA: solute molar volume at its normal boiling point.
        VB: solvent molar volume at its normal boiling point.
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
        >>> vb = q(8.816478555304741e-5, "m**3/mol")
        >>> eta = q(9.163064908813372e-4, "Pa*s")
        >>> r = siddiqi_lucas_diffusivity("organic", va, vb, q(298.15, "K"), eta)
        >>> round(r.d.magnitude, 16)
        1.9844756155084e-09
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "VA": input_to_si(spec, "VA", VA),
        "VB": input_to_si(spec, "VB", VB),
        "T": input_to_si(spec, "T", T),
        "eta": input_to_si(spec, "eta", eta),
    }
    apply_checks(checks.on_input, values.get, warnings)

    va_cm3 = values["VA"] * 1.0e6
    vb_cm3 = values["VB"] * 1.0e6
    eta_cp = max(0.01, values["eta"] * 1000.0)

    if form == "aqueous":
        d_cm2s = 2.98e-7 * math.pow(eta_cp, -1.026) * math.pow(va_cm3, -0.5473) * values["T"]
    else:
        d_cm2s = (
            9.89e-8
            * math.pow(eta_cp, -0.907)
            * math.pow(va_cm3, -0.45)
            * math.pow(vb_cm3, 0.265)
            * values["T"]
        )
    d = d_cm2s * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return SiddiqiLucasDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
