"""``eos.costald_molar_volume`` - the saturated liquid molar volume from the COSTALD
equation.

Spec: ``specs/calcs/eos/costald_molar_volume.toml``, which records the two
dimensionless volume functions and the characteristic-volume back-calculation.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CostaldMolarVolumeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.costald_molar_volume"


def _vr0(tr: float) -> float:
    tau = 1.0 - tr
    if tau <= 0.0:
        return 1.0
    t13 = math.pow(tau, 1.0 / 3.0)
    return 1.0 - 1.52816 * t13 + 1.43907 * t13 * t13 - 0.81446 * tau + 0.190454 * t13 * tau


def _vrdelta(tr: float) -> float:
    return (-0.296123 + 0.386914 * tr - 0.0427258 * tr * tr - 0.0480645 * tr**3) / (tr - 1.00001)


def costald_molar_volume(
    omega: float, Tc: Q, Vc: Q, M: Q, rho_normal: Q, T: Q
) -> CostaldMolarVolumeResult:
    """The saturated liquid molar volume of a pure component, from the COSTALD equation.

    ``Vc`` and ``rho_normal`` feed the characteristic volume: it is back-calculated
    from the normal liquid density at 288.71 K, falling back to ``Vc`` when the
    density is absent or the standard reduced temperature reaches 0.9.

    Args:
        omega: acentric factor. Dimensionless.
        Tc: critical temperature.
        Vc: critical volume, the characteristic-volume fallback.
        M: molar mass.
        rho_normal: normal liquid density at 288.71 K; zero when unknown.
        T: absolute temperature, below the critical temperature.

    Returns:
        The saturated liquid molar volume, in ``m**3/mol``.

    Raises:
        OutOfRangeError: if ``Tc``, ``Vc`` or ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> tc = q(507.6, "K")
        >>> vc = q(3.7e-4, "m**3/mol")
        >>> m = q(0.086177, "kg/mol")
        >>> rho = q(664.0, "kg/m**3")
        >>> r = costald_molar_volume(0.3013, tc, vc, m, rho, q(298.15, "K"))
        >>> round(r.v.magnitude, 18)
        0.000131500647538388
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "omega": omega,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Vc": input_to_si(spec, "Vc", Vc),
        "M": input_to_si(spec, "M", M),
        "rho_normal": input_to_si(spec, "rho_normal", rho_normal),
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    tc = values["Tc"]
    v_star = values["Vc"]
    if values["rho_normal"] > 0.0:
        tr_std = 288.71 / tc
        if tr_std < 0.9:
            factor = _vr0(tr_std) * (1.0 - omega * _vrdelta(tr_std))
            if abs(factor) > 1e-15:
                estimate = values["M"] / values["rho_normal"] / factor
                if estimate > 0.0:
                    v_star = estimate

    tr = values["T"] / tc
    v = v_star * _vr0(tr) * (1.0 - omega * _vrdelta(tr))

    apply_checks(checks.derived, lambda name: v if name == "v" else None, warnings)

    return CostaldMolarVolumeResult(v=from_si(v, "m**3/mol"), warnings=tuple(warnings))
