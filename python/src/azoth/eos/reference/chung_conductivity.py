"""``eos.chung_conductivity`` - the gas thermal conductivity from the Chung
correlation.

Spec: ``specs/calcs/eos/chung_conductivity.toml``, which records the dilute-gas
viscosity, the Eucken correction and the gas constant.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ChungConductivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.chung_conductivity"

#: NeqSim's gas constant, not the exact SI value.
_R = 8.3144621


def chung_conductivity(
    Cv0: Q,
    M: Q,
    omega: float,
    Tc: Q,
    Vc: Q,
    dipole: float,
    kappa: float,
    T: Q,
) -> ChungConductivityResult:
    """The gas thermal conductivity of a pure component, from Chung's correlation.

    ``Cv0`` is the ideal-gas heat capacity at constant volume at ``T``; ``dipole``
    is in debye and ``kappa`` the dimensionless viscosity correction factor. The
    viscosity inside the correlation is the *dilute-gas* Chung viscosity, not
    ``eos.chung_viscosity``'s dense-gas form.

    Args:
        Cv0: ideal-gas heat capacity at constant volume at ``T``; ``Cp0 - R``.
        M: molar mass.
        omega: acentric factor.
        Tc: critical temperature.
        Vc: critical volume.
        dipole: dipole moment in debye; zero for a nonpolar component.
        kappa: the viscosity correction factor.
        T: absolute temperature.

    Returns:
        The gas thermal conductivity, in ``W/(m*K)``.

    Raises:
        OutOfRangeError: if ``Tc``, ``Vc`` or ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> cv = q(27.544151394, "J/(mol*K)")
        >>> m = q(0.016043, "kg/mol")
        >>> tc = q(190.56, "K")
        >>> vc = q(9.9e-5, "m**3/mol")
        >>> r = chung_conductivity(cv, m, 0.0115, tc, vc, 0.0, 0.0, q(300.0, "K"))
        >>> round(r.k.magnitude, 16)
        0.03384023073765177
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "Cv0": input_to_si(spec, "Cv0", Cv0),
        "M": input_to_si(spec, "M", M),
        "omega": omega,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Vc": input_to_si(spec, "Vc", Vc),
        "dipole": dipole,
        "kappa": kappa,
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    vc_cm3 = values["Vc"] * 1.0e6
    m_g = values["M"] * 1000.0

    rel_visc = 131.3 * dipole / math.sqrt(vc_cm3 * values["Tc"])
    rel4 = rel_visc**4
    fc = 1.0 - 0.2756 * omega + 0.059035 * rel4 + kappa

    t_star = 1.2593 * values["T"] / values["Tc"]
    omega_v = (
        1.16145 * math.pow(t_star, -0.14874)
        + 0.52487 * math.exp(-0.77320 * t_star)
        + 2.16178 * math.exp(-2.43787 * t_star)
    )

    eta0 = (
        40.785 * fc * math.sqrt(m_g * values["T"]) / (math.pow(vc_cm3, 2.0 / 3.0) * omega_v)
    ) * 1.0e-7

    alpha = values["Cv0"] / _R - 1.5
    beta = 0.7862 - 0.7109 * omega + 1.3168 * omega * omega
    z = 2.0 + 10.5 * (values["T"] / values["Tc"]) ** 2
    psi = 1.0 + alpha * (
        (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
        / (0.6366 + beta * z + 1.061 * alpha * beta)
    )

    k = 3.75 * _R * psi * eta0 / values["M"]

    apply_checks(checks.derived, lambda name: k if name == "k" else None, warnings)

    return ChungConductivityResult(k=from_si(k, "W/(m*K)"), warnings=tuple(warnings))
