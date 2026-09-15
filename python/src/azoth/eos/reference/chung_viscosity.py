"""``eos.chung_viscosity`` - the gas viscosity from the Chung correlation.

Spec: ``specs/calcs/eos/chung_viscosity.toml``, which records the ten coefficient
rows, the collision integral, and the dense-gas correction.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ChungViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.chung_viscosity"

#: The ten rows of Chung's high-pressure coefficient table (TPoLG Table 9-5).
_CHUNG_HP = (
    (6.324, 50.412, -51.680, 1189.0),
    (1.210e-3, -1.154e-3, -6.257e-3, 0.03728),
    (5.283, 254.209, -168.48, 3898.0),
    (6.623, 38.096, -8.464, 31.42),
    (19.745, 7.630, -14.354, 31.53),
    (-1.9, -12.537, 4.985, -18.15),
    (24.275, 3.450, -11.291, 69.35),
    (0.7972, 1.117, 0.01235, -4.117),
    (-0.2382, 0.06770, -0.8163, 4.025),
    (0.06863, 0.3479, 0.5926, -0.727),
)


def chung_viscosity(
    omega: float,
    Tc: Q,
    Vc: Q,
    M: Q,
    dipole: float,
    kappa: float,
    T: Q,
    V: Q,
) -> ChungViscosityResult:
    """The gas dynamic viscosity of a pure component, from Chung's correlation.

    ``dipole`` is in debye and ``kappa`` the dimensionless viscosity correction
    factor; ``V`` is the gas's molar volume, which the dense-gas correction
    ``y = Vc/(6*V)`` consumes.

    Args:
        omega: acentric factor. Dimensionless.
        Tc: critical temperature.
        Vc: critical volume.
        M: molar mass.
        dipole: dipole moment in debye; zero for a nonpolar component.
        kappa: the viscosity correction factor.
        T: absolute temperature.
        V: molar volume of the gas at ``T``.

    Returns:
        The gas dynamic viscosity, in ``Pa*s``.

    Raises:
        OutOfRangeError: if ``Tc``, ``Vc``, ``T`` or ``V`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> tc = q(190.56, "K")
        >>> vc = q(9.9e-5, "m**3/mol")
        >>> m = q(0.016043, "kg/mol")
        >>> v = q(2.4409707154781444e-3, "m**3/mol")
        >>> r = chung_viscosity(0.0115, tc, vc, m, 0.0, 0.0, q(300.0, "K"), v)
        >>> round(r.mu.magnitude, 16)
        1.1240094154301e-05
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "omega": omega,
        "Tc": input_to_si(spec, "Tc", Tc),
        "Vc": input_to_si(spec, "Vc", Vc),
        "M": input_to_si(spec, "M", M),
        "dipole": dipole,
        "kappa": kappa,
        "T": input_to_si(spec, "T", T),
        "V": input_to_si(spec, "V", V),
    }

    apply_checks(checks.on_input, values.get, warnings)

    tc = values["Tc"]
    vc_cm3 = values["Vc"] * 1.0e6
    m_g = values["M"] * 1000.0

    rel_visc = 131.3 * dipole / math.sqrt(vc_cm3 * tc)
    rel4 = rel_visc**4
    fc = 1.0 - 0.2756 * omega + 0.059035 * rel4 + kappa
    e = [row[0] + row[1] * omega + row[2] * rel4 + row[3] * kappa for row in _CHUNG_HP]

    temp_var = 1.2593 * values["T"] / tc
    omega_visc = (
        1.16145 * math.pow(temp_var, -0.14874)
        + 0.52487 * math.exp(-0.77320 * temp_var)
        + 2.16178 * math.exp(-2.43787 * temp_var)
    )

    chungy = values["Vc"] / (6.0 * values["V"])
    g1 = (1.0 - 0.5 * chungy) / (1.0 - chungy) ** 3
    g2 = (
        e[0] * ((1.0 - math.exp(-e[3] * chungy)) / chungy)
        + e[1] * g1 * math.exp(e[4] * chungy)
        + e[2] * g1
    ) / (e[0] * e[3] + e[1] + e[2])
    viskstarstar = (
        e[6]
        * chungy
        * chungy
        * g2
        * math.exp(e[7] + e[8] / temp_var + e[9] * math.pow(temp_var, -2.0))
    )
    viskstar = math.sqrt(temp_var) / omega_visc * (fc * (1.0 / g2 + e[5] * chungy)) + viskstarstar

    mu_micropoise = viskstar * 36.344 * math.sqrt(m_g * tc) / math.pow(vc_cm3, 2.0 / 3.0)
    mu = mu_micropoise * 1.0e-7

    apply_checks(checks.derived, lambda name: mu if name == "mu" else None, warnings)

    return ChungViscosityResult(mu=from_si(mu, "Pa*s"), warnings=tuple(warnings))
