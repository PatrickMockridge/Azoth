"""``eos.friction_theory_conductivity`` - NeqSim's ``FrictionTheoryConductivityMethod``.

A Chung dilute-gas conductivity mixed by Mason-Saxena, plus a dense term that converts
the phase's viscosity excess with the Eucken factor.

Spec: ``specs/models/eos/friction_theory_conductivity.toml``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FrictionTheoryConductivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.chung_conductivity import _R

CALC_ID = "eos.friction_theory_conductivity"


def _chung_dilute(
    cv0: float,
    m_kg: float,
    omega: float,
    tc: float,
    vc: float,
    dipole: float,
    temperature: float,
) -> tuple[float, float]:
    """Chung's dilute-gas conductivity and viscosity, as ``eos.chung_conductivity`` has them."""
    vc_cm3 = vc * 1.0e6
    m_g = m_kg * 1000.0

    rel_visc = 131.3 * dipole / math.sqrt(vc_cm3 * tc)
    fc = 1.0 - 0.2756 * omega + 0.059035 * rel_visc**4

    t_star = 1.2593 * temperature / tc
    omega_v = (
        1.16145 * math.pow(t_star, -0.14874)
        + 0.52487 * math.exp(-0.77320 * t_star)
        + 2.16178 * math.exp(-2.43787 * t_star)
    )
    eta0 = (
        40.785
        * fc
        * math.sqrt(m_g * temperature)
        / (math.pow(vc_cm3, 2.0 / 3.0) * omega_v)
        * 1.0e-7
    )

    alpha = cv0 / _R - 1.5
    beta = 0.7862 - 0.7109 * omega + 1.3168 * omega * omega
    z = 2.0 + 10.5 * (temperature / tc) ** 2
    psi = 1.0 + alpha * (
        (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
        / (0.6366 + beta * z + 1.061 * alpha * beta)
    )
    return 3.75 * _R * psi * eta0 / m_kg, eta0


def friction_theory_conductivity(
    cv0: Sequence[Q],
    tc: Sequence[Q],
    critical_volume: Sequence[Q],
    omega: Sequence[float],
    dipole: Sequence[float],
    molar_mass: Sequence[Q],
    x: Sequence[float],
    T: Q,
    viscosity: Q,
) -> FrictionTheoryConductivityResult:
    """The dense-fluid thermal conductivity from a phase's own viscosity.

    Args:
        cv0: each component's ideal-gas heat capacity at constant volume; ``Cp0 - R``.
        tc: each component's critical temperature.
        critical_volume: each component's critical volume.
        omega: each component's acentric factor.
        dipole: each component's dipole moment, in debye.
        molar_mass: each component's molar mass.
        x: the phase's mole fractions.
        T: absolute temperature.
        viscosity: the phase's total dynamic viscosity, from whatever model it carries.

    Returns:
        The conductivity and its two terms, in ``W/(m*K)``.

    Raises:
        InvalidInputError: if the per-component vectors disagree in length.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = friction_theory_conductivity(
        ...     [q(27.544151394, "J/(mol*K)")], [q(190.56, "K")], [q(9.9e-5, "m**3/mol")],
        ...     [0.0115], [0.0], [q(0.016043, "kg/mol")], [1.0], q(300.0, "K"),
        ...     q(1.1403910448332615e-5, "Pa*s"))
        >>> abs(r.conductivity.magnitude - 0.034385915611225566) < 1e-14
        True
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(x)
    if (
        len(cv0) != count
        or len(tc) != count
        or len(critical_volume) != count
        or len(omega) != count
        or len(dipole) != count
        or len(molar_mass) != count
    ):
        raise InvalidInputError(
            "x",
            f"{count} mole fraction(s) against {len(cv0)} Cv0, {len(tc)} Tc, "
            f"{len(critical_volume)} Vc, {len(omega)} omega, {len(dipole)} dipole and "
            f"{len(molar_mass)} molar mass value(s): every one is one entry per component",
        )

    cv0_values = [input_to_si(spec, "cv0", value) for value in cv0]
    tc_values = [input_to_si(spec, "tc", value) for value in tc]
    vc_values = [input_to_si(spec, "critical_volume", value) for value in critical_volume]
    mass_values = [input_to_si(spec, "molar_mass", value) for value in molar_mass]
    temperature = input_to_si(spec, "T", T)
    viscosity_value = input_to_si(spec, "viscosity", viscosity)
    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "viscosity": viscosity_value,
            "critical_volume": min(vc_values) if vc_values else None,
            "molar_mass": min(mass_values) if mass_values else None,
            "x": min(float(value) for value in x) if x else None,
        }.get,
        warnings,
    )

    pure_conductivity: list[float] = []
    pure_viscosity: list[float] = []
    for i in range(count):
        conductivity, eta0 = _chung_dilute(
            cv0_values[i],
            mass_values[i],
            float(omega[i]),
            tc_values[i],
            vc_values[i],
            float(dipole[i]),
            temperature,
        )
        pure_conductivity.append(max(conductivity, 1.0e-50))
        pure_viscosity.append(eta0)

    # **The class's two mass exponents disagree**: `M_i/M_j` in the conductivity rule and
    # `M_j/M_i` in the viscosity rule. Reproduced, not repaired.
    dilute_conductivity = 0.0
    for i in range(count):
        denominator = 0.0
        for j in range(count):
            ratio = mass_values[i] / mass_values[j]
            aij = (
                1.0 + math.sqrt(pure_conductivity[i] / pure_conductivity[j]) * ratio**0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + ratio))
            denominator += x[j] * aij
        dilute_conductivity += x[i] * pure_conductivity[i] / denominator

    dilute_viscosity = 0.0
    for i in range(count):
        denominator = 0.0
        for j in range(count):
            ratio = mass_values[i] / mass_values[j]
            phi = (
                1.0
                + math.sqrt(pure_viscosity[i] / pure_viscosity[j])
                * (mass_values[j] / mass_values[i]) ** 0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + ratio))
            denominator += x[j] * phi
        dilute_viscosity += x[i] * pure_viscosity[i] / denominator

    cv0_mix = sum(x[i] * cv0_values[i] for i in range(count))
    eucken_factor = cv0_mix + 1.25 * _R
    molar_mass_mix = sum(x[i] * mass_values[i] for i in range(count))

    excess = viscosity_value - dilute_viscosity
    dense_conductivity = (
        excess * eucken_factor / molar_mass_mix
        if excess > 0.0 and molar_mass_mix > 1.0e-30
        else 0.0
    )

    conductivity = max(dilute_conductivity + dense_conductivity, 1.0e-10)

    apply_checks(
        checks.derived, lambda name: conductivity if name == "conductivity" else None, warnings
    )

    return FrictionTheoryConductivityResult(
        conductivity=from_si(conductivity, "W/(m*K)"),
        pure_component_conductivity=tuple(from_si(value, "W/(m*K)") for value in pure_conductivity),
        dilute_conductivity=from_si(dilute_conductivity, "W/(m*K)"),
        dense_conductivity=from_si(dense_conductivity, "W/(m*K)"),
        warnings=tuple(warnings),
    )
