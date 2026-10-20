"""``eos.friction_theory_viscosity`` - the friction-theory viscosity of
Quinones-Cisneros and Firoozabadi, over a Chung baseline and Wilke mixing.

Spec: ``specs/models/eos/friction_theory_viscosity.toml``, which records
NeqSim's ``FrictionTheoryViscosityMethod.calcViscosity`` and its helpers.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FrictionTheoryViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.friction_theory_viscosity"

#: NeqSim's gas constant, ``ThermodynamicConstantsInterface.R``.
R = 8.3144621
#: The floor ``MIN_VISCOSITY_PA_S`` puts under the answer.
MIN_VISCOSITY_PA_S = 1.0e-6
#: Micropoise to Pa*s.
MICRO_POISE_TO_PA_S = 1.0e-7

#: The PR and SRK reduced-constant sets, as ``setPengRobinsonConstants`` and
#: ``setSrkConstants`` write them.
_PR = {
    "kapac": -0.140464,
    "kaprc": 1.19902e-2,
    "kaprrc": 8.55115e-4,
    "kapa": [
        [-4.89197e-2, 0.270572, -4.48111e-2],
        [0.270572, -1.10473e-4, 4.08972e-5],
        [-4.48111e-2, 4.08972e-5, -5.79765e-9],
    ],
    "kapr": [
        [-0.357875, 0.637572, -7.9024e-2],
        [0.637572, -6.02128e-5, 3.72408e-5],
        [-7.9024e-2, 3.72408e-5, -5.65610e-9],
    ],
    "kaprr": 1.37290e-8,
}
_SRK = {
    "kapac": -0.165302,
    "kaprc": 6.99574e-3,
    "kaprrc": 1.26358e-3,
    "kapa": [
        [-0.114804, 0.246622, -3.94638e-2],
        [0.246622, -1.15648e-4, 4.18863e-5],
        [-3.94638e-2, 4.18863e-5, -5.91999e-9],
    ],
    "kapr": [
        [-0.315903, 0.566713, -7.29995e-2],
        [0.566713, -1.0086e-4, 5.17459e-5],
        [-7.29995e-2, 5.17459e-5, -5.68708e-9],
    ],
    "kaprr": 1.35994e-8,
}


def _kapr(phi: float, gamma: float, c: dict[str, Any]) -> float:
    matrix = c["kapr"]
    return float(
        c["kaprc"]
        + matrix[0][0] * (gamma - 1.0)
        + (matrix[1][0] + matrix[1][1] * phi) * math.expm1(gamma - 1.0)
        + (matrix[2][0] + matrix[2][1] * phi + matrix[2][2] * phi * phi)
        * math.expm1(2.0 * gamma - 2.0)
    )


def _kaprr(phi: float, gamma: float, c: dict[str, Any]) -> float:
    return float(c["kaprrc"] + c["kaprr"] * phi * math.expm1(2.0 * gamma) * (gamma - 1.0) ** 2)


def _kapa(phi: float, gamma: float, c: dict[str, Any]) -> float:
    matrix = c["kapa"]
    return float(
        c["kapac"]
        + matrix[0][0] * (gamma - 1.0)
        + (matrix[1][0] + matrix[1][1] * phi) * math.expm1(gamma - 1.0)
        + (matrix[2][0] + matrix[2][1] * phi + matrix[2][2] * phi * phi)
        * math.expm1(2.0 * gamma - 2.0)
    )


def friction_theory_viscosity(
    tc: Sequence[Q],
    pc: Sequence[Q],
    critical_volume: Sequence[Q],
    omega: Sequence[float],
    molar_mass: Sequence[Q],
    critical_viscosity: Sequence[Q],
    viscosity_correction_factor: Sequence[float],
    dipole: Sequence[float],
    x: Sequence[float],
    T: Q,
    p_repulsive: Q,
    p_attractive: Q,
    peng_robinson: bool,
) -> FrictionTheoryViscosityResult:
    """The friction-theory viscosity of a mixture at a state.

    The pure-component baseline is the class's own Chung form, the dilute mixture is
    Wilke's rule, and the excess friction is
    ``kappa_r P_rep + kappa_a P_att + kappa_rr P_rep^2``. ``p_repulsive`` and
    ``p_attractive`` are the caller's EOS values; this id does not build one.

    Args:
        tc: each component's critical temperature.
        pc: each component's critical pressure.
        critical_volume: each component's critical molar volume.
        omega: each component's acentric factor.
        molar_mass: each component's molar mass.
        critical_viscosity: each component's critical viscosity.
        viscosity_correction_factor: NeqSim's ``VISCFACT``, the Chung ``Fc`` term.
        dipole: each component's dipole moment, in debye.
        x: the phase's mole fractions.
        T: absolute temperature.
        p_repulsive: the EOS's repulsive pressure.
        p_attractive: the EOS's attractive pressure.
        peng_robinson: true for ``PhasePrEos``, false for the SRK constant set.

    Returns:
        The mixture viscosity, each component's pure baseline, and the two terms.

    Raises:
        InvalidInputError: if the vectors differ in length.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = friction_theory_viscosity(
        ...     [q(190.56, "K"), q(425.12, "K")],
        ...     [q(45.99, "bar"), q(37.96, "bar")],
        ...     [q(99.0, "cm**3/mol"), q(255.0, "cm**3/mol")],
        ...     [0.0115, 0.2002],
        ...     [q(0.016043, "kg/mol"), q(0.058123, "kg/mol")],
        ...     [q(1.5293e-5, "Pa*s"), q(2.5768e-5, "Pa*s")],
        ...     [0.0, 0.0], [0.0, 0.0],
        ...     [0.09542082642782022, 0.9045791735721797],
        ...     q(300.0, "K"), q(829.1753495002592, "Pa"),
        ...     q(-809.1753495002592, "Pa"), False)
        >>> round(r.mu.magnitude, 16)
        0.0001383523806148
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(x)
    if (
        len(tc) != count
        or len(pc) != count
        or len(critical_volume) != count
        or len(omega) != count
        or len(molar_mass) != count
        or len(critical_viscosity) != count
        or len(viscosity_correction_factor) != count
        or len(dipole) != count
    ):
        raise InvalidInputError(
            "x",
            f"{count} mole fraction(s) against {len(tc)} critical temperature(s), "
            f"{len(pc)} critical pressure(s), {len(critical_volume)} critical volume(s), "
            f"{len(omega)} acentric factor(s), {len(molar_mass)} molar mass(es), "
            f"{len(critical_viscosity)} critical viscosity(ies), "
            f"{len(viscosity_correction_factor)} correction factor(s) and {len(dipole)} "
            f"dipole moment(s): every one is one entry per component",
        )

    temperatures = [input_to_si(spec, "tc", value) for value in tc]
    pressures = [input_to_si(spec, "pc", value) for value in pc]
    volumes = [input_to_si(spec, "critical_volume", value) for value in critical_volume]
    masses = [input_to_si(spec, "molar_mass", value) for value in molar_mass]
    viscosities = [input_to_si(spec, "critical_viscosity", value) for value in critical_viscosity]
    temperature = input_to_si(spec, "T", T)
    p_rep = input_to_si(spec, "p_repulsive", p_repulsive)
    p_att = input_to_si(spec, "p_attractive", p_attractive)

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "tc": min(temperatures) if temperatures else None,
            "pc": min(pressures) if pressures else None,
            "critical_volume": min(volumes) if volumes else None,
            "molar_mass": min(masses) if masses else None,
            "critical_viscosity": min(viscosities) if viscosities else None,
            "x": min(float(value) for value in x) if x else None,
            "p_repulsive": p_rep,
        }.get,
        warnings,
    )

    pure: list[float] = []
    for i in range(count):
        vc_cm3 = volumes[i] * 1.0e6
        relative_dipole = 131.3 * float(dipole[i]) / math.sqrt(vc_cm3 * temperatures[i])
        fc = (
            1.0
            - 0.2756 * float(omega[i])
            + 0.059035 * relative_dipole**4
            + float(viscosity_correction_factor[i])
        )
        t_star = 1.2593 * temperature / temperatures[i]
        var_last = -6.435e-4 * t_star**0.14874 * math.sin(18.0323 * t_star**-0.76830 - 7.27371)
        omega_visc = (
            1.16145 / t_star**0.14874
            + 0.52487 / math.exp(0.77320 * t_star)
            + 2.16178 / math.exp(2.43787 * t_star)
            + var_last
        )
        pure.append(
            40.785
            * math.sqrt(masses[i] * 1000.0 * temperature)
            / (vc_cm3 ** (2.0 / 3.0) * omega_visc)
            * fc
        )

    molar_mass_pow = [mass**0.3 for mass in masses]
    mass_sum = sum(float(x[i]) / molar_mass_pow[i] for i in range(count))
    visk0 = 0.0
    for i in range(count):
        denominator = 0.0
        for j in range(count):
            phi = (
                1.0 + math.sqrt(pure[i] / pure[j]) * (masses[j] / masses[i]) ** 0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + masses[i] / masses[j]))
            denominator += float(x[j]) * phi
        visk0 += float(x[i]) * pure[i] / denominator

    c = _PR if peng_robinson else _SRK
    kaprmx = 0.0
    kapamx = 0.0
    kaprrmx = 0.0
    for i in range(count):
        nci = viscosities[i] * 1.0e7
        pc_bar = pressures[i] / 1.0e5
        phi = 1.0e1 * R * temperatures[i] / pc_bar
        gamma = temperatures[i] / temperature
        zi = float(x[i]) / (molar_mass_pow[i] * mass_sum)
        kaprmx += zi * nci * _kapr(phi, gamma, c) / pc_bar
        kapamx += zi * nci * _kapa(phi, gamma, c) / pc_bar
        kaprrmx += zi * nci * _kaprr(phi, gamma, c) / (pc_bar * pc_bar)
    visk1 = kaprmx * p_rep + kapamx * p_att + kaprrmx * p_rep * p_rep

    mu = max((visk0 + visk1) * MICRO_POISE_TO_PA_S, MIN_VISCOSITY_PA_S)

    apply_checks(checks.derived, lambda name: mu if name == "mu" else None, warnings)

    return FrictionTheoryViscosityResult(
        mu=from_si(mu, "Pa*s"),
        pure_component_viscosity=tuple(
            from_si(value * MICRO_POISE_TO_PA_S, "Pa*s") for value in pure
        ),
        dilute_viscosity=from_si(visk0 * MICRO_POISE_TO_PA_S, "Pa*s"),
        friction_viscosity=from_si(visk1 * MICRO_POISE_TO_PA_S, "Pa*s"),
        warnings=tuple(warnings),
    )
