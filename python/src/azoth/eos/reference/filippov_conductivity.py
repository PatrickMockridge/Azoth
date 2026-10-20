"""``eos.filippov_conductivity`` - a liquid's thermal conductivity from the Filippov
mixing rule, with the optional Missenard pressure correction.

Spec: ``specs/models/eos/filippov_conductivity.toml``, which records
NeqSim's ``FilippovConductivityMethod.calcConductivity`` and its two helpers.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FilippovConductivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.filippov_conductivity"

#: The floor NeqSim puts under a pure component's polynomial, in W/(m K).
MIN_PURE_CONDUCTIVITY = 1.0e-10


def filippov_conductivity(
    liquid_conductivity: Sequence[Sequence[float]],
    molar_mass: Sequence[Q],
    z: Sequence[float],
    tc: Sequence[Q],
    pc: Sequence[Q],
    T: Q,
    P: Q,
    filippov_coefficient: float,
    use_pressure_correction: bool,
) -> FilippovConductivityResult:
    """A liquid mixture's thermal conductivity from the components' polynomials.

    The pure conductivities are the LIQCOND polynomials at ``T``, each floored at
    ``1e-10`` W/(m K), mixed by **mass** fraction:
    ``k = sum_i w_i k_i - c sum_{i<j} w_i w_j |k_i - k_j|``. The Missenard factor
    ``0.98 + 0.0079 Pr sqrt(Tr)`` is applied where ``use_pressure_correction`` is set
    and the state is inside its window (``Tr < 0.95`` and ``Pr >= 0.01``); it is never
    below one.

    Args:
        liquid_conductivity: one row of three coefficients (``LIQCOND1``-``LIQCOND3``)
            per component.
        molar_mass: each component's molar mass.
        z: the phase's mole fractions.
        tc: each component's critical temperature.
        pc: each component's critical pressure.
        T: absolute temperature.
        P: absolute pressure.
        filippov_coefficient: the interaction coefficient, ``0.72`` in the class.
        use_pressure_correction: whether Missenard's correction is applied.

    Returns:
        The mixture's thermal conductivity, the Missenard factor actually applied, and
        each component's floored polynomial value.

    Raises:
        InvalidInputError: if the vectors differ in length.
        OutOfRangeError: if an input is outside its declared range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = filippov_conductivity(
        ...     [[0.251502, 0.0005238919, -3.82111e-6], [-0.384, 0.00525, -6.37e-6]],
        ...     [q(0.04401, "kg/mol"), q(0.018015, "kg/mol")],
        ...     [0.0006691762234084198, 0.9993308237765915],
        ...     [q(304.19, "K"), q(647.3, "K")],
        ...     [q(73.815, "bar"), q(220.89, "bar")],
        ...     q(313.15, "K"), q(50.0, "bar"), 0.72, True)
        >>> round(r.k.magnitude, 13)
        0.6337077403925
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(z)
    if (
        len(liquid_conductivity) != count
        or len(molar_mass) != count
        or len(tc) != count
        or len(pc) != count
    ):
        raise InvalidInputError(
            "liquid_conductivity",
            f"{len(liquid_conductivity)} coefficient row(s), {len(molar_mass)} molar "
            f"mass(es), {len(tc)} critical temperature(s), {len(pc)} critical pressure(s) "
            f"and {count} mole fraction(s): the mean is over the components, so all five "
            f"are one entry each",
        )

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    molar_masses = [input_to_si(spec, "molar_mass", value) for value in molar_mass]
    masses = [float(z[i]) * molar_masses[i] for i in range(count)]
    total_mass = sum(masses)
    critical_temperatures = [input_to_si(spec, "tc", value) for value in tc]
    critical_pressures = [input_to_si(spec, "pc", value) for value in pc]

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "P": pressure,
            "z": min(float(value) for value in z) if z else None,
            "molar_mass": min(molar_masses) if molar_masses else None,
            "tc": min(critical_temperatures) if critical_temperatures else None,
            "pc": min(critical_pressures) if critical_pressures else None,
        }.get,
        warnings,
    )

    pure: list[float] = []
    for row in liquid_conductivity:
        value = float(row[0]) + float(row[1]) * temperature + float(row[2]) * temperature**2
        pure.append(max(value, MIN_PURE_CONDUCTIVITY))

    weights = [mass / total_mass if total_mass > 0.0 else 0.0 for mass in masses]
    linear = sum(w * k for w, k in zip(weights, pure, strict=True))
    correction = 0.0
    for i in range(count):
        for j in range(i + 1, count):
            correction += weights[i] * weights[j] * abs(pure[i] - pure[j])
    k = max(linear - filippov_coefficient * correction, MIN_PURE_CONDUCTIVITY)

    tc_mix = sum(float(z[i]) * critical_temperatures[i] for i in range(count))
    pc_mix = sum(float(z[i]) * critical_pressures[i] for i in range(count))
    applied = 1.0
    if use_pressure_correction and tc_mix > 0.0 and pc_mix > 0.0:
        reduced_temperature = temperature / tc_mix
        reduced_pressure = pressure / pc_mix
        if reduced_temperature < 0.95 and reduced_pressure >= 0.01:
            applied = max(0.98 + 0.0079 * reduced_pressure * reduced_temperature**0.5, 1.0)
    k *= applied

    apply_checks(checks.derived, lambda name: k if name == "k" else None, warnings)

    return FilippovConductivityResult(
        k=from_si(k, "W/(m*K)"),
        pressure_correction=applied,
        pure_component_conductivity=tuple(from_si(value, "W/(m*K)") for value in pure),
        warnings=tuple(warnings),
    )
