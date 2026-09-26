"""``eos.liquid_conductivity_polynom`` - a liquid's thermal conductivity from the LIQCOND
polynomial.

Spec: ``specs/models/eos/liquid_conductivity_polynom.toml``, which records the mass-fraction mean,
the polynomial and the floor a negative one takes.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LiquidConductivityPolynomResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.liquid_conductivity_polynom"

#: The floor NeqSim puts under a pure component's polynomial, in W/(m K).
MIN_PURE_CONDUCTIVITY = 1.0e-10


def liquid_conductivity_polynom(
    liquid_conductivity: Sequence[Sequence[float]],
    molar_mass: Sequence[Q],
    z: Sequence[float],
    T: Q,
) -> LiquidConductivityPolynomResult:
    """A liquid mixture's thermal conductivity, from its components' LIQCOND polynomials.

    **A mass-fraction mean**, which is ``calcConductivity``: each component's polynomial is
    evaluated at ``T``, floored at ``1e-10`` W/(m K), and weighted by its mass fraction - so
    ``molar_mass`` and ``z`` are both read even though the polynomial itself is a function of
    temperature alone.

    **This is the aqueous phase's conductivity.** A gas or a hydrocarbon liquid takes PFCT
    (``azoth.eos.thermal_conductivity``) and the two answer different numbers at the same state.

    Args:
        liquid_conductivity: one row of three coefficients (``LIQCOND1``-``LIQCOND3``) per
            component.
        molar_mass: each component's molar mass.
        z: the phase's mole fractions.
        T: absolute temperature.

    Returns:
        The mixture's thermal conductivity, in ``W/(m*K)``.

    Raises:
        InvalidInputError: if the vectors differ in length.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = liquid_conductivity_polynom(
        ...     [[0.251502, 0.0005238919, -3.82111e-6], [-0.384, 0.00525, -6.37e-6]],
        ...     [q(0.04401, "kg/mol"), q(0.018015, "kg/mol")],
        ...     [0.0006691762234084198, 0.9993308237765915], q(313.15, "K"))
        >>> round(r.k.magnitude, 13)
        0.6344057039895
    """
    from azoth._models_gen import model

    spec: dict[str, Any] = model(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    count = len(z)
    if len(liquid_conductivity) != count or len(molar_mass) != count:
        raise InvalidInputError(
            "liquid_conductivity",
            f"{len(liquid_conductivity)} row(s) of coefficients, {len(molar_mass)} molar "
            f"mass(es) and {count} mole fraction(s): the mean is over the components, so the "
            f"three are one entry each",
        )

    temperature = input_to_si(spec, "T", T)
    masses = [float(z[i]) * input_to_si(spec, "molar_mass", molar_mass[i]) for i in range(count)]
    total_mass = sum(masses)

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            # The class reads a mole fraction without checking it; a negative one would take
            # mass away from the mean.
            "z": min(float(value) for value in z) if z else None,
        }.get,
        warnings,
    )

    k = 0.0
    for i, row in enumerate(liquid_conductivity):
        weight = masses[i] / total_mass if total_mass > 0.0 else 0.0
        pure = float(row[0]) + float(row[1]) * temperature + float(row[2]) * temperature**2
        k += weight * max(pure, MIN_PURE_CONDUCTIVITY)

    apply_checks(checks.derived, lambda name: k if name == "k" else None, warnings)

    return LiquidConductivityPolynomResult(k=from_si(k, "W/(m*K)"), warnings=tuple(warnings))
