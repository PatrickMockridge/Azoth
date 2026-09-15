"""The temperature iteration ``eos.bubble_temperature`` and ``eos.dew_temperature`` share.

Specs: ``specs/models/eos/bubble_temperature.toml`` and ``specs/models/eos/dew_temperature.toml``

The mirror image of :mod:`azoth.eos.reference._phase_boundary`: *given one phase's
composition, at what temperature does the other phase appear at a fixed pressure?* The
same successive substitution on the K-values, with the outer loop moving temperature by
a Newton step on ``S = sum_i x_i K_i`` (or ``sum_i x_i / K_i`` for a dew point) whose
slope is the central-difference of the fugacity coefficients with respect to
temperature.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    normalise,
    phase_state,
    reduced_parameters,
    wilson_saturation_pressures,
)
from azoth.eos.reference._phase_boundary import (
    TRIVIAL_TOLERANCE,
    VAPOUR,
    _check_composition,
)


def _dfugdt(
    mixture: Mixture, pressure: float, x: list[float], *, liquid: bool, temperature: float
) -> list[float]:
    """``d ln phi_i / dT`` for a phase at a composition and root, by central difference.

    The temperature analogue of the K-value substitution's ingredient: the Newton step
    on ``S`` needs ``dK_i/dT = K_i (d ln phi_i^liquid/dT - d ln phi_i^vapour/dT)``. Taken
    from the fugacity coefficients themselves rather than an analytic derivative, for the
    same reason the flash-property solver takes its slopes by central difference.
    """
    delta = max(1.0e-4 * temperature, 1.0e-6)
    above = phase_state(
        reduced_parameters(mixture, temperature + delta, pressure),
        mixture.kij,
        list(x),
        liquid=liquid,
    )
    below = phase_state(
        reduced_parameters(mixture, max(temperature - delta, 1.0), pressure),
        mixture.kij,
        list(x),
        liquid=liquid,
    )
    return [(a - b) / (2.0 * delta) for a, b in zip(above.ln_phi, below.ln_phi, strict=True)]


def phase_boundary_temperature(
    mixture: Mixture,
    pressure: float,
    held: list[float],
    incipient: str,
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """The temperature at which the incipient phase appears, at a fixed pressure.

    Args:
        mixture: the components and their interaction parameters.
        pressure: absolute pressure in pascal.
        held: the composition of the phase that is present, in mole fractions.
        incipient: :data:`VAPOUR` for a bubble point, :data:`LIQUID` for a dew point.
        algorithm: the spec's algorithm block, which carries the tolerance and cap.

    Returns:
        A dict of ``temperature``, ``incipient``, ``k``, ``z_held``, ``z_incipient``,
        ``iterations``, ``residual`` and ``warnings``.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``held`` is the
            wrong length, has a negative entry, or does not sum to one.
        OutOfRangeError: if the mixture has no boundary at this pressure.
        SolverNotConvergedError: if the iteration hits its cap.
    """
    n = len(mixture)
    which = "bubble" if incipient == VAPOUR else "dew"
    if n < 2:
        raise InvalidInputError(
            "components",
            f"a {which} point needs two phases with different compositions, and one "
            f"component cannot have them. A pure component's is its saturation "
            f"pressure, which `eos.pure_saturation` computes.",
        )
    _check_composition(held, n, "held")

    temperature = max(algorithm.get("initial_temperature", 300.0), 50.0)

    # Seed the incipient phase *away* from the held one, as the pressure search does:
    # the two equal would make the first step trivial.
    psat = wilson_saturation_pressures(mixture.components, temperature)
    k_wilson = [psat[i] / pressure for i in range(n)]
    other = normalise(
        [held[i] * k_wilson[i] if incipient == VAPOUR else held[i] / k_wilson[i] for i in range(n)]
    )

    warnings: list[Warning] = []
    iterations = 0
    residual = math.nan

    for step in range(1, algorithm["max_iterations"] + 1):
        iterations = step

        reduced = reduced_parameters(mixture, temperature, pressure)
        warnings.extend(reduced.warnings)

        if incipient == VAPOUR:
            liquid, vapour = held, other
        else:
            liquid, vapour = other, held
        liquid_state = phase_state(reduced, mixture.kij, list(liquid), liquid=True)
        vapour_state = phase_state(reduced, mixture.kij, list(vapour), liquid=False)

        k = [
            math.exp(lp - lv)
            for lp, lv in zip(liquid_state.ln_phi, vapour_state.ln_phi, strict=True)
        ]

        # Checked before the update, and on the K-values rather than on the residual:
        # `S - 1` cancels, so it is small here long before the K-values are near 1.
        if all(abs(math.log(value)) < TRIVIAL_TOLERANCE for value in k):
            raise OutOfRangeError(
                "min_t_over_tc",
                min(temperature / c.Tc.to_base_units().magnitude for c in mixture.components),
                f"the two phases converged onto the same composition at every temperature "
                f"(max |ln K| fell below {TRIVIAL_TOLERANCE:e}), so this mixture has no "
                f"{which} point at this pressure. Unlike the flash, there is no vapour "
                f"fraction to report as absent: the temperature itself is what was being "
                f"solved for and it is not determined.",
            )

        # The amount of incipient phase the held phase would produce.
        if incipient == VAPOUR:
            s = sum(held[i] * k[i] for i in range(n))
            other = [held[i] * k[i] / s for i in range(n)]
        else:
            s = sum(held[i] / k[i] for i in range(n))
            other = [held[i] / k[i] / s for i in range(n)]
        residual = abs(s - 1.0)

        if residual <= algorithm["tolerance"]:
            if incipient == VAPOUR:
                z_held, z_incipient = liquid_state.z, vapour_state.z
            else:
                z_held, z_incipient = vapour_state.z, liquid_state.z
            return {
                "temperature": temperature,
                "incipient": other,
                "k": k,
                "z_held": z_held,
                "z_incipient": z_incipient,
                "iterations": iterations,
                "residual": residual,
                "warnings": warnings,
            }

        # The Newton step on temperature, with the derivative of `S` assembled from the
        # per-component fugacity temperature derivatives.
        if incipient == VAPOUR:
            d_liquid = _dfugdt(mixture, pressure, list(held), liquid=True, temperature=temperature)
            d_vapour = _dfugdt(
                mixture, pressure, list(other), liquid=False, temperature=temperature
            )
            dsdt = sum(held[i] * k[i] * (d_liquid[i] - d_vapour[i]) for i in range(n))
        else:
            d_liquid = _dfugdt(mixture, pressure, list(other), liquid=True, temperature=temperature)
            d_vapour = _dfugdt(mixture, pressure, list(held), liquid=False, temperature=temperature)
            dsdt = -sum(held[i] / k[i] * (d_liquid[i] - d_vapour[i]) for i in range(n))

        if not math.isfinite(dsdt) or dsdt == 0.0:
            temperature += 5.0 if s > 1.0 else -5.0
        else:
            temperature -= (s - 1.0) / dsdt
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])

    raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])
