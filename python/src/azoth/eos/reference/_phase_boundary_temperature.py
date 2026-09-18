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
    phase_derivatives,
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
    """``d ln phi_i / dT`` for a phase at a composition and root.

    The temperature analogue of the K-value substitution's ingredient: the Newton step
    on ``S`` needs ``dK_i/dT = K_i (d ln phi_i^liquid/dT - d ln phi_i^vapour/dT)``.
    Analytic, from :func:`azoth.eos.reference._mixture_state.phase_derivatives`, which is
    the surface the second-order flashes are written in.

    It was a central difference of :func:`phase_state` at ``h = 1e-4 T`` until the
    analytic form existed, and the difference was never the problem: it is six
    ``phase_state`` evaluations per iteration, and it agrees with this one to about
    ``1e-9`` - which is the difference quotient's own truncation error, not this
    derivative's accuracy.
    """
    reduced = reduced_parameters(mixture, temperature, pressure)
    state = phase_state(reduced, mixture.kij, list(x), liquid=liquid)
    return phase_derivatives(
        reduced,
        mixture.kij,
        list(x),
        state.z,
        temperature=temperature,
        pressure=pressure,
    ).d_ln_phi_dt


def phase_boundary_temperature(
    mixture: Mixture,
    pressure: float,
    held: list[float],
    incipient: str,
    algorithm: dict[str, Any],
    capillary_pressure: float = 0.0,
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

        # The Kelvin shift, when the boundary is drawn inside a pore rather than on a flat
        # surface: ``K_cap = K exp(-Vm_L dP_cap / (R T))``, and the liquid's molar volume at the
        # cubic's own root is ``z R T / P``, so the whole exponent is
        #
        #     Vm_L dP_cap / (R T) = z_liquid dP_cap / P
        #
        # which needs neither a gas constant nor a molar-volume lookup. NeqSim writes it the long
        # way, from ``getMolarVolume("m3/mol")`` with a ``1e-4`` fallback and a ``vmL > 0.01``
        # guard; taken from the root it cannot fail either way.
        if capillary_pressure:
            kelvin = math.exp(-capillary_pressure * liquid_state.z / pressure)
            k = [value * kelvin for value in k]

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
