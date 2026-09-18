"""The pressure iteration ``eos.bubble_pressure`` and ``eos.dew_pressure`` share.

Specs: ``specs/models/eos/bubble_pressure.toml`` and ``specs/models/eos/dew_pressure.toml``

Both models ask the same question with the phases exchanged: *given one phase's
composition, at what pressure does the other phase appear?* Their equations differ,
``sum_i x_i K_i = 1`` against ``sum_i y_i / K_i = 1``, and everything else - the
initialisation, the update, the tolerance and the guard against the trivial solution
- is one piece of code run in two directions.
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

#: The ``max_i |ln K_i|`` below which the two phases have merged.
#:
#: Not tight, and deliberately. The guard is on the K-values because ``S - 1`` is a
#: weighted sum whose terms cancel, so it can sit inside the convergence tolerance
#: while the K-values are still far from one; and the threshold has to clear the
#: lowest value a genuine boundary reaches at any step, not the value a degenerate
#: state ends at. The measurements behind the constant, and why it is ``1e-2`` here
#: and ``1e-8`` in the flash, are in ``specs/models/eos/bubble_pressure.toml``.
TRIVIAL_TOLERANCE = 1.0e-02

#: The largest change in the incipient composition at which the answer is a settled state.
#:
#: The pressure and the composition are two halves of one answer and both have to converge.
#: ``S = 1`` fixes the pressure; the K-values fix the composition, and the composition is one
#: update behind whatever ``S`` last read. A Newton on ``P`` converges in half the steps the
#: fixed-point update took, and so leaves the composition half as settled.
COMPOSITION_TOL = 1.0e-12

#: The largest change in ``ln P`` one Newton step may make, ``ln(10)`` - a decade.
#:
#: The fixed-point update could not run away and a Newton can: ``P * S`` moves the pressure by
#: the residual itself, bounded by the state it was measured at, where a Newton divides by a
#: derivative that goes to zero where the locus turns. Measured before this cap: a methane-rich
#: vapour with no dew point at any pressure drove ``B = b P/RT`` to ``2.9e31`` and tripped a
#: debug assertion in ``pr_z_factor``, where the fixed-point update walked the pressure up until
#: the trivial-solution guard refused it.
MAX_LN_STEP = 2.302585092994046


def _dfugdp(
    mixture: Mixture, temperature: float, x: list[float], *, liquid: bool, pressure: float
) -> list[float]:
    """``d ln phi_i / dP`` for a phase at a composition and root.

    The pressure analogue of ``_phase_boundary_temperature``'s ``_dfugdt``: the Newton step on
    ``S`` needs ``dK_i/dP = K_i (d ln phi_i^liquid/dP - d ln phi_i^vapour/dP)``. Analytic, from
    ``phase_derivatives`` - the surface P6 item 2 ported.

    The temperature side has used the analytic form since that surface existed and this side
    never did: it moved the pressure by ``P * S``, the fixed-point update NeqSim's ``...Der``
    classes replace with exactly this derivative.
    """
    reduced = reduced_parameters(mixture, temperature, pressure)
    state = phase_state(reduced, mixture.kij, x, liquid=liquid)
    return phase_derivatives(
        reduced, mixture.kij, x, state.z, temperature=temperature, pressure=pressure
    ).d_ln_phi_dp


#: Which phase appears at the boundary being sought.
VAPOUR = "vapour"
LIQUID = "liquid"


def phase_boundary_pressure(
    mixture: Mixture,
    temperature: float,
    held: list[float],
    incipient: str,
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """The pressure at which the incipient phase appears, at a fixed temperature.

    Args:
        mixture: the components and their interaction parameters.
        temperature: absolute temperature in kelvin.
        held: the composition of the phase that is present, in mole fractions.
        incipient: :data:`VAPOUR` for a bubble point, :data:`LIQUID` for a dew point.
        algorithm: the spec's algorithm block, which carries the tolerance and cap.

    Returns:
        A dict of ``pressure``, ``incipient``, ``k``, ``z_held``, ``z_incipient``,
        ``iterations``, ``residual`` and ``warnings``.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``held`` is the
            wrong length, has a negative entry, or does not sum to one.
        OutOfRangeError: if the mixture has no boundary at this temperature.
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

    psat = wilson_saturation_pressures(mixture.components, temperature)

    # Raoult's law, which is exact in the ideal limit and a fine starting guess
    # everywhere else.
    if incipient == VAPOUR:
        pressure = sum(held[i] * psat[i] for i in range(n))
    else:
        pressure = 1.0 / sum(held[i] / psat[i] for i in range(n))
    if not math.isfinite(pressure) or pressure <= 0.0:
        raise InvalidInputError(
            "held",
            "the saturation-pressure estimate for this composition is not a positive "
            "finite number, so there is no pressure to start the search from",
        )

    # The incipient phase starts *away* from the held one. Beginning with the two
    # equal would make the first step evaluate both phases at one composition, where
    # the K-values are 1 unless the cubic has two roots to tell them apart - which is
    # the trivial solution wearing a converged face.
    k_wilson = [psat[i] / pressure for i in range(n)]
    other = normalise(
        [held[i] * k_wilson[i] if incipient == VAPOUR else held[i] / k_wilson[i] for i in range(n)]
    )

    warnings: list[Warning] = []
    iterations = 0
    residual = math.nan
    previous_other = list(other)

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
                f"the two phases converged onto the same composition at every pressure "
                f"(max |ln K| fell below {TRIVIAL_TOLERANCE:e}), so this mixture has no "
                f"{which} point at this temperature - it is at or above its critical "
                f"condition. Unlike the flash, there is no vapour fraction to report as "
                f"absent: the pressure itself is what was being solved for and it is "
                f"not determined.",
            )

        # The amount of incipient phase the held phase would produce.
        previous_other = list(other)
        if incipient == VAPOUR:
            s = sum(held[i] * k[i] for i in range(n))
            other = [held[i] * k[i] / s for i in range(n)]
        else:
            s = sum(held[i] / k[i] for i in range(n))
            other = [held[i] / k[i] / s for i in range(n)]
        residual = abs(s - 1.0)
        # The composition's own convergence: the largest move of any mole fraction.
        #
        # The rust kernel carries the reasoning. `S = 1` fixes the pressure; the K-values fix
        # the composition, and the composition is one update behind whatever `S` last read. A
        # Newton on `P` converges in half the steps the fixed-point update took, and so leaves
        # the composition half as settled - measured, the same pressure to `1e-12` with the
        # vapour fraction differing by `2.3e-8`.
        settled = max(abs(a - b) for a, b in zip(other, previous_other, strict=True))

        if residual <= algorithm["tolerance"] and settled <= COMPOSITION_TOL:
            if incipient == VAPOUR:
                z_held, z_incipient = liquid_state.z, vapour_state.z
            else:
                z_held, z_incipient = vapour_state.z, liquid_state.z
            return {
                "pressure": pressure,
                "incipient": other,
                "k": k,
                "z_held": z_held,
                "z_incipient": z_incipient,
                "iterations": iterations,
                "residual": residual,
                "warnings": warnings,
            }

        # **The Newton step on `ln P`**, replacing the fixed-point `P * S`. `S` is `sum x_i K_i`
        # for a bubble point and `sum y_i / K_i` for a dew point, so
        #
        #     d S / d ln P = P * sum_i (held_i K_i^{+-1}) * (d ln phi_i^L/dP - d ln phi_i^V/dP)
        #
        # with the sign of the `K` power carrying the dew point's reciprocal. The fixed-point
        # update is a contraction only near the answer; a Newton more than squares the error.
        if incipient == VAPOUR:
            liquid, vapour = held, other
        else:
            liquid, vapour = other, held
        d_liquid = _dfugdp(mixture, temperature, list(liquid), liquid=True, pressure=pressure)
        d_vapour = _dfugdp(mixture, temperature, list(vapour), liquid=False, pressure=pressure)
        dsdp = 0.0
        for i in range(n):
            if incipient == VAPOUR:
                dsdp += held[i] * k[i] * (d_liquid[i] - d_vapour[i])
            else:
                dsdp -= held[i] / k[i] * (d_liquid[i] - d_vapour[i])
        dsdx = dsdp * pressure
        if math.isfinite(dsdx) and dsdx != 0.0:
            ln_step = max(-MAX_LN_STEP, min(MAX_LN_STEP, -(s - 1.0) / dsdx))
            pressure = math.exp(math.log(pressure) + ln_step)
        else:
            # The fixed-point update as the fallback, for a slope the derivative surface cannot
            # supply.
            pressure = pressure * s if incipient == VAPOUR else pressure / s
        if not math.isfinite(pressure) or pressure <= 0.0:
            raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])

    raise SolverNotConvergedError(iterations, residual, algorithm["tolerance"])


def _check_composition(values: list[float], n: int, field: str) -> None:
    """A composition must be one entry per component, in ``[0, 1]``, summing to one."""
    if len(values) != n:
        raise InvalidInputError(
            field, f"a composition for {n} components has {len(values)} entries"
        )
    for i, value in enumerate(values):
        if value < 0.0:
            raise InvalidInputError(
                field, f"{field}[{i}] is {value} but a mole fraction cannot be negative"
            )
    if abs(sum(values) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            field,
            f"the composition sums to {sum(values)}, not to one. Renormalising it here "
            f"would make a caller's error invisible in every number downstream, so it "
            f"is refused instead",
        )
