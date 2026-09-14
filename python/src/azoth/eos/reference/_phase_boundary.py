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
                "pressure": pressure,
                "incipient": other,
                "k": k,
                "z_held": z_held,
                "z_incipient": z_incipient,
                "iterations": iterations,
                "residual": residual,
                "warnings": warnings,
            }

        # Too much incipient phase means move away from it: for a bubble point `S` is
        # how much vapour the liquid would give off, so `S > 1` is too much and `P * S`
        # raises the pressure; for a dew point `S` is how much liquid the vapour would
        # condense, so `S < 1` is too little and `P / S` raises it.
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
