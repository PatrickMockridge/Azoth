"""Inverting a molar property for temperature, at a fixed pressure.

The machinery behind ``eos.ph_flash`` and ``eos.ps_flash``, which differ in exactly one
thing: whether the property being inverted is the enthalpy or the entropy. Everything
else - how the state at a trial temperature is assembled, why the bracket is a scan,
why the answer is narrowed on the temperature rather than on the property, why the
warnings are deduplicated - is the same for both, and it is subtle enough that writing
it twice would invite the two copies to disagree.

The precedent is :mod:`azoth.eos.reference._phase_boundary`, shared by
``eos.bubble_pressure`` and ``eos.dew_pressure`` because "the guard against the trivial
solution has to be written once". The same argument applies here, and it is stronger
for the single-phase branch below: a feed that is entirely one phase has no vapour
fraction, and the flash's value for it is an *extrapolation* rather than a number
nobody should use.

# What the two models share, and what they do not

They share the bracket, the bisection, the branch on the phase, and the warning
handling. They do not share the property itself: an enthalpy and an entropy are
different functions, each comes from ``eos.molar_enthalpy_entropy``, and each has its
own units. So the caller names which one it wants and this module sums it the same way
either time.
"""

from __future__ import annotations

from typing import Any, Literal

from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.units import Q, from_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel, molar_enthalpy_entropy
from azoth.eos.reference.pt_flash import pt_flash

#: Which molar property is being inverted.
Property = Literal["h", "s"]

#: The unit each property is reported in, for the conversion at the boundary.
UNIT_OF: dict[str, str] = {"h": "J/mol", "s": "J/(mol*K)"}

TWO_PHASE = "two_phase"
ALL_LIQUID = "all_liquid"


def property_at(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
    which: Property,
) -> tuple[float, dict[str, Any]]:
    """The molar property of a mixture at a temperature and pressure, and its split.

    The composition the flash settles on is the equilibrium one, so this is the
    property of the *feed* at that state - which is what makes it comparable with a
    duty or a change a caller supplied.

    Returns:
        ``(the property in its base unit, the flash's own result as a mapping)``.
    """
    flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), z)
    phase = str(flash.phase)

    if phase == TWO_PHASE:
        beta = flash.beta
        if beta is None:
            # Unreachable: the flash reports a vapour fraction in exactly the two-phase
            # case. Checked rather than asserted away because a `float | None`
            # multiplied into a property is the kind of thing that becomes a
            # `TypeError` at a caller's feet rather than here.
            raise InvalidInputError("phase", "a two-phase flash reported no vapour fraction")

        liquid = _phase_property(
            mixture, ideal_gas, t_si, p_si, list(flash.x), flash.z_liquid, which
        )
        vapour = _phase_property(
            mixture, ideal_gas, t_si, p_si, list(flash.y), flash.z_vapour, which
        )
        value = (1.0 - beta) * liquid + beta * vapour
    else:
        # One phase, so the whole feed is in it and its own root describes it. `beta`
        # is deliberately not used - see the module docstring.
        root = flash.z_liquid if phase == ALL_LIQUID else flash.z_vapour
        value = _phase_property(mixture, ideal_gas, t_si, p_si, z, root, which)

    return value, {"phase": phase, "flash": flash, "warnings": list(flash.warnings)}


def _phase_property(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    composition: list[float],
    root: float,
    which: Property,
) -> float:
    """One phase's molar property, in its base unit."""
    result = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        from_si(t_si, "K"),
        from_si(p_si, "Pa"),
        composition,
        root,
    )
    quantity: Q = getattr(result, which)
    return float(quantity.to_base_units().magnitude)


def solve_temperature(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    which: Property,
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert ``X(T, P) = target`` for ``T`` by bracketing, then bisecting.

    Raises:
        SolverNotConvergedError: if the scan finds no sign change - the requested
            property is outside the range this model's bracket covers - or if the
            bisection reaches its cap. The two are reported apart because they mean
            different things: one is a state this model cannot represent, the other is
            an iteration that did not settle.
    """
    bracket = algorithm["bracket"]
    lower = float(bracket["lower"])
    upper = float(bracket["upper"])
    steps = int(bracket["steps"])

    warnings: list[Warning] = []
    lo, hi = bracket_by_scan(mixture, ideal_gas, p_si, target, z, which, lower, upper, steps)
    lo_value, _ = property_at(mixture, ideal_gas, lo, p_si, z, which)

    tolerance = float(algorithm["tolerance"])
    max_iterations = int(algorithm["max_iterations"])
    iterations = 0
    mid = lo
    mid_value = lo_value
    mid_state: dict[str, Any] = {}
    while iterations < max_iterations:
        iterations += 1
        mid = 0.5 * (lo + hi)
        mid_value, mid_state = property_at(mixture, ideal_gas, mid, p_si, z, which)
        warnings.extend(mid_state["warnings"])

        # The bracket is narrowed on the *temperature*, not on the property. Relative
        # convergence on a property that crosses zero is ill-conditioned - an enthalpy
        # is genuinely zero at some temperature for a datum that puts it there - while
        # the temperature interval is well behaved and is what the answer is.
        if (hi - lo) <= tolerance * mid:
            break

        if (mid_value - target) * (lo_value - target) <= 0.0:
            hi = mid
        else:
            lo, lo_value = mid, mid_value
    else:
        # The cap is reached without the bracket meeting its tolerance. The residual
        # reported is the one a caller can act on.
        raise SolverNotConvergedError(
            iterations, abs(mid_value - target) / max(abs(target), 1.0), tolerance
        )

    return {
        "T": mid,
        "residual": abs(mid_value - target) / max(abs(target), 1.0),
        "state": mid_state,
        "iterations": iterations,
        "warnings": warnings,
    }


def bracket_by_scan(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    which: Property,
    lower: float,
    upper: float,
    steps: int,
) -> tuple[float, float]:
    """The narrowest interval of the scan that contains the requested property.

    Scanned from the bottom up and returning the *first* sign change, so the interval
    is determined by the bracket alone rather than by where a search happened to start
    - which is what lets the two implementations agree on the iteration count.

    **A temperature with no state is skipped rather than fatal.** Some ``(T, P)`` pairs
    on the scan have no admissible liquid root - the cubic's smallest root falls below
    the mixture's ``B``, so ``ln(Z - B)`` is the logarithm of a negative number - and
    that is not confined to the ends of the range: on methane/n-butane at 15 bar it
    happens at 160 K and nowhere else between 100 K and 400 K. Such a point has no
    property, so it cannot bracket anything. Aborting the search on one - which this did
    until the process layer needed a valve at 15 bar - made ``eos.ph_flash`` unusable at
    ordinary states, and the failure looked like a caller's mistake.

    **Only an out-of-range state is skipped.** An :class:`InvalidInputError` - a
    composition that is not a composition, a vector of the wrong length - does not depend
    on the temperature, so it is the caller's error at every point and it propagates
    immediately. Skipping those too would turn "your ``z`` is wrong" into "the solver did
    not converge", which reports a failure of the search where the search was never the
    problem.

    Raises:
        SolverNotConvergedError: if no sign change is found between two temperatures that
            both have states. The residual is infinite when nothing on the scan was
            evaluable at all, which is the honest reading rather than a number.
        InvalidInputError: from :func:`property_at`, for arguments that are the caller's
            error at every temperature.
    """
    previous: tuple[float, float] | None = None
    for index in range(steps + 1):
        t = lower + (upper - lower) * index / steps
        try:
            value, _ = property_at(mixture, ideal_gas, t, p_si, z, which)
        except OutOfRangeError:
            continue
        if previous is not None and (value - target) * (previous[1] - target) <= 0.0:
            return previous[0], t
        previous = (t, value)

    raise SolverNotConvergedError(
        steps,
        float("inf") if previous is None else abs(previous[1] - target) / max(abs(target), 1.0),
        0.0,
    )


def distinct_warnings(warnings: list[Warning]) -> tuple[Warning, ...]:
    """One of each distinct warning, in first-seen order.

    The search evaluates the flash thousands of times, so the same caveat arrives
    thousands of times. Returning them all would make the result depend on how many
    iterations the search took, which is a property of the algorithm rather than of the
    state the caller asked about.
    """
    seen: dict[tuple[str, str | None, str], Warning] = {}
    for warning in warnings:
        seen.setdefault((warning.code.value, warning.field, warning.message), warning)
    return tuple(seen.values())
