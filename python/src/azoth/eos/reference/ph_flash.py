"""``eos.ph_flash`` - the temperature a mixture reaches at a given pressure and enthalpy.

Spec: ``specs/models/eos/ph_flash.yaml``

# What this model is, and why it is the one the process layer needs

Every unit operation that adds or removes energy - a heater, a cooler, a compressor, a
valve - knows the pressure it leaves a stream at and the duty it put in, and does not
know the temperature that results. That is this model. It is the second of the two
that the port needs before a flowsheet can run at all, the first being
:mod:`azoth.eos.reference.pt_flash`.

# The procedure is an outer solve over two things that already exist

The enthalpy of a mixture at a pressure is a **strictly increasing** function of
temperature, which is what makes a bisection well posed, and it is assembled from
parts this library already has:

* the phase split at a trial temperature, from ``eos.pt_flash``;
* each phase's enthalpy, from ``eos.molar_enthalpy_entropy``.

    H(T) = (1 - beta) * H_liquid(x, Z_l) + beta * H_vapour(y, Z_v)

A single-phase state is the same expression with the whole of it on one side, and the
root that describes that phase. It is deliberately *not* expressed through ``beta``
there: the flash reports a beta outside ``[0, 1]`` for a single-phase feed - it is the
extrapolated split, not a physical one - and multiplying an enthalpy by 1.888 would be
a wrong answer shaped exactly like a right one.

# Why the bracket is a scan and not an expansion

The search needs an interval containing the answer before bisection can start. An
adaptive expansion from a guess is faster and has no hard limits, and it is also a
second thing for the two implementations to agree about: two expansions that stop at
different points take different numbers of iterations, and this project requires the
iteration counts to *match*. A fixed scan over a stated range is the same in both
languages by construction, and its one cost - a range that does not cover every
possible answer - is stated in the spec rather than hidden.

# Warnings are deduplicated, and the reason is not tidiness

The search evaluates the flash thousands of times, and the flash warns. Carrying every
warning from every trial temperature would return a tuple thousands of entries long
describing states nobody asked about, and - worse - it would make the result depend on
how many iterations the search happened to take. What a caller needs is the set of
distinct caveats that apply to the *answer*, so that is what is returned.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PhFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel, molar_enthalpy_entropy
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.ph_flash"

#: The phase the flash reports when the feed is entirely one phase. Named here
#: because the split below branches on it, and a bare string would be the same
#: literal written twice in two languages.
TWO_PHASE = "two_phase"
ALL_LIQUID = "all_liquid"
ALL_VAPOUR = "all_vapour"


def enthalpy_at(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
) -> tuple[float, dict[str, Any]]:
    """The molar enthalpy of a mixture at a temperature and pressure, and its split.

    The composition the flash settles on is the equilibrium one, so this is the
    enthalpy of the *feed* at that state - which is what makes it comparable with a
    duty a caller supplied.

    Returns:
        ``(enthalpy in J/mol, the flash's own result as a mapping)``.
    """
    flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), z)
    phase = str(flash.phase)

    if phase == TWO_PHASE:
        beta = flash.beta
        if beta is None:
            # Unreachable: the flash reports a vapour fraction in exactly the
            # two-phase case. Checked rather than asserted away because a `float |
            # None` multiplied into an enthalpy is the kind of thing that becomes a
            # `TypeError` at a caller's feet rather than here.
            raise InvalidInputError(
                "phase", "a two-phase flash reported no vapour fraction"
            )
        liquid = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            from_si(t_si, "K"),
            from_si(p_si, "Pa"),
            list(flash.x),
            flash.z_liquid,
        )
        vapour = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            from_si(t_si, "K"),
            from_si(p_si, "Pa"),
            list(flash.y),
            flash.z_vapour,
        )
        h = (1.0 - beta) * liquid.h.to_base_units().magnitude
        h += beta * vapour.h.to_base_units().magnitude
    else:
        # One phase, so the whole feed is in it and its own root describes it. beta is
        # deliberately not used - see the module docstring.
        root = flash.z_liquid if phase == ALL_LIQUID else flash.z_vapour
        one = molar_enthalpy_entropy(
            mixture, ideal_gas, from_si(t_si, "K"), from_si(p_si, "Pa"), z, root
        )
        h = one.h.to_base_units().magnitude

    return h, {"phase": phase, "flash": flash, "warnings": list(flash.warnings)}


def solve_temperature(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    h_si: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert ``H(T, P) = h_si`` for ``T`` by bracketing, then bisecting.

    Raises:
        SolverNotConvergedError: if the scan finds no sign change - the requested
            enthalpy is outside the range this model's bracket covers - or if the
            bisection reaches its cap. The two are reported apart because they mean
            different things: one is a state this model cannot represent, the other
            is an iteration that did not settle.
    """
    bracket = algorithm["bracket"]
    lower = float(bracket["lower"])
    upper = float(bracket["upper"])
    steps = int(bracket["steps"])

    warnings: list[Warning] = []
    lo, hi = _bracket_by_scan(mixture, ideal_gas, p_si, h_si, z, lower, upper, steps)
    lo_h, lo_state = enthalpy_at(mixture, ideal_gas, lo, p_si, z)

    tolerance = float(algorithm["tolerance"])
    max_iterations = int(algorithm["max_iterations"])
    iterations = 0
    mid_state = lo_state
    mid_h = lo_h
    while iterations < max_iterations:
        iterations += 1
        mid = 0.5 * (lo + hi)
        mid_h, mid_state = enthalpy_at(mixture, ideal_gas, mid, p_si, z)
        warnings.extend(mid_state["warnings"])

        # The bracket is narrowed on the *temperature*, not on the enthalpy. Relative
        # convergence on an enthalpy that crosses zero is ill-conditioned - it is
        # genuinely zero at some temperature for a datum that puts it there - while
        # the temperature interval is well behaved and is what the answer is.
        if (hi - lo) <= tolerance * mid:
            lo, hi = mid, mid
            break

        if (mid_h - h_si) * (lo_h - h_si) <= 0.0:
            hi = mid
        else:
            lo, lo_h = mid, mid_h
    else:
        # The cap is reached without the bracket meeting its tolerance. The residual
        # reported is the one the caller can act on - how far the enthalpy is from the
        # one asked for - rather than an iteration count.
        raise SolverNotConvergedError(
            iterations, abs(mid_h - h_si) / max(abs(h_si), 1.0), tolerance
        )

    return {
        "T": mid,
        "residual": abs(mid_h - h_si) / max(abs(h_si), 1.0),
        "state": mid_state,
        "iterations": iterations,
        "warnings": warnings,
    }


def _bracket_by_scan(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    h_si: float,
    z: list[float],
    lower: float,
    upper: float,
    steps: int,
) -> tuple[float, float]:
    """The narrowest interval of the scan that contains the requested enthalpy.

    Scanned from the bottom up and returning the *first* sign change, so the interval
    is determined by the bracket alone rather than by where a search happened to
    start - which is what lets the two implementations agree on the iteration count.

    Raises:
        SolverNotConvergedError: if no sign change is found on the scan.
    """
    previous_t = lower
    previous_h, _ = enthalpy_at(mixture, ideal_gas, previous_t, p_si, z)
    for index in range(1, steps + 1):
        t = lower + (upper - lower) * index / steps
        h, _ = enthalpy_at(mixture, ideal_gas, t, p_si, z)
        if (h - h_si) * (previous_h - h_si) <= 0.0:
            return previous_t, t
        previous_t, previous_h = t, h

    raise SolverNotConvergedError(
        steps,
        abs(previous_h - h_si) / max(abs(h_si), 1.0),
        0.0,
    )


def distinct_warnings(warnings: list[Warning]) -> tuple[Warning, ...]:
    """One of each distinct warning, in first-seen order.

    The search evaluates the flash thousands of times, so the same caveat arrives
    thousands of times. Returning them all would make the result depend on how many
    iterations the search took, which is a property of the algorithm rather than of
    the state the caller asked about.
    """
    seen: dict[tuple[str, str | None, str], Warning] = {}
    for warning in warnings:
        seen.setdefault((warning.code.value, warning.field, warning.message), warning)
    return tuple(seen.values())


def ph_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    P: Q,
    H: Q,
    z: list[float],
) -> PhFlashResult:
    """The temperature at which a mixture has a given molar enthalpy at a pressure.

    ``H`` is a *difference* from the datum ``ideal_gas`` carries, not an absolute
    quantity: two calls with different reference values are not comparable, and their
    difference is a plausible number rather than an error.

    Raises:
        InvalidInputError: if the composition or any ideal-gas vector is the wrong
            length, or if ``z`` is not a composition.
        OutOfRangeError: if ``P`` is not positive, or a range check on the answer
            fails.
        SolverNotConvergedError: if no temperature on the bracket covers the requested
            enthalpy, or if the bisection reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    h_si = input_to_si(spec, "H", H)
    apply_checks(checks.on_input, {"P": p_si, "H": h_si}.get, warnings)

    solved = solve_temperature(mixture, ideal_gas, p_si, h_si, list(z), spec["algorithm"])
    warnings.extend(solved["warnings"])
    state = solved["state"]
    flash = state["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["T"] if name == "T" else None,
        warnings,
    )

    return PhFlashResult(
        T=from_si(solved["T"], "K"),
        beta=flash.beta,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        phase=flash.phase,
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=solved["iterations"],
        residual=solved["residual"],
        warnings=distinct_warnings(warnings),
    )
