"""Inverting a molar property for temperature, at a fixed pressure.

The machinery behind ``eos.ph_flash`` and ``eos.ps_flash``, which differ in two things:
which property is being inverted, and the variable the iteration runs in. Everything
else - how the state at a trial temperature is assembled, the damping, the step clamp,
what happens when a trial temperature cannot be evaluated, why the warnings are
deduplicated - is the same for both, and it is subtle enough that writing it twice would
invite the two copies to disagree.

The precedent is :mod:`azoth.eos.reference._phase_boundary`, shared by
``eos.bubble_pressure`` and ``eos.dew_pressure`` because "the guard against the trivial
solution has to be written once".

The iteration is upstream's: ``thermodynamicoperations/flashops/PSFlash.java`` and
``PHflash.java``, NeqSim 3.20.0. Both are quasi-Newton in the temperature - entropy in
``T``, enthalpy in ``1/T`` - damped by a factor that halves whenever the residual grows,
and neither is fatal when a trial temperature cannot be evaluated.
"""

from __future__ import annotations

import math
from typing import Any, Literal

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.units import Q, from_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel, molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import pt_flash

#: Which molar property is being inverted.
Property = Literal["h", "s", "v", "u"]

#: The unit each property is reported in, for the conversion at the boundary.
UNIT_OF: dict[str, str] = {"h": "J/mol", "s": "J/(mol*K)", "v": "m**3/mol", "u": "J/mol"}

TWO_PHASE = "two_phase"
ALL_LIQUID = "all_liquid"

#: The largest temperature step one iteration may take, in kelvin.
#:
#: Upstream clamps to ten in both solvers. It is what keeps a Newton step taken far from
#: the root - where the derivative is a poor local model - from throwing the iterate into
#: a region the cubic cannot describe.
MAX_STEP = 10.0

#: How much of a Newton step the first iteration takes, before any damping.
INITIAL_FACTOR = 0.8

#: Below this the damping stops halving: a smaller factor cannot make progress and only
#: drives the iteration into its cap.
MIN_FACTOR = 0.1

#: The relative term in the entropy solver's tolerance, upstream's
#: ``RELATIVE_ENTROPY_FLASH_TOLERANCE``. The tolerance is
#: ``max(algorithm tolerance, |target| * this)``.
RELATIVE_ENTROPY_TOLERANCE = 1.0e-10

#: A residual this small that has stopped improving is accepted rather than driven on.
STAGNANT_RESIDUAL = 1.0e-4

#: Consecutive non-improving iterations at a low residual before the answer is accepted.
STAGNANT_LIMIT = 5

#: How many times the enthalpy solver may halve the gap back towards a temperature that
#: worked, before the failure is reported.
RETRY_LIMIT = 15


def property_at(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
    which: Property,
) -> tuple[float, dict[str, Any]]:
    """The molar property of a mixture at a temperature and pressure, and its split.

    The composition the flash settles on is the equilibrium one, so this is the property
    of the *feed* at that state - which is what makes it comparable with a duty or a
    change a caller supplied.

    Returns:
        ``(the property in its base unit, the flash's own result as a mapping)``.
    """
    value, _, _, flash = _evaluate(mixture, ideal_gas, t_si, p_si, z, which)
    return value, {"phase": str(flash.phase), "flash": flash, "warnings": list(flash.warnings)}


def _evaluate(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
    which: Property,
) -> tuple[float, float, float, Any]:
    """The property, the heat capacity and the volume at a state, on the branch picked.

    Carries the heat capacity and the volume alongside the value because the Newton step
    needs a derivative: ``cp`` for a temperature step and the volume for a pressure step.
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

        liquid, liquid_v = _phase_state(
            mixture, ideal_gas, t_si, p_si, list(flash.x), flash.z_liquid
        )
        vapour, vapour_v = _phase_state(
            mixture, ideal_gas, t_si, p_si, list(flash.y), flash.z_vapour
        )
        value = (1.0 - beta) * _property_of(liquid, liquid_v, which, p_si) + beta * _property_of(
            vapour, vapour_v, which, p_si
        )
        cp = (1.0 - beta) * float(liquid.cp.to_base_units().magnitude) + beta * float(
            vapour.cp.to_base_units().magnitude
        )
        volume = (1.0 - beta) * liquid_v + beta * vapour_v
    else:
        # One phase, so the whole feed is in it and *its* root describes it. `beta` is
        # deliberately not used - see the module docstring - and neither are
        # `z_liquid`/`z_vapour`, which belong to the extrapolated phase compositions:
        # where the flash reports a negative flash its phase composition is not a state,
        # and its cubic root is a different number from the feed's.
        reduced = reduced_parameters(mixture, t_si, p_si)
        root = phase_state(reduced, mixture.kij, z, liquid=(phase == ALL_LIQUID)).z
        state, volume = _phase_state(mixture, ideal_gas, t_si, p_si, z, root)
        value = _property_of(state, volume, which, p_si)
        cp = float(state.cp.to_base_units().magnitude)

    return value, cp, volume, flash


def _phase_state(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    composition: list[float],
    root: float,
) -> tuple[Any, float]:
    """One phase's state, in its base unit, and its molar volume."""
    state = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        from_si(t_si, "K"),
        from_si(p_si, "Pa"),
        composition,
        root,
    )
    return state, root * MOLAR_GAS_CONSTANT * t_si / p_si


def _property_of(state: Any, volume: float, which: Property, p_si: float) -> float:
    """The property a state carries, in its base unit."""
    if which == "v":
        return volume
    if which == "u":
        return float(state.h.to_base_units().magnitude) - p_si * volume
    quantity: Q = getattr(state, which)
    return float(quantity.to_base_units().magnitude)


def _temperature_dependent(error: Exception) -> bool:
    """Whether a failed trial is a property of *that* temperature rather than the call.

    An :class:`InvalidInputError` - a composition that is not a composition, a vector of
    the wrong length - does not depend on the temperature, so it is the caller's error at
    every point and must not be absorbed by the recovery. Everything else a trial
    temperature can raise is a property of that temperature, and the solvers recover.
    """
    return not isinstance(error, InvalidInputError)


def _slope(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    z: list[float],
    which: Property,
    temperature: float,
    fallback: float,
) -> float:
    """The derivative of the property with respect to temperature, at constant pressure.

    Taken from the property itself rather than from the heat capacity, because across a
    phase boundary the two are not the same quantity. ``cp`` is the phase-fraction-
    weighted heat capacity of the two phases; the equilibrium ``dS/dT`` at constant
    pressure also carries the latent heat of the split changing with temperature, and
    that term is the **larger** of the two. Measured on methane/n-butane at 5 bar,
    ``dS/dT`` reaches 1.3 J/(mol*K**2) where ``cp/T`` is 0.5. Stepping on the smaller one
    makes the iteration a fixed point with a gain of about three, which oscillates
    between two temperatures instead of converging.

    A central difference costs two evaluations and is exact to second order in ``delta``,
    chosen relative to the temperature so the step is the same fraction of the state in
    either implementation. Where either side cannot be evaluated the heat capacity's
    value is used instead, which is the frozen-composition slope.
    """
    delta = max(1.0e-4 * temperature, 1.0e-6)
    try:
        above, _, _, _ = _evaluate(mixture, ideal_gas, temperature + delta, p_si, z, which)
        below, _, _, _ = _evaluate(
            mixture, ideal_gas, max(temperature - delta, 1.0), p_si, z, which
        )
    except Exception:
        return fallback
    return (above - below) / (2.0 * delta)


def _slope_pressure(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    z: list[float],
    which: Property,
    pressure: float,
    fallback: float,
) -> float:
    """The derivative of the property with respect to pressure, at constant temperature.

    The pressure-side twin of :func:`_slope`, a central difference over pressure. The
    fallback is a crude single-phase estimate: ``dV/dP ~ -V/P``, ``dH/dP ~ V`` and, by
    the Maxwell relation ``dS/dP = -dV/dT``, ``dS/dP ~ -V/T``.
    """
    delta = max(1.0e-4 * pressure, 1.0)
    try:
        above, _, _, _ = _evaluate(mixture, ideal_gas, t_si, pressure + delta, z, which)
        below, _, _, _ = _evaluate(mixture, ideal_gas, t_si, max(pressure - delta, 1.0), z, which)
    except Exception:
        return fallback
    return (above - below) / (2.0 * delta)


def _start_temperature(algorithm: dict[str, Any]) -> float:
    """The starting temperature, from the spec, floored where upstream floors it."""
    return max(float(algorithm.get("initial_temperature", 300.0)), 50.0)


def solve_temperature(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    which: Property,
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert ``X(T, P) = target`` for ``T``.

    ``"h"`` runs upstream's ``PHflash.solveQ`` and ``"s"`` its ``PSFlash.solveQ``. They
    share the damping, the step clamp and the recovery, and differ in the variable the
    step is taken in and in the residual that is tested.

    Returns:
        A mapping with ``T``, ``residual``, ``state``, ``iterations`` and ``warnings``.

    Raises:
        SolverNotConvergedError: if the iteration reaches its cap without the residual
            falling below the tolerance.
        InvalidInputError: from :func:`property_at`, for arguments that are the caller's
            error at every temperature.
    """
    if which == "h":
        return _solve_enthalpy(mixture, ideal_gas, p_si, target, z, algorithm)
    if which == "s":
        return _solve_entropy(mixture, ideal_gas, p_si, target, z, algorithm)
    if which == "v":
        return _solve_volume(mixture, ideal_gas, p_si, target, z, algorithm)
    raise InvalidInputError(
        "property", "internal energy is inverted over pressure, not temperature"
    )


def _solve_entropy(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert the entropy, by upstream's ``PSFlash.solveQ``."""
    tolerance = max(float(algorithm["tolerance"]), abs(target) * RELATIVE_ENTROPY_TOLERANCE)
    stagnation = min(STAGNANT_RESIDUAL, tolerance)
    cap = int(algorithm["max_iterations"])

    temperature = _start_temperature(algorithm)
    value, cp, _, flash = _evaluate(mixture, ideal_gas, temperature, p_si, z, "s")
    warnings: list[Warning] = list(flash.warnings)

    # Upstream's initial values, kept because the damping rule reads them: a first
    # iteration counts as an improvement on ``1.0e10`` and so takes a half step, where
    # seeding the errors at infinity would make the relaxation below evaluate to zero
    # and the iteration would never move.
    iterations = 1
    error = 1.0
    error_old = 1.0e10
    factor = INITIAL_FACTOR
    correct_factor = True
    stagnant = 0

    while True:
        if error > error_old and factor > MIN_FACTOR and correct_factor:
            factor *= 0.5
        elif error < error_old and correct_factor:
            factor = 1.0
        iterations += 1

        # The residual ``target - S`` falls as ``S`` rises, so its derivative is the
        # property's slope negated.
        residual = target - value
        derivative = -_slope(mixture, ideal_gas, p_si, z, "s", temperature, cp / temperature)
        if math.isfinite(derivative) and derivative != 0.0:
            candidate = temperature - factor * residual / derivative

            if not math.isfinite(candidate):
                candidate = temperature + 1.0
                correct_factor = False
            elif candidate < 0.0:
                candidate = abs(temperature - MAX_STEP)
                correct_factor = False
            elif abs(temperature - candidate) > MAX_STEP:
                candidate = temperature - math.copysign(MAX_STEP, temperature - candidate)
                correct_factor = False
            else:
                correct_factor = True

            try:
                value, cp, _, flash = _evaluate(mixture, ideal_gas, candidate, p_si, z, "s")
            except Exception as error_raised:
                if not _temperature_dependent(error_raised):
                    raise
                # The step is undone and the damping halved, which is upstream's
                # response. The state that was good stays good, so nothing is lost but
                # the progress this step would have made, and the residual is read from
                # the state that remains.
                factor *= 0.5
                correct_factor = False
            else:
                temperature = candidate
                warnings.extend(flash.warnings)

        error_old = error
        error = abs(target - value)
        if iterations > 3 and abs(error - error_old) <= tolerance and error <= stagnation:
            stagnant += 1
        else:
            stagnant = 0

        if (error + error_old) <= tolerance and iterations >= 3:
            break
        if stagnant >= STAGNANT_LIMIT or iterations >= cap:
            break

    if error > tolerance:
        raise SolverNotConvergedError(iterations, error, tolerance)

    return {
        "T": temperature,
        "residual": error,
        "state": {"phase": str(flash.phase), "flash": flash, "warnings": list(flash.warnings)},
        "iterations": iterations,
        "warnings": distinct_warnings(warnings),
    }


#: The relative term in the volume solver's tolerance, upstream's ``PVflash`` test
#: ``|V - Vspec|/Vspec < 1e-9``.
RELATIVE_VOLUME_TOLERANCE = 1.0e-9


def _solve_volume(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert the volume, by upstream's ``PVflash.solveQ``, in temperature."""
    tolerance = max(float(algorithm["tolerance"]), abs(target) * RELATIVE_VOLUME_TOLERANCE)
    stagnation = min(STAGNANT_RESIDUAL, tolerance)
    cap = int(algorithm["max_iterations"])

    temperature = _start_temperature(algorithm)
    value, _, volume, flash = _evaluate(mixture, ideal_gas, temperature, p_si, z, "v")
    warnings: list[Warning] = list(flash.warnings)

    iterations = 1
    error = 1.0
    error_old = 1.0e10
    factor = INITIAL_FACTOR
    correct_factor = True
    stagnant = 0

    while True:
        if error > error_old and factor > MIN_FACTOR and correct_factor:
            factor *= 0.5
        elif error < error_old and correct_factor:
            factor = 1.0
        iterations += 1

        residual = target - value
        derivative = -_slope(mixture, ideal_gas, p_si, z, "v", temperature, volume / temperature)
        if math.isfinite(derivative) and derivative != 0.0:
            candidate = temperature - factor * residual / derivative

            if not math.isfinite(candidate):
                candidate = temperature + 1.0
                correct_factor = False
            elif candidate < 0.0:
                candidate = abs(temperature - MAX_STEP)
                correct_factor = False
            elif abs(temperature - candidate) > MAX_STEP:
                candidate = temperature - math.copysign(MAX_STEP, temperature - candidate)
                correct_factor = False
            else:
                correct_factor = True

            try:
                value, _, volume, flash = _evaluate(mixture, ideal_gas, candidate, p_si, z, "v")
            except Exception as error_raised:
                if not _temperature_dependent(error_raised):
                    raise
                factor *= 0.5
                correct_factor = False
            else:
                temperature = candidate
                warnings.extend(flash.warnings)

        error_old = error
        error = abs(target - value)
        if iterations > 3 and abs(error - error_old) <= tolerance and error <= stagnation:
            stagnant += 1
        else:
            stagnant = 0

        if (error + error_old) <= tolerance and iterations >= 3:
            break
        if stagnant >= STAGNANT_LIMIT or iterations >= cap:
            break

    if error > tolerance:
        raise SolverNotConvergedError(iterations, error, tolerance)

    return {
        "T": temperature,
        "residual": error,
        "state": {"phase": str(flash.phase), "flash": flash, "warnings": list(flash.warnings)},
        "iterations": iterations,
        "warnings": distinct_warnings(warnings),
    }


def _solve_enthalpy(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    target: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert the enthalpy, by upstream's ``PHflash.solveQ``, in reciprocal temperature.

    The variable is ``1/T`` rather than ``T`` because an enthalpy against temperature is
    close to linear in the reciprocal, which makes the Newton step a good model over a
    much wider range.
    """
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])
    # The residual is relative, so a target of zero has no scale to be relative to.
    # Upstream divides by ``|Hspec|`` unconditionally; refusing is the honest alternative
    # to returning whatever the division produces.
    scale = abs(target)
    if scale == 0.0:
        raise InvalidInputError(
            "H",
            "an enthalpy of exactly zero has no scale for the solver's relative residual "
            "to be measured against, so the inversion would divide by it",
        )

    temperature = _start_temperature(algorithm)
    value, cp, _, flash = _evaluate(mixture, ideal_gas, temperature, p_si, z, "h")
    warnings: list[Warning] = list(flash.warnings)

    # Upstream's initial values, kept because the damping rule reads them: a first
    # iteration counts as an improvement on ``1.0e10`` and so takes a half step.
    iterations = 1
    error = 1.0
    error_old = 1.0e10
    factor = INITIAL_FACTOR
    correct_factor = True
    retries = 0
    # The bracket upstream keeps from the sign of the residual. It starts open, because
    # the first iteration has seen one temperature and nothing bounds it from the other
    # side yet.
    min_temperature = 0.0
    max_temperature = 1.0e10

    while True:
        if error > error_old and factor > MIN_FACTOR and correct_factor:
            factor *= 0.5
        elif error < error_old and correct_factor:
            factor = iterations / (iterations + 1.0)
        iterations += 1

        # The step is taken in ``1/T``, where the residual ``(H - target)/scale`` has
        # derivative ``-T**2 * dH/dT / scale`` - so the update is applied to the
        # reciprocal and the clamps below are applied to the temperature it converts
        # back to.
        residual = (value - target) / scale
        derivative = (
            -temperature
            * temperature
            * _slope(mixture, ideal_gas, p_si, z, "h", temperature, cp)
            / scale
        )
        if math.isfinite(derivative) and derivative != 0.0:
            reciprocal = 1.0 / temperature - factor * residual / derivative
            candidate = math.inf if reciprocal == 0.0 else 1.0 / reciprocal

            if not math.isfinite(candidate):
                candidate = temperature + 1.0
                correct_factor = False
            elif candidate < 0.0:
                candidate = abs(temperature + MAX_STEP)
                correct_factor = False
            elif abs(temperature - candidate) > MAX_STEP:
                candidate = temperature - math.copysign(MAX_STEP, temperature - candidate)
                correct_factor = False
            else:
                correct_factor = True
            candidate = min(max(candidate, min_temperature + 0.1), max_temperature - 0.1)

            # A trial temperature the inner flash cannot settle is backed off towards
            # the last one that worked rather than reported. Upstream's comment is the
            # reason: a trial temperature can land where the cubic has no valid root, and
            # aborting the whole inversion there would make the model fail wherever its
            # path crossed such a region.
            accepted = False
            while True:
                try:
                    value, cp, _, flash = _evaluate(mixture, ideal_gas, candidate, p_si, z, "h")
                    accepted = True
                    break
                except Exception as error_raised:
                    if not _temperature_dependent(error_raised):
                        raise
                    retries += 1
                    if retries > RETRY_LIMIT:
                        raise SolverNotConvergedError(iterations, error, tolerance) from None
                    candidate = 0.5 * (candidate + temperature)
                    if abs(candidate - temperature) < 1.0e-09:
                        break

            if accepted:
                temperature = candidate
                warnings.extend(flash.warnings)

                # The bracket tightens from the sign of the residual at the temperature
                # just evaluated, which is how a step that overshoots is detected.
                if residual > 0.0 and temperature > max_temperature:
                    max_temperature = temperature
                elif residual < 0.0 and temperature < min_temperature:
                    min_temperature = temperature

        error_old = error
        error = abs((value - target) / scale)

        if (error + error_old) <= tolerance and iterations >= 3:
            break
        if iterations >= cap:
            break

    if error > tolerance:
        raise SolverNotConvergedError(iterations, error, tolerance)

    return {
        "T": temperature,
        "residual": error,
        "state": {"phase": str(flash.phase), "flash": flash, "warnings": list(flash.warnings)},
        "iterations": iterations,
        "warnings": distinct_warnings(warnings),
    }


def solve_pressure(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    target: float,
    z: list[float],
    which: Property,
    algorithm: dict[str, Any],
    start_pressure: float,
) -> dict[str, Any]:
    """Invert ``X(T, P) = target`` for ``P``, at a fixed temperature.

    Newton over pressure with a central-difference slope, the same damping, step clamp
    and recovery as :func:`solve_temperature`. Upstream's ``TVflash``, ``THflash``,
    ``TSflash`` and ``TUflash`` all reduce to this; they differ only in which property is
    inverted.

    Returns:
        A mapping with ``T``, ``P``, ``residual``, ``state``, ``iterations`` and
        ``warnings``.

    Raises:
        SolverNotConvergedError: if the iteration reaches its cap without the residual
            falling below the tolerance.
    """
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])
    # The residual is relative to the target's magnitude, so a value like a molar volume
    # - whose answer differs from the target only in the last bits - is not compared by a
    # subtraction that cancels to machine noise.
    scale = max(abs(target), 1.0)

    pressure = max(start_pressure, 1.0)
    value, _, volume, flash = _evaluate(mixture, ideal_gas, t_si, pressure, z, which)
    warnings: list[Warning] = list(flash.warnings)

    # Upstream's initial values, kept because the damping rule reads them: a first
    # iteration counts as an improvement on ``1.0e10`` and so takes a half step.
    iterations = 1
    error = 1.0
    error_old = 1.0e10
    factor = INITIAL_FACTOR
    correct_factor = True
    retries = 0

    while True:
        if error > error_old and factor > MIN_FACTOR and correct_factor:
            factor *= 0.5
        elif error < error_old and correct_factor:
            factor = 1.0
        iterations += 1

        residual = (value - target) / scale
        if which == "v":
            fallback = -volume / pressure
        elif which == "s":
            fallback = -volume / t_si
        else:
            fallback = volume
        derivative = _slope_pressure(mixture, ideal_gas, t_si, z, which, pressure, fallback) / scale
        if math.isfinite(derivative) and derivative != 0.0:
            candidate = pressure - factor * residual / derivative

            if not math.isfinite(candidate):
                candidate = pressure * 1.1
                correct_factor = False
            elif candidate <= 0.0:
                candidate = pressure / 2.0
                correct_factor = False
            elif abs(pressure - candidate) > 0.5 * pressure:
                candidate = pressure - math.copysign(0.5 * pressure, pressure - candidate)
                correct_factor = False
            else:
                correct_factor = True

            accepted = False
            while True:
                try:
                    value, _, volume, flash = _evaluate(
                        mixture, ideal_gas, t_si, candidate, z, which
                    )
                    accepted = True
                    break
                except Exception as error_raised:
                    if not _temperature_dependent(error_raised):
                        raise
                    retries += 1
                    if retries > RETRY_LIMIT:
                        raise SolverNotConvergedError(iterations, error, tolerance) from None
                    candidate = 0.5 * (candidate + pressure)
                    if abs(candidate - pressure) < 1.0:
                        break

            if accepted:
                pressure = candidate
                warnings.extend(flash.warnings)

        error_old = error
        error = abs((value - target) / scale)

        if (error + error_old) <= tolerance and iterations >= 3:
            break
        if iterations >= cap:
            break

    if error > tolerance:
        raise SolverNotConvergedError(iterations, error, tolerance)

    return {
        "T": t_si,
        "P": pressure,
        "residual": error,
        "state": {"phase": str(flash.phase), "flash": flash, "warnings": list(flash.warnings)},
        "iterations": iterations,
        "warnings": distinct_warnings(warnings),
    }


def distinct_warnings(warnings: list[Warning]) -> tuple[Warning, ...]:
    """One of each distinct warning, in first-seen order.

    The search evaluates the flash many times, so the same caveat arrives many times.
    Returning them all would make the result depend on how many iterations the search
    took, which is a property of the algorithm rather than of the state the caller asked
    about.
    """
    seen: dict[tuple[str, str | None, str], Warning] = {}
    for warning in warnings:
        seen.setdefault((warning.code.value, warning.field, warning.message), warning)
    return tuple(seen.values())
