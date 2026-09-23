"""The reactive PH flash's outer loop, the Python twin of
``crates/azoth-reactions/src/reactive_ph_flash.rs``.

Given a pressure and a total enthalpy it looks for the temperature: an outer loop on ``T``
wraps the reactive TP flash, and each pass asks the inner flash what the fluid's enthalpy is
at that temperature and steps towards the one that was specified.

**The loop is not the one the class's docstring describes.** It says a Newton on ``1/T``
"following Michelsen 1987", and the file carries a ``computeEnthalpyDerivative`` computing
exactly ``dQ/d(1/T)`` - **which nothing calls**. ``solveEnthalpySpec`` is a secant on ``T``
with a bisection fallback bracketed on ``[50, 5000] K``, and ``1/T`` appears nowhere but the
comments. This port follows the code.

**The enthalpy is thermochemical**, because NeqSim's process-stream enthalpy is a sensible one
that excludes the formation enthalpies and a reactive calculation's composition moves. The
caller supplies it - and the class's own is **twice** the fluid's, because the two phase
objects its constructor leaves each report the whole fluid's moles, so both halves of the
comparison carry the factor. It cancels in the residual, which is a ratio.

The inner flash is the caller's too: a pass needs its pass count, its thermochemical enthalpy
and its heat capacity, and nothing here knows about equations of state or reactions.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import NamedTuple

#: The outer loop's pass cap, from ``MAX_OUTER_ITER``.
MAX_OUTER_ITERATIONS = 200

#: The tolerance on the normalised enthalpy residual, from ``TOL``.
TOL = 1.0e-8

#: The largest temperature step a pass may take, from ``MAX_T_STEP``.
MAX_TEMPERATURE_STEP = 50.0

#: The lowest temperature the loop will try, from ``T_MIN``.
T_MIN = 50.0

#: The highest, from ``T_MAX``.
T_MAX = 5000.0

#: The bracket width below which the loop declares convergence.
BRACKET_TOLERANCE = 1.0e-6


class PhState(NamedTuple):
    """What one pass of the inner flash answers with."""

    #: The flash's own passes at this temperature, which the loop reports summed.
    iterations: int
    #: The state's **thermochemical** enthalpy.
    thermochemical_enthalpy: float
    #: The state's heat capacity at constant pressure.
    cp: float


class PhFlashOutcome(NamedTuple):
    """What the PH flash answers with."""

    #: The temperature the loop stopped at.
    temperature: float
    #: ``isConverged``. **Also true where the *bracket* closed rather than the residual.**
    converged: bool
    #: ``getOuterIterations``.
    outer_iterations: int
    #: ``getTotalInnerIterations``: every inner flash's passes, summed.
    total_inner_iterations: int


def reactive_ph_flash(
    initial_temperature: float,
    thermochemical_enthalpy_spec: float,
    inner: Callable[[float], PhState],
) -> PhFlashOutcome:
    """``solveEnthalpySpec``: the secant-with-bisection loop on ``T``.

    ``thermochemical_enthalpy_spec`` is the specification **with the formation inventory
    already added**, which is what the class's constructor builds.
    """
    absolute_spec = abs(thermochemical_enthalpy_spec)
    if absolute_spec < 1.0e-10:
        # The class's own guard: a near-zero enthalpy would divide by nothing.
        absolute_spec = 1.0

    def residual(state: PhState) -> float:
        return (state.thermochemical_enthalpy - thermochemical_enthalpy_spec) / absolute_spec

    state = inner(initial_temperature)
    total_inner_iterations = state.iterations
    temperature = initial_temperature
    error = residual(state)

    if abs(error) < TOL:
        return PhFlashOutcome(temperature, True, 0, total_inner_iterations)

    # The bracket, tracked on the sign of the residual: `H(T)` rises with `T`, so a positive
    # residual means the temperature is too high.
    low = T_MIN
    high = T_MAX
    if error > 0.0:
        high = temperature
    else:
        low = temperature

    previous_temperature = temperature
    previous_error = math.nan
    converged = False
    outer_iterations = 0

    for iteration in range(MAX_OUTER_ITERATIONS):
        outer_iterations = iteration + 1
        has_bracket = low > T_MIN and high < T_MAX

        def newton(temperature: float, state: PhState = state) -> float:
            cp = 100.0 if abs(state.cp) < 1.0e-20 else state.cp
            return temperature - (state.thermochemical_enthalpy - thermochemical_enthalpy_spec) / cp

        if (
            iteration == 0
            or math.isnan(previous_error)
            or abs(temperature - previous_temperature) < 1.0e-12
        ):
            # The first pass, and any pass that did not move: a heat-capacity Newton step.
            following = newton(temperature)
        else:
            # The secant step, which is what captures `dH/dT` including the reaction's own
            # contribution - the heat the shifting equilibrium absorbs or releases.
            slope = (error - previous_error) / (temperature - previous_temperature)
            following = temperature - error / slope if abs(slope) > 1.0e-30 else newton(temperature)

        # A step larger than the cap is trimmed to it, and a step that leaves the bracket is
        # replaced by the bracket's midpoint - the guaranteed-convergence fallback.
        if abs(following - temperature) > MAX_TEMPERATURE_STEP:
            following = temperature + math.copysign(MAX_TEMPERATURE_STEP, following - temperature)
        if has_bracket and (following <= low or following >= high):
            following = 0.5 * (low + high)
        following = max(T_MIN, min(T_MAX, following))

        # A step that went nowhere: bisect if there is a bracket, else move a kelvin the way
        # the residual points.
        if abs(following - temperature) < 1.0e-12:
            following = (
                0.5 * (low + high) if has_bracket else temperature + (-1.0 if error > 0.0 else 1.0)
            )

        previous_temperature = temperature
        previous_error = error

        state = inner(following)
        total_inner_iterations += state.iterations
        temperature = following
        error = residual(state)

        if error > 0.0 and temperature < high:
            high = temperature
        elif error < 0.0 and temperature > low:
            low = temperature

        if abs(error) < TOL:
            converged = True
            break
        has_bracket = low > T_MIN and high < T_MAX
        if has_bracket and (high - low) < BRACKET_TOLERANCE:
            converged = True
            break

    return PhFlashOutcome(temperature, converged, outer_iterations, total_inner_iterations)
