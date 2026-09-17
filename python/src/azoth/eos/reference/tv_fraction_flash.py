"""``eos.tv_fraction_flash`` - the pressure at which a gas volume fraction is a given value.

Spec: ``specs/models/eos/tv_fraction_flash.toml``. NeqSim's ``TVfractionFlash``.

A *volume* fraction, not a mole fraction: the gas phase's share of the mixture's volume.
It is what ASTM D6377's vapour pressure is defined on, where the four-to-one ratio is
volumes and not moles, and at a fixed temperature it rises monotonically as the pressure
falls - so a specified fraction is a pressure.

The residual's derivative is a central difference of the residual itself rather than the
quotient rule upstream writes. That rule mixes a volume-*corrected* residual with
uncorrected derivatives of the mixture and the gas phase separately; differencing the
quantity the iteration actually drives is self-consistent.
"""

from __future__ import annotations

import math
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, TvFractionFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.tv_fraction_flash"

#: How far the pressure is walked down when the feed is single phase, and how many times.
WALK_FACTOR = 0.9
WALK_LIMIT = 20

#: The damping's start and its floor, upstream's adaptive parameters.
DAMPING_START = 100.0
DAMPING_FLOOR = 20.0

#: The largest pressure step one iteration may take, and the fraction of the current
#: pressure beyond which it is capped anyway.
#:
#: **Upstream's ten is ten *bar*** - NeqSim holds pressures in bar and this crate holds
#: them in pascals - so the constant is written in the unit the cap is expressed in and
#: converted once, here.
MAX_STEP_BAR = 10.0
MAX_STEP_PA = MAX_STEP_BAR * 1.0e5
MAX_STEP_FRACTION = 0.5

#: The fewest steps the iteration takes whatever the residual does.
MINIMUM_STEPS = 6


def _volume_fraction(mixture: Mixture, t: float, p: float, feed: list[float]) -> tuple[float, Any]:
    """The gas phase's share of the mixture's volume, at a pressure.

    The phase volumes are the cubic's own, `z R T / P`. Upstream's default applies the
    Peneloux translation as well, `v - sum_i x_i c_i`; this does not, because nothing
    wires a translation onto a databank-built mixture on either side.
    """
    flash = pt_flash(mixture, T=from_si(t, "K"), P=from_si(p, "Pa"), z=feed)
    if flash.beta is None:
        # A single-phase feed has no gas to take a share of. Upstream reads the phase
        # fraction as one or zero there, and so does this: the residual is then a
        # constant, which is what tells the iteration the outlet is single phase.
        return (0.0 if flash.phase is Phase.ALL_LIQUID else 1.0), flash

    reduced = reduced_parameters(mixture, t, p)
    x = list(flash.x)
    y = list(flash.y)
    v_liquid = phase_state(reduced, mixture.kij, x, liquid=True).z * MOLAR_GAS_CONSTANT * t / p
    v_vapour = phase_state(reduced, mixture.kij, y, liquid=False).z * MOLAR_GAS_CONSTANT * t / p
    total = (1.0 - flash.beta) * v_liquid + flash.beta * v_vapour
    return flash.beta * v_vapour / total, flash


def tv_fraction_flash(
    mixture: Mixture, T: Q, fraction: float, P: Q, z: list[float]
) -> TvFractionFlashResult:
    """The pressure at which a feed's gas volume fraction at a temperature is ``fraction``.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is not positive, or the feed has no gas phase.
        InvalidInputError: if ``fraction`` is not strictly inside ``(0, 1)``.
        SolverNotConvergedError: if the iteration reaches its cap.

    Example:
        >>> import azoth
        >>> from azoth.eos import components
        >>> q = azoth.ureg.Quantity
        >>> mix, _ = components.mixture_of(["methane", "n-butane"])
        >>> r = tv_fraction_flash(mix, q(330.0, "K"), 0.5, q(2.5e6, "Pa"), [0.6, 0.4])
        >>> round(r.P.to("bar").magnitude, 3)
        107.506
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    start = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": start}.get, warnings)
    if fraction <= 0.0 or fraction >= 1.0:
        raise InvalidInputError(
            "fraction",
            f"a volume fraction of {fraction} is not inside (0, 1). At zero the whole "
            f"mixture is liquid and at one it is all gas, and neither has a two-phase "
            f"pressure to find",
        )

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    feed = list(z)

    pressure = start
    found, flash = _volume_fraction(mixture, t_si, pressure, feed)

    # The preamble: a feed that is single phase at the starting pressure has no gas to
    # take a share of, so the pressure is walked down until one appears.
    attempts = 0
    while flash.beta is None and attempts < WALK_LIMIT:
        pressure *= WALK_FACTOR
        attempts += 1
        found, flash = _volume_fraction(mixture, t_si, pressure, feed)
    if flash.beta is None:
        raise OutOfRangeError(
            "fraction",
            fraction,
            f"the feed is single phase at every pressure from the {start} Pa asked for "
            f"down to {pressure} Pa, so there is no volume fraction to solve for",
        )

    error = 100.0
    error_old = error
    damping = DAMPING_START
    iterations = 0

    for step in range(1, int(algorithm["max_iterations"]) + 1):
        iterations = step
        # The derivative is a central difference of the residual the iteration drives,
        # taken at a step relative to the pressure so it is the same fraction of the
        # state in either implementation.
        h = 1.0e-4 * pressure
        above = _volume_fraction(mixture, t_si, pressure + h, feed)[0]
        below = _volume_fraction(mixture, t_si, max(pressure - h, 1.0), feed)[0]
        slope = (above - below) / (2.0 * h)
        if not math.isfinite(slope) or abs(slope) < 1.0e-30:
            raise SolverNotConvergedError(iterations, error, tolerance)

        if step > 3 and error < error_old * 0.9:
            damping = max(damping * 0.9, DAMPING_FLOOR)
        step_fraction = step / (step + damping)
        candidate = pressure - step_fraction * (found - fraction) / slope
        if candidate <= 0.0:
            candidate = pressure * 0.9
        cap = min(MAX_STEP_PA, MAX_STEP_FRACTION * abs(pressure))
        if abs(candidate - pressure) > cap:
            candidate = pressure + math.copysign(cap, candidate - pressure)

        moved = abs(candidate - pressure)
        pressure = candidate
        found, flash = _volume_fraction(mixture, t_si, pressure, feed)
        error_old = error
        error = abs(found - fraction)

        if error < 1.0e-8 and step > 3:
            break
        done = error <= tolerance and moved <= 1.0e-6
        if (done and step >= MINIMUM_STEPS) or step >= int(algorithm["max_iterations"]):
            break

    if error > 1.0e-4:
        raise SolverNotConvergedError(iterations, error, tolerance)
    if error > tolerance:
        warnings.append(
            Warning(
                code=WarningCode.SOLVER_NOT_CONVERGED,
                message=(
                    f"the volume fraction at the answer is {found}, against the {fraction} "
                    f"asked for. The iteration left by its cap rather than by its "
                    f"tolerance, and upstream accepts this to 1e-4 where its own loop "
                    f"tests 1e-6"
                ),
            )
        )

    warnings.extend(flash.warnings)
    return TvFractionFlashResult(
        P=from_si(pressure, "Pa"),
        T=from_si(t_si, "K"),
        beta=flash.beta,
        volume_fraction=found,
        phase=flash.phase,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=iterations,
        residual=error,
        warnings=tuple(warnings),
    )
