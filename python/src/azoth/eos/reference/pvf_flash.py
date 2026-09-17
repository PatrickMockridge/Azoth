"""``eos.pvf_flash`` - the temperature at which a feed's vapour fraction is a given value.

Spec: ``specs/models/eos/pvf_flash.toml``. NeqSim's ``PVFflash``.

At a fixed pressure the vapour fraction rises monotonically through the two-phase
region, from zero at the bubble point to one at the dew point, so a specified fraction
is a temperature. The search is Illinois' - a bracketing regula falsi - and it needs no
derivative of the vapour fraction, which is what makes it robust where a Newton on a
*flashed* quantity is not: the fraction is a maximum of zero and one at the ends of the
region, so its slope vanishes exactly where the iteration has to cross out of one phase
into the other.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, PvfFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.pvf_flash"

#: Upstream's `bracketAttempts` loop: ten kelvin a side, twenty times, floored at 50 K
#: and capped at 2000 K. It is what lets a feed whose two-phase region lies well away
#: from its own temperature be reached rather than reported as unanswerable.
WIDEN_STEP = 10.0
WIDEN_LIMIT = 20
COLDEST = 50.0
HOTTEST = 2000.0


def _beta_at(mixture: Mixture, p: float, feed: list[float], t: float) -> tuple[float, Any]:
    """The vapour fraction at a temperature, which is what the search brackets on.

    A single-phase feed answers one or zero rather than an absence: the search is *for*
    a fraction, and a trial outside the two-phase region has to return something ordered
    against the specification for the bracket to move.
    """
    flash = pt_flash(mixture, T=from_si(t, "K"), P=from_si(p, "Pa"), z=feed)
    if flash.phase is Phase.TWO_PHASE:
        beta = flash.beta if flash.beta is not None else 0.0
    elif flash.phase is Phase.ALL_LIQUID:
        beta = 0.0
    else:
        beta = 1.0
    return beta, flash


def pvf_flash(
    mixture: Mixture, P: Q, beta: float, temperature: Q, z: list[float]
) -> PvfFlashResult:
    """The temperature at which a feed's vapour fraction at a pressure is ``beta``.

    ``temperature`` centres the bracket - upstream searches the feed's own temperature
    span - and is not an initial guess at the answer.

    Raises:
        OutOfRangeError: if ``P`` is not positive, or ``beta`` outside ``(0, 1)``.
        InvalidInputError: if ``beta`` is exactly an endpoint, which is the bubble or
            dew point and belongs to ``eos.bubble_temperature`` / ``eos.dew_temperature``.
        SolverNotConvergedError: if the search brackets nothing, or reaches its cap.

    Example:
        >>> import azoth
        >>> from azoth.eos import components
        >>> q = azoth.ureg.Quantity
        >>> mix, _ = components.mixture_of(["methane", "n-butane"])
        >>> r = pvf_flash(mix, q(2.5e6, "Pa"), 0.8422055475803881, q(330.0, "K"), [0.6, 0.4])
        >>> round(r.T.to("K").magnitude, 6)
        330.0
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    feed_temperature = input_to_si(spec, "temperature", temperature)
    # Checked before the range checks, because the spec's own bound on `beta` is the
    # same interval and would answer first with a message that names a bound rather
    # than the two models whose job the endpoints are.
    if beta <= 0.0 or beta >= 1.0:
        which = "bubble" if beta <= 0.0 else "dew"
        raise InvalidInputError(
            "beta",
            f"a vapour fraction of {beta} is the {which} point, and those are "
            f"`eos.bubble_temperature` and `eos.dew_temperature` - calculations with "
            f"their own procedures. This model solves the interior of the two-phase "
            f"region, which is the part they do not",
        )
    apply_checks(checks.on_input, {"P": p_si, "beta": beta}.get, warnings)

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    bracket = algorithm["bracket"]
    centre = feed_temperature
    cold = float(bracket["lower"]) * centre
    hot = float(bracket["upper"]) * centre

    cold_beta = _beta_at(mixture, p_si, list(z), cold)[0]
    hot_beta = _beta_at(mixture, p_si, list(z), hot)[0]
    # Widen until the specification is inside, which is what makes the search a
    # bracketing one rather than a descent from a guess.
    attempts = 0
    while cold_beta > beta and attempts < WIDEN_LIMIT:
        cold = max(cold - WIDEN_STEP, COLDEST)
        cold_beta = _beta_at(mixture, p_si, list(z), cold)[0]
        attempts += 1
    attempts = 0
    while hot_beta < beta and attempts < WIDEN_LIMIT:
        hot = min(hot + WIDEN_STEP, HOTTEST)
        hot_beta = _beta_at(mixture, p_si, list(z), hot)[0]
        attempts += 1
    if not (cold_beta <= beta <= hot_beta):
        raise SolverNotConvergedError(0, beta - min(cold_beta, hot_beta), tolerance)

    # Illinois' method: regula falsi, with the end that has not moved halved before the
    # next step. Plain regula falsi converges from one side only and crawls when the
    # function is curved; halving the stale end restores the bisection's guarantee
    # without giving up the secant's speed.
    t_a, t_b = cold, hot
    f_a, f_b = cold_beta - beta, hot_beta - beta
    iterations = 0
    residual = float("nan")
    answer = float("nan")

    for step in range(1, int(algorithm["max_iterations"]) + 1):
        iterations = step
        t_c = min(max(t_a - f_a * (t_b - t_a) / (f_b - f_a), COLDEST), HOTTEST)
        f_c = _beta_at(mixture, p_si, list(z), t_c)[0] - beta
        residual = abs(f_c)
        if residual < tolerance:
            answer = t_c
            break
        if f_c * f_b < 0.0:
            t_a, f_a = t_c, f_c
        else:
            t_b, f_b = t_c, f_c
            f_a *= 0.5
        if abs(t_b - t_a) < 1.0e-10:
            answer = t_c
            break
    if answer != answer:  # NaN
        raise SolverNotConvergedError(iterations, residual, tolerance)

    found, flash = _beta_at(mixture, p_si, list(z), answer)
    warnings.extend(flash.warnings)
    return PvfFlashResult(
        T=from_si(answer, "K"),
        beta=found,
        phase=flash.phase,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=iterations,
        residual=abs(found - beta),
        warnings=tuple(warnings),
    )
