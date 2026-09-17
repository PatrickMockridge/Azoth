"""``eos.pv_reflux_flash`` - the temperature at which a phase ratio is a given value.

Spec: ``specs/models/eos/pv_reflux_flash.toml``. NeqSim's ``PVrefluxflash``.

A distillation column's condenser fixes how much liquid it returns against how much it
takes off, and a reboiler fixes the same ratio the other way round. Both are a ratio of
the two phase amounts, and at a fixed pressure that ratio moves with the temperature, so
a specified ratio is a temperature.

The derivative is a secant over the last two iterates rather than a slope of the flashed
ratio, which is what makes the first step a *probe*: there is no derivative until there
are two iterates, and inventing one from the first residual alone would take a step sized
by the temperature itself.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PtFlashResult, PvRefluxFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.pv_reflux_flash"

#: Upstream's two-kelvin cap. It is what keeps a secant taken from two nearly-equal
#: residuals from throwing the iterate across the two-phase region.
MAX_STEP = 2.0

#: How far the first iteration probes, in kelvin, before any derivative exists.
PROBE = 0.1

#: Which phase the reflux ratio is of.
#:
#: A condenser asks for the vapour's - the liquid it returns over the vapour it takes -
#: and a reboiler for the liquid's. The two ratios are reciprocals, so the choice is not
#: cosmetic: it decides which end of the column's operating line the answer sits on.
VAPOUR = "vapour"
LIQUID = "liquid"


def _ratio(beta: float | None, phase: str) -> float:
    """``1 / beta_phase - 1``, which is the other phase's amount over this one's."""
    # A single-phase feed has no ratio: one of the amounts is zero, and upstream reads
    # the phase fraction as exactly one there rather than refusing. `1/1 - 1 = 0`.
    if beta is None:
        fraction = 1.0
    elif phase == VAPOUR:
        fraction = beta
    else:
        fraction = 1.0 - beta
    return 1.0 / fraction - 1.0


def pv_reflux_flash(
    mixture: Mixture, P: Q, reflux: float, phase: str, temperature: Q, z: list[float]
) -> PvRefluxFlashResult:
    """The temperature at which a phase ratio at a pressure is ``reflux``.

    Raises:
        OutOfRangeError: if ``P`` or the starting temperature is not positive.
        SolverNotConvergedError: if the iteration reaches its cap without the ratio
            meeting the one asked for.

    Example:
        >>> import azoth
        >>> from azoth.eos import components
        >>> q = azoth.ureg.Quantity
        >>> mix, _ = components.mixture_of(["methane", "n-butane"])
        >>> r = pv_reflux_flash(mix, q(2.5e6, "Pa"), 0.1873586001338361, "vapour",
        ...                     q(330.0, "K"), [0.6, 0.4])
        >>> round(r.T.to("K").magnitude, 3)
        330.0
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    start = input_to_si(spec, "temperature", temperature)
    apply_checks(checks.on_input, {"P": p_si, "temperature": start}.get, warnings)

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])

    def at(t: float) -> PtFlashResult:
        return pt_flash(mixture, T=from_si(t, "K"), P=from_si(p_si, "Pa"), z=list(z))

    # Upstream's own arrangement of the three iterates: `t_old` holds the temperature the
    # flash on hand was taken at, `t_older` the one before it, and the secant is over the
    # two. The shift happens before `t_old` is read, which is what makes the first
    # iteration's difference `T_start - 0` - and why the first step is a probe rather
    # than a secant.
    current = start
    t_old = 0.0
    t_older = 0.0
    f_old = 0.0
    f = 0.0
    iterations = 0
    flash = at(current)

    for step in range(1, int(algorithm["max_iterations"]) + 1):
        iterations = step
        f_old = f
        t_older = t_old
        t_old = current
        f = reflux - _ratio(flash.beta, phase)

        measured = t_old - t_older
        slope = (f - f_old) / measured if abs(measured) > 1.0e-12 else 0.0

        if abs(f) <= tolerance:
            break
        if step >= int(algorithm["max_iterations"]):
            raise SolverNotConvergedError(iterations, abs(f), tolerance)

        if step < 2 or slope == 0.0 or not math.isfinite(slope):
            # No derivative yet, or a degenerate one: step a fixed tenth of a kelvin
            # *away* from the residual's sign, so the second iterate brackets it.
            probe = PROBE if f > 0.0 else (-PROBE if f < 0.0 else 0.0)
            current = t_old + probe
        else:
            full = min(max(f / slope, -MAX_STEP), MAX_STEP)
            damping = min(1.0, 0.4 + 0.06 * step)
            current = t_old - full * damping
        flash = at(current)

    warnings.extend(flash.warnings)
    residual = abs(reflux - _ratio(flash.beta, phase))
    return PvRefluxFlashResult(
        T=from_si(current, "K"),
        beta=flash.beta,
        phase=flash.phase,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
