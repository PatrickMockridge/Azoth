"""``eos.water_dew_point_temperature``: the Python reference.

The mirror of ``crates/azoth-eos/src/water_dew_point_temperature.rs``. NeqSim's
``WaterDewPointTemperatureFlash``: the second phase is *prescribed* - phase 1 set to ``1e-30``
on every component but one ``water``-or-``MEG`` at ``1.0`` - and the temperature walks by
``100 K`` per unit residual to the fugacity equality.

**Every ``water`` and ``MEG`` component is walked whatever its mole fraction, and the largest
temperature reached is the answer.** NeqSim keeps the maximum over the loop it runs; a
zero-mole target is refused here rather than walked, and a fluid naming neither is refused
rather than answered with the ``0.0`` K the class writes.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WaterDewPointTemperatureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.water_dew_point_temperature"

#: The residual the fixed point stops on, in mole fraction. NeqSim's own ``1e-7``.
TOLERANCE = 1.0e-7

#: NeqSim's own ``maxNumberOfIterations``.
MAXIMUM_ITERATIONS = 10000

#: The temperature step per unit of residual, in kelvin: NeqSim's own ``100.0 * funk``.
GAIN = 100.0

#: The mole fraction every non-target component of the incipient phase is given.
TRACE = 1.0e-30

#: The two components the class seeds an incipient phase from. NeqSim matches them with
#: ``equals("water") || equals("MEG")``; this library's table spells the second ``meg``, so
#: the match is case-insensitive and the substances are the same two.
TARGETS = ("water", "meg")


def _incipient(index: int, count: int) -> list[float]:
    """The incipient phase's composition: every component at ``TRACE``, the target at one."""
    x = [TRACE] * count
    x[index] = 1.0
    return x


def water_dew_point_temperature(
    components: list[str],
    z: list[float],
    P: Q,
    eos: str = "srk",
) -> WaterDewPointTemperatureResult:
    """The temperature at which a gas condenses water.

    Raises:
        InvalidInputError: if ``z`` does not match ``components``, if no component is named
            ``water`` or ``MEG``, if one of them carries a zero mole fraction, or if it carries
            no melting-point temperature for the fixed point to start from.
        OutOfRangeError: from a trial state.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_si = float(input_to_si(spec, "P", P))
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError(
            "z",
            f"{len(components)} components and {len(z)} mole fractions",
        )

    from azoth.eos import components as databank

    mixture, _ = databank.mixture_of(components, eos=eos)
    lowered = [name.strip().lower() for name in components]

    best: tuple[float, float, int, bool] | None = None
    for index, name in enumerate(lowered):
        if name not in TARGETS:
            continue
        if not math.isfinite(z[index]) or z[index] <= 0.0:
            raise InvalidInputError(
                "z",
                f"{components[index]!r} is the phase this model seeds from and its mole "
                f"fraction is {z[index]}, so there is no dew point to find - NeqSim walks it "
                f"anyway and answers a number that is not one",
            )
        seed = mixture.components[index].melting_point_temperature
        if not seed > 0.0:
            raise InvalidInputError(
                "components",
                f"{components[index]!r} carries no melting-point temperature, which is the "
                f"temperature the fixed point starts from; a card's substance has none, and "
                f"seeding at zero kelvin is not a state this library represents",
            )
        x = _incipient(index, len(components))
        temperature = float(seed)
        residual = 0.0
        iterations = 0
        converged = False
        while iterations < MAXIMUM_ITERATIONS:
            iterations += 1
            reduced = reduced_parameters(mixture, temperature, p_si)
            kij = mixture.effective_kij(temperature)
            vapour = phase_state(reduced, kij, list(z), liquid=False)
            aqueous = phase_state(reduced, kij, x, liquid=True)
            # `beta_0` is one: the incipient phase carries no moles, so phase 0 is the whole
            # of the fluid. The ratio is taken in logs because the two coefficients span
            # orders of magnitude at a water mole fraction of `3.5e-5`.
            residual = z[index] - math.exp(aqueous.ln_phi[index] - vapour.ln_phi[index])
            temperature += GAIN * residual
            if abs(residual) < TOLERANCE:
                converged = True
                break
        if best is None or temperature > best[0]:
            best = (temperature, residual, iterations, converged)

    if best is None:
        raise InvalidInputError(
            "components",
            "the fluid names neither 'water' nor 'MEG', which are the two components this "
            "model seeds an incipient phase from. NeqSim writes a zero temperature there; "
            "this refuses.",
        )

    temperature, residual, iterations, converged = best
    return WaterDewPointTemperatureResult(
        temperature=from_si(temperature, "K"),
        residual=residual,
        iterations=iterations,
        converged=converged,
        warnings=tuple(warnings),
    )
