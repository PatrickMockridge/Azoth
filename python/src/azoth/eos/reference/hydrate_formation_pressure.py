"""``eos.hydrate_formation_pressure``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_formation_pressure.rs``:

```text
solve  f_w^hydrate(T, P, fugacities) / f_w^fluid(T, P) - 1 = 0   for P
```

The same equilibrium :mod:`azoth.eos.reference.hydrate_formation_temperature` reads the
other way round, with the same kernel underneath and the same divergence stated: NeqSim
iterates ``P <- P (f_w^hydrate/f_w^fluid)`` and this brackets and bisects, so the two agree
on the root and differ on where to stop - measured at 288.15 K, NeqSim's answer is
``1.03e-5`` relative low with a residual of ``-1.15e-6`` at it where this bisects to
``1.65e-9``.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import (
    HydrateFormationPressureResult,
    HydrateStructure,
    PtFlashResult,
)
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import _hydrate
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.hydrate_formation_pressure"


def _reference_water_fugacity(eos: str, t: float, p: float) -> float:
    """The fugacity of pure water at a state on the cubic the fluid runs, in Pa."""
    from azoth.eos import components as databank

    water, _ = databank.mixture_of(["water"], eos=eos)
    reduced = reduced_parameters(water, t, p)
    state = phase_state(reduced, water.kij, [1.0], liquid=True)
    return math.exp(state.ln_phi[0]) * p


def _guest_fugacities(flash: PtFlashResult, z: list[float], p: float) -> list[float]:
    """The guests' fugacities at a flash, in Pa, from the phase ``setFug`` reads."""
    composition = list(flash.y) if flash.phase == "two_phase" else list(z)
    coefficients = flash.ln_phi_liquid if flash.phase == "all_liquid" else flash.ln_phi_vapour
    return [
        fraction * math.exp(coefficient) * p
        for fraction, coefficient in zip(composition, coefficients, strict=True)
    ]


def hydrate_formation_pressure(
    components: list[str],
    T: Q,
    z: list[float],
    eos: str = "srk",
    hydrate_model: str = _hydrate.PVTSIM,
) -> HydrateFormationPressureResult:
    """The pressure at which a fluid's hydrate appears, at a temperature.

    **The components cross by name**, as for the formation temperature, because the hydrate's
    guest tables are keyed by name and a :class:`~azoth.eos.mixture.Component` carries none.

    Args:
        components: the substances, by name.
        T: absolute temperature; the pressure is what is solved for.
        z: the overall mole fractions. Checked rather than renormalised.
        eos: the cubic the fluid runs, which is also the one the hydrate's reference water
            phase is built from.

    Raises:
        InvalidInputError: if the fluid has no water, or nothing in it occupies a cage.
        OutOfRangeError: if ``T`` is not positive or a trial's cavity sum has no value.
        SolverNotConvergedError: if the scan finds no sign change or the bisection hits its
            cap.
    """
    from azoth import _models_gen
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = input_to_si(spec, "T", T)
    apply_checks(checks.on_input, {"T": t_si}.get, warnings)

    hydration = _hydrate.hydration_for(components)
    if hydration.water_index is None:
        raise InvalidInputError("components", "the hydration carries no water index")
    water_index = hydration.water_index
    mixture, _ = databank.mixture_of(components, eos=eos)

    algorithm = spec["algorithm"]
    bracket = algorithm["bracket"]
    tolerance = float(algorithm["tolerance"])
    steps = max(int(bracket["steps"]), 2)

    def residual(p: float) -> float:
        flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p, "Pa"), z)
        fugacities = _guest_fugacities(flash, z, p)
        reference = _reference_water_fugacity(eos, t_si, p)
        _structure, coefficient = _hydrate.stable_structure(
            hydration.guests, fugacities, hydrate_model, t_si, p, reference
        )
        return coefficient * p / fugacities[water_index] - 1.0

    # **The scan steps in the logarithm of the pressure**, because that is the variable a
    # hydrate curve is drawn on; a linear scan would put every point of interest in its first
    # two steps.
    lower = float(bracket["lower"])
    upper = float(bracket["upper"])
    iterations = 0
    previous: tuple[float, float] | None = None
    found: tuple[float, float] | None = None
    last: tuple[float, float] | None = None
    for step in range(steps):
        p = lower * (upper / lower) ** (step / (steps - 1))
        value = residual(p)
        iterations += 1
        if previous is not None and previous[1] * value <= 0.0:
            found = (previous[0], p)
            break
        previous = (p, value)
        last = (p, value)

    if found is None:
        raise SolverNotConvergedError(
            iterations=iterations,
            residual=math.nan if last is None else last[1],
            tolerance=tolerance,
        )

    low, high = found
    low_value = residual(low)
    iterations += 1
    for _ in range(int(algorithm["max_iterations"])):
        mid = 0.5 * (low + high)
        value = residual(mid)
        iterations += 1
        if abs(value) < tolerance or high - low < 1.0e-3:
            flash = pt_flash(mixture, from_si(t_si, "K"), from_si(mid, "Pa"), z)
            fugacities = _guest_fugacities(flash, z, mid)
            reference = _reference_water_fugacity(eos, t_si, mid)
            structure, _coefficient = _hydrate.stable_structure(
                hydration.guests, fugacities, hydrate_model, t_si, mid, reference
            )
            return HydrateFormationPressureResult(
                pressure=from_si(mid, "Pa"),
                structure=(
                    HydrateStructure.STRUCTURE_I
                    if structure == 0
                    else HydrateStructure.STRUCTURE_II
                ),
                iterations=iterations,
                residual=value,
                warnings=tuple(warnings),
            )
        if low_value * value <= 0.0:
            high = mid
        else:
            low, low_value = mid, value

    raise SolverNotConvergedError(
        iterations=iterations,
        residual=0.5 * (low + high),
        tolerance=tolerance,
    )


__all__ = ["hydrate_formation_pressure"]
