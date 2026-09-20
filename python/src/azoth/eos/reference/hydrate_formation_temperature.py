"""``eos.hydrate_formation_temperature``: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate_formation_temperature.rs``:

```text
solve  f_w^hydrate(T, P, fugacities) / f_w^fluid(T, P) - 1 = 0   for T
```

Each trial is a full flash of the fluid, because the guests' fugacities come from the gas it
produces - what NeqSim's ``setFug`` copies in before every occupancy evaluation - and the
hydrate's water fugacity is then built from those by :mod:`azoth.eos.reference._hydrate`.

The answer is settled by **water**: the hydrate and the fluid meet where water's fugacity is
the same in both. The structure is not an input but an output of the same comparison, because
the lower of the two coefficients is the stable hydrate.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import (
    HydrateFormationTemperatureResult,
    HydrateStructure,
    PtFlashResult,
)
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import _hydrate
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.hydrate_formation_temperature"


def _reference_water_fugacity(eos: str, t: float, p: float) -> float:
    """The fugacity of pure water at a state on the cubic the fluid runs, in Pa.

    NeqSim's ``refPhase`` is a one-component phase of the *host's* class - an SRK phase for
    an SRK fluid, not ``ComponentWater`` - and this builds the same thing: the same cubic,
    water alone, on the liquid root.
    """
    from azoth.eos import components as databank

    water, _ = databank.mixture_of(["water"], eos=eos)
    reduced = reduced_parameters(water, t, p)
    state = phase_state(reduced, water.kij, [1.0], liquid=True)
    return math.exp(state.ln_phi[0]) * p


def _guest_fugacities(flash: PtFlashResult, z: list[float], p: float) -> list[float]:
    """The guests' fugacities at a flash, in Pa, from the phase ``setFug`` reads.

    The **vapour** phase's, which is what NeqSim's ``setFug`` copies in. A feed the flash
    finds single-phase has no vapour to read, so its own composition is its phase
    composition, with the coefficients of whichever phase the flash says it is.
    """
    two_phase = getattr(flash, "phase", "") == "two_phase"
    composition = list(flash.y) if two_phase else list(z)
    coefficients = flash.ln_phi_liquid if flash.phase == "all_liquid" else flash.ln_phi_vapour
    return [
        fraction * math.exp(coefficient) * p
        for fraction, coefficient in zip(composition, coefficients, strict=True)
    ]


def hydrate_formation_temperature(
    components: list[str], P: Q, z: list[float], eos: str = "srk"
) -> HydrateFormationTemperatureResult:
    """The temperature at which a fluid's hydrate appears, at a pressure.

    **The components cross by name**, as they do for the Fürst phase, because the hydrate's
    guest tables are keyed by name and a :class:`~azoth.eos.mixture.Component` carries none.

    Args:
        components: the substances, by name.
        P: absolute pressure; the temperature is what is solved for.
        z: the overall mole fractions. Checked rather than renormalised.
        eos: the cubic the fluid runs, which is also the one the hydrate's reference water
            phase is built from.

    Raises:
        InvalidInputError: if the fluid has no water, or nothing in it occupies a cage.
        OutOfRangeError: if ``P`` is not positive or a trial's cavity sum has no value.
        SolverNotConvergedError: if the scan finds no sign change or the bisection hits its
            cap.
    """
    from azoth import _models_gen
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": p_si}.get, warnings)

    hydration = _hydrate.hydration_for(components)
    if hydration.water_index is None:
        raise InvalidInputError("components", "the hydration carries no water index")
    water_index = hydration.water_index
    mixture, _ = databank.mixture_of(components, eos=eos)

    algorithm = spec["algorithm"]
    bracket = algorithm["bracket"]
    tolerance = float(algorithm["tolerance"])
    steps = max(int(bracket["steps"]), 2)

    def residual(t: float) -> float:
        flash = pt_flash(mixture, from_si(t, "K"), from_si(p_si, "Pa"), z)
        fugacities = _guest_fugacities(flash, z, p_si)
        reference = _reference_water_fugacity(eos, t, p_si)
        _structure, coefficient = _hydrate.stable_structure(
            hydration.guests, fugacities, t, p_si, reference
        )
        return coefficient * p_si / fugacities[water_index] - 1.0

    iterations = 0
    previous: tuple[float, float] | None = None
    found: tuple[float, float] | None = None
    lower = float(bracket["lower"])
    upper = float(bracket["upper"])
    for step in range(steps):
        t = lower + (upper - lower) * (step / (steps - 1))
        value = residual(t)
        iterations += 1
        if previous is not None and previous[1] * value <= 0.0:
            found = (previous[0], t)
            break
        previous = (t, value)

    if found is None:
        raise SolverNotConvergedError(
            iterations=iterations,
            residual=math.nan if previous is None else previous[1],
            tolerance=tolerance,
        )

    low, high = found
    low_value = residual(low)
    iterations += 1
    for _ in range(int(algorithm["max_iterations"])):
        mid = 0.5 * (low + high)
        value = residual(mid)
        iterations += 1
        if abs(value) < tolerance or high - low < 1.0e-9:
            flash = pt_flash(mixture, from_si(mid, "K"), from_si(p_si, "Pa"), z)
            fugacities = _guest_fugacities(flash, z, p_si)
            reference = _reference_water_fugacity(eos, mid, p_si)
            structure, _coefficient = _hydrate.stable_structure(
                hydration.guests, fugacities, mid, p_si, reference
            )
            return HydrateFormationTemperatureResult(
                temperature=from_si(mid, "K"),
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


__all__ = ["hydrate_formation_temperature"]
