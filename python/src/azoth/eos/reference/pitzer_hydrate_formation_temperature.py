"""``eos.pitzer_hydrate_formation_temperature``: the Python reference.

The mirror of ``crates/azoth-eos/src/pitzer_hydrate_formation_temperature.rs``: NeqSim's
``PitzerHydrateFlash`` temperature solve. Each trial is the hybrid EoS/GE flash on a
``SystemPitzer``, and the hydrate's water fugacity coefficient is
``ComponentHydratePitzer``'s: the PVTsim exponent against the pure-water Antoine reference.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateStructure, PitzerHydrateFormationTemperatureResult
from azoth.core.units import Q, from_si, input_to_si, to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import _hydrate
from azoth.eos.reference.antoine_vapor_pressure import saturation_pressure
from azoth.eos.reference.hybrid_eos_ge_flash import hybrid_eos_ge_flash
from azoth.eos.reference.iapws_henry_law import ALIASES, TABLE

MODEL_ID = "eos.pitzer_hydrate_formation_temperature"

TEMPERATURE_STEP = 5.0
MINIMUM_TEMPERATURE = 273.15
MAXIMUM_TEMPERATURE = 323.15
MAXIMUM_PRESSURE = 1000.0e5
WATER_ACTIVITY_TOLERANCE = 1.0e-8


def _dataset(name: str | None) -> str:
    if name in (None, "auto"):
        return "auto"
    if name == "phreeqc_co2_chloride":
        return "phreeqc_co2_chloride"
    raise InvalidInputError("pitzer_dataset", f"`{name}` is not one of auto, phreeqc_co2_chloride")


def _structure(name: str | None) -> int | None:
    if name in (None, "auto"):
        return None
    if name == "si":
        return 0
    if name == "sii":
        return 1
    raise InvalidInputError("structure", f"`{name}` is not one of auto, si, sii")


def _henry_floor(
    components: list[str], moles: list[float], guests: list[_hydrate.HydrateGuest]
) -> float:
    floor = MINIMUM_TEMPERATURE
    for index, name in enumerate(components):
        if moles[index] <= 0.0 or not guests[index].former:
            continue
        gas = ALIASES.get(name.strip().lower())
        if gas is not None:
            floor = max(floor, TABLE[gas][3])
    return floor


def _trial(
    components: list[str],
    cubic: str,
    moles: list[float],
    guests: list[_hydrate.HydrateGuest],
    water_index: int,
    t: float,
    p: float,
    requested: int | None,
    dataset: str,
    warnings: list[Warning],
) -> tuple[float, int, float]:
    """One trial: the hybrid flash, the guests' aqueous fugacities, the hydrate coefficient."""
    flash = hybrid_eos_ge_flash(
        components,
        cubic,
        from_si(t, "K"),
        from_si(p, "Pa"),
        moles,
        dataset=dataset,
    )
    aqueous = len(flash.x) - 1
    composition = flash.x[aqueous]
    ln_phi = flash.ln_phi[aqueous]
    fugacity = [composition[i] * math.exp(ln_phi[i]) * p for i in range(len(components))]
    f_water = fugacity[water_index]
    if not f_water > 0.0:
        raise OutOfRangeError(
            "components",
            f_water,
            "the aqueous phase carries no water, so there is no solvent fugacity to match",
        )
    guest_fugacities = [fugacity[i] if guest.former else 0.0 for i, guest in enumerate(guests)]
    p_sat = saturation_pressure(_components.entry("water"), t, warnings)
    structure, coefficient = _hydrate.pitzer_structure(
        guests, guest_fugacities, t, p, p_sat, requested
    )
    water_activity = f_water / p_sat
    if not (water_activity > 0.0 and water_activity <= 1.0 + WATER_ACTIVITY_TOLERANCE):
        raise OutOfRangeError(
            "water_activity",
            water_activity,
            "the Pitzer aqueous phase's water activity is outside (0, 1]; downstream NeqSim "
            "refuses it as a non-physical solvent reference",
        )
    return math.log(coefficient * p / f_water), structure, water_activity


def pitzer_hydrate_formation_temperature(
    components: list[str],
    cubic: str,
    moles: list[Q],
    P: Q,
    structure: str | None = None,
    minimum_temperature: Q | None = None,
    maximum_temperature: Q | None = None,
    pitzer_dataset: str | None = None,
) -> PitzerHydrateFormationTemperatureResult:
    """The temperature at which a Pitzer brine's hydrate appears, at a pressure.

    Raises:
        InvalidInputError: for an unknown structure or dataset name, or a feed with no water
            or no hydrate former.
        OutOfRangeError: for a pressure or temperature outside the class's window, or a
            water activity outside ``(0, 1]``.
        SolverNotConvergedError: if the scan finds no sign change or the bisection hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    minimum = (
        None
        if minimum_temperature is None
        else input_to_si(spec, "minimum_temperature", minimum_temperature)
    )
    maximum = (
        None
        if maximum_temperature is None
        else input_to_si(spec, "maximum_temperature", maximum_temperature)
    )
    apply_checks(
        checks.on_input,
        {"P": p_si, "minimum_temperature": minimum, "maximum_temperature": maximum}.get,
        warnings,
    )
    requested = _structure(structure)
    dataset = _dataset(pitzer_dataset)
    if not (p_si > 0.0 and p_si <= MAXIMUM_PRESSURE):
        raise OutOfRangeError(
            "P", p_si, "the Pitzer hydrate calculation requires 0 < P <= 1000 bara"
        )

    walk = [to_si(value, "mol", "moles") for value in moles]
    if len(walk) != len(components):
        raise InvalidInputError(
            "moles", f"{len(components)} component(s) and {len(walk)} mole number(s)"
        )
    hydration = _hydrate.hydration_for(components)
    floor = _henry_floor(components, walk, hydration.guests)
    lower = max(
        MINIMUM_TEMPERATURE
        if minimum is None
        else min(max(minimum, MINIMUM_TEMPERATURE), MAXIMUM_TEMPERATURE),
        floor,
    )
    upper = (
        MAXIMUM_TEMPERATURE
        if maximum is None
        else min(max(maximum, MINIMUM_TEMPERATURE), MAXIMUM_TEMPERATURE)
    )
    if upper < floor:
        raise InvalidInputError(
            "maximum_temperature",
            f"the temperature bounds do not overlap the guests' Henry range, whose floor "
            f"is {floor} K",
        )

    tolerance = float(spec["algorithm"]["tolerance"])
    cap = int(spec["algorithm"]["max_iterations"])
    water_index = hydration.water_index
    if water_index is None:
        raise InvalidInputError(
            "components", "no water: a hydrate is water's, so a fluid without it has none"
        )

    def evaluate(t: float) -> tuple[float, int, float]:
        return _trial(
            components,
            cubic,
            walk,
            hydration.guests,
            water_index,
            t,
            p_si,
            requested,
            dataset,
            warnings,
        )

    iterations = 1
    last = evaluate(lower)
    if abs(last[0]) <= tolerance:
        return _finish(last, lower, iterations, warnings)
    current = lower
    current_residual = last[0]
    while current < upper:
        nxt = min(current + TEMPERATURE_STEP, upper)
        value = evaluate(nxt)
        iterations += 1
        if abs(value[0]) <= tolerance:
            return _finish(value, nxt, iterations, warnings)
        if math.copysign(1.0, current_residual) != math.copysign(1.0, value[0]):
            lo, hi = current, nxt
            lo_residual = current_residual
            for _ in range(cap):
                middle = 0.5 * (lo + hi)
                middle_value = evaluate(middle)
                iterations += 1
                if abs(middle_value[0]) <= tolerance:
                    return _finish(middle_value, middle, iterations, warnings)
                if math.copysign(1.0, lo_residual) != math.copysign(1.0, middle_value[0]):
                    hi = middle
                else:
                    lo = middle
                    lo_residual = middle_value[0]
            raise SolverNotConvergedError(iterations, 0.5 * (lo + hi), tolerance)
        current = nxt
        current_residual = value[0]
        last = value
    raise SolverNotConvergedError(iterations, last[0], tolerance)


def _finish(
    value: tuple[float, int, float],
    temperature: float,
    iterations: int,
    warnings: list[Warning],
) -> PitzerHydrateFormationTemperatureResult:
    residual, structure, water_activity = value
    return PitzerHydrateFormationTemperatureResult(
        temperature=from_si(temperature, "K"),
        structure=HydrateStructure.STRUCTURE_I if structure == 0 else HydrateStructure.STRUCTURE_II,
        water_activity=water_activity,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
