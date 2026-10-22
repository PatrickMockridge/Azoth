"""``eos.pitzer_hydrate_formation_pressure``: the Python reference.

The mirror of ``crates/azoth-eos/src/pitzer_hydrate_formation_pressure.rs``: NeqSim's
``PitzerHydrateFlash`` pressure solve, the same trial as the temperature model with the
pressure scanned at a fixed temperature.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateStructure, PitzerHydrateFormationPressureResult
from azoth.core.units import Q, from_si, input_to_si, to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import _hydrate
from azoth.eos.reference.pitzer_hydrate_formation_temperature import (
    MAXIMUM_TEMPERATURE,
    MINIMUM_TEMPERATURE,
    _dataset,
    _henry_floor,
    _structure,
    _trial,
)

MODEL_ID = "eos.pitzer_hydrate_formation_pressure"

MINIMUM_PRESSURE = 1.0e5
MAXIMUM_PRESSURE = 1000.0e5
PRESSURE_STEP = 0.5e5
PRESSURE_RATIO = 1.3


def pitzer_hydrate_formation_pressure(
    components: list[str],
    cubic: str,
    moles: list[Q],
    T: Q,
    structure: str | None = None,
    minimum_pressure: Q | None = None,
    maximum_pressure: Q | None = None,
    pitzer_dataset: str | None = None,
) -> PitzerHydrateFormationPressureResult:
    """The pressure at which a Pitzer brine's hydrate appears, at a temperature.

    Raises:
        InvalidInputError: for an unknown structure or dataset name, or a feed with no water
            or no hydrate former.
        OutOfRangeError: for a temperature outside 273.15-323.15 K or below the guests' Henry
            floor, a pressure bound that is not positive, or a water activity outside ``(0, 1]``.
        SolverNotConvergedError: if the scan finds no sign change or the bisection hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    minimum = (
        None
        if minimum_pressure is None
        else input_to_si(spec, "minimum_pressure", minimum_pressure)
    )
    maximum = (
        None
        if maximum_pressure is None
        else input_to_si(spec, "maximum_pressure", maximum_pressure)
    )
    apply_checks(
        checks.on_input,
        {"T": t_si, "minimum_pressure": minimum, "maximum_pressure": maximum}.get,
        warnings,
    )
    requested = _structure(structure)
    dataset = _dataset(pitzer_dataset)
    if not (MINIMUM_TEMPERATURE <= t_si <= MAXIMUM_TEMPERATURE):
        raise OutOfRangeError("T", t_si, "the Pitzer hydrate calculation requires 273.15-323.15 K")

    walk = [to_si(value, "mol", "moles") for value in moles]
    if len(walk) != len(components):
        raise InvalidInputError(
            "moles", f"{len(components)} component(s) and {len(walk)} mole number(s)"
        )
    hydration = _hydrate.hydration_for(components)
    floor = _henry_floor(components, walk, hydration.guests)
    if t_si < floor:
        raise OutOfRangeError(
            "T", t_si, f"the temperature is below the guests' Henry reference limit of {floor} K"
        )
    lower = min(
        max(MINIMUM_PRESSURE if minimum is None else minimum, MINIMUM_PRESSURE), MAXIMUM_PRESSURE
    )
    upper = (
        MAXIMUM_PRESSURE
        if maximum is None
        else min(max(maximum, MINIMUM_PRESSURE), MAXIMUM_PRESSURE)
    )

    tolerance = float(spec["algorithm"]["tolerance"])
    cap = int(spec["algorithm"]["max_iterations"])
    water_index = hydration.water_index
    if water_index is None:
        raise InvalidInputError(
            "components", "no water: a hydrate is water's, so a fluid without it has none"
        )

    def evaluate(p: float) -> tuple[float, int, float]:
        return _trial(
            components,
            cubic,
            walk,
            hydration.guests,
            water_index,
            t_si,
            p,
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
        nxt = min(max(current * PRESSURE_RATIO, current + PRESSURE_STEP), upper)
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
    pressure: float,
    iterations: int,
    warnings: list[Warning],
) -> PitzerHydrateFormationPressureResult:
    residual, structure, water_activity = value
    return PitzerHydrateFormationPressureResult(
        pressure=from_si(pressure, "Pa"),
        structure=HydrateStructure.STRUCTURE_I if structure == 0 else HydrateStructure.STRUCTURE_II,
        water_activity=water_activity,
        iterations=iterations,
        residual=residual,
        warnings=tuple(warnings),
    )
