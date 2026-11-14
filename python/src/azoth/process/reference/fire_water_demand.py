"""``process.fire_water_demand`` - the fire-water demand.

As NeqSim's ``FireWaterDemandCalculator`` computes it. Spec:
``specs/models/process/fire_water_demand.toml``

The Python twin of ``crates/azoth-process/src/models/fire_water_demand.rs``: the area and
object demand, the simultaneous factor, the hourly rate and the water and foam volumes.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FireWaterDemandResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.fire_water_demand"

#: The class's NORSOK process-area rate, L/(min m2).
_DEFAULT_APPLICATION_RATE_LPM_M2 = 10.0
#: The class's default duration, s.
_DEFAULT_DURATION_S = 1800.0
#: The class's default foam concentrate percent.
_DEFAULT_FOAM_CONCENTRATE_PERCENT = 0.0
#: The class's default simultaneous area factor.
_DEFAULT_SIMULTANEOUS_AREA_FACTOR = 1.0


class _Solution(NamedTuple):
    """Every number the calculation produces, in SI, before the units are attached."""

    area_demand: float
    object_demand: float
    total_demand: float
    total_demand_m3_per_hour: float
    water_volume: float
    foam_concentrate_volume: float
    ratio: float | None


def _require_finite(name: str, value: float) -> None:
    """Refuse a value the arithmetic cannot carry to a published number."""
    if not math.isfinite(value):
        raise InvalidInputError(
            name, "a non-finite value is not a number this calculation publishes"
        )


def _spec() -> dict[str, Any]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def _evaluate(
    protected_area: Q,
    application_rate: float | None = None,
    duration: Q | None = None,
    foam_concentrate_percent: float | None = None,
    simultaneous_area_factor: float | None = None,
    protected_objects: list[dict[str, object]] | None = None,
) -> tuple[_Solution, list[Warning]]:
    """The class's arithmetic over SI values; shared with the layer dumper."""
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    area = input_to_si(spec, "protected_area", protected_area)
    stated_rate = None if application_rate is None else float(application_rate)
    stated_duration = None if duration is None else input_to_si(spec, "duration", duration)
    stated_foam = None if foam_concentrate_percent is None else float(foam_concentrate_percent)
    stated_factor = (
        None if simultaneous_area_factor is None else float(simultaneous_area_factor)
    )
    rows = input_records_to_si(spec, "protected_objects", protected_objects or [])
    apply_checks(
        checks.on_input,
        {
            "protected_area": area,
            "application_rate": stated_rate,
            "duration": stated_duration,
            "foam_concentrate_percent": stated_foam,
            "simultaneous_area_factor": stated_factor,
        }.get,
        warnings,
    )

    _require_finite("protected_area", area)
    rate = _DEFAULT_APPLICATION_RATE_LPM_M2 if stated_rate is None else stated_rate
    _require_finite("application_rate", rate)
    duration_s = _DEFAULT_DURATION_S if stated_duration is None else stated_duration
    _require_finite("duration", duration_s)
    foam = _DEFAULT_FOAM_CONCENTRATE_PERCENT if stated_foam is None else stated_foam
    _require_finite("foam_concentrate_percent", foam)
    factor = _DEFAULT_SIMULTANEOUS_AREA_FACTOR if stated_factor is None else stated_factor
    _require_finite("simultaneous_area_factor", factor)

    object_demand = 0.0
    for row in rows:
        surface = float(row["surface_area"])
        row_rate = float(row["application_rate"])
        _require_finite("protected_objects.surface_area", surface)
        _require_finite("protected_objects.application_rate", row_rate)
        if surface < 0.0 or row_rate < 0.0:
            raise InvalidInputError(
                "protected_objects",
                "a negative surface area or application rate is not an object the demand can "
                "cover; the class's `addObject` accepts it and publishes a negative demand",
            )
        object_demand += surface * row_rate

    area_demand = area * rate
    total = (area_demand + object_demand) * factor
    m3_per_hour = total * 60.0 / 1000.0
    water = total * duration_s / 60000.0
    foam_volume = water * foam / 100.0
    ratio = None if area_demand <= 0.0 else object_demand / area_demand

    # The class works in L/min and m3/h; the solution is in SI base magnitudes, as every
    # reference's is, and `from_si` re-expresses them in the spec's own units.
    return (
        _Solution(
            area_demand=area_demand / 60000.0,
            object_demand=object_demand / 60000.0,
            total_demand=total / 60000.0,
            total_demand_m3_per_hour=m3_per_hour / 3600.0,
            water_volume=water,
            foam_concentrate_volume=foam_volume,
            ratio=ratio,
        ),
        warnings,
    )


def fire_water_demand(
    protected_area: Q,
    application_rate: float | None = None,
    duration: Q | None = None,
    foam_concentrate_percent: float | None = None,
    simultaneous_area_factor: float | None = None,
    protected_objects: list[dict[str, object]] | None = None,
) -> FireWaterDemandResult:
    """The area and object demand, the simultaneous factor and the volumes.

    Args:
        protected_area: the protected deck or ground area.
        application_rate: the general-area rate in ``L/(min m2)``; absent is ``10``.
        duration: the fire-water duration; absent is ``1800 s``, the class's 30 minutes.
        foam_concentrate_percent: the foam concentrate percent in ``[0, 10]``; absent is 0.
        simultaneous_area_factor: the factor applied to the summed demand; absent is 1.
        protected_objects: the dedicated objects, each ``{surface_area, application_rate}``.

    Returns:
        The area and object demand, the total, the hourly rate and the water and foam
        volumes.

    Raises:
        InvalidInputError: for a negative area or rate, a non-positive duration, a foam
            percent outside ``[0, 10]``, a factor below one, a negative row area or rate,
            or a non-finite input.
    """
    solution, warnings = _evaluate(
        protected_area,
        application_rate,
        duration,
        foam_concentrate_percent,
        simultaneous_area_factor,
        protected_objects,
    )
    return FireWaterDemandResult(
        area_demand=from_si(solution.area_demand, "L/min"),
        object_demand=from_si(solution.object_demand, "L/min"),
        total_demand=from_si(solution.total_demand, "L/min"),
        total_demand_m3_per_hour=from_si(solution.total_demand_m3_per_hour, "m**3/h"),
        water_volume=from_si(solution.water_volume, "m**3"),
        foam_concentrate_volume=from_si(solution.foam_concentrate_volume, "m**3"),
        object_to_area_demand_ratio=solution.ratio,
        warnings=tuple(warnings),
    )
