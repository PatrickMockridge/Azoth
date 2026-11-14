"""``process.vce`` - the TNO Multi-Energy overpressure.

As NeqSim's ``VCEModel`` computes it. Spec: ``specs/models/process/vce.toml``

The Python twin of ``crates/azoth-process/src/models/vce.rs``: the Sachs-scaled
distance, the strength-class plateau and power-law fit, and its bisection inverse.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import VceResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.vce"

#: The ten strength-class plateaus at `Rbar <= 1`, dimensionless.
_PLATEAU = (0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 5.00, 10.0)
#: The bisection's cap.
_BISECTION_STEPS = 80
#: The bisection's near end, m.
_SEARCH_MIN = 0.1
#: The bisection's far end, m.
_SEARCH_MAX = 50000.0
#: The class's default ambient pressure, Pa.
_DEFAULT_AMBIENT_PRESSURE = 101325.0


def _combustion_energy(mass: float, heat_of_combustion: float) -> float:
    """`m dHc`, J."""
    return mass * heat_of_combustion


def _scaled_distance(distance: float, energy: float, ambient_pressure: float) -> float:
    """`d / cbrt(E / p0)`, the Sachs-scaled distance; infinite at zero energy."""
    scale = math.cbrt(energy / ambient_pressure)
    return math.inf if scale == 0.0 else distance / scale


def _dimensionless_overpressure(
    distance: float, energy: float, ambient_pressure: float, strength_class: int
) -> float:
    """The class's plateau or power-law dimensionless overpressure."""
    rbar = _scaled_distance(distance, energy, ambient_pressure)
    if rbar <= 0.0:
        return math.inf
    plateau = _PLATEAU[strength_class - 1]
    if rbar <= 1.0:
        return plateau
    beta = 1.0 if strength_class <= 5 else 1.2
    return plateau * math.pow(rbar, -beta)


def _overpressure(
    distance: float, energy: float, ambient_pressure: float, strength_class: int
) -> float:
    """The side-on overpressure at a distance, Pa."""
    return (
        _dimensionless_overpressure(distance, energy, ambient_pressure, strength_class)
        * ambient_pressure
    )


def _distance_to_overpressure(
    energy: float, ambient_pressure: float, strength_class: int, target: float
) -> float | None:
    """The class's bisection, absent when the cloud never reaches the target."""
    lo = _SEARCH_MIN
    hi = _SEARCH_MAX
    if _overpressure(lo, energy, ambient_pressure, strength_class) < target:
        return None
    for _ in range(_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        if _overpressure(mid, energy, ambient_pressure, strength_class) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _require_finite(name: str, value: float) -> None:
    """Refuse a value the arithmetic cannot carry to a published number."""
    if not math.isfinite(value):
        raise InvalidInputError(
            name, "a non-finite value is not a number this calculation publishes"
        )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def vce(
    cloud_fuel_mass: Q,
    heat_of_combustion: float,
    strength_class: float,
    ambient_pressure: Q | None = None,
    distance: Q | None = None,
    target_overpressure: Q | None = None,
) -> VceResult:
    """The TNO Multi-Energy overpressure at one distance, or its inverse.

    Args:
        cloud_fuel_mass: the mass of fuel in the flammable cloud.
        heat_of_combustion: the lower heating value in J/kg; declared dimensionless.
        strength_class: the TNO strength class, a whole number from 1 to 10.
        ambient_pressure: the pressure the scaled distance reads; absent is ``101325`` Pa.
        distance: the distance from the cloud centre at which the overpressure is asked.
        target_overpressure: the side-on overpressure whose distance is solved for.

    Returns:
        The combustion energy, and with a distance the scaled distance, the
        dimensionless overpressure and the overpressure, and with a target the contour
        distance when the cloud reaches it.

    Raises:
        InvalidInputError: for a value the class's constructor refuses, a non-positive
            distance or target, a non-integral strength class, or a non-finite input.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mass = input_to_si(spec, "cloud_fuel_mass", cloud_fuel_mass)
    stated_ambient = (
        None
        if ambient_pressure is None
        else input_to_si(spec, "ambient_pressure", ambient_pressure)
    )
    stated_distance = None if distance is None else input_to_si(spec, "distance", distance)
    stated_target = (
        None
        if target_overpressure is None
        else input_to_si(spec, "target_overpressure", target_overpressure)
    )
    apply_checks(
        checks.on_input,
        {
            "cloud_fuel_mass": mass,
            "heat_of_combustion": float(heat_of_combustion),
            "strength_class": float(strength_class),
            "ambient_pressure": stated_ambient,
            "distance": stated_distance,
            "target_overpressure": stated_target,
        }.get,
        warnings,
    )

    _require_finite("cloud_fuel_mass", mass)
    _require_finite("heat_of_combustion", float(heat_of_combustion))
    _require_finite("strength_class", float(strength_class))
    ambient = _DEFAULT_AMBIENT_PRESSURE if stated_ambient is None else stated_ambient
    _require_finite("ambient_pressure", ambient)
    if stated_distance is not None:
        _require_finite("distance", stated_distance)
    if stated_target is not None:
        _require_finite("target_overpressure", stated_target)
    if float(strength_class) != int(strength_class):
        raise InvalidInputError(
            "strength_class", "the TNO strength class is a whole number from 1 to 10"
        )

    strength = int(strength_class)
    energy = _combustion_energy(mass, float(heat_of_combustion))
    contour = (
        None
        if stated_target is None
        else _distance_to_overpressure(energy, ambient, strength, stated_target)
    )

    return VceResult(
        combustion_energy=from_si(energy, "J"),
        scaled_distance=(
            None if stated_distance is None else _scaled_distance(stated_distance, energy, ambient)
        ),
        dimensionless_overpressure=(
            None
            if stated_distance is None
            else _dimensionless_overpressure(stated_distance, energy, ambient, strength)
        ),
        overpressure=(
            None
            if stated_distance is None
            else from_si(_overpressure(stated_distance, energy, ambient, strength), "Pa")
        ),
        distance_to_overpressure=None if contour is None else from_si(contour, "m"),
        warnings=tuple(warnings),
    )
