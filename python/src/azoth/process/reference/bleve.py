"""``process.bleve`` - the CCPS/TNO fireball.

As NeqSim's ``BLEVECalculator`` computes it. Spec: ``specs/models/process/bleve.toml``

The Python twin of ``crates/azoth-process/src/models/bleve.rs``: the fireball's diameter,
duration and surface emissive power, the capped point-source view factor and its
bisection. The transmissivity is ``process.jet_fire``'s private helper.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BleveResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference.jet_fire import _wayne_transmissivity

MODEL_ID = "process.bleve"

#: The diameter correlation's coefficient, m per kg^(1/3).
_DIAMETER_COEFFICIENT = 5.8
#: The duration correlation's coefficient, s per kg^(1/3).
_DURATION_COEFFICIENT = 0.45
#: The bisection's cap.
_BISECTION_STEPS = 60
#: The bisection's far end, m.
_SEARCH_MAX = 5000.0
#: The near bound's clearance beyond the fireball radius, m.
_NEAR_CLEARANCE = 0.1
#: The class's default relative humidity.
_DEFAULT_HUMIDITY = 0.5


def _fireball_diameter(mass: float) -> float:
    """`5.8 M^(1/3)`, m."""
    return _DIAMETER_COEFFICIENT * math.cbrt(mass)


def _fireball_duration(mass: float) -> float:
    """`0.45 M^(1/3)`, s."""
    return _DURATION_COEFFICIENT * math.cbrt(mass)


def _surface_emissive_power(
    mass: float, heat_of_combustion: float, radiative_fraction: float
) -> float:
    """`eta M dHc / (pi D^2 t)`, W/m2."""
    diameter = _fireball_diameter(mass)
    duration = _fireball_duration(mass)
    return (
        radiative_fraction * mass * heat_of_combustion / (math.pi * diameter * diameter * duration)
    )


def _incident_heat_flux(
    mass: float,
    heat_of_combustion: float,
    radiative_fraction: float,
    relative_humidity: float,
    distance: float,
) -> float:
    """The incident flux at a distance, W/m2."""
    if distance <= 0.0:
        return math.inf
    diameter = _fireball_diameter(mass)
    view = min(1.0, diameter * diameter / (4.0 * distance * distance))
    return (
        _surface_emissive_power(mass, heat_of_combustion, radiative_fraction)
        * view
        * _wayne_transmissivity(distance, relative_humidity)
    )


def _distance_to_flux(
    mass: float,
    heat_of_combustion: float,
    radiative_fraction: float,
    relative_humidity: float,
    target: float,
) -> float | None:
    """The class's bisection, absent when the fireball never reaches the target."""
    lo = _fireball_diameter(mass) / 2.0 + _NEAR_CLEARANCE
    hi = _SEARCH_MAX
    if (
        _incident_heat_flux(mass, heat_of_combustion, radiative_fraction, relative_humidity, lo)
        < target
    ):
        return None
    for _ in range(_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        if (
            _incident_heat_flux(
                mass, heat_of_combustion, radiative_fraction, relative_humidity, mid
            )
            > target
        ):
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


def bleve(
    fuel_mass: Q,
    heat_of_combustion: float,
    radiative_fraction: float,
    relative_humidity: float | None = None,
    distance: Q | None = None,
    target_flux: Q | None = None,
) -> BleveResult:
    """The fireball radiation at one distance, or the distance to a target flux.

    Args:
        fuel_mass: the mass of vaporising fuel.
        heat_of_combustion: the lower heating value in J/kg; declared dimensionless.
        radiative_fraction: the fraction of the combustion energy radiated.
        relative_humidity: the ambient humidity in ``(0, 1]``; absent is ``0.5``.
        distance: the distance from the fireball centre at which the flux is asked.
        target_flux: the flux whose distance is solved for.

    Returns:
        The fireball's diameter, duration and SEP, and with a distance the transmissivity
        and flux, and with a target the contour distance when the fireball reaches it.

    Raises:
        InvalidInputError: for a value the class's constructor refuses, a radiative
            fraction outside ``[0, 1]``, a non-positive distance or target, a humidity
            outside ``(0, 1]``, or a non-finite input.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mass = input_to_si(spec, "fuel_mass", fuel_mass)
    stated_distance = None if distance is None else input_to_si(spec, "distance", distance)
    stated_target = None if target_flux is None else input_to_si(spec, "target_flux", target_flux)
    humidity_input = None if relative_humidity is None else float(relative_humidity)
    apply_checks(
        checks.on_input,
        {
            "fuel_mass": mass,
            "heat_of_combustion": float(heat_of_combustion),
            "radiative_fraction": float(radiative_fraction),
            "relative_humidity": humidity_input,
            "distance": stated_distance,
            "target_flux": stated_target,
        }.get,
        warnings,
    )

    _require_finite("fuel_mass", mass)
    _require_finite("heat_of_combustion", float(heat_of_combustion))
    _require_finite("radiative_fraction", float(radiative_fraction))
    humidity = _DEFAULT_HUMIDITY if humidity_input is None else humidity_input
    _require_finite("relative_humidity", humidity)
    if stated_distance is not None:
        _require_finite("distance", stated_distance)
    if stated_target is not None:
        _require_finite("target_flux", stated_target)

    contour = (
        None
        if stated_target is None
        else _distance_to_flux(
            mass, float(heat_of_combustion), float(radiative_fraction), humidity, stated_target
        )
    )

    return BleveResult(
        fireball_diameter=from_si(_fireball_diameter(mass), "m"),
        fireball_duration=from_si(_fireball_duration(mass), "s"),
        surface_emissive_power=from_si(
            _surface_emissive_power(mass, float(heat_of_combustion), float(radiative_fraction)),
            "W/m**2",
        ),
        transmissivity=(
            None if stated_distance is None else _wayne_transmissivity(stated_distance, humidity)
        ),
        incident_heat_flux=(
            None
            if stated_distance is None
            else from_si(
                _incident_heat_flux(
                    mass,
                    float(heat_of_combustion),
                    float(radiative_fraction),
                    humidity,
                    stated_distance,
                ),
                "W/m**2",
            )
        ),
        distance_to_flux=None if contour is None else from_si(contour, "m"),
        warnings=tuple(warnings),
    )
