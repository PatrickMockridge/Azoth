"""``process.jet_fire`` - the API 521 point-source jet fire.

As NeqSim's ``JetFireModel`` computes it. Spec: ``specs/models/process/jet_fire.toml``

The Python twin of ``crates/azoth-process/src/models/jet_fire.rs``: the radiated power,
Wayne's transmissivity, the incident flux and its fixed-point inverse. ``process.bleve``
imports the private transmissivity helper.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import JetFireResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.jet_fire"

#: The fixed point's iteration cap, the class's own.
_INVERSE_STEPS = 25
#: The fixed point's stopping width, m.
_INVERSE_WIDTH = 0.01
#: The class's default relative humidity.
_DEFAULT_HUMIDITY = 0.5


def _wayne_transmissivity(distance: float, relative_humidity: float) -> float:
    """Wayne's simplified atmospheric transmissivity; shared with ``process.bleve``."""
    water = relative_humidity * 1700.0
    xpw = water * max(distance, 1.0)
    log = math.log10(xpw)
    return min(1.0, max(0.0, 1.006 - 0.01171 * log - 0.02368 * log * log))


def _distance_to_flux(power: float, relative_humidity: float, target: float) -> float:
    """The class's fixed-point inverse of the point-source flux."""
    denominator = 4.0 * math.pi * target
    radius = math.sqrt(power / denominator)
    for _ in range(_INVERSE_STEPS):
        tau = _wayne_transmissivity(radius, relative_humidity)
        following = math.sqrt(power * tau / denominator)
        if abs(following - radius) < _INVERSE_WIDTH:
            radius = following
            break
        radius = following
    return radius


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


def jet_fire(
    mass_flow: Q,
    heat_of_combustion: float,
    radiative_fraction: float,
    relative_humidity: float | None = None,
    distance: Q | None = None,
    target_flux: Q | None = None,
) -> JetFireResult:
    """The point-source jet-fire radiation at one distance, or the distance to a target flux.

    Args:
        mass_flow: the fuel mass flow burning at the jet.
        heat_of_combustion: the lower heating value in J/kg; declared dimensionless.
        radiative_fraction: the fraction of the combustion heat radiated, in ``[0, 1]``.
        relative_humidity: the ambient humidity in ``(0, 1]``; absent is ``0.5``.
        distance: the distance from the flame centre at which the flux is asked.
        target_flux: the flux whose distance is solved for.

    Returns:
        The radiated power, and with a distance the transmissivity and flux, and with a
        target the fixed-point distance.

    Raises:
        InvalidInputError: for a value the class's constructor refuses, a non-positive
            distance or target, a humidity outside ``(0, 1]``, or a non-finite input.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    mass = input_to_si(spec, "mass_flow", mass_flow)
    stated_distance = None if distance is None else input_to_si(spec, "distance", distance)
    stated_target = None if target_flux is None else input_to_si(spec, "target_flux", target_flux)
    humidity_input = None if relative_humidity is None else float(relative_humidity)
    apply_checks(
        checks.on_input,
        {
            "mass_flow": mass,
            "heat_of_combustion": float(heat_of_combustion),
            "radiative_fraction": float(radiative_fraction),
            "relative_humidity": humidity_input,
            "distance": stated_distance,
            "target_flux": stated_target,
        }.get,
        warnings,
    )

    _require_finite("mass_flow", mass)
    _require_finite("heat_of_combustion", float(heat_of_combustion))
    _require_finite("radiative_fraction", float(radiative_fraction))
    humidity = _DEFAULT_HUMIDITY if humidity_input is None else humidity_input
    _require_finite("relative_humidity", humidity)
    if stated_distance is not None:
        _require_finite("distance", stated_distance)
    if stated_target is not None:
        _require_finite("target_flux", stated_target)

    power = float(radiative_fraction) * mass * float(heat_of_combustion)

    return JetFireResult(
        total_radiative_power=from_si(power, "W"),
        transmissivity=(
            None if stated_distance is None else _wayne_transmissivity(stated_distance, humidity)
        ),
        incident_heat_flux=(
            None
            if stated_distance is None
            else from_si(
                power
                * _wayne_transmissivity(stated_distance, humidity)
                / (4.0 * math.pi * stated_distance * stated_distance),
                "W/m**2",
            )
        ),
        distance_to_flux=(
            None
            if stated_target is None
            else from_si(_distance_to_flux(power, humidity, stated_target), "m")
        ),
        warnings=tuple(warnings),
    )
