"""``process.pool_fire`` - the API 521 / Mudan pool fire.

As NeqSim's ``PoolFireModel`` computes it. Spec: ``specs/models/process/pool_fire.toml``

The Python twin of ``crates/azoth-process/src/models/pool_fire.rs``: Thomas's flame
height, Mudan's smoke-corrected SEP, the vertical-cylinder view factor and its bisection.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PoolFireResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.pool_fire"

#: Thomas's gravitational constant, m/s2.
_GRAVITY = 9.81
#: Mudan's clean-flame emissive power, W/m2.
_SEP_MAX = 140_000.0
#: Mudan's smoke-blocked emissive power, W/m2.
_SEP_SMOKE = 20_000.0
#: Mudan's smoke-decay constant, 1/m.
_SEP_DECAY = 0.12
#: The class's default air density, kg/m3.
_DEFAULT_AIR_DENSITY = 1.2
#: The bisection's cap.
_BISECTION_STEPS = 60
#: The bisection's far end, m.
_SEARCH_MAX = 1000.0
#: The near bound's clearance beyond the pool radius, m.
_NEAR_CLEARANCE = 0.1
#: The view factor's clearance inside the pool radius, m.
_VIEW_CLEARANCE = 0.001


def _total_burn_rate(diameter: float, burn_rate: float) -> float:
    """`m'' pi D^2 / 4`, kg/s."""
    return burn_rate * 0.25 * math.pi * diameter * diameter


def _flame_height(diameter: float, burn_rate: float, air_density: float) -> float:
    """Thomas's no-wind flame height, m."""
    term = burn_rate / (air_density * math.sqrt(_GRAVITY * diameter))
    return diameter * 42.0 * math.pow(term, 0.61)


def _surface_emissive_power(diameter: float) -> float:
    """Mudan's smoke-corrected SEP, W/m2."""
    clean = math.exp(-_SEP_DECAY * diameter)
    return _SEP_MAX * clean + _SEP_SMOKE * (1.0 - clean)


def _view_factor(diameter: float, burn_rate: float, air_density: float, distance: float) -> float:
    """Mudan's vertical-cylinder view factor at a ground-level receiver."""
    radius = diameter / 2.0
    height = _flame_height(diameter, burn_rate, air_density) / radius
    s = max(distance, radius + _VIEW_CLEARANCE) / radius
    a = (height * height + s * s + 1.0) / (2.0 * s)
    f = (1.0 / math.pi) * (
        math.atan(height / math.sqrt(s * s - 1.0))
        - (a - (s * s + 1.0) / (2.0 * s)) / math.sqrt(a * a - 1.0)
    )
    return min(1.0, max(0.0, f))


def _incident_heat_flux(
    diameter: float, burn_rate: float, air_density: float, distance: float
) -> float:
    """The incident flux at a distance, W/m2."""
    return _surface_emissive_power(diameter) * _view_factor(
        diameter, burn_rate, air_density, distance
    )


def _distance_to_flux(
    diameter: float, burn_rate: float, air_density: float, target: float
) -> float | None:
    """The class's bisection, absent when the fire never reaches the target."""
    lo = diameter / 2.0 + _NEAR_CLEARANCE
    hi = _SEARCH_MAX
    if _incident_heat_flux(diameter, burn_rate, air_density, lo) < target:
        return None
    for _ in range(_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        if _incident_heat_flux(diameter, burn_rate, air_density, mid) > target:
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


def pool_fire(
    pool_diameter: Q,
    burn_rate: Q,
    air_density: float | None = None,
    distance: Q | None = None,
    target_flux: Q | None = None,
) -> PoolFireResult:
    """The pool-fire radiation at one distance, or the distance to a target flux.

    Args:
        pool_diameter: the equivalent circular pool diameter.
        burn_rate: the mass-burning rate per unit area.
        air_density: the ambient air density; absent is ``1.2`` kg/m3.
        distance: the horizontal distance from the pool centre.
        target_flux: the flux whose distance is solved for.

    Returns:
        The burn rate, flame height and SEP, and with a distance the view factor and
        flux, and with a target the contour distance when the fire reaches it.

    Raises:
        InvalidInputError: for a value the class's constructor refuses, a non-positive
            target, or a non-finite input.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    diameter = input_to_si(spec, "pool_diameter", pool_diameter)
    burning = input_to_si(spec, "burn_rate", burn_rate)
    stated_distance = None if distance is None else input_to_si(spec, "distance", distance)
    stated_target = None if target_flux is None else input_to_si(spec, "target_flux", target_flux)
    air_input = None if air_density is None else float(air_density)
    apply_checks(
        checks.on_input,
        {
            "pool_diameter": diameter,
            "burn_rate": burning,
            "air_density": air_input,
            "distance": stated_distance,
            "target_flux": stated_target,
        }.get,
        warnings,
    )

    _require_finite("pool_diameter", diameter)
    _require_finite("burn_rate", burning)
    air = _DEFAULT_AIR_DENSITY if air_input is None else air_input
    _require_finite("air_density", air)
    if stated_distance is not None:
        _require_finite("distance", stated_distance)
    if stated_target is not None:
        _require_finite("target_flux", stated_target)

    contour = (
        None if stated_target is None else _distance_to_flux(diameter, burning, air, stated_target)
    )

    return PoolFireResult(
        total_burn_rate=from_si(_total_burn_rate(diameter, burning), "kg/s"),
        flame_height=from_si(_flame_height(diameter, burning, air), "m"),
        surface_emissive_power=from_si(_surface_emissive_power(diameter), "W/m**2"),
        view_factor=(
            None
            if stated_distance is None
            else _view_factor(diameter, burning, air, stated_distance)
        ),
        incident_heat_flux=(
            None
            if stated_distance is None
            else from_si(_incident_heat_flux(diameter, burning, air, stated_distance), "W/m**2")
        ),
        distance_to_flux=None if contour is None else from_si(contour, "m"),
        warnings=tuple(warnings),
    )
