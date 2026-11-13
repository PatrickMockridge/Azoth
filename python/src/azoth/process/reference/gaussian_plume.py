"""``process.gaussian_plume`` - the Briggs plume, as NeqSim's ``GaussianPlume`` computes it.

Spec: ``specs/models/process/gaussian_plume.toml``

The Python twin of ``crates/azoth-process/src/models/gaussian_plume.rs``: the rural/urban
Briggs sigma tables, the two-image concentration, and the log scan with its bisection.
``gas_dispersion`` imports the private coefficient and concentration helpers, so the
analyzer's Gaussian arm is the same arithmetic as this id.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import GaussianPlumeResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.gaussian_plume"

#: The scan's logarithmic sample count above the first point.
_SCAN_STEPS = 200
#: The scan's lower end, m.
_SCAN_MIN_DISTANCE = 1.0
#: The scan's upper end, and its answer when the profile never falls below the target, m.
_SCAN_MAX_DISTANCE = 100_000.0
#: The bisection cap.
_BISECTION_STEPS = 60
#: The bisection's stopping width, m.
_BISECTION_WIDTH = 0.5


def _stable(word: str) -> str:
    """The class's stability class, validated against its six."""
    if word not in ("A", "B", "C", "D", "E", "F"):
        raise InvalidInputError("stability", f"a Pasquill class is one of A-F, and {word!r} is not")
    return word


def _terrain(word: str) -> str:
    """The terrain category, validated against the class's two."""
    if word not in ("rural", "urban"):
        raise InvalidInputError(
            "terrain", f"a terrain category is `rural` or `urban`, and {word!r} is not"
        )
    return word


def _sigma_y(distance: float, stability: str, terrain: str) -> float:
    """``GaussianPlume.sigmaY``."""
    a, b = {
        ("A", "rural"): (0.22, 0.0001),
        ("A", "urban"): (0.32, 0.0004),
        ("B", "rural"): (0.16, 0.0001),
        ("B", "urban"): (0.32, 0.0004),
        ("C", "rural"): (0.11, 0.0001),
        ("C", "urban"): (0.22, 0.0004),
        ("D", "rural"): (0.08, 0.0001),
        ("D", "urban"): (0.16, 0.0004),
        ("E", "rural"): (0.06, 0.0001),
        ("E", "urban"): (0.11, 0.0004),
        ("F", "rural"): (0.04, 0.0001),
        ("F", "urban"): (0.11, 0.0004),
    }[(stability, terrain)]
    return a * distance / math.pow(1.0 + b * distance, 0.5)


def _sigma_z(distance: float, stability: str, terrain: str) -> float:
    """``GaussianPlume.sigmaZ``."""
    a, b, c = {
        ("A", "rural"): (0.20, 0.0, 1.0),
        ("A", "urban"): (0.24, 0.001, 0.5),
        ("B", "rural"): (0.12, 0.0, 1.0),
        ("B", "urban"): (0.24, 0.001, 0.5),
        ("C", "rural"): (0.08, 0.0002, 0.5),
        ("C", "urban"): (0.20, 0.0, 1.0),
        ("D", "rural"): (0.06, 0.0015, 0.5),
        ("D", "urban"): (0.14, 0.0003, 0.5),
        ("E", "rural"): (0.03, 0.0003, 1.0),
        ("E", "urban"): (0.08, 0.00015, 0.5),
        ("F", "rural"): (0.016, 0.0003, 1.0),
        ("F", "urban"): (0.08, 0.00015, 0.5),
    }[(stability, terrain)]
    return a * distance / math.pow(1.0 + b * distance, c)


def sigma_y(distance: float, stability: str, terrain: str) -> float:
    """The horizontal coefficient at a distance, m; the analyzer's Gaussian arm uses this."""
    return _sigma_y(distance, _stable(stability), _terrain(terrain))


def sigma_z(distance: float, stability: str, terrain: str) -> float:
    """The vertical coefficient at a distance, m."""
    return _sigma_z(distance, _stable(stability), _terrain(terrain))


def concentration(
    emission_rate: float,
    stack_height: float,
    wind_speed: float,
    stability: str,
    terrain: str,
    x: float,
    y: float,
    z: float,
) -> float:
    """The two-image Gaussian concentration, kg/m3."""
    if x <= 0.0 or emission_rate == 0.0:
        return 0.0
    sy = _sigma_y(x, _stable(stability), _terrain(terrain))
    sz = _sigma_z(x, _stable(stability), _terrain(terrain))
    h = stack_height
    term1 = emission_rate / (2.0 * math.pi * wind_speed * sy * sz)
    y_term = math.exp(-0.5 * (y * y) / (sy * sy))
    z_term1 = math.exp(-0.5 * (z - h) ** 2 / (sz * sz))
    z_term2 = math.exp(-0.5 * (z + h) ** 2 / (sz * sz))
    return term1 * y_term * (z_term1 + z_term2)


def distance_to_threshold(
    emission_rate: float,
    stack_height: float,
    wind_speed: float,
    stability: str,
    terrain: str,
    threshold: float,
) -> float | None:
    """The class's log scan and bisection; ``None`` when the profile never reaches it."""

    def centreline(x: float) -> float:
        return concentration(
            emission_rate, stack_height, wind_speed, stability, terrain, x, 0.0, 0.0
        )

    previous_distance = _SCAN_MIN_DISTANCE
    previous_concentration = centreline(previous_distance)
    above: float | None = previous_distance if previous_concentration >= threshold else None
    below: float | None = None
    for step in range(1, _SCAN_STEPS + 1):
        fraction = step / _SCAN_STEPS
        distance = _SCAN_MIN_DISTANCE * math.pow(_SCAN_MAX_DISTANCE / _SCAN_MIN_DISTANCE, fraction)
        value = centreline(distance)
        if value >= threshold:
            above = distance
        elif previous_concentration >= threshold:
            above = previous_distance
            below = distance
        previous_distance = distance
        previous_concentration = value

    if above is None:
        return None
    if below is None:
        return _SCAN_MAX_DISTANCE
    lo = above
    hi = below
    for _ in range(_BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        if centreline(mid) > threshold:
            lo = mid
        else:
            hi = mid
        if hi - lo < _BISECTION_WIDTH:
            break
    return 0.5 * (lo + hi)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)


def gaussian_plume(
    emission_rate: Q,
    effective_stack_height: Q,
    wind_speed: Q,
    stability: str,
    terrain: str,
    distance: Q,
    crosswind: Q,
    height: Q,
    threshold: Q,
) -> GaussianPlumeResult:
    """The concentration of a continuous Gaussian plume at one point.

    Args:
        emission_rate: the source strength.
        effective_stack_height: the effective release height; the class applies no plume rise.
        wind_speed: the wind speed at the 10 m reference.
        stability: the Pasquill class, ``A``-``F``.
        terrain: ``rural`` or ``urban``.
        distance: the downwind distance at which the concentration is asked.
        crosswind: the receptor's crosswind offset from the centreline.
        height: the receptor's height above the ground.
        threshold: the concentration whose downwind distance is asked.

    Returns:
        The two coefficients, the concentration, the ground centreline concentration and the
        distance to the threshold, absent when the scan never reaches it.

    Raises:
        InvalidInputError: for a value the class's constructor refuses, a non-positive
            threshold, a distance the spec's range closes, or an unknown stability or terrain
            word.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    emission = input_to_si(spec, "emission_rate", emission_rate)
    stack = input_to_si(spec, "effective_stack_height", effective_stack_height)
    wind = input_to_si(spec, "wind_speed", wind_speed)
    x = input_to_si(spec, "distance", distance)
    y = input_to_si(spec, "crosswind", crosswind)
    z = input_to_si(spec, "height", height)
    limit = input_to_si(spec, "threshold", threshold)

    apply_checks(
        checks.on_input,
        {
            "emission_rate": emission,
            "wind_speed": wind,
            "distance": x,
            "threshold": limit,
        }.get,
        warnings,
    )
    stability_word = _stable(stability)
    terrain_word = _terrain(terrain)

    distance_value = distance_to_threshold(
        emission, stack, wind, stability_word, terrain_word, limit
    )

    return GaussianPlumeResult(
        sigma_y=from_si(_sigma_y(x, stability_word, terrain_word), "m"),
        sigma_z=from_si(_sigma_z(x, stability_word, terrain_word), "m"),
        concentration=from_si(
            concentration(emission, stack, wind, stability_word, terrain_word, x, y, z),
            "kg/m**3",
        ),
        centerline_ground_concentration=from_si(
            concentration(emission, stack, wind, stability_word, terrain_word, x, 0.0, 0.0),
            "kg/m**3",
        ),
        distance_to_threshold=(None if distance_value is None else from_si(distance_value, "m")),
        warnings=tuple(warnings),
    )
