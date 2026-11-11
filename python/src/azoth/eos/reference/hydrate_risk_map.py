"""``eos.hydrate_risk_map`` - NeqSim's ``HydrateRiskMapper``.

The hydrate formation temperature at each profile point, its subcooling against the
actual temperature, and the risk band the two decide.

Spec: ``specs/models/eos/hydrate_risk_map.toml``.
"""

from __future__ import annotations

from typing import Any, cast

from azoth.core.errors import AzothError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrateRiskLevel, HydrateRiskMapResult, HydrateRiskPoint
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import entry
from azoth.eos.reference.hydrate_formation_temperature import hydrate_formation_temperature

CALC_ID = "eos.hydrate_risk_map"

DEFAULT_HIGH_RISK_SUBCOOLING = 3.0
DEFAULT_MEDIUM_RISK_SUBCOOLING = 6.0


def _rank(risk: HydrateRiskLevel) -> int:
    return {
        HydrateRiskLevel.CRITICAL: 0,
        HydrateRiskLevel.HIGH: 1,
        HydrateRiskLevel.MEDIUM: 2,
        HydrateRiskLevel.LOW: 3,
        HydrateRiskLevel.UNKNOWN: 4,
    }[risk]


def hydrate_risk_map(
    components: list[str],
    z: list[float],
    distances: list[Q],
    pressures: list[Q],
    temperatures: list[Q],
    eos: str | None = None,
    hydrate_model: str | None = None,
    high_risk_subcooling: Q | None = None,
    medium_risk_subcooling: Q | None = None,
) -> HydrateRiskMapResult:
    """The hydrate risk along a profile.

    Args:
        components: the fluid's substances, by name.
        z: the overall mole fractions.
        distances: each profile point's distance from the inlet.
        pressures: each point's pressure.
        temperatures: each point's actual temperature.
        eos: the cubic the hydrate reference water phase is built from; absent takes srk.
        hydrate_model: which fitted component model fills the cages; absent takes pvtsim.
        high_risk_subcooling: the high-risk threshold; absent takes 3 K.
        medium_risk_subcooling: the medium-risk threshold; absent takes 6 K.

    Returns:
        One row per point, the worst band, the minimum subcooling and the counts.

    Raises:
        InvalidInputError: if the point vectors differ in length, or a name is unknown.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = hydrate_risk_map(
        ...     ["methane", "ethane", "propane", "water"], [0.82, 0.06, 0.03, 0.09],
        ...     [q(0.0, "m")], [q(1.0e7, "Pa")], [q(300.0, "K")], "srk", "pvtsim")
        >>> r.total_points
        1
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if not (len(distances) == len(pressures) == len(temperatures)):
        from azoth.core.errors import InvalidInputError

        raise InvalidInputError(
            "distances",
            "the profile needs one distance, pressure and temperature per point; got "
            f"{len(distances)}, {len(pressures)} and {len(temperatures)}",
        )
    for name in components:
        entry(name)

    values: dict[str, float] = {
        "high_risk_subcooling": (
            from_si(DEFAULT_HIGH_RISK_SUBCOOLING, "K").to("K").magnitude
            if high_risk_subcooling is None
            else input_to_si(spec, "high_risk_subcooling", high_risk_subcooling)
        ),
        "medium_risk_subcooling": (
            from_si(DEFAULT_MEDIUM_RISK_SUBCOOLING, "K").to("K").magnitude
            if medium_risk_subcooling is None
            else input_to_si(spec, "medium_risk_subcooling", medium_risk_subcooling)
        ),
    }
    apply_checks(checks.on_input, values.get, warnings)

    cubic = eos or "srk"
    hydrate_component = hydrate_model or "pvtsim"

    rows: list[HydrateRiskPoint] = []
    critical_count = 0
    unknown_count = 0
    minimum_subcooling: float | None = None
    worst = HydrateRiskLevel.LOW

    for distance, pressure, temperature in zip(distances, pressures, temperatures, strict=True):
        pressure_si = input_to_si(spec, "pressures", pressure)
        actual = input_to_si(spec, "temperatures", temperature)
        try:
            result = hydrate_formation_temperature(
                components, from_si(pressure_si, "Pa"), z, cubic, hydrate_component
            )
            hydrate_temperature = result.temperature.to("K").magnitude
            subcooling = actual - hydrate_temperature
            if subcooling < 0.0:
                risk = HydrateRiskLevel.CRITICAL
            elif subcooling < values["high_risk_subcooling"]:
                risk = HydrateRiskLevel.HIGH
            elif subcooling < values["medium_risk_subcooling"]:
                risk = HydrateRiskLevel.MEDIUM
            else:
                risk = HydrateRiskLevel.LOW
        except AzothError:
            hydrate_temperature = None
            subcooling = None
            risk = HydrateRiskLevel.UNKNOWN

        if risk is HydrateRiskLevel.CRITICAL:
            critical_count += 1
        if risk is HydrateRiskLevel.UNKNOWN:
            unknown_count += 1
        else:
            known_subcooling = cast(float, subcooling)
            if minimum_subcooling is None or known_subcooling < minimum_subcooling:
                minimum_subcooling = known_subcooling
            if _rank(risk) < _rank(worst):
                worst = risk

        rows.append(
            HydrateRiskPoint(
                distance=from_si(input_to_si(spec, "distances", distance), "m"),
                pressure=from_si(pressure_si, "Pa"),
                actual_temperature=from_si(actual, "K"),
                hydrate_temperature=(
                    None if hydrate_temperature is None else from_si(hydrate_temperature, "K")
                ),
                subcooling=None if subcooling is None else from_si(subcooling, "K"),
                risk_level=risk.value,
            )
        )

    overall = (
        HydrateRiskLevel.UNKNOWN
        if unknown_count > 0 and worst is not HydrateRiskLevel.CRITICAL
        else worst
    )

    return HydrateRiskMapResult(
        points=tuple(rows),
        overall_risk=overall,
        minimum_subcooling=(
            None if minimum_subcooling is None else from_si(minimum_subcooling, "K")
        ),
        critical_point_count=critical_count,
        unknown_point_count=unknown_count,
        total_points=len(rows),
        all_points_evaluated=unknown_count == 0,
        warnings=tuple(warnings),
    )
