"""``process.pfp_demand`` - the passive fire protection demand.

As NeqSim's ``PfpDemandCalculator`` computes it. Spec: ``specs/models/process/pfp_demand.toml``

The Python twin of ``crates/azoth-process/src/models/pfp_demand.rs``: the lumped bare-steel
time to the critical temperature, the strict PFP test, the conduction-limited thickness and
the H/J rating.
"""

from __future__ import annotations

import math
from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PfpDemandResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "process.pfp_demand"

#: Carbon steel density, kg/m3.
_STEEL_DENSITY = 7850.0
#: Carbon steel specific heat, J/(kg K).
_STEEL_CP = 490.0
#: The class's defaults.
_DEFAULT_ABSORPTIVITY = 0.9
_DEFAULT_INITIAL_TEMPERATURE_K = 288.15
_DEFAULT_CRITICAL_TEMPERATURE_K = 673.15
_DEFAULT_FIRE_TEMPERATURE_K = 1173.15
_DEFAULT_PFP_CONDUCTIVITY = 0.2
#: The rating's 60-minute boundary, s; exactly 3600 s is short.
_LONG_DURATION_S = 3600.0
_LONG_DURATION_EPSILON_S = 1.0e-6
#: The floor on the fire-to-steel driving temperature, K.
_MIN_DRIVING_TEMPERATURE_K = 1.0


class _Solution(NamedTuple):
    """Every number the calculation produces, in SI, before the units are attached."""

    mass_per_area: float
    bare: float
    survival: float
    required: bool
    thickness: float
    rating: str
    q_abs: float
    q_allow: float | None
    steel_average: float | None
    driving: float | None


def _fire_type(value: str | None) -> str:
    """The stated fire type, the class's ``pool`` default when absent."""
    if value is None or value == "pool":
        return "pool"
    if value == "jet":
        return "jet"
    raise InvalidInputError("fire_type", f"the fire type is `pool` or `jet`, and {value!r} is not")


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


def _evaluate(
    incident_heat_flux: Q,
    steel_thickness: Q,
    required_survival_time: Q,
    absorptivity: float | None = None,
    initial_temperature: Q | None = None,
    critical_temperature: Q | None = None,
    fire_temperature: Q | None = None,
    pfp_conductivity: Q | None = None,
    fire_type: str | None = None,
) -> tuple[_Solution, list[Warning]]:
    """The class's arithmetic over SI values; shared with the layer dumper."""
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    flux = input_to_si(spec, "incident_heat_flux", incident_heat_flux)
    thickness = input_to_si(spec, "steel_thickness", steel_thickness)
    survival = input_to_si(spec, "required_survival_time", required_survival_time)
    stated_alpha = None if absorptivity is None else float(absorptivity)
    stated_initial = (
        None
        if initial_temperature is None
        else input_to_si(spec, "initial_temperature", initial_temperature)
    )
    stated_critical = (
        None
        if critical_temperature is None
        else input_to_si(spec, "critical_temperature", critical_temperature)
    )
    stated_fire = (
        None
        if fire_temperature is None
        else input_to_si(spec, "fire_temperature", fire_temperature)
    )
    stated_conductivity = (
        None
        if pfp_conductivity is None
        else input_to_si(spec, "pfp_conductivity", pfp_conductivity)
    )
    apply_checks(
        checks.on_input,
        {
            "incident_heat_flux": flux,
            "steel_thickness": thickness,
            "required_survival_time": survival,
            "absorptivity": stated_alpha,
            "initial_temperature": stated_initial,
            "critical_temperature": stated_critical,
            "fire_temperature": stated_fire,
            "pfp_conductivity": stated_conductivity,
        }.get,
        warnings,
    )

    _require_finite("incident_heat_flux", flux)
    _require_finite("steel_thickness", thickness)
    _require_finite("required_survival_time", survival)
    alpha = _DEFAULT_ABSORPTIVITY if stated_alpha is None else stated_alpha
    _require_finite("absorptivity", alpha)
    initial = _DEFAULT_INITIAL_TEMPERATURE_K if stated_initial is None else stated_initial
    _require_finite("initial_temperature", initial)
    critical = _DEFAULT_CRITICAL_TEMPERATURE_K if stated_critical is None else stated_critical
    _require_finite("critical_temperature", critical)
    fire = _DEFAULT_FIRE_TEMPERATURE_K if stated_fire is None else stated_fire
    _require_finite("fire_temperature", fire)
    conductivity = _DEFAULT_PFP_CONDUCTIVITY if stated_conductivity is None else stated_conductivity
    _require_finite("pfp_conductivity", conductivity)
    family = _fire_type(fire_type)

    if critical <= initial:
        raise InvalidInputError(
            "critical_temperature",
            "the critical temperature must exceed the initial one; at or below it the class "
            "publishes a zero or negative bare-steel time with an infinite or negative "
            "thickness",
        )

    mass_per_area = _STEEL_DENSITY * thickness
    delta_t = critical - initial
    q_abs = alpha * flux
    bare = mass_per_area * _STEEL_CP * delta_t / q_abs
    required = bare < survival
    if required:
        q_allow_value = mass_per_area * _STEEL_CP * delta_t / survival
        steel_average_value = 0.5 * (initial + critical)
        driving_value = max(_MIN_DRIVING_TEMPERATURE_K, fire - steel_average_value)
        required_thickness = conductivity * driving_value / q_allow_value
        q_allow: float | None = q_allow_value
        steel_average: float | None = steel_average_value
        driving: float | None = driving_value
    else:
        q_allow = None
        steel_average = None
        driving = None
        required_thickness = 0.0
    if not required:
        rating = "none"
    elif survival > _LONG_DURATION_S + _LONG_DURATION_EPSILON_S:
        rating = "j120" if family == "jet" else "h120"
    else:
        rating = "j60" if family == "jet" else "h60"

    return (
        _Solution(
            mass_per_area=mass_per_area,
            bare=bare,
            survival=survival,
            required=required,
            thickness=required_thickness,
            rating=rating,
            q_abs=q_abs,
            q_allow=q_allow,
            steel_average=steel_average,
            driving=driving,
        ),
        warnings,
    )


def pfp_demand(
    incident_heat_flux: Q,
    steel_thickness: Q,
    required_survival_time: Q,
    absorptivity: float | None = None,
    initial_temperature: Q | None = None,
    critical_temperature: Q | None = None,
    fire_temperature: Q | None = None,
    pfp_conductivity: Q | None = None,
    fire_type: str | None = None,
) -> PfpDemandResult:
    """The bare-steel time, the strict verdict, the thickness and the rating.

    Args:
        incident_heat_flux: the incident heat flux on the exposed steel.
        steel_thickness: the exposed steel wall thickness.
        required_survival_time: the time the wall must survive the fire.
        absorptivity: the steel surface absorptivity in ``(0, 1]``; absent is ``0.9``.
        initial_temperature: the wall's initial temperature; absent is ``288.15 K``.
        critical_temperature: the critical steel temperature; absent is ``673.15 K``.
        fire_temperature: the effective fire temperature; absent is ``1173.15 K``.
        pfp_conductivity: the protection's conductivity; absent is ``0.2 W/(m K)``.
        fire_type: ``pool`` (H-class) or ``jet`` (J-class); absent is ``pool``.

    Returns:
        The steel mass per area, the bare time, the echoed survival time, the verdict, the
        thickness in metres and millimetres, and the rating.

    Raises:
        InvalidInputError: for a non-positive flux, thickness or survival time, an
            absorptivity outside ``(0, 1]``, a non-positive conductivity or temperature, a
            critical temperature at or below the initial one, an unknown fire type, or a
            non-finite input.
    """
    solution, warnings = _evaluate(
        incident_heat_flux,
        steel_thickness,
        required_survival_time,
        absorptivity,
        initial_temperature,
        critical_temperature,
        fire_temperature,
        pfp_conductivity,
        fire_type,
    )
    return PfpDemandResult(
        steel_mass_per_area=solution.mass_per_area,
        bare_steel_time_to_critical=from_si(solution.bare, "s"),
        required_survival_time=from_si(solution.survival, "s"),
        pfp_required=solution.required,
        required_pfp_thickness=from_si(solution.thickness, "m"),
        required_pfp_thickness_mm=from_si(solution.thickness, "mm"),
        rating=solution.rating,
        warnings=tuple(warnings),
    )
