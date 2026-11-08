"""``thermal.zone_rating`` - the shell-and-tube rating run zone by zone.

Spec: ``specs/calcs/thermal/zone_rating.toml``

The Python twin of ``crates/azoth-thermal/src/zone_rating.rs``, and both mirror
``ThermalDesignCalculator.calculateZones``. Each zone overrides the eight transport properties
and reruns the rating's permissive arithmetic; mass flows, fouling, the deposit, the method and
every geometry stay the base request's. **The permissive path's unanswered states are refused**:
a zone whose transport properties are not positive finite is an error, not a zero or NaN
coefficient.
"""

from __future__ import annotations

import math
from typing import Any

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import ZoneRatingRow
from azoth.core.result import ZoneRatingResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.thermal.reference.shell_and_tube_rating import (
    BELL_DELAWARE,
    KERN,
    _evaluate_si,
    _validate_strict,
)

CALC_ID = "thermal.zone_rating"

#: The base's own transport properties. No zone reads them; they are the class's defaults and
#: exist only so the shared strict validation sees a well-formed request.
_BASE_DEFAULTS = {
    "tube_density": 800.0,
    "tube_viscosity": 5e-4,
    "tube_cp": 2500.0,
    "tube_conductivity": 0.13,
    "shell_density": 1.5,
    "shell_viscosity": 1.5e-5,
    "shell_cp": 1100.0,
    "shell_conductivity": 0.025,
}


def zone_rating(
    tube_outside_diameter: Q,
    tube_inside_diameter: Q,
    tube_length: Q,
    tube_count: float,
    tube_passes: float,
    tube_pitch: Q,
    triangular_pitch: bool,
    shell_inside_diameter: Q,
    baffle_spacing: Q,
    baffle_count: float,
    baffle_cut: float,
    tube_wall_conductivity: Q,
    tube_mass_flow: Q,
    shell_mass_flow: Q,
    fouling_tube: float,
    fouling_shell: float,
    shell_side_method: str,
    tube_fouling_layer_thickness: Q | None = None,
    tube_fouling_layer_conductivity: Q | None = None,
    shell_viscosity_wall: Q | None = None,
    has_sealing: bool | None = None,
    sealing_pairs: float | None = None,
    tube_to_baffle_clearance: Q | None = None,
    shell_to_baffle_clearance: Q | None = None,
    bypass_area: Q | None = None,
    *,
    zones: list[dict[str, Any]],
) -> ZoneRatingResult:
    """Rates every zone in the caller's order.

    Raises:
        InvalidInputError: for a base state the strict rating refuses and for a zone whose
            stated transport properties are not positive finite, or whose duty or LMTD is not
            finite.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    raw: dict[str, float | None] = {
        "tube_outside_diameter": input_to_si(spec, "tube_outside_diameter", tube_outside_diameter),
        "tube_inside_diameter": input_to_si(spec, "tube_inside_diameter", tube_inside_diameter),
        "tube_length": input_to_si(spec, "tube_length", tube_length),
        "tube_count": tube_count,
        "tube_passes": tube_passes,
        "tube_pitch": input_to_si(spec, "tube_pitch", tube_pitch),
        "shell_inside_diameter": input_to_si(spec, "shell_inside_diameter", shell_inside_diameter),
        "baffle_spacing": input_to_si(spec, "baffle_spacing", baffle_spacing),
        "baffle_count": baffle_count,
        "baffle_cut": baffle_cut,
        "tube_wall_conductivity": input_to_si(
            spec, "tube_wall_conductivity", tube_wall_conductivity
        ),
        "tube_mass_flow": input_to_si(spec, "tube_mass_flow", tube_mass_flow),
        "shell_mass_flow": input_to_si(spec, "shell_mass_flow", shell_mass_flow),
        "fouling_tube": fouling_tube,
        "fouling_shell": fouling_shell,
        "tube_fouling_layer_thickness": None
        if tube_fouling_layer_thickness is None
        else input_to_si(spec, "tube_fouling_layer_thickness", tube_fouling_layer_thickness),
        "tube_fouling_layer_conductivity": None
        if tube_fouling_layer_conductivity is None
        else input_to_si(spec, "tube_fouling_layer_conductivity", tube_fouling_layer_conductivity),
        "shell_viscosity_wall": None
        if shell_viscosity_wall is None
        else input_to_si(spec, "shell_viscosity_wall", shell_viscosity_wall),
        "sealing_pairs": sealing_pairs,
        "tube_to_baffle_clearance": None
        if tube_to_baffle_clearance is None
        else input_to_si(spec, "tube_to_baffle_clearance", tube_to_baffle_clearance),
        "shell_to_baffle_clearance": None
        if shell_to_baffle_clearance is None
        else input_to_si(spec, "shell_to_baffle_clearance", shell_to_baffle_clearance),
        "bypass_area": None
        if bypass_area is None
        else input_to_si(spec, "bypass_area", bypass_area),
        **_BASE_DEFAULTS,
    }
    apply_checks(checks.on_input, raw.get, warnings)

    rows: list[dict[str, float | None]] = input_records_to_si(spec, "zones", zones)
    _validate_strict(raw, shell_side_method, has_sealing, triangular_pitch)

    rated: list[ZoneRatingRow] = []
    for row in rows:
        name = row["zone_name"]
        properties = (
            "tube_density",
            "tube_viscosity",
            "tube_cp",
            "tube_conductivity",
            "shell_density",
            "shell_viscosity",
            "shell_cp",
            "shell_conductivity",
        )
        for field in properties:
            value = row[field]
            if value is None or not math.isfinite(value) or value <= 0.0:
                raise InvalidInputError(
                    field,
                    f"the zone {name} states {value}, which is not positive finite",
                )
        duty_fraction = row["duty_fraction"]
        total_duty = row["total_duty"]
        lmtd = row["lmtd"]
        if (
            duty_fraction is None
            or total_duty is None
            or lmtd is None
            or not math.isfinite(duty_fraction)
            or not math.isfinite(total_duty)
            or not math.isfinite(lmtd)
        ):
            raise InvalidInputError(
                "zones", f"the zone {name} states a duty or LMTD that is not finite"
            )

        request = dict(raw)
        request.update({field: row[field] for field in properties})
        values = _evaluate_si(request, shell_side_method, has_sealing, triangular_pitch)
        overall_u = values["overall_u"]
        area = (
            duty_fraction * total_duty / (overall_u * lmtd)
            if overall_u > 0.0 and lmtd > 0.0
            else 0.0
        )
        rated.append(
            ZoneRatingRow(
                zone_name=str(name),
                duty_fraction=duty_fraction,
                tube_side_heat_transfer_coefficient=from_si(values["tube_htc"], "W/(m**2*K)"),
                shell_side_heat_transfer_coefficient=from_si(values["shell_htc"], "W/(m**2*K)"),
                overall_heat_transfer_coefficient=from_si(overall_u, "W/(m**2*K)"),
                lmtd=from_si(lmtd, "K"),
                required_area=from_si(area, "m**2"),
            )
        )

    return ZoneRatingResult(
        zone_count=float(len(rated)),
        zones=tuple(rated),
        warnings=tuple(warnings),
    )


__all__ = ["BELL_DELAWARE", "CALC_ID", "KERN", "zone_rating"]
