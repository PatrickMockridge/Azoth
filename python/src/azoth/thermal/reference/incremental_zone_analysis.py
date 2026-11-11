"""``thermal.incremental_zone_analysis`` - the HTRI-style incremental rating.

Spec: ``specs/calcs/thermal/incremental_zone_analysis.toml``

The Python twin of ``crates/azoth-thermal/src/incremental_zone_analysis.rs``, and both mirror
``IncrementalZoneAnalysis``. The class is a single pass over the caller's zones: no subdivision,
no iteration and no convergence test. **The tube side dispatches on the regime** - Shah for
``condensing``, Gungor-Winterton for ``evaporating``, the single-phase path for everything else,
``two_phase`` included - and the shell side is always Kern and ignores its regime. A zero or
non-positive flow is refused here; the class answers it with a coefficient.
"""

from __future__ import annotations

import math
import sys
from typing import Any

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import IncrementalZoneRow
from azoth.core.result import IncrementalZoneAnalysisResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si
from azoth.core.warnings import Warning
from azoth.thermal.reference.boiling_heat_transfer import average_htc as boiling_average_htc
from azoth.thermal.reference.shah_condensation import (
    average_htc as shah_average_htc,
)
from azoth.thermal.reference.shah_condensation import (
    liquid_only_htc,
)
from azoth.thermal.reference.shell_and_tube_rating import (
    _crossflow_area,
    _kern_shell_side_htc,
    _kern_shell_side_pressure_drop,
    _shell_equivalent_diameter,
    darcy_friction,
)
from azoth.thermal.reference.two_phase_pressure_drop import friedel_pressure_drop

CALC_ID = "thermal.incremental_zone_analysis"

VAPOR = "vapor"
LIQUID = "liquid"
CONDENSING = "condensing"
EVAPORATING = "evaporating"
TWO_PHASE = "two_phase"

_REGIME_SPELLINGS = {
    VAPOR: "VAPOR",
    LIQUID: "LIQUID",
    CONDENSING: "CONDENSING",
    EVAPORATING: "EVAPORATING",
    TWO_PHASE: "TWO_PHASE",
}

LIQUID_DENSITY_DEFAULT = 800.0
VAPOR_DENSITY_DEFAULT = 5.0
LIQUID_VISCOSITY_DEFAULT = 5e-4
VAPOR_VISCOSITY_DEFAULT = 1e-5
LIQUID_CP_DEFAULT = 2500.0
LIQUID_CONDUCTIVITY_DEFAULT = 0.13
QUALITY_IN_DEFAULT = 0.9
QUALITY_OUT_DEFAULT = 0.7
REDUCED_PRESSURE_DEFAULT = 0.1
SURFACE_TENSION_DEFAULT = 0.01
HEAT_OF_VAPORIZATION_DEFAULT = 200000.0


def _positive(field: str, value: float | None) -> float:
    if value is None or not math.isfinite(value) or value <= 0.0:
        raise InvalidInputError(field, f"{value} is not positive finite")
    return value


def _whole(value: float, field: str) -> int:
    if not math.isfinite(value) or value <= 0.0 or value != math.floor(value):
        raise InvalidInputError(field, f"{value} is not a positive whole count")
    return int(value)


def _single_phase_htc(
    mass_flux: float, viscosity: float, cp: float, k: float, tube_id: float
) -> float:
    if viscosity <= 0.0 or k <= 0.0 or tube_id <= 0.0:
        return 0.0
    re = mass_flux * tube_id / viscosity
    pr = cp * viscosity / k
    if re > 3000.0:
        f = darcy_friction(re)
        nu = (
            (f / 8.0)
            * (re - 1000.0)
            * pr
            / (1.0 + 12.7 * math.sqrt(f / 8.0) * (math.pow(pr, 2.0 / 3.0) - 1.0))
        )
        nu = max(nu, 0.0)
    elif re > 2300.0:
        f = darcy_friction(3000.0)
        turbulent = (
            (f / 8.0)
            * 2000.0
            * pr
            / (1.0 + 12.7 * math.sqrt(f / 8.0) * (math.pow(pr, 2.0 / 3.0) - 1.0))
        )
        fraction = (re - 2300.0) / 700.0
        nu = 3.66 + fraction * (turbulent - 3.66)
    else:
        nu = 3.66
    return nu * k / tube_id


def _zone_lmtd(row: dict[str, Any]) -> float:
    dt1: float = float(row["hot_inlet_temp"]) - float(row["cold_outlet_temp"])
    dt2: float = float(row["hot_outlet_temp"]) - float(row["cold_inlet_temp"])
    if dt1 <= 0.0 or dt2 <= 0.0:
        return max(abs(dt1), abs(dt2))
    if abs(dt1 - dt2) < 0.01:
        return (dt1 + dt2) / 2.0
    return (dt1 - dt2) / math.log(dt1 / dt2)


def incremental_zone_analysis(
    tube_inside_diameter: Q,
    tube_outside_diameter: Q,
    tube_length: Q,
    tube_count: float,
    tube_passes: float,
    shell_inside_diameter: Q,
    baffle_spacing: Q,
    tube_pitch: Q,
    triangular_pitch: bool,
    tube_wall_conductivity: Q,
    fouling_tube: float,
    fouling_shell: float,
    zones: list[dict[str, Any]],
) -> IncrementalZoneAnalysisResult:
    """Rates every zone in the caller's order.

    Raises:
        InvalidInputError: for a non-positive geometry, a non-integral count, or a zone whose
            regime-read properties are not positive finite.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    raw = {
        "tube_inside_diameter": input_to_si(spec, "tube_inside_diameter", tube_inside_diameter),
        "tube_outside_diameter": input_to_si(spec, "tube_outside_diameter", tube_outside_diameter),
        "tube_length": input_to_si(spec, "tube_length", tube_length),
        "tube_count": tube_count,
        "tube_passes": tube_passes,
        "shell_inside_diameter": input_to_si(spec, "shell_inside_diameter", shell_inside_diameter),
        "baffle_spacing": input_to_si(spec, "baffle_spacing", baffle_spacing),
        "tube_pitch": input_to_si(spec, "tube_pitch", tube_pitch),
        "tube_wall_conductivity": input_to_si(
            spec, "tube_wall_conductivity", tube_wall_conductivity
        ),
        "fouling_tube": fouling_tube,
        "fouling_shell": fouling_shell,
    }
    apply_checks(checks.on_input, raw.get, warnings)

    positives = [
        raw["tube_inside_diameter"],
        raw["tube_outside_diameter"],
        raw["tube_length"],
        raw["shell_inside_diameter"],
        raw["baffle_spacing"],
        raw["tube_pitch"],
        raw["tube_wall_conductivity"],
    ]
    if any(value is None or not math.isfinite(value) or value <= 0.0 for value in positives):
        raise InvalidInputError(
            "geometry", "the incremental analysis requires positive finite geometry"
        )
    if (
        not math.isfinite(fouling_tube)
        or fouling_tube < 0.0
        or not math.isfinite(fouling_shell)
        or fouling_shell < 0.0
    ):
        raise InvalidInputError(
            "fouling", "the incremental analysis requires finite nonnegative fouling resistances"
        )

    tube_count_i = _whole(tube_count, "tube_count")
    tube_passes_i = _whole(tube_passes, "tube_passes")
    rows: list[dict[str, Any]] = input_records_to_si(spec, "zones", zones)
    n_zones = len(rows)

    area_per_tube = math.pi / 4.0 * raw["tube_inside_diameter"] ** 2
    tubes_per_pass = tube_count_i // tube_passes_i
    flow_area = tubes_per_pass * area_per_tube
    if flow_area <= 0.0:
        raise InvalidInputError("tube_count", "the tube flow area is not positive")

    total_area = 0.0
    total_duty = 0.0
    total_tube_dp = 0.0
    total_shell_dp = 0.0
    total_u_times_area = 0.0
    minimum_approach = sys.float_info.max
    out: list[IncrementalZoneRow] = []

    for index, row in enumerate(rows):
        tube_regime = str(row["tube_phase_regime"])
        shell_regime = str(row["shell_phase_regime"])
        if tube_regime not in _REGIME_SPELLINGS:
            raise InvalidInputError("tube_phase_regime", f"{tube_regime} is not a declared regime")
        if shell_regime not in _REGIME_SPELLINGS:
            raise InvalidInputError(
                "shell_phase_regime", f"{shell_regime} is not a declared regime"
            )
        _positive("tube_mass_flow", row["tube_mass_flow"])
        _positive("shell_mass_flow", row["shell_mass_flow"])
        _positive("shell_density", row["shell_density"])
        _positive("shell_viscosity", row["shell_viscosity"])
        _positive("shell_cp", row["shell_cp"])
        _positive("shell_conductivity", row["shell_conductivity"])
        for field in (
            "duty",
            "hot_inlet_temp",
            "hot_outlet_temp",
            "cold_inlet_temp",
            "cold_outlet_temp",
        ):
            value = row[field]
            if value is None or not math.isfinite(value):
                raise InvalidInputError(
                    "zones", f"the zone {row['zone_name']} states a non-finite {field}"
                )

        liquid_density = row["tube_liquid_density"]
        liquid_density = LIQUID_DENSITY_DEFAULT if liquid_density is None else liquid_density
        liquid_viscosity = row["tube_liquid_viscosity"]
        liquid_viscosity = (
            LIQUID_VISCOSITY_DEFAULT if liquid_viscosity is None else liquid_viscosity
        )
        liquid_cp = row["tube_liquid_cp"]
        liquid_cp = LIQUID_CP_DEFAULT if liquid_cp is None else liquid_cp
        liquid_k = row["tube_liquid_conductivity"]
        liquid_k = LIQUID_CONDUCTIVITY_DEFAULT if liquid_k is None else liquid_k
        vapor_density = row["tube_vapor_density"]
        vapor_density = VAPOR_DENSITY_DEFAULT if vapor_density is None else vapor_density
        vapor_viscosity = row["tube_vapor_viscosity"]
        vapor_viscosity = VAPOR_VISCOSITY_DEFAULT if vapor_viscosity is None else vapor_viscosity
        quality_in = row["tube_quality_in"]
        quality_in = QUALITY_IN_DEFAULT if quality_in is None else quality_in
        quality_out = row["tube_quality_out"]
        quality_out = QUALITY_OUT_DEFAULT if quality_out is None else quality_out
        reduced_pressure = row["tube_reduced_pressure"]
        reduced_pressure = (
            REDUCED_PRESSURE_DEFAULT if reduced_pressure is None else reduced_pressure
        )
        surface_tension = row["tube_surface_tension"]
        surface_tension = SURFACE_TENSION_DEFAULT if surface_tension is None else surface_tension
        latent = row["tube_heat_of_vaporization"]
        latent = HEAT_OF_VAPORIZATION_DEFAULT if latent is None else latent

        mass_flux = row["tube_mass_flow"] / flow_area
        if tube_regime in (CONDENSING, EVAPORATING):
            _positive("tube_liquid_density", liquid_density)
            _positive("tube_liquid_viscosity", liquid_viscosity)
            _positive("tube_liquid_cp", liquid_cp)
            _positive("tube_liquid_conductivity", liquid_k)
            if tube_regime == CONDENSING:
                _positive("tube_reduced_pressure", reduced_pressure)
            else:
                _positive("tube_vapor_density", vapor_density)
                _positive("tube_vapor_viscosity", vapor_viscosity)
                _positive("tube_heat_of_vaporization", latent)
        else:
            _positive("tube_density", row["tube_density"])
            _positive("tube_viscosity", row["tube_viscosity"])
            _positive("tube_cp", row["tube_cp"])
            _positive("tube_conductivity", row["tube_conductivity"])

        if tube_regime == CONDENSING:
            h_lo = liquid_only_htc(
                mass_flux,
                raw["tube_inside_diameter"],
                liquid_density,
                liquid_viscosity,
                liquid_cp,
                liquid_k,
            )
            h_tube = shah_average_htc(h_lo, reduced_pressure, quality_in, quality_out, 10)
        elif tube_regime == EVAPORATING:
            heat_flux = (
                row["duty"]
                / (
                    tube_count_i
                    * math.pi
                    * raw["tube_inside_diameter"]
                    * raw["tube_length"]
                    / n_zones
                )
                if tube_count_i > 0 and raw["tube_inside_diameter"] > 0 and n_zones > 0
                else 0.0
            )
            h_tube = boiling_average_htc(
                mass_flux,
                raw["tube_inside_diameter"],
                liquid_density,
                vapor_density,
                liquid_viscosity,
                liquid_cp,
                liquid_k,
                heat_flux,
                latent,
                quality_in,
                quality_out,
                10,
            )
        else:
            h_tube = _single_phase_htc(
                mass_flux,
                row["tube_viscosity"],
                row["tube_cp"],
                row["tube_conductivity"],
                raw["tube_inside_diameter"],
            )

        equivalent = _shell_equivalent_diameter(
            raw["tube_outside_diameter"], raw["tube_pitch"], triangular_pitch
        )
        crossflow = _crossflow_area(
            raw["shell_inside_diameter"],
            raw["baffle_spacing"],
            raw["tube_outside_diameter"],
            raw["tube_pitch"],
        )
        shell_flux = row["shell_mass_flow"] / crossflow if crossflow > 0.0 else 0.0
        h_shell = (
            _kern_shell_side_htc(
                shell_flux,
                equivalent,
                row["shell_viscosity"],
                row["shell_cp"],
                row["shell_conductivity"],
                row["shell_viscosity"],
            )
            if equivalent > 0.0
            else 0.0
        )

        if h_tube <= 0.0 or h_shell <= 0.0:
            overall_u = 0.0
        else:
            ratio = raw["tube_outside_diameter"] / raw["tube_inside_diameter"]
            total_resistance = (
                1.0 / h_shell
                + fouling_shell
                + raw["tube_outside_diameter"]
                * math.log(ratio)
                / (2.0 * raw["tube_wall_conductivity"])
                + ratio * fouling_tube
                + ratio / h_tube
            )
            overall_u = 1.0 / total_resistance if total_resistance > 0.0 else 0.0

        lmtd = _zone_lmtd(row)
        area = row["duty"] / (overall_u * lmtd) if overall_u > 0.0 and lmtd > 0.0 else 0.0

        zone_length = raw["tube_length"] / max(n_zones, 1)
        if tube_regime in (CONDENSING, EVAPORATING):
            tube_dp = friedel_pressure_drop(
                mass_flux,
                (quality_in + quality_out) / 2.0,
                raw["tube_inside_diameter"],
                zone_length,
                liquid_density,
                vapor_density,
                liquid_viscosity,
                vapor_viscosity,
                surface_tension,
            )
        elif row["tube_viscosity"] <= 0.0 or row["tube_density"] <= 0.0 or mass_flux <= 0.0:
            tube_dp = 0.0
        else:
            velocity = mass_flux / row["tube_density"]
            re = mass_flux * raw["tube_inside_diameter"] / row["tube_viscosity"]
            friction = darcy_friction(re)
            tube_dp = (
                friction
                * zone_length
                / raw["tube_inside_diameter"]
                * row["tube_density"]
                * velocity
                * velocity
                / 2.0
            )

        zone_baffles = max(1, round(raw["tube_length"] / raw["baffle_spacing"] / max(n_zones, 1)))
        shell_dp = (
            _kern_shell_side_pressure_drop(
                shell_flux,
                equivalent,
                raw["shell_inside_diameter"],
                zone_baffles,
                row["shell_density"],
                row["shell_viscosity"],
                row["shell_viscosity"],
            )
            if equivalent > 0.0
            else 0.0
        )

        total_area += area
        total_duty += row["duty"]
        total_tube_dp += tube_dp
        total_shell_dp += shell_dp
        total_u_times_area += overall_u * area
        minimum_approach = min(
            minimum_approach,
            min(
                abs(row["hot_inlet_temp"] - row["cold_outlet_temp"]),
                abs(row["hot_outlet_temp"] - row["cold_inlet_temp"]),
            ),
        )

        out.append(
            IncrementalZoneRow(
                zone_index=float(index + 1),
                zone_name=str(row["zone_name"]),
                tube_phase_regime=_REGIME_SPELLINGS[tube_regime],
                shell_phase_regime=_REGIME_SPELLINGS[shell_regime],
                tube_side_heat_transfer_coefficient=from_si(h_tube, "W/(m**2*K)"),
                shell_side_heat_transfer_coefficient=from_si(h_shell, "W/(m**2*K)"),
                overall_heat_transfer_coefficient=from_si(overall_u, "W/(m**2*K)"),
                lmtd=from_si(lmtd, "K"),
                duty=from_si(row["duty"], "W"),
                required_area=from_si(area, "m**2"),
                tube_side_pressure_drop=from_si(tube_dp, "Pa"),
                shell_side_pressure_drop=from_si(shell_dp, "Pa"),
                hot_inlet_temp=from_si(row["hot_inlet_temp"], "K"),
                hot_outlet_temp=from_si(row["hot_outlet_temp"], "K"),
                cold_inlet_temp=from_si(row["cold_inlet_temp"], "K"),
                cold_outlet_temp=from_si(row["cold_outlet_temp"], "K"),
            )
        )

    weighted = total_u_times_area / total_area if total_area > 0.0 else 0.0
    return IncrementalZoneAnalysisResult(
        number_of_zones=float(n_zones),
        total_required_area=from_si(total_area, "m**2"),
        total_duty=from_si(total_duty, "W"),
        weighted_average_u=from_si(weighted, "W/(m**2*K)"),
        total_tube_side_pressure_drop=from_si(total_tube_dp, "Pa"),
        total_shell_side_pressure_drop=from_si(total_shell_dp, "Pa"),
        minimum_approach_temperature=from_si(minimum_approach, "K"),
        zones=tuple(out),
        warnings=tuple(warnings),
    )


__all__ = ["CALC_ID", "incremental_zone_analysis"]
