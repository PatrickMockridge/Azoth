"""``thermal.boiling_heat_transfer`` - Chen and Gungor-Winterton flow boiling.

Spec: ``specs/calcs/thermal/boiling_heat_transfer.toml``

The Python twin of ``crates/azoth-thermal/src/boiling_heat_transfer.rs``, and both mirror
``BoilingHeatTransfer`` composed as ``IncrementalZoneAnalysis`` uses it: the average
Gungor-Winterton coefficient over a quality range, with the Chen local, the two factors and
the horizontal correction beside it.

**The boiling quality clamp is ``[0.001, 0.999]``**, unlike the condensation correlation's
``[0, 1]``. The horizontal correction applies below ``Fr_l = 0.05`` and the average uses
the uncorrected coefficient, as the class does.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BoilingHeatTransferResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.boiling_heat_transfer"
MIN_QUALITY = 0.001
MAX_QUALITY = 0.999
GRAVITY = 9.81


def martinelli_parameter(
    x: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
) -> float:
    """The Lockhart-Martinelli parameter, answering 1 at an endpoint."""
    if (
        x <= 0.0
        or x >= 1.0
        or vapor_density <= 0.0
        or liquid_density <= 0.0
        or vapor_viscosity <= 0.0
    ):
        return 1.0
    return (
        math.pow((1.0 - x) / x, 0.9)
        * math.pow(vapor_density / liquid_density, 0.5)
        * math.pow(liquid_viscosity / vapor_viscosity, 0.1)
    )


def chen_enhancement_factor(x_tt: float) -> float:
    """The Chen enhancement factor."""
    if x_tt <= 0.0:
        return 1.0
    inverse = 1.0 / x_tt
    if inverse <= 0.1:
        return 1.0
    return 2.35 * math.pow(0.213 + inverse, 0.736)


def chen_suppression_factor(re_two_phase: float) -> float:
    """The Chen suppression factor."""
    if re_two_phase <= 0.0:
        return 1.0
    return 1.0 / (1.0 + 2.53e-6 * math.pow(re_two_phase, 1.17))


def chen_htc(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
    liquid_cp: float,
    liquid_conductivity: float,
    surface_tension: float,
    heat_of_vaporization: float,
    wall_superheat: float,
    saturation_pressure_difference: float,
) -> float:
    """The Chen local coefficient."""
    if liquid_viscosity <= 0.0 or liquid_conductivity <= 0.0 or tube_id <= 0.0 or mass_flux <= 0.0:
        return 0.0
    x = min(max(vapor_quality, MIN_QUALITY), MAX_QUALITY)
    re_liquid = mass_flux * (1.0 - x) * tube_id / liquid_viscosity
    pr_liquid = liquid_cp * liquid_viscosity / liquid_conductivity
    h_liquid = (
        0.023
        * math.pow(max(re_liquid, 1.0), 0.8)
        * math.pow(pr_liquid, 0.4)
        * liquid_conductivity
        / tube_id
    )
    x_tt = martinelli_parameter(x, liquid_density, vapor_density, liquid_viscosity, vapor_viscosity)
    enhancement = chen_enhancement_factor(x_tt)

    nucleate = 0.0
    if (
        wall_superheat > 0.0
        and saturation_pressure_difference > 0.0
        and surface_tension > 0.0
        and heat_of_vaporization > 0.0
    ):
        numerator = (
            math.pow(liquid_conductivity, 0.79)
            * math.pow(liquid_cp, 0.45)
            * math.pow(liquid_density, 0.49)
        )
        denominator = (
            math.pow(surface_tension, 0.5)
            * math.pow(liquid_viscosity, 0.29)
            * math.pow(heat_of_vaporization, 0.24)
            * math.pow(vapor_density, 0.24)
        )
        if denominator > 0.0:
            nucleate = (
                0.00122
                * (numerator / denominator)
                * math.pow(wall_superheat, 0.24)
                * math.pow(saturation_pressure_difference, 0.75)
            )
    re_two_phase = re_liquid * math.pow(enhancement, 1.25)
    suppression = chen_suppression_factor(re_two_phase)
    return h_liquid * enhancement + nucleate * suppression


def gungor_winterton_htc(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    liquid_cp: float,
    liquid_conductivity: float,
    heat_flux: float,
    heat_of_vaporization: float,
) -> float:
    """The Gungor-Winterton local coefficient."""
    if liquid_viscosity <= 0.0 or liquid_conductivity <= 0.0 or tube_id <= 0.0 or mass_flux <= 0.0:
        return 0.0
    x = min(max(vapor_quality, MIN_QUALITY), MAX_QUALITY)
    re_liquid = mass_flux * (1.0 - x) * tube_id / liquid_viscosity
    pr_liquid = liquid_cp * liquid_viscosity / liquid_conductivity
    h_liquid = (
        0.023
        * math.pow(max(re_liquid, 1.0), 0.8)
        * math.pow(pr_liquid, 0.4)
        * liquid_conductivity
        / tube_id
    )
    boiling_number = (
        heat_flux / (mass_flux * heat_of_vaporization)
        if heat_of_vaporization > 0.0 and mass_flux > 0.0
        else 0.0
    )
    convective = 0.0
    if x < MAX_QUALITY and vapor_density > 0.0:
        convective = (
            1.12 * math.pow(x / (1.0 - x), 0.75) * math.pow(liquid_density / vapor_density, 0.41)
        )
    nucleate = 3000.0 * math.pow(max(boiling_number, 0.0), 0.86)
    return h_liquid * (1.0 + nucleate + convective)


def gungor_winterton_corrected_htc(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    liquid_cp: float,
    liquid_conductivity: float,
    heat_flux: float,
    heat_of_vaporization: float,
    horizontal: bool,
) -> float:
    """The Gungor-Winterton coefficient with the horizontal stratified correction."""
    h = gungor_winterton_htc(
        mass_flux,
        vapor_quality,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        liquid_cp,
        liquid_conductivity,
        heat_flux,
        heat_of_vaporization,
    )
    if horizontal and liquid_density > 0.0:
        x = min(max(vapor_quality, MIN_QUALITY), MAX_QUALITY)
        liquid_flux = mass_flux * (1.0 - x)
        liquid_velocity = liquid_flux / liquid_density
        froude = liquid_velocity * liquid_velocity / (GRAVITY * tube_id)
        if froude < 0.05:
            h *= math.pow(froude, 0.1 - 2.0 * froude)
    return h


def average_htc(
    mass_flux: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    liquid_cp: float,
    liquid_conductivity: float,
    heat_flux: float,
    heat_of_vaporization: float,
    quality_in: float,
    quality_out: float,
    intervals: int,
) -> float:
    """Simpson's-rule average of the uncorrected Gungor-Winterton coefficient."""
    if mass_flux <= 0.0 or liquid_viscosity <= 0.0 or tube_id <= 0.0:
        return 0.0
    x_in = min(max(quality_in, MIN_QUALITY), MAX_QUALITY)
    x_out = min(max(quality_out, MIN_QUALITY), MAX_QUALITY)
    if abs(x_in - x_out) < 1e-10:
        return gungor_winterton_htc(
            mass_flux,
            x_in,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            liquid_cp,
            liquid_conductivity,
            heat_flux,
            heat_of_vaporization,
        )
    n = max(4, intervals)
    if n % 2 != 0:
        n += 1
    dx = (x_out - x_in) / n
    total = gungor_winterton_htc(
        mass_flux,
        x_in,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        liquid_cp,
        liquid_conductivity,
        heat_flux,
        heat_of_vaporization,
    ) + gungor_winterton_htc(
        mass_flux,
        x_out,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        liquid_cp,
        liquid_conductivity,
        heat_flux,
        heat_of_vaporization,
    )
    for i in range(1, n):
        x = x_in + i * dx
        h = gungor_winterton_htc(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            liquid_cp,
            liquid_conductivity,
            heat_flux,
            heat_of_vaporization,
        )
        total += (2.0 if i % 2 == 0 else 4.0) * h
    return abs(total * dx / 3.0 / (x_out - x_in))


def boiling_heat_transfer(
    mass_flux: Q,
    tube_inner_diameter: Q,
    liquid_density: Q,
    vapor_density: Q,
    liquid_viscosity: Q,
    vapor_viscosity: Q,
    liquid_cp: Q,
    liquid_conductivity: Q,
    surface_tension: Q,
    heat_of_vaporization: float,
    wall_superheat: Q,
    saturation_pressure_difference: Q,
    heat_flux: Q,
    vapor_quality: float,
    quality_in: float,
    quality_out: float,
    intervals: float | None = None,
    horizontal: bool | None = None,
) -> BoilingHeatTransferResult:
    """The composed boiling calculation.

    Raises:
        InvalidInputError: for a non-integral interval count.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    values = {
        "mass_flux": input_to_si(spec, "mass_flux", mass_flux),
        "tube_inner_diameter": input_to_si(spec, "tube_inner_diameter", tube_inner_diameter),
        "liquid_density": input_to_si(spec, "liquid_density", liquid_density),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
        "vapor_viscosity": input_to_si(spec, "vapor_viscosity", vapor_viscosity),
        "liquid_cp": input_to_si(spec, "liquid_cp", liquid_cp),
        "liquid_conductivity": input_to_si(spec, "liquid_conductivity", liquid_conductivity),
        "surface_tension": input_to_si(spec, "surface_tension", surface_tension),
        "heat_of_vaporization": heat_of_vaporization,
        "wall_superheat": input_to_si(spec, "wall_superheat", wall_superheat),
        "saturation_pressure_difference": input_to_si(
            spec, "saturation_pressure_difference", saturation_pressure_difference
        ),
        "heat_flux": input_to_si(spec, "heat_flux", heat_flux),
    }
    apply_checks(checks.on_input, values.get, warnings)

    if intervals is None:
        n = 20
    else:
        if not math.isfinite(intervals) or intervals != math.floor(intervals):
            raise InvalidInputError("intervals", f"{intervals} is not a whole interval count")
        n = int(intervals)

    horizontal_value = bool(horizontal)
    # The published factors are the class's public helpers at the stated quality, while the
    # Chen local clamps the quality first; at a clamp the two differ.
    x_tt = martinelli_parameter(
        vapor_quality,
        values["liquid_density"],
        values["vapor_density"],
        values["liquid_viscosity"],
        values["vapor_viscosity"],
    )
    enhancement = chen_enhancement_factor(x_tt)
    x = min(max(vapor_quality, MIN_QUALITY), MAX_QUALITY)
    re_liquid = (
        values["mass_flux"] * (1.0 - x) * values["tube_inner_diameter"] / values["liquid_viscosity"]
    )
    return BoilingHeatTransferResult(
        chen_htc=from_si(
            chen_htc(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["vapor_viscosity"],
                values["liquid_cp"],
                values["liquid_conductivity"],
                values["surface_tension"],
                values["heat_of_vaporization"],
                values["wall_superheat"],
                values["saturation_pressure_difference"],
            ),
            "W/(m**2*K)",
        ),
        gungor_winterton_htc=from_si(
            gungor_winterton_htc(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["liquid_cp"],
                values["liquid_conductivity"],
                values["heat_flux"],
                values["heat_of_vaporization"],
            ),
            "W/(m**2*K)",
        ),
        gungor_winterton_corrected_htc=from_si(
            gungor_winterton_corrected_htc(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["liquid_cp"],
                values["liquid_conductivity"],
                values["heat_flux"],
                values["heat_of_vaporization"],
                horizontal_value,
            ),
            "W/(m**2*K)",
        ),
        average_htc=from_si(
            average_htc(
                values["mass_flux"],
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["liquid_cp"],
                values["liquid_conductivity"],
                values["heat_flux"],
                values["heat_of_vaporization"],
                quality_in,
                quality_out,
                n,
            ),
            "W/(m**2*K)",
        ),
        martinelli_parameter=x_tt,
        chen_enhancement_factor=enhancement,
        chen_suppression_factor=chen_suppression_factor(re_liquid * math.pow(enhancement, 1.25)),
        warnings=tuple(warnings),
    )
