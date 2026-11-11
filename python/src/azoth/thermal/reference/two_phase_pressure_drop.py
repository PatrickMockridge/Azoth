"""``thermal.two_phase_pressure_drop`` - Friedel and Muller-Steinhagen-Heck.

Spec: ``specs/calcs/thermal/two_phase_pressure_drop.toml``

The Python twin of ``crates/azoth-thermal/src/two_phase_pressure_drop.rs``, and both mirror
``TwoPhasePressureDrop``: the Friedel frictional gradient and its total and quality-averaged
drops, the Muller-Steinhagen-Heck gradient, and the gravitational and acceleration terms of a
homogeneous two-phase mixture.

**At either quality endpoint the two frictional gradients agree**: at zero both reduce to the
all-liquid Fanning term, at one to the all-vapour term. **A zero surface tension is not a
refusal**: the class answers its ``We = 1e6`` sentinel. And a falling quality answers a
negative acceleration term, a pressure recovery rather than a drop.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TwoPhasePressureDropResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.two_phase_pressure_drop"
GRAVITY = 9.81


def fanning_friction(reynolds: float) -> float:
    """The smooth-tube Fanning factor: 16/Re below 2300, Blasius above."""
    if reynolds <= 0.0:
        return 0.0
    if reynolds < 2300.0:
        return 16.0 / reynolds
    return 0.079 / math.pow(reynolds, 0.25)


def friedel_gradient(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
    surface_tension: float,
) -> float:
    """The Friedel two-phase frictional gradient."""
    if (
        mass_flux <= 0.0
        or tube_id <= 0.0
        or liquid_density <= 0.0
        or vapor_density <= 0.0
        or liquid_viscosity <= 0.0
        or vapor_viscosity <= 0.0
    ):
        return 0.0
    x = min(max(vapor_quality, 0.0), 1.0)
    re_lo = mass_flux * tube_id / liquid_viscosity
    f_lo = fanning_friction(re_lo)
    dpdz_lo = 4.0 * f_lo * mass_flux * mass_flux / (2.0 * liquid_density * tube_id)
    if x < 1e-10:
        return dpdz_lo
    re_vo = mass_flux * tube_id / vapor_viscosity
    f_vo = fanning_friction(re_vo)
    e = math.pow(1.0 - x, 2) + x * x * (liquid_density * f_vo) / (vapor_density * f_lo)
    f = math.pow(x, 0.78) * math.pow(1.0 - x, 0.224)
    h = (
        math.pow(liquid_density / vapor_density, 0.91)
        * math.pow(vapor_viscosity / liquid_viscosity, 0.19)
        * math.pow(1.0 - vapor_viscosity / liquid_viscosity, 0.7)
    )
    rho_h = 1.0 / (x / vapor_density + (1.0 - x) / liquid_density)
    fr = mass_flux * mass_flux / (rho_h * rho_h * GRAVITY * tube_id)
    we = (
        mass_flux * mass_flux * tube_id / (rho_h * surface_tension)
        if surface_tension > 0.0
        else 1e6
    )
    fr_we = math.pow(max(fr, 1e-10), 0.045) * math.pow(max(we, 1e-10), 0.035)
    phi_lo_sq = e + 3.24 * f * h / fr_we
    return phi_lo_sq * dpdz_lo


def friedel_pressure_drop(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    tube_length: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
    surface_tension: float,
) -> float:
    """The Friedel drop over a tube length at a constant quality."""
    return (
        friedel_gradient(
            mass_flux,
            vapor_quality,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        )
        * tube_length
    )


def friedel_average_pressure_drop(
    mass_flux: float,
    quality_in: float,
    quality_out: float,
    tube_id: float,
    tube_length: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
    surface_tension: float,
    intervals: int,
) -> float:
    """Simpson's-rule average of the Friedel gradient over a range, times the length."""
    if tube_length <= 0.0:
        return 0.0
    x_in = min(max(quality_in, 0.0), 1.0)
    x_out = min(max(quality_out, 0.0), 1.0)
    if abs(x_in - x_out) < 1e-10:
        return friedel_pressure_drop(
            mass_flux,
            x_in,
            tube_id,
            tube_length,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        )
    n = max(4, intervals)
    if n % 2 != 0:
        n += 1
    dx = (x_out - x_in) / n
    total = friedel_gradient(
        mass_flux,
        x_in,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        vapor_viscosity,
        surface_tension,
    ) + friedel_gradient(
        mass_flux,
        x_out,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        vapor_viscosity,
        surface_tension,
    )
    for i in range(1, n):
        x = x_in + i * dx
        gradient = friedel_gradient(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        )
        total += (2.0 if i % 2 == 0 else 4.0) * gradient
    return abs(total * dx / 3.0 / (x_out - x_in)) * tube_length


def muller_steinhagen_heck_gradient(
    mass_flux: float,
    vapor_quality: float,
    tube_id: float,
    liquid_density: float,
    vapor_density: float,
    liquid_viscosity: float,
    vapor_viscosity: float,
) -> float:
    """The Muller-Steinhagen-Heck two-phase frictional gradient."""
    if (
        mass_flux <= 0.0
        or tube_id <= 0.0
        or liquid_density <= 0.0
        or vapor_density <= 0.0
        or liquid_viscosity <= 0.0
        or vapor_viscosity <= 0.0
    ):
        return 0.0
    x = min(max(vapor_quality, 0.0), 1.0)
    re_lo = mass_flux * tube_id / liquid_viscosity
    f_lo = fanning_friction(re_lo)
    a = 4.0 * f_lo * mass_flux * mass_flux / (2.0 * liquid_density * tube_id)
    re_vo = mass_flux * tube_id / vapor_viscosity
    f_vo = fanning_friction(re_vo)
    b = 4.0 * f_vo * mass_flux * mass_flux / (2.0 * vapor_density * tube_id)
    g_factor = a + 2.0 * (b - a) * x
    return g_factor * math.pow(1.0 - x, 1.0 / 3.0) + b * x * x * x


def homogeneous_void_fraction(x: float, liquid_density: float, vapor_density: float) -> float:
    """The homogeneous void fraction, with the class's endpoint pins."""
    if x < 1e-10:
        return 0.0
    if x > 1.0 - 1e-10:
        return 1.0
    return 1.0 / (1.0 + (1.0 - x) / x * vapor_density / liquid_density)


def gravitational_gradient(
    vapor_quality: float, liquid_density: float, vapor_density: float
) -> float:
    """The gravitational gradient of the homogeneous mixture."""
    if liquid_density <= 0.0 or vapor_density <= 0.0:
        return 0.0
    x = min(max(vapor_quality, 0.0), 1.0)
    alpha = homogeneous_void_fraction(x, liquid_density, vapor_density)
    rho_tp = liquid_density * (1.0 - alpha) + vapor_density * alpha
    return rho_tp * GRAVITY


def momentum_flux(quality: float, liquid_density: float, vapor_density: float) -> float:
    """The momentum-flux term at a quality, with the class's endpoint pins."""
    x = min(max(quality, 0.0), 1.0)
    if x < 1e-10:
        return 1.0 / liquid_density
    if x > 1.0 - 1e-10:
        return 1.0 / vapor_density
    alpha = homogeneous_void_fraction(x, liquid_density, vapor_density)
    liquid_term = math.pow(1.0 - x, 2) / ((1.0 - alpha) * liquid_density)
    vapor_term = x * x / (alpha * vapor_density)
    return liquid_term + vapor_term


def acceleration_pressure_drop(
    mass_flux: float,
    quality_in: float,
    quality_out: float,
    liquid_density: float,
    vapor_density: float,
) -> float:
    """The acceleration pressure drop over a quality change."""
    if mass_flux <= 0.0 or liquid_density <= 0.0 or vapor_density <= 0.0:
        return 0.0
    flux_out = momentum_flux(quality_out, liquid_density, vapor_density)
    flux_in = momentum_flux(quality_in, liquid_density, vapor_density)
    return mass_flux * mass_flux * (flux_out - flux_in)


def two_phase_pressure_drop(
    mass_flux: Q,
    vapor_quality: float,
    tube_inner_diameter: Q,
    tube_length: Q,
    liquid_density: Q,
    vapor_density: Q,
    liquid_viscosity: Q,
    vapor_viscosity: Q,
    surface_tension: Q,
    quality_in: float,
    quality_out: float,
    intervals: float | None = None,
) -> TwoPhasePressureDropResult:
    """The composed two-phase pressure-drop calculation.

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
        "tube_length": input_to_si(spec, "tube_length", tube_length),
        "liquid_density": input_to_si(spec, "liquid_density", liquid_density),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
        "vapor_viscosity": input_to_si(spec, "vapor_viscosity", vapor_viscosity),
        "surface_tension": input_to_si(spec, "surface_tension", surface_tension),
    }
    apply_checks(checks.on_input, values.get, warnings)

    if intervals is None:
        n = 20
    else:
        if not math.isfinite(intervals) or intervals != math.floor(intervals):
            raise InvalidInputError("intervals", f"{intervals} is not a whole interval count")
        n = int(intervals)

    return TwoPhasePressureDropResult(
        friedel_gradient=from_si(
            friedel_gradient(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["vapor_viscosity"],
                values["surface_tension"],
            ),
            "Pa/m",
        ),
        friedel_pressure_drop=from_si(
            friedel_pressure_drop(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["tube_length"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["vapor_viscosity"],
                values["surface_tension"],
            ),
            "Pa",
        ),
        friedel_average_pressure_drop=from_si(
            friedel_average_pressure_drop(
                values["mass_flux"],
                quality_in,
                quality_out,
                values["tube_inner_diameter"],
                values["tube_length"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["vapor_viscosity"],
                values["surface_tension"],
                n,
            ),
            "Pa",
        ),
        muller_steinhagen_heck_gradient=from_si(
            muller_steinhagen_heck_gradient(
                values["mass_flux"],
                vapor_quality,
                values["tube_inner_diameter"],
                values["liquid_density"],
                values["vapor_density"],
                values["liquid_viscosity"],
                values["vapor_viscosity"],
            ),
            "Pa/m",
        ),
        gravitational_gradient=from_si(
            gravitational_gradient(
                vapor_quality, values["liquid_density"], values["vapor_density"]
            ),
            "Pa/m",
        ),
        acceleration_pressure_drop=from_si(
            acceleration_pressure_drop(
                values["mass_flux"],
                quality_in,
                quality_out,
                values["liquid_density"],
                values["vapor_density"],
            ),
            "Pa",
        ),
        warnings=tuple(warnings),
    )
