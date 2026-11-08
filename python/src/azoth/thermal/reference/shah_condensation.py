"""``thermal.shah_condensation`` - Shah's in-tube condensation correlation.

Spec: ``specs/calcs/thermal/shah_condensation.toml``

The Python twin of ``crates/azoth-thermal/src/shah_condensation.rs``, and both mirror
``ShahCondensation`` composed as ``IncrementalZoneAnalysis`` uses it. The liquid-only
Dittus-Boelter coefficient feeds the local, average and vertical entries.

**Two behaviours not to tidy.** At a quality of one both Shah terms vanish and the local
coefficient is zero; a quality at or below zero clamps to zero and answers the liquid-only
coefficient. And a liquid density not above the vapour makes the vertical ``Jg``
denominator non-positive and falls back to the shear-dominated local coefficient.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ShahCondensationResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.shah_condensation"
MIN_REDUCED_PRESSURE = 0.001
MAX_REDUCED_PRESSURE = 0.44
GRAVITY = 9.81


def liquid_only_htc(
    mass_flux: float,
    tube_id: float,
    _liquid_density: float,
    liquid_viscosity: float,
    liquid_cp: float,
    liquid_conductivity: float,
) -> float:
    """The all-liquid Dittus-Boelter coefficient."""
    if liquid_viscosity <= 0.0 or liquid_conductivity <= 0.0 or tube_id <= 0.0 or mass_flux <= 0.0:
        return 0.0
    re_lo = mass_flux * tube_id / liquid_viscosity
    pr_l = liquid_cp * liquid_viscosity / liquid_conductivity
    nu_lo = 0.023 * math.pow(re_lo, 0.8) * math.pow(pr_l, 0.4)
    return nu_lo * liquid_conductivity / tube_id


def local_htc(h_lo: float, vapor_quality: float, reduced_pressure: float) -> float:
    """The local Shah coefficient at a quality and reduced pressure."""
    if h_lo <= 0.0 or reduced_pressure <= 0.0:
        return 0.0
    x = min(max(vapor_quality, 0.0), 1.0)
    pr = min(max(reduced_pressure, MIN_REDUCED_PRESSURE), MAX_REDUCED_PRESSURE)
    term1 = math.pow(1.0 - x, 0.8)
    term2 = 3.8 * math.pow(x, 0.76) * math.pow(1.0 - x, 0.04) / math.pow(pr, 0.38)
    return h_lo * (term1 + term2)


def average_htc(
    h_lo: float,
    reduced_pressure: float,
    quality_in: float,
    quality_out: float,
    intervals: int,
) -> float:
    """Simpson's-rule average over a quality range."""
    if h_lo <= 0.0 or reduced_pressure <= 0.0:
        return 0.0
    x_in = min(max(quality_in, 0.0), 1.0)
    x_out = min(max(quality_out, 0.0), 1.0)
    if abs(x_in - x_out) < 1e-10:
        return local_htc(h_lo, x_in, reduced_pressure)
    n = max(4, intervals)
    if n % 2 != 0:
        n += 1
    dx = (x_out - x_in) / n
    total = local_htc(h_lo, x_in, reduced_pressure) + local_htc(h_lo, x_out, reduced_pressure)
    for i in range(1, n):
        x = x_in + i * dx
        h = local_htc(h_lo, x, reduced_pressure)
        total += (2.0 if i % 2 == 0 else 4.0) * h
    # The class's expression is `sum*|dx|/3/|x_out-x_in|*|x_out-x_in|`; the two range
    # factors cancel, and the literal form is carried so the rounding is the same.
    return total * abs(dx) / 3.0 / abs(x_out - x_in) * abs(x_out - x_in)


def vertical_htc(
    h_lo: float,
    vapor_quality: float,
    reduced_pressure: float,
    mass_flux: float,
    tube_id: float,
    vapor_density: float,
    liquid_density: float,
) -> float:
    """The vertical-tube coefficient, with the class's gravity branch."""
    if h_lo <= 0.0 or reduced_pressure <= 0.0 or tube_id <= 0.0:
        return 0.0
    x = min(max(vapor_quality, 0.0), 1.0)
    argument = GRAVITY * tube_id * vapor_density * (liquid_density - vapor_density)
    h_shah = local_htc(h_lo, x, reduced_pressure)
    # A non-positive argument is the class's sentinel for `Jg` (its `sqrt` answers NaN and
    # the guard is `> 0`), which then takes the shear-dominated branch.
    if argument <= 0.0:
        return h_shah
    jg = x * mass_flux / math.sqrt(argument)
    if jg >= 0.98:
        return h_shah
    pr = max(reduced_pressure, MIN_REDUCED_PRESSURE)
    h_nu = 1.32 * math.pow(pr, -1.0 / 3.0) * h_lo
    h_grav = h_nu * math.pow(1.0 - x, 0.8)
    return max(h_shah, h_grav)


def is_in_valid_range(reduced_pressure: float, re_liquid_only: float, mass_flux: float) -> bool:
    """The class's validated-range predicate, over the stated values."""
    if reduced_pressure < MIN_REDUCED_PRESSURE or reduced_pressure > MAX_REDUCED_PRESSURE:
        return False
    if re_liquid_only < 350.0:
        return False
    return 4.0 <= mass_flux <= 820.0


def shah_condensation(
    mass_flux: Q,
    tube_inner_diameter: Q,
    liquid_density: Q,
    liquid_viscosity: Q,
    liquid_cp: Q,
    liquid_conductivity: Q,
    vapor_density: Q,
    reduced_pressure: float,
    vapor_quality: float,
    quality_in: float,
    quality_out: float,
    intervals: float | None = None,
) -> ShahCondensationResult:
    """The composed condensation calculation.

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
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
        "liquid_cp": input_to_si(spec, "liquid_cp", liquid_cp),
        "liquid_conductivity": input_to_si(spec, "liquid_conductivity", liquid_conductivity),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "reduced_pressure": reduced_pressure,
    }
    apply_checks(checks.on_input, values.get, warnings)

    if intervals is None:
        n = 20
    else:
        if not math.isfinite(intervals) or intervals != math.floor(intervals):
            raise InvalidInputError("intervals", f"{intervals} is not a whole interval count")
        n = int(intervals)

    h_lo = liquid_only_htc(
        values["mass_flux"],
        values["tube_inner_diameter"],
        values["liquid_density"],
        values["liquid_viscosity"],
        values["liquid_cp"],
        values["liquid_conductivity"],
    )
    re_lo = values["mass_flux"] * values["tube_inner_diameter"] / values["liquid_viscosity"]
    return ShahCondensationResult(
        liquid_only_htc=from_si(h_lo, "W/(m**2*K)"),
        liquid_only_reynolds=re_lo,
        local_htc=from_si(local_htc(h_lo, vapor_quality, reduced_pressure), "W/(m**2*K)"),
        average_htc=from_si(
            average_htc(h_lo, reduced_pressure, quality_in, quality_out, n),
            "W/(m**2*K)",
        ),
        vertical_htc=from_si(
            vertical_htc(
                h_lo,
                vapor_quality,
                reduced_pressure,
                values["mass_flux"],
                values["tube_inner_diameter"],
                values["vapor_density"],
                values["liquid_density"],
            ),
            "W/(m**2*K)",
        ),
        in_valid_range=is_in_valid_range(reduced_pressure, re_lo, values["mass_flux"]),
        warnings=tuple(warnings),
    )
