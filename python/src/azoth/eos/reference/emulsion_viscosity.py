"""``eos.emulsion_viscosity`` - NeqSim's ``EmulsionViscosityCalculator``.

The Arirachakaran inversion point and one of six relative-viscosity models, with the
class's tightness and demulsifier corrections.

Spec: ``specs/calcs/eos/emulsion_viscosity.toml``.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EmulsionType, EmulsionViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.emulsion_viscosity"

MODELS = {"einstein", "taylor", "brinkman", "pal_rhodes", "woelflin", "richardson"}


def _relative_viscosity(
    model: str,
    phi: float,
    continuous: float,
    dispersed: float,
    max_packing_fraction: float,
    tightness_factor: float,
) -> float:
    if model == "einstein":
        return 1.0 + 2.5 * phi
    if model == "taylor":
        return 1.0 + 2.5 * phi * (dispersed + 0.4 * continuous) / (dispersed + continuous)
    if model == "brinkman":
        return math.pow(1.0 - phi, -2.5)
    if model == "pal_rhodes":
        if phi >= max_packing_fraction:
            return 1.0e6
        return math.pow(1.0 - phi / max_packing_fraction, -2.5 * max_packing_fraction)
    if model == "woelflin":
        k = 2.5 if tightness_factor < 0.33 else 4.0 if tightness_factor < 0.67 else 6.0
        return math.exp(k * phi)
    return math.exp((2.5 + 5.0 * tightness_factor) * phi)


def emulsion_viscosity(
    oil_viscosity: Q,
    water_viscosity: Q,
    oil_density: Q,
    water_density: Q,
    water_cut: float,
    max_packing_fraction: float,
    model: str,
    tightness_factor: float,
    demulsifier_present: bool,
    demulsifier_efficiency: float,
) -> EmulsionViscosityResult:
    """The emulsion viscosity, its inversion point and the continuous phase.

    Args:
        oil_viscosity: the oil phase's dynamic viscosity.
        water_viscosity: the water phase's dynamic viscosity.
        oil_density: the oil phase's density.
        water_density: the water phase's density.
        water_cut: the water volume fraction.
        max_packing_fraction: the maximum dispersed packing fraction.
        model: the relative-viscosity model.
        tightness_factor: the emulsion tightness.
        demulsifier_present: whether a demulsifier is dosed.
        demulsifier_efficiency: the demulsifier efficiency.

    Returns:
        The effective and relative viscosities, the inversion point and the phase type.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if the model name is unknown.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = emulsion_viscosity(
        ...     q(0.01, "Pa*s"), q(0.0005, "Pa*s"), q(850.0, "kg/m**3"),
        ...     q(1020.0, "kg/m**3"), 0.30, 0.74, "pal_rhodes", 0.5, False, 0.0)
        >>> r.emulsion_type
        <EmulsionType.WATER_IN_OIL: 'water_in_oil'>
    """
    from azoth._registry_gen import spec as _spec_for

    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "oil_viscosity": oil_viscosity,
            "water_viscosity": water_viscosity,
            "oil_density": oil_density,
            "water_density": water_density,
        }.items()
    }
    values["water_cut"] = water_cut
    values["max_packing_fraction"] = max_packing_fraction
    apply_checks(checks.on_input, values.get, warnings)

    if model not in MODELS:
        raise InvalidInputError("model", f"`{model}` is not a declared model")

    viscosity_ratio_phases = values["oil_viscosity"] / values["water_viscosity"]
    density_ratio = values["water_density"] / values["oil_density"]
    inversion = 1.0 / (1.0 + math.pow(viscosity_ratio_phases, 0.25) * math.pow(density_ratio, 0.5))
    inversion = min(0.90, max(0.20, inversion))
    inversion = min(0.95, max(0.15, inversion + (tightness_factor - 0.5) * 0.10))

    is_inverted = water_cut >= inversion
    if is_inverted:
        continuous = values["water_viscosity"]
        dispersed = values["oil_viscosity"]
        phi = 1.0 - water_cut
    else:
        continuous = values["oil_viscosity"]
        dispersed = values["water_viscosity"]
        phi = water_cut

    if phi <= 0.0:
        relative = 1.0
    elif phi >= 1.0:
        relative = dispersed / continuous
    else:
        relative = _relative_viscosity(
            model, phi, continuous, dispersed, max_packing_fraction, tightness_factor
        )

    effective = continuous * relative
    if tightness_factor > 0.5:
        effective *= 1.0 + (tightness_factor - 0.5) * 2.0
    if demulsifier_present and demulsifier_efficiency > 0.0:
        effective *= 1.0 - demulsifier_efficiency * 0.5

    return EmulsionViscosityResult(
        effective_viscosity=from_si(effective, "Pa*s"),
        relative_viscosity=relative,
        viscosity_ratio=effective / continuous,
        inversion_water_cut=inversion,
        emulsion_type=EmulsionType.OIL_IN_WATER if is_inverted else EmulsionType.WATER_IN_OIL,
        is_inverted=is_inverted,
        warnings=tuple(warnings),
    )
