"""``eos.scale_kinetics`` - NeqSim's ``ScaleKinetics``.

A classical-nucleation induction time, a parabolic surface-growth rate and a
mass-transport rate over the supersaturation a saturation index states, with the
rate-controlling regime reported.

Spec: ``specs/calcs/eos/scale_kinetics.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ScaleKineticsResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.scale_kinetics"

#: The class's year and calcite molar volume, in seconds and cubic metres per mole.
SECONDS_PER_YEAR = 3.1536e7
CALCITE_MOLAR_VOLUME = 37.0e-6


def scale_kinetics(
    saturation_index: float,
    nucleation_constant: Q,
    nucleation_shape_factor: float,
    surface_rate_constant: Q,
    surface_reaction_order: float,
    mass_transfer_coefficient: Q,
    bulk_concentration: Q,
    equilibrium_concentration: Q,
) -> ScaleKineticsResult:
    """The screening kinetics of a supersaturated solution.

    Args:
        saturation_index: the thermodynamic index, ``log10(IAP/Ksp)``.
        nucleation_constant: the nucleation pre-exponential constant.
        nucleation_shape_factor: the nucleation exponential shape factor ``B``.
        surface_rate_constant: the surface-reaction rate constant.
        surface_reaction_order: the reaction order in ``S - 1``.
        mass_transfer_coefficient: the mass-transfer coefficient; zero means no data.
        bulk_concentration: the bulk concentration of the limiting scaling ion.
        equilibrium_concentration: its equilibrium concentration.

    Returns:
        The induction time, both growth rates, the effective rate and the regime.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = scale_kinetics(
        ...     1.0, q(3.6, "s"), 3.0, q(1.5854895991882293e-10, "m/s"), 2.0,
        ...     q(1.0e-5, "m/s"), q(0.0, "mol/m**3"), q(0.0, "mol/m**3"))
        >>> r.limiting_regime
        'REACTION'
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "nucleation_constant": input_to_si(spec, "nucleation_constant", nucleation_constant),
        "surface_rate_constant": input_to_si(spec, "surface_rate_constant", surface_rate_constant),
        "mass_transfer_coefficient": input_to_si(
            spec, "mass_transfer_coefficient", mass_transfer_coefficient
        ),
        "bulk_concentration": input_to_si(spec, "bulk_concentration", bulk_concentration),
        "equilibrium_concentration": input_to_si(
            spec, "equilibrium_concentration", equilibrium_concentration
        ),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # The setters' clamps: a negative duration, rate or concentration is zero.
    constant = max(0.0, values["nucleation_constant"])
    shape = max(0.0, nucleation_shape_factor)
    rate = max(0.0, values["surface_rate_constant"])
    coefficient = max(0.0, values["mass_transfer_coefficient"])
    bulk = max(0.0, values["bulk_concentration"])
    equilibrium = max(0.0, values["equilibrium_concentration"])

    supersaturation_ratio = math.pow(10.0, saturation_index)

    if supersaturation_ratio <= 1.0:
        induction_time = math.inf
        surface_growth = 0.0
        transport_growth = 0.0
        effective_growth = 0.0
        regime = "NONE"
    else:
        ln_s = math.log(supersaturation_ratio)
        induction_hours = constant / 3600.0 * math.exp(shape / (ln_s * ln_s))

        rate_mm_yr = rate * 1000.0 * SECONDS_PER_YEAR
        surface_mm_yr = rate_mm_yr * math.pow(supersaturation_ratio - 1.0, surface_reaction_order)

        driving_mol_m3 = max(0.0, bulk - equilibrium)
        transport_m_s = coefficient * driving_mol_m3 * CALCITE_MOLAR_VOLUME
        transport_mm_yr = transport_m_s * 1000.0 * SECONDS_PER_YEAR

        if transport_mm_yr <= 0.0 or surface_mm_yr <= transport_mm_yr:
            effective_mm_yr, regime = surface_mm_yr, "REACTION"
        else:
            effective_mm_yr, regime = transport_mm_yr, "TRANSPORT"

        induction_time = induction_hours * 3600.0
        surface_growth = surface_mm_yr * 1e-3 / SECONDS_PER_YEAR
        transport_growth = transport_mm_yr * 1e-3 / SECONDS_PER_YEAR
        effective_growth = effective_mm_yr * 1e-3 / SECONDS_PER_YEAR

    return ScaleKineticsResult(
        supersaturation_ratio=supersaturation_ratio,
        induction_time=from_si(induction_time, "s"),
        surface_limited_growth=from_si(surface_growth, "m/s"),
        transport_limited_growth=from_si(transport_growth, "m/s"),
        effective_growth_rate=from_si(effective_growth, "m/s"),
        limiting_regime=regime,
        warnings=tuple(warnings),
    )
