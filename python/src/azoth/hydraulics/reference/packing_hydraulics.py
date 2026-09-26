"""``hydraulics.packing_hydraulics`` - a packed bed's flooding, load, pressure drop and film
coefficients.

Spec: ``specs/calcs/hydraulics/packing_hydraulics.toml``

The Python twin of ``crates/azoth-hydraulics/src/packing_hydraulics.rs``, and both mirror
``PackingHydraulicsCalculator``. **Three of the four correlations are Onda's** - the wetted area
and the two film coefficients - and the fourth is Leva's pressure drop; the flooding velocity is
the Eckert GPDC fit, and the HETP is a two-resistance combination of the film coefficients with a
stripping factor of one.

# The packing table

The library registers 22 built-ins and *then* loads ``designdata/Packing.csv``, so a row whose
normalized name collides replaces the earlier registration: ``Pall-Ring-50`` resolves to the
file's **plastic** row (111.1 m**2/m**3, void 0.919, factor 180) rather than to the built-in metal
one (120.0, 0.96, 66). Both are in :mod:`azoth.hydraulics.packing`.

# The HETP is clamped, and on these states the clamp is what answers

The two-resistance combination gives transfer units that can be hundreds of metres tall - the
captured absorber states measure `336` m and `203` m - and the class reports the empirical
estimate's two-fold band instead: `1.44` m and `1.0` m. A port that published the HTU sum would be
answering a different question.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PackingHydraulicsResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.hydraulics.packing import PackingSpecification, packing_or_default

CALC_ID = "hydraulics.packing_hydraulics"

G = 9.81
WATER_VISCOSITY = 0.001
ECKERT_C0, ECKERT_C1, ECKERT_C2 = -1.668, -1.085, -0.297
MIN_FLOW_PARAMETER, MAX_FLOW_PARAMETER = 0.005, 5.0
ONDA_GAS_CONSTANT, ONDA_LIQUID_CONSTANT, ONDA_WETTING_CONSTANT = 5.23, 0.0051, -1.45
MIN_WETTED_RATIO = 0.2
MIN_WETTING_RATE_STRUCTURED, MIN_WETTING_RATE_RANDOM = 2.5e-5, 5.0e-5
DESIGN_MIN_FLOOD, DESIGN_MAX_FLOOD = 40.0, 80.0
MIN_HETP_M = 0.1


def packing_hydraulics(
    packing: str,
    column_diameter: Q,
    packed_height: Q,
    vapor_mass_flow: Q,
    liquid_mass_flow: Q,
    vapor_density: Q,
    liquid_density: Q,
    vapor_viscosity: Q,
    liquid_viscosity: Q,
    surface_tension: Q,
    vapor_diffusivity: Q,
    liquid_diffusivity: Q,
    hydraulic_capacity_factor: float,
) -> PackingHydraulicsResult:
    """A packed bed's hydraulics at a stated state.

    Args:
        packing: the packing's name or an alias; a name nothing matches is ``Pall-Ring-50``.
        column_diameter: the column's internal diameter.
        packed_height: the packed height.
        vapor_mass_flow: the vapour's mass flow.
        liquid_mass_flow: the liquid's mass flow.
        vapor_density: the vapour's mass density.
        liquid_density: the liquid's mass density.
        vapor_viscosity: the vapour's dynamic viscosity.
        liquid_viscosity: the liquid's dynamic viscosity.
        surface_tension: the interfacial surface tension.
        vapor_diffusivity: the vapour's diffusivity.
        liquid_diffusivity: the liquid's diffusivity.
        hydraulic_capacity_factor: the packing's relative hydraulic capacity.

    Returns:
        The bed's rates, its two film coefficients, the HTU ladder and the HETP, with the
        wetting and design verdicts.

    Raises:
        OutOfRangeError: if the diameter, either mass flow, either density or either viscosity is
            not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = packing_hydraulics(
        ...     "Pall-Ring-50", q(1.0, "m"), q(5.0, "m"), q(0.35, "kg/s"), q(3.5, "kg/s"),
        ...     q(45.0, "kg/m**3"), q(990.0, "kg/m**3"), q(1.8e-5, "Pa*s"), q(6.5e-4, "Pa*s"),
        ...     q(0.072, "N/m"), q(3.3e-7, "m**2/s"), q(1.9e-9, "m**2/s"), 1.0)
        >>> round(r.wetted_area, 10)
        47.6320858192
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "column_diameter": input_to_si(spec, "column_diameter", column_diameter),
        "packed_height": input_to_si(spec, "packed_height", packed_height),
        "vapor_mass_flow": input_to_si(spec, "vapor_mass_flow", vapor_mass_flow),
        "liquid_mass_flow": input_to_si(spec, "liquid_mass_flow", liquid_mass_flow),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "liquid_density": input_to_si(spec, "liquid_density", liquid_density),
        "vapor_viscosity": input_to_si(spec, "vapor_viscosity", vapor_viscosity),
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
        "surface_tension": input_to_si(spec, "surface_tension", surface_tension),
        "vapor_diffusivity": input_to_si(spec, "vapor_diffusivity", vapor_diffusivity),
        "liquid_diffusivity": input_to_si(spec, "liquid_diffusivity", liquid_diffusivity),
    }
    apply_checks(
        checks.on_input,
        {**values, "hydraulic_capacity_factor": hydraulic_capacity_factor}.get,
        warnings,
    )

    resolved = packing_or_default(packing)
    area = math.pi / 4.0 * values["column_diameter"] ** 2

    # ---- The flooding velocity: the Eckert fit, solved for the velocity.
    flow_parameter = 0.0
    if values["vapor_mass_flow"] > 0.0:
        flow_parameter = (values["liquid_mass_flow"] / values["vapor_mass_flow"]) * math.sqrt(
            values["vapor_density"] / values["liquid_density"]
        )
    flow_parameter = min(max(flow_parameter, MIN_FLOW_PARAMETER), MAX_FLOW_PARAMETER)
    x = math.log10(flow_parameter)
    y_flood = math.pow(10.0, ECKERT_C0 + ECKERT_C1 * x + ECKERT_C2 * x * x)
    rho_factor = values["vapor_density"] / (
        G * (values["liquid_density"] - values["vapor_density"])
    )
    mu_factor = math.pow(values["liquid_viscosity"] / WATER_VISCOSITY, 0.1)
    flooding_velocity = (
        math.sqrt(max(y_flood / (resolved.packing_factor * rho_factor * mu_factor), 0.0))
        * hydraulic_capacity_factor
    )

    vapor_velocity = (values["vapor_mass_flow"] / values["vapor_density"]) / area
    f_factor = vapor_velocity * math.sqrt(values["vapor_density"])
    liquid_velocity = (values["liquid_mass_flow"] / values["liquid_density"]) / area
    percent_flood = (vapor_velocity / flooding_velocity) * 100.0 if flooding_velocity > 0.0 else 0.0

    # ---- Leva's pressure drop.
    liquid_loading = values["liquid_mass_flow"] / area
    dry = (
        resolved.packing_factor
        * values["vapor_density"]
        * vapor_velocity
        * vapor_velocity
        / resolved.void_fraction**3
        * 0.01
    )
    leva_c = 0.010 if resolved.category.lower() == "structured" else 0.015
    pressure_drop_per_meter = max(dry * math.pow(10.0, leva_c * liquid_loading), 0.0)

    # ---- Onda's wetted area and the two film coefficients.
    a = resolved.specific_surface_area
    re_l = liquid_velocity * values["liquid_density"] / (values["liquid_viscosity"] * a)
    fr_l = liquid_velocity * liquid_velocity * a / G
    we_l = (
        liquid_velocity
        * liquid_velocity
        * values["liquid_density"]
        / (values["surface_tension"] * a)
    )
    sigma_ratio = resolved.critical_surface_tension / values["surface_tension"]
    exponent = (
        ONDA_WETTING_CONSTANT
        * math.pow(sigma_ratio, 0.75)
        * math.pow(re_l, 0.1)
        * math.pow(fr_l, -0.05)
        * math.pow(we_l, 0.2)
    )
    wetted_area = min(max(1.0 - math.exp(exponent), MIN_WETTED_RATIO), 1.0) * a

    diameter = effective_packing_diameter(resolved)
    re_v = (values["vapor_mass_flow"] / area) / (a * values["vapor_viscosity"])
    sc_v = values["vapor_viscosity"] / (values["vapor_density"] * values["vapor_diffusivity"])
    k_ga = (
        ONDA_GAS_CONSTANT
        * (a * values["vapor_diffusivity"] / (diameter * diameter))
        * math.pow(re_v, 0.7)
        * math.pow(sc_v, 1.0 / 3.0)
    )
    re_l2 = (values["liquid_mass_flow"] / area) / (wetted_area * values["liquid_viscosity"])
    sc_l = values["liquid_viscosity"] / (values["liquid_density"] * values["liquid_diffusivity"])
    scale = math.pow(values["liquid_density"] / (values["liquid_viscosity"] * G), 1.0 / 3.0)
    k_la = (
        ONDA_LIQUID_CONSTANT
        * math.pow(re_l2, 2.0 / 3.0)
        * math.pow(sc_l, -0.5)
        * math.pow(a * diameter, 0.4)
        / scale
        * wetted_area
    )

    # ---- The HETP's three rungs.
    minimum_wetting_rate = (
        MIN_WETTING_RATE_STRUCTURED
        if resolved.category.lower() == "structured"
        else MIN_WETTING_RATE_RANDOM
    )
    wetting_rate = (values["liquid_mass_flow"] / values["liquid_density"]) / (area * a)
    wetting_ok = wetting_rate >= minimum_wetting_rate

    estimated = estimate_hetp(resolved, values["column_diameter"])
    if (
        k_ga <= 0.0
        or k_la <= 0.0
        or values["vapor_mass_flow"] <= 0.0
        or values["liquid_mass_flow"] <= 0.0
    ):
        htu_g = htu_l = htu_og = 0.0
        hetp = estimated
    else:
        htu_g = (values["vapor_mass_flow"] / area) / (k_ga * values["vapor_density"])
        htu_l = (values["liquid_mass_flow"] / area) / k_la
        htu_og = htu_g + htu_l
        hetp = min(max(htu_og, MIN_HETP_M), estimated * 2.0)
    theoretical_stages = values["packed_height"] / hetp if hetp > 0.0 else 0.0
    design_ok = wetting_ok and DESIGN_MIN_FLOOD <= percent_flood <= DESIGN_MAX_FLOOD

    apply_checks(
        checks.derived,
        lambda name: (
            percent_flood
            if name == "percent_flood"
            else (wetted_area if name == "wetted_area" else None)
        ),
        warnings,
    )

    return PackingHydraulicsResult(
        packing_name=resolved.name,
        packing_category=resolved.category,
        specific_surface_area=resolved.specific_surface_area,
        void_fraction=resolved.void_fraction,
        packing_factor=resolved.packing_factor,
        flooding_velocity=flooding_velocity,
        vapor_velocity=vapor_velocity,
        liquid_velocity=liquid_velocity,
        f_factor=f_factor,
        percent_flood=percent_flood,
        pressure_drop_per_meter=from_si(pressure_drop_per_meter, "Pa"),
        total_pressure_drop=from_si(pressure_drop_per_meter * values["packed_height"], "Pa"),
        wetted_area=wetted_area,
        k_ga=k_ga,
        k_la=k_la,
        htu_g=htu_g,
        htu_l=htu_l,
        htu_og=htu_og,
        hetp=from_si(hetp, "m"),
        theoretical_stages=theoretical_stages,
        wetting_rate=wetting_rate,
        minimum_wetting_rate=minimum_wetting_rate,
        wetting_ok=wetting_ok,
        design_ok=design_ok,
        warnings=tuple(warnings),
    )


def estimate_hetp(packing: PackingSpecification, column_diameter: float) -> float:
    """``estimateHETP``: the empirical band the two-resistance HETP is clamped to."""
    if packing.category.lower() == "structured":
        return max(100.0 / packing.specific_surface_area + 0.1, 0.15)
    estimate = 0.12 + 0.012 * packing.nominal_size_mm
    if column_diameter > 1.0:
        estimate *= 1.0 + 0.1 * (column_diameter - 1.0)
    return max(estimate, 0.15)


def effective_packing_diameter(packing: PackingSpecification) -> float:
    """``getEffectivePackingDiameter``: the nominal size, or ``4 eps / a`` where there is none."""
    if packing.nominal_size_mm > 0.0:
        return packing.nominal_size_mm / 1000.0
    if packing.specific_surface_area > 0.0 and packing.void_fraction > 0.0:
        return max(4.0 * packing.void_fraction / packing.specific_surface_area, 1.0e-4)
    return 0.025
