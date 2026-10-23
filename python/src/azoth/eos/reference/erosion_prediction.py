"""``eos.erosion_prediction`` - NeqSim's ``ErosionPredictionCalculator``.

The API RP 14E erosional velocity and the DNV RP O501 sand-erosion rate, with the
geometry, material and particle branches the class carries, then the cumulative wall
loss and the risk band.

Spec: ``specs/calcs/eos/erosion_prediction.toml``.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ErosionPredictionResult, ErosionRiskLevel
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.erosion_prediction"

LB_FT3_PER_KG_M3 = 0.062428
SECONDS_PER_YEAR = 31_557_600.0
REFERENCE_PARTICLE_DIAMETER = 0.00025
REFERENCE_PARTICLE_DENSITY = 2650.0

GEOMETRIES = {"straight", "elbow", "tee", "blind_tee", "reducer", "choke", "weld"}
MATERIALS = {
    "carbon_steel",
    "duplex_steel",
    "super_duplex",
    "13cr",
    "inconel",
    "titanium",
    "tungsten_carbide",
}


def _geometry_factor(geometry: str, bend_radius_ratio: float) -> float:
    if geometry == "straight":
        return 0.0
    if geometry == "elbow":
        if bend_radius_ratio <= 1.0:
            return 1.0
        if bend_radius_ratio <= 1.5:
            return 0.8
        if bend_radius_ratio <= 3.0:
            return 0.5
        return 0.3
    return {"tee": 0.6, "blind_tee": 0.3, "reducer": 0.4, "choke": 2.0, "weld": 1.2}[geometry]


def _impact_angle(geometry: str, bend_radius_ratio: float) -> float:
    if geometry == "elbow":
        return max(45.0 / bend_radius_ratio, 5.0)
    return {
        "straight": 1.0,
        "tee": 45.0,
        "blind_tee": 90.0,
        "reducer": 15.0,
        "choke": 30.0,
        "weld": 20.0,
    }[geometry]


def _eroded_area(geometry: str, diameter: float) -> float:
    radius = diameter / 2.0
    if geometry == "straight":
        return math.pi * diameter
    if geometry == "elbow":
        return math.pi * diameter * diameter * 0.1
    if geometry in ("tee", "blind_tee"):
        return math.pi * radius * radius * 0.5
    if geometry == "reducer":
        return math.pi * radius * radius * 0.3
    if geometry == "choke":
        return math.pi * radius * radius * 0.2
    return math.pi * diameter * 0.005


def _erosion_constant(material: str) -> float:
    if material in ("inconel", "titanium"):
        return 1.0e-9
    if material == "tungsten_carbide":
        return 5.0e-10
    return 2.0e-9


def _velocity_exponent(material: str) -> float:
    return 2.3 if material == "tungsten_carbide" else 2.6


def _material_density(material: str) -> float:
    if material == "titanium":
        return 4500.0
    if material == "inconel":
        return 8440.0
    return 7850.0


def erosion_prediction(
    mixture_density: Q,
    mixture_velocity: Q,
    pipe_diameter: Q,
    wall_thickness: Q,
    sand_rate: Q,
    sand_particle_diameter: Q,
    sand_particle_density: Q,
    mixture_viscosity: Q,
    geometry: str,
    bend_radius_ratio: float,
    pipe_material: str,
    api_c_factor: float,
    corrosion_allowance: Q,
    design_life: Q,
) -> ErosionPredictionResult:
    """The API RP 14E erosional velocity and the DNV RP O501 erosion rate.

    Args:
        mixture_density: the flowing mixture's density.
        mixture_velocity: the mixture velocity.
        pipe_diameter: the pipe's internal diameter.
        wall_thickness: the clean wall thickness.
        sand_rate: the produced sand rate.
        sand_particle_diameter: the sand particle diameter.
        sand_particle_density: the sand particle density.
        mixture_viscosity: the mixture's dynamic viscosity.
        geometry: the fitting the erosion is assessed at.
        bend_radius_ratio: the bend radius over the diameter.
        pipe_material: the wall material.
        api_c_factor: the API RP 14E C factor.
        corrosion_allowance: the allowance the remaining wall subtracts.
        design_life: the exposure time the cumulative loss is scaled by.

    Returns:
        The erosional velocity, ratio, rate, cumulative loss, remaining wall and risk band.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.
        InvalidInputError: if a geometry or material name is unknown.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = erosion_prediction(
        ...     q(100.0, "kg/m**3"), q(10.0, "m/s"), q(0.1524, "m"), q(0.01, "m"),
        ...     q(0.0, "kg/(m**2*s)"), q(0.00025, "m"), q(2650.0, "kg/m**3"),
        ...     q(0.001, "Pa*s"), "elbow", 1.5, "carbon_steel", 100.0,
        ...     q(0.003, "m"), q(788_940_000.0, "s"))
        >>> r.risk_level
        <ErosionRiskLevel.MEDIUM: 'medium'>
    """
    from azoth._registry_gen import spec as _spec_for

    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "mixture_density": mixture_density,
            "mixture_velocity": mixture_velocity,
            "pipe_diameter": pipe_diameter,
            "wall_thickness": wall_thickness,
            "sand_rate": sand_rate,
            "sand_particle_diameter": sand_particle_diameter,
            "sand_particle_density": sand_particle_density,
            "mixture_viscosity": mixture_viscosity,
            "corrosion_allowance": corrosion_allowance,
            "design_life": design_life,
        }.items()
    }
    values["bend_radius_ratio"] = bend_radius_ratio
    values["api_c_factor"] = api_c_factor
    apply_checks(checks.on_input, values.get, warnings)

    if geometry not in GEOMETRIES:
        raise InvalidInputError("geometry", f"`{geometry}` is not a declared geometry")
    if pipe_material not in MATERIALS:
        raise InvalidInputError("pipe_material", f"`{pipe_material}` is not a declared material")

    diameter = values["pipe_diameter"]
    velocity = values["mixture_velocity"]
    erosional = (
        api_c_factor / math.sqrt(values["mixture_density"] * LB_FT3_PER_KG_M3) * 0.3048
        if values["mixture_density"] > 0.0
        else 0.0
    )
    ratio = velocity / erosional if erosional > 0.0 else 0.0
    within_api = ratio <= 1.0

    particle_diameter = values["sand_particle_diameter"]
    stokes = (
        values["sand_particle_density"]
        * particle_diameter
        * particle_diameter
        * velocity
        / (18.0 * values["mixture_viscosity"] * diameter)
    )
    particle_ratio = 0.9 if stokes > 10.0 else 0.8 if stokes > 1.0 else 0.7
    particle_velocity = velocity * particle_ratio
    sand_rate_kg_yr = values["sand_rate"] * SECONDS_PER_YEAR

    rate_mm_yr: float
    if values["sand_rate"] > 0.0 and geometry != "straight":
        impact = math.sin(math.radians(_impact_angle(geometry, bend_radius_ratio)))
        impact_function = impact * (1.0 + 0.6 * (1.0 - impact))
        erosion_kg_kg = (
            _erosion_constant(pipe_material)
            * impact_function
            * particle_velocity ** _velocity_exponent(pipe_material)
            * (particle_diameter / REFERENCE_PARTICLE_DIAMETER) ** 0.2
            * _geometry_factor(geometry, bend_radius_ratio)
            * (values["sand_particle_density"] / REFERENCE_PARTICLE_DENSITY)
        )
        area = _eroded_area(geometry, diameter)
        rate_mm_yr = (
            erosion_kg_kg * sand_rate_kg_yr / _material_density(pipe_material) / area * 1000.0
            if area > 0.0
            else 0.0
        )
    elif values["sand_rate"] > 0.0:
        rate_mm_yr = (
            _erosion_constant(pipe_material)
            * particle_velocity ** _velocity_exponent(pipe_material)
            * sand_rate_kg_yr
            * 0.01
        )
    else:
        rate_mm_yr = 0.0

    cumulative_mm = rate_mm_yr * values["design_life"] / SECONDS_PER_YEAR
    remaining_mm = (
        values["wall_thickness"] * 1000.0
        - cumulative_mm
        - values["corrosion_allowance"] * 1000.0
    )

    if rate_mm_yr > 1.0 or ratio > 1.5:
        risk = ErosionRiskLevel.CRITICAL
    elif rate_mm_yr > 0.5 or ratio > 1.0:
        risk = ErosionRiskLevel.HIGH
    elif rate_mm_yr > 0.1 or ratio > 0.8:
        risk = ErosionRiskLevel.MEDIUM
    else:
        risk = ErosionRiskLevel.LOW

    return ErosionPredictionResult(
        erosional_velocity=from_si(erosional, "m/s"),
        velocity_ratio=ratio,
        erosion_rate=from_si(rate_mm_yr * 1.0e-3 / SECONDS_PER_YEAR, "m/s"),
        cumulative_erosion=from_si(cumulative_mm * 1.0e-3, "m"),
        remaining_wall_thickness=from_si(remaining_mm * 1.0e-3, "m"),
        within_api_limits=within_api,
        within_erosion_limits=remaining_mm > 0.0,
        risk_level=risk,
        warnings=tuple(warnings),
    )
