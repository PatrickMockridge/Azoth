"""Heat transfer calculations.

The first namespace beside :mod:`azoth.hydraulics`. Its existence is the point:
everything here - the spec, the generated range checks, the test cases, the
documentation - flows from the same ``specs/calcs/`` tree and the same tooling,
with no hydraulics-shaped assumption left in the pipeline.

* :func:`conduction_plane_wall` - steady conduction through a slab

# Which implementation answers

As in the hydraulics namespace, each function dispatches to the Rust extension
when it is built and to :mod:`azoth.thermal.reference` otherwise. Both are always
reachable - see :func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from typing import Any

from azoth._dispatch import resolve
from azoth.core.result import (
    BoilingHeatTransferResult,
    ConductionPlaneWallResult,
    FoulingModelResult,
    LmtdCorrectionFactorResult,
    ShahCondensationResult,
    ShellAndTubeRatingResult,
    ShellAndTubeVibrationResult,
    TemaStandardResult,
    TubeInsertModelResult,
    TwoPhasePressureDropResult,
    ZoneRatingResult,
)
from azoth.core.units import Q

__all__ = [
    "boiling_heat_transfer",
    "conduction_plane_wall",
    "fouling_model",
    "lmtd_correction_factor",
    "shah_condensation",
    "shell_and_tube_rating",
    "shell_and_tube_vibration",
    "tema_standard",
    "tube_insert_model",
    "two_phase_pressure_drop",
    "zone_rating",
]

_BOILING_HEAT_TRANSFER = "thermal.boiling_heat_transfer"
_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"
_FOULING_MODEL = "thermal.fouling_model"
_LMTD_CORRECTION_FACTOR = "thermal.lmtd_correction_factor"
_SHAH_CONDENSATION = "thermal.shah_condensation"
_SHELL_AND_TUBE_RATING = "thermal.shell_and_tube_rating"
_SHELL_AND_TUBE_VIBRATION = "thermal.shell_and_tube_vibration"
_TEMA_STANDARD = "thermal.tema_standard"
_TUBE_INSERT_MODEL = "thermal.tube_insert_model"
_TWO_PHASE_PRESSURE_DROP = "thermal.two_phase_pressure_drop"
_ZONE_RATING = "thermal.zone_rating"


def conduction_plane_wall(k: Q, A: Q, dT: Q, L: Q) -> ConductionPlaneWallResult:
    """Steady heat flow through a plane wall.

    ``dT`` is a temperature *difference*, not an absolute temperature: pass
    ``Q(30, "delta_degC")`` or ``Q(30, "K")``. ``q`` follows the sign of ``dT``.

    Raises:
        OutOfRangeError: if ``k``, ``A`` or ``L`` is not positive.

    See :func:`azoth.thermal.reference.conduction_plane_wall`.
    """
    return resolve(_CONDUCTION_PLANE_WALL)(k=k, A=A, dT=dT, L=L)  # type: ignore[no-any-return]


def shell_and_tube_rating(
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
    tube_density: Q,
    tube_viscosity: Q,
    tube_cp: Q,
    tube_conductivity: Q,
    tube_mass_flow: Q,
    shell_density: Q,
    shell_viscosity: Q,
    shell_cp: Q,
    shell_conductivity: Q,
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
) -> ShellAndTubeRatingResult:
    """A shell-and-tube exchanger's thermal-hydraulic rating.

    ``fouling_tube`` and ``fouling_shell`` are in m**2*K/W, declared dimensionless
    because no vocabulary unit expresses the area-specific resistance.

    Raises:
        InvalidInputError: for a state the class's ``calculateStrict`` refuses.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.shell_and_tube_rating`.
    """
    return resolve(_SHELL_AND_TUBE_RATING)(  # type: ignore[no-any-return]
        tube_outside_diameter=tube_outside_diameter,
        tube_inside_diameter=tube_inside_diameter,
        tube_length=tube_length,
        tube_count=tube_count,
        tube_passes=tube_passes,
        tube_pitch=tube_pitch,
        triangular_pitch=triangular_pitch,
        shell_inside_diameter=shell_inside_diameter,
        baffle_spacing=baffle_spacing,
        baffle_count=baffle_count,
        baffle_cut=baffle_cut,
        tube_wall_conductivity=tube_wall_conductivity,
        tube_density=tube_density,
        tube_viscosity=tube_viscosity,
        tube_cp=tube_cp,
        tube_conductivity=tube_conductivity,
        tube_mass_flow=tube_mass_flow,
        shell_density=shell_density,
        shell_viscosity=shell_viscosity,
        shell_cp=shell_cp,
        shell_conductivity=shell_conductivity,
        shell_mass_flow=shell_mass_flow,
        fouling_tube=fouling_tube,
        fouling_shell=fouling_shell,
        shell_side_method=shell_side_method,
        tube_fouling_layer_thickness=tube_fouling_layer_thickness,
        tube_fouling_layer_conductivity=tube_fouling_layer_conductivity,
        shell_viscosity_wall=shell_viscosity_wall,
        has_sealing=has_sealing,
        sealing_pairs=sealing_pairs,
        tube_to_baffle_clearance=tube_to_baffle_clearance,
        shell_to_baffle_clearance=shell_to_baffle_clearance,
        bypass_area=bypass_area,
    )


def shell_and_tube_vibration(
    tube_outside_diameter: Q,
    tube_inside_diameter: Q,
    unsupported_span: Q,
    tube_pitch: Q,
    tube_material_e: Q,
    tube_density: Q,
    crossflow_velocity: Q,
    shell_fluid_density: Q,
    tube_fluid_density: Q,
    shell_inside_diameter: Q,
    sonic_velocity: Q,
    damping_ratio: float,
    triangular_pitch: bool,
) -> ShellAndTubeVibrationResult:
    """A shell-and-tube exchanger's tube-vibration screen.

    Raises:
        InvalidInputError: for an inside diameter not smaller than the outside one, a
            negative damping ratio or velocity.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.shell_and_tube_vibration`.
    """
    return resolve(_SHELL_AND_TUBE_VIBRATION)(  # type: ignore[no-any-return]
        tube_outside_diameter=tube_outside_diameter,
        tube_inside_diameter=tube_inside_diameter,
        unsupported_span=unsupported_span,
        tube_pitch=tube_pitch,
        tube_material_e=tube_material_e,
        tube_density=tube_density,
        crossflow_velocity=crossflow_velocity,
        shell_fluid_density=shell_fluid_density,
        tube_fluid_density=tube_fluid_density,
        shell_inside_diameter=shell_inside_diameter,
        sonic_velocity=sonic_velocity,
        damping_ratio=damping_ratio,
        triangular_pitch=triangular_pitch,
    )


def lmtd_correction_factor(
    hot_inlet_temperature: Q,
    hot_outlet_temperature: Q,
    cold_inlet_temperature: Q,
    cold_outlet_temperature: Q,
    shell_passes: float | None = None,
) -> LmtdCorrectionFactorResult:
    """The LMTD correction factor for a shell-and-tube exchanger.

    Raises:
        InvalidInputError: for a non-integral shell-pass count.
        OutOfRangeError: for a temperature outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.lmtd_correction_factor`.
    """
    return resolve(_LMTD_CORRECTION_FACTOR)(  # type: ignore[no-any-return]
        hot_inlet_temperature=hot_inlet_temperature,
        hot_outlet_temperature=hot_outlet_temperature,
        cold_inlet_temperature=cold_inlet_temperature,
        cold_outlet_temperature=cold_outlet_temperature,
        shell_passes=shell_passes,
    )


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
    """Shah's in-tube condensation correlation, composed as the zone analysis uses it.

    Raises:
        InvalidInputError: for a non-integral interval count.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.shah_condensation`.
    """
    return resolve(_SHAH_CONDENSATION)(  # type: ignore[no-any-return]
        mass_flux=mass_flux,
        tube_inner_diameter=tube_inner_diameter,
        liquid_density=liquid_density,
        liquid_viscosity=liquid_viscosity,
        liquid_cp=liquid_cp,
        liquid_conductivity=liquid_conductivity,
        vapor_density=vapor_density,
        reduced_pressure=reduced_pressure,
        vapor_quality=vapor_quality,
        quality_in=quality_in,
        quality_out=quality_out,
        intervals=intervals,
    )


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
    """Chen and Gungor-Winterton flow boiling, composed as the zone analysis uses it.

    Raises:
        InvalidInputError: for a non-integral interval count.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.boiling_heat_transfer`.
    """
    return resolve(_BOILING_HEAT_TRANSFER)(  # type: ignore[no-any-return]
        mass_flux=mass_flux,
        tube_inner_diameter=tube_inner_diameter,
        liquid_density=liquid_density,
        vapor_density=vapor_density,
        liquid_viscosity=liquid_viscosity,
        vapor_viscosity=vapor_viscosity,
        liquid_cp=liquid_cp,
        liquid_conductivity=liquid_conductivity,
        surface_tension=surface_tension,
        heat_of_vaporization=heat_of_vaporization,
        wall_superheat=wall_superheat,
        saturation_pressure_difference=saturation_pressure_difference,
        heat_flux=heat_flux,
        vapor_quality=vapor_quality,
        quality_in=quality_in,
        quality_out=quality_out,
        intervals=intervals,
        horizontal=horizontal,
    )


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
    """Friedel and Muller-Steinhagen-Heck two-phase pressure drop.

    Raises:
        InvalidInputError: for a non-integral interval count.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.two_phase_pressure_drop`.
    """
    return resolve(_TWO_PHASE_PRESSURE_DROP)(  # type: ignore[no-any-return]
        mass_flux=mass_flux,
        vapor_quality=vapor_quality,
        tube_inner_diameter=tube_inner_diameter,
        tube_length=tube_length,
        liquid_density=liquid_density,
        vapor_density=vapor_density,
        liquid_viscosity=liquid_viscosity,
        vapor_viscosity=vapor_viscosity,
        surface_tension=surface_tension,
        quality_in=quality_in,
        quality_out=quality_out,
        intervals=intervals,
    )


def fouling_model(
    model: str,
    velocity: Q,
    fluid_density: Q,
    fluid_viscosity: Q,
    wall_temperature: Q,
    tube_inner_diameter: Q,
    alpha: float,
    beta: float,
    activation_energy: Q,
    gamma: float,
    fixed_fouling_resistance: float,
    rf_max: float,
    time_constant: Q,
    time_hours: Q,
    target_fouling_resistance: float,
    threshold_rate: float,
    euler_steps: float | None = None,
    euler_step_hours: Q | None = None,
) -> FoulingModelResult:
    """Fixed, Kern-Seaton and Ebert-Panchal heat-exchanger fouling.

    ``alpha``, ``beta``, ``gamma``, the two resistances, the target and the threshold rate
    are dimensionless: the vocabulary carries no inverse-pressure-time or area-specific
    resistance unit, and their SI units are in the spec's descriptions.

    Raises:
        InvalidInputError: for a non-integral or negative Euler step count.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.fouling_model`.
    """
    return resolve(_FOULING_MODEL)(  # type: ignore[no-any-return]
        model=model,
        velocity=velocity,
        fluid_density=fluid_density,
        fluid_viscosity=fluid_viscosity,
        wall_temperature=wall_temperature,
        tube_inner_diameter=tube_inner_diameter,
        alpha=alpha,
        beta=beta,
        activation_energy=activation_energy,
        gamma=gamma,
        fixed_fouling_resistance=fixed_fouling_resistance,
        rf_max=rf_max,
        time_constant=time_constant,
        time_hours=time_hours,
        target_fouling_resistance=target_fouling_resistance,
        threshold_rate=threshold_rate,
        euler_steps=euler_steps,
        euler_step_hours=euler_step_hours,
    )


def tube_insert_model(
    insert: str,
    reynolds_number: float,
    prandtl_number: float,
    plain_heat_transfer_coefficient: Q,
    plain_pressure_drop: Q,
    twist_ratio: float | None = None,
    tape_thickness: Q | None = None,
    tape_width_fraction: float | None = None,
    wire_diameter: Q | None = None,
    wire_coil_pitch: Q | None = None,
    matrix_density: float | None = None,
    roughness_ratio: float | None = None,
    helix_angle: float | None = None,
) -> TubeInsertModelResult:
    """Twisted tape, wire matrix and coiled wire enhancement.

    Every insert parameter is optional and carries the class's default when absent;
    ``helix_angle`` is in degrees and ``roughness_ratio`` and the fractions are dimensionless.

    Raises:
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.tube_insert_model`.
    """
    return resolve(_TUBE_INSERT_MODEL)(  # type: ignore[no-any-return]
        insert=insert,
        reynolds_number=reynolds_number,
        prandtl_number=prandtl_number,
        plain_heat_transfer_coefficient=plain_heat_transfer_coefficient,
        plain_pressure_drop=plain_pressure_drop,
        twist_ratio=twist_ratio,
        tape_thickness=tape_thickness,
        tape_width_fraction=tape_width_fraction,
        wire_diameter=wire_diameter,
        wire_coil_pitch=wire_coil_pitch,
        matrix_density=matrix_density,
        roughness_ratio=roughness_ratio,
        helix_angle=helix_angle,
    )


def tema_standard(
    front_head: str,
    shell: str,
    rear_head: str,
    tema_class: str,
    tube_size: str,
    pitch_pattern: str,
    baffle_type: str,
    shell_inside_diameter: Q,
    tube_pitch: Q,
    tube_passes: float,
    needs_cleanable: bool,
    large_temperature_difference: bool,
    high_pressure: bool,
    hazardous: bool,
    tube_material: str | None = None,
) -> TemaStandardResult:
    """TEMA configuration, class limits and layout rules.

    The head and class letters are the class's own; ``tube_material`` is matched by
    upper-cased substring for the unsupported span and is optional.

    Raises:
        InvalidInputError: for a non-integral pass count or an undeclared enum value.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.tema_standard`.
    """
    return resolve(_TEMA_STANDARD)(  # type: ignore[no-any-return]
        front_head=front_head,
        shell=shell,
        rear_head=rear_head,
        tema_class=tema_class,
        tube_size=tube_size,
        pitch_pattern=pitch_pattern,
        baffle_type=baffle_type,
        shell_inside_diameter=shell_inside_diameter,
        tube_pitch=tube_pitch,
        tube_passes=tube_passes,
        needs_cleanable=needs_cleanable,
        large_temperature_difference=large_temperature_difference,
        high_pressure=high_pressure,
        hazardous=hazardous,
        tube_material=tube_material,
    )


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
    """The shell-and-tube rating run zone by zone.

    Each zone states its own eight transport properties and its own LMTD; the mass flows,
    fouling, deposit, method and geometry are the base request's. The zone rows are
    mappings in the spec's declared units.

    Raises:
        InvalidInputError: for a base state the strict rating refuses, and for a zone whose
            transport properties, duty or LMTD are not finite or not positive.
        OutOfRangeError: for an input outside the spec's declared bounds.

    See :func:`azoth.thermal.reference.zone_rating`.
    """
    return resolve(_ZONE_RATING)(  # type: ignore[no-any-return]
        tube_outside_diameter=tube_outside_diameter,
        tube_inside_diameter=tube_inside_diameter,
        tube_length=tube_length,
        tube_count=tube_count,
        tube_passes=tube_passes,
        tube_pitch=tube_pitch,
        triangular_pitch=triangular_pitch,
        shell_inside_diameter=shell_inside_diameter,
        baffle_spacing=baffle_spacing,
        baffle_count=baffle_count,
        baffle_cut=baffle_cut,
        tube_wall_conductivity=tube_wall_conductivity,
        tube_mass_flow=tube_mass_flow,
        shell_mass_flow=shell_mass_flow,
        fouling_tube=fouling_tube,
        fouling_shell=fouling_shell,
        shell_side_method=shell_side_method,
        tube_fouling_layer_thickness=tube_fouling_layer_thickness,
        tube_fouling_layer_conductivity=tube_fouling_layer_conductivity,
        shell_viscosity_wall=shell_viscosity_wall,
        has_sealing=has_sealing,
        sealing_pairs=sealing_pairs,
        tube_to_baffle_clearance=tube_to_baffle_clearance,
        shell_to_baffle_clearance=shell_to_baffle_clearance,
        bypass_area=bypass_area,
        zones=zones,
    )
