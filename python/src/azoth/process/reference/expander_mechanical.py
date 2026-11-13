"""``ExpanderMechanicalDesign.calcDesign`` - the API 617 surface.

The twin of ``crates/azoth-process/src/rotating_design/expander.rs``.

The class exposes no setter: every quantity is the expander's run state and the class's
constants. The arithmetic is transcribed in the order ``calcDesign`` runs it, and the
dataclass converts the class's millimetres, bar and degrees Celsius at the boundary.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.eos import components as _components
from azoth.process.reference import capacity

#: The material limit for the wheel tip speed, m/s.
MAX_TIP_SPEED = 450.0
#: The casing design pressure's margin.
DESIGN_PRESSURE_MARGIN = 1.10
#: The casing design temperature's margin, degC.
DESIGN_TEMPERATURE_MARGIN = 30.0
#: Steel's density, kg/m3.
STEEL_DENSITY = 7850.0
#: The radial inflow turbine's velocity ratio.
VELOCITY_RATIO = 0.7
#: The casing's allowable stress, MPa.
CASING_ALLOWABLE_STRESS_MPA = 137.9
#: The shaft's allowable shear stress, Pa.
SHAFT_ALLOWABLE_SHEAR_PA = 50.0e6


@dataclass(frozen=True, slots=True)
class ExpanderMechanicalDesign:
    """``calcDesign``'s stored result, in SI."""

    expander_type: str
    load_type: str
    number_of_stages: int
    wheel_diameter: float
    tip_speed: float
    shaft_diameter: float
    rated_speed: float
    max_continuous_speed: float
    trip_speed: float
    first_critical_speed: float
    casing_pressure: float
    casing_temperature: float
    casing_wall_thickness: float
    bearing_span: float
    wheel_weight: float
    shaft_weight: float
    generator_weight: float
    seal_type: str
    bearing_type: str
    inner_diameter: float
    outer_diameter: float
    tantan_length: float
    weight_vessel: float
    weight_internals: float
    weight_piping: float
    weight_electro_instrument: float
    weight_structural_steel: float
    weight_total: float
    module_length: float
    module_width: float
    module_height: float


def design(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_t: float,
    inlet_p: float,
    outlet_t: float,
    power: float,
) -> ExpanderMechanicalDesign:
    """``calcDesign`` on the expander's run state, in SI.

    ``outlet_t`` is the run's outlet temperature and ``power`` its shaft power in watts
    (negative for an expansion). The class stores the run's isentropic efficiency but no
    sizing arithmetic reads it.

    Example:
        >>> design(["methane"], 1.0, [1.0], 320.0, 60e5, 286.0, -1.1e3).expander_type
        'radial_inflow'
    """
    mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))

    mass_flow_kg_hr = inlet_n * molar_mass * 3600.0
    # The feed's root, not a negative flash's trial root: NeqSim's `getDensity` reads the
    # system's, which is what the capacity family's `from_feed` selects.
    volume_flow_m3_hr = (
        capacity.system_volumetric_flow(
            components, inlet_t, inlet_p, inlet_z, inlet_n, from_feed=True
        )
        * 3600.0
    )
    inlet_pressure_bara = inlet_p / 1.0e5
    inlet_temperature_c = inlet_t - 273.15
    outlet_temperature_c = outlet_t - 273.15

    recovered_power_kw = abs(power / 1000.0)
    enthalpy_drop_kj_kg = (
        recovered_power_kw * 3600.0 / mass_flow_kg_hr if mass_flow_kg_hr > 0.0 else 0.0
    )

    casing_pressure_bara = inlet_pressure_bara * DESIGN_PRESSURE_MARGIN
    casing_temperature_c = (
        max(inlet_temperature_c, outlet_temperature_c) + DESIGN_TEMPERATURE_MARGIN
    )

    expander_type, load_type, number_of_stages = _select_expander_type(
        volume_flow_m3_hr, enthalpy_drop_kj_kg, recovered_power_kw
    )
    wheel = _wheel_sizing(volume_flow_m3_hr, enthalpy_drop_kj_kg, number_of_stages)
    shaft_diameter = _shaft_sizing(wheel.rated_speed, recovered_power_kw)
    rotor = _rotor_dynamics(
        wheel.wheel_diameter, shaft_diameter, wheel.rated_speed, casing_pressure_bara
    )
    casing = _casing_design(wheel.wheel_diameter, casing_pressure_bara, casing_temperature_c)
    generator_weight = _generator_sizing(recovered_power_kw, load_type)
    weights = _weights(
        wheel.wheel_diameter,
        shaft_diameter,
        rotor.bearing_span,
        casing,
        generator_weight,
        load_type,
        rotor.seal_type,
        recovered_power_kw,
    )
    module_length, module_width, module_height = _module_dimensions(
        rotor.bearing_span, casing.outer_diameter, load_type, recovered_power_kw
    )

    return ExpanderMechanicalDesign(
        expander_type=expander_type,
        load_type=load_type,
        number_of_stages=number_of_stages,
        wheel_diameter=wheel.wheel_diameter / 1000.0,
        tip_speed=wheel.tip_speed,
        shaft_diameter=shaft_diameter / 1000.0,
        rated_speed=wheel.rated_speed,
        max_continuous_speed=wheel.max_continuous_speed,
        trip_speed=wheel.trip_speed,
        first_critical_speed=rotor.first_critical_speed,
        casing_pressure=casing_pressure_bara * 1.0e5,
        casing_temperature=casing_temperature_c + 273.15,
        casing_wall_thickness=casing.wall_thickness / 1000.0,
        bearing_span=rotor.bearing_span / 1000.0,
        wheel_weight=weights.wheel,
        shaft_weight=weights.shaft,
        generator_weight=generator_weight,
        seal_type=rotor.seal_type,
        bearing_type=rotor.bearing_type,
        inner_diameter=casing.inner_diameter / 1000.0,
        outer_diameter=casing.outer_diameter / 1000.0,
        tantan_length=rotor.bearing_span * 2.0 / 1000.0,
        weight_vessel=weights.vessel,
        weight_internals=weights.internals,
        weight_piping=weights.piping,
        weight_electro_instrument=weights.electro_instrument,
        weight_structural_steel=weights.structural_steel,
        weight_total=weights.total,
        module_length=module_length,
        module_width=module_width,
        module_height=module_height,
    )


def _select_expander_type(
    volume_flow_m3_hr: float, enthalpy_drop_kj_kg: float, recovered_power_kw: float
) -> tuple[str, str, int]:
    """``selectExpanderType``: axial above the volume-flow threshold, generator above 100 kW."""
    if volume_flow_m3_hr > 100_000.0:
        expander_type = "axial"
        number_of_stages = max(0, math.ceil(enthalpy_drop_kj_kg / 50.0))
    else:
        expander_type = "radial_inflow"
        number_of_stages = 1
    load_type = "generator" if recovered_power_kw > 100.0 else "brake"
    return expander_type, load_type, number_of_stages


@dataclass(frozen=True, slots=True)
class _Wheel:
    wheel_diameter: float
    tip_speed: float
    rated_speed: float
    max_continuous_speed: float
    trip_speed: float


def _wheel_sizing(volume_flow_m3_hr: float, enthalpy_drop_kj_kg: float, stages: int) -> _Wheel:
    """``calculateWheelSizing``, in the class's millimetres and rpm."""
    if enthalpy_drop_kj_kg <= 0.0 or volume_flow_m3_hr <= 0.0:
        return _Wheel(300.0, 300.0, 20_000.0, 21_000.0, 23_100.0)

    spouting_velocity = math.sqrt(2.0 * enthalpy_drop_kj_kg * 1000.0)
    tip_speed = min(VELOCITY_RATIO * spouting_velocity, MAX_TIP_SPEED)
    volume_flow_m3_s = volume_flow_m3_hr / 3600.0

    enthalpy_per_stage = enthalpy_drop_kj_kg / max(1, stages)
    rated_speed = 70.0 * enthalpy_per_stage**0.75 / math.sqrt(volume_flow_m3_s * 60.0)
    rated_speed = min(max(rated_speed, 3000.0), 100_000.0)

    wheel_diameter = tip_speed * 60.0 / (math.pi * rated_speed) * 1000.0
    wheel_diameter = min(max(wheel_diameter, 100.0), 1000.0)
    rated_speed = tip_speed * 60.0 / (math.pi * wheel_diameter / 1000.0)

    return _Wheel(
        wheel_diameter,
        tip_speed,
        rated_speed,
        rated_speed * 1.05,
        rated_speed * 1.10,
    )


def _shaft_sizing(rated_speed: float, power_kw: float) -> float:
    """``calculateShaftSizing``, in millimetres."""
    if rated_speed <= 0.0 or power_kw <= 0.0:
        return 50.0
    torque_knm = power_kw * 60.0 / (2.0 * math.pi * rated_speed)
    torque_nm = torque_knm * 1000.0
    diameter_m = float((16.0 * torque_nm / (math.pi * SHAFT_ALLOWABLE_SHEAR_PA)) ** (1.0 / 3.0))
    diameter_mm = diameter_m * 1000.0 * 1.5
    return max(math.ceil(diameter_mm / 5.0) * 5.0, 30.0)


@dataclass(frozen=True, slots=True)
class _Rotor:
    bearing_span: float
    first_critical_speed: float
    seal_type: str
    bearing_type: str


def _rotor_dynamics(
    wheel_diameter: float, shaft_diameter: float, rated_speed: float, casing_pressure_bara: float
) -> _Rotor:
    """``calculateRotorDynamics``: the span, the Rayleigh-Ritz speed and the two bands."""
    bearing_span = wheel_diameter * 2.5
    d = shaft_diameter / 1000.0
    length = bearing_span / 1000.0
    mass_per_length = STEEL_DENSITY * math.pi * d * d / 4.0
    first_critical_speed = 946.0 * math.sqrt(d**4 / (mass_per_length * length**3))
    first_critical_speed = max(first_critical_speed, rated_speed * 0.5)

    if rated_speed > 50_000.0:
        bearing_type = "active_magnetic"
    elif rated_speed > 20_000.0:
        bearing_type = "tilting_pad"
    else:
        bearing_type = "sleeve"
    seal_type = "dry_gas" if casing_pressure_bara > 50.0 else "labyrinth"
    return _Rotor(bearing_span, first_critical_speed, seal_type, bearing_type)


@dataclass(frozen=True, slots=True)
class _Casing:
    inner_diameter: float
    outer_diameter: float
    wall_thickness: float


def _casing_design(
    wheel_diameter: float, casing_pressure_bara: float, casing_temperature_c: float
) -> _Casing:
    """``calculateCasingDesign``, in millimetres."""
    inner_diameter = wheel_diameter * 1.3
    pressure_mpa = casing_pressure_bara * 0.1
    allowable_stress = CASING_ALLOWABLE_STRESS_MPA
    if casing_temperature_c > 200.0:
        allowable_stress *= 0.85
    elif casing_temperature_c < -50.0:
        allowable_stress *= 0.9
    joint_efficiency = 0.85
    wall_thickness = (
        pressure_mpa
        * inner_diameter
        / (2.0 * allowable_stress * joint_efficiency - 0.6 * pressure_mpa)
        + 3.0
    )
    wall_thickness = max(wall_thickness, 10.0)
    return _Casing(
        inner_diameter,
        inner_diameter + 2.0 * wall_thickness,
        wall_thickness,
    )


def _generator_sizing(recovered_power_kw: float, load_type: str) -> float:
    """``calculateGeneratorSizing``, in kilograms."""
    if load_type != "generator":
        return 0.0
    if recovered_power_kw <= 100.0:
        return 100.0 + recovered_power_kw * 3.0
    if recovered_power_kw <= 1000.0:
        return 400.0 + recovered_power_kw * 2.0
    return 2400.0 + recovered_power_kw


@dataclass(frozen=True, slots=True)
class _Weights:
    wheel: float
    shaft: float
    vessel: float
    internals: float
    piping: float
    electro_instrument: float
    structural_steel: float
    total: float


def _weights(
    wheel_diameter: float,
    shaft_diameter: float,
    bearing_span: float,
    casing: _Casing,
    generator_weight: float,
    load_type: str,
    seal_type: str,
    recovered_power_kw: float,
) -> _Weights:
    """``calculateWeights``, in kilograms."""
    wheel_volume = math.pi * (wheel_diameter / 2000.0) ** 2 * (wheel_diameter / 3000.0)
    wheel = wheel_volume * STEEL_DENSITY * 0.7

    shaft_length = bearing_span * 1.5
    shaft_volume = math.pi * (shaft_diameter / 2000.0) ** 2 * (shaft_length / 1000.0)
    shaft = shaft_volume * STEEL_DENSITY

    casing_length = bearing_span * 2.0
    casing_volume = (
        math.pi
        / 4.0
        * ((casing.outer_diameter / 1000.0) ** 2 - (casing.inner_diameter / 1000.0) ** 2)
        * (casing_length / 1000.0)
    )
    vessel = casing_volume * STEEL_DENSITY * 1.5

    internals = wheel + shaft + wheel * 0.3
    equipment = vessel + internals + generator_weight
    lube_system = 200.0 + recovered_power_kw * 0.2
    seal_system = 500.0 if seal_type == "dry_gas" else 100.0
    baseplate = equipment * 0.2
    piping = equipment * 0.15
    electrical = recovered_power_kw * 0.3 if load_type == "generator" else 50.0
    structural = equipment * 0.1

    return _Weights(
        wheel=wheel,
        shaft=shaft,
        vessel=vessel,
        internals=internals,
        piping=piping + lube_system + seal_system,
        electro_instrument=electrical,
        structural_steel=structural + baseplate,
        total=equipment + lube_system + seal_system + baseplate + piping + electrical + structural,
    )


def _module_dimensions(
    bearing_span: float, outer_diameter: float, load_type: str, recovered_power_kw: float
) -> tuple[float, float, float]:
    """``calculateModuleDimensions``, in metres."""
    expander_length = bearing_span / 1000.0 * 2.5
    generator_length = (
        1.0 + math.sqrt(recovered_power_kw) * 0.03 if load_type == "generator" else 0.5
    )
    return (
        expander_length + generator_length + 2.0,
        max(2.0, outer_diameter / 1000.0 * 2.0 + 1.0),
        max(1.5, outer_diameter / 1000.0 + 1.0),
    )
