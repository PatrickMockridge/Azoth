"""``TurboExpanderCompressorMechanicalDesign.calcDesign`` - the coupled machine's sizing.

The twin of ``crates/azoth-process/src/rotating_design/turbo_expander_compressor.rs``.

The design reads both ends: the expander's states size its wheel and casing, the
compressor's suction and discharge size its impeller and casing, and one shaft takes the
two powers. The class's auxiliary surfaces - shear pins, the seal-gas envelope, the thrust
balance and the anti-surge report - are configured fields and are not published.

**The equipment's design field is null until ``initTECMechanicalDesign()``**, and nothing in
``run`` constructs it; the probe calls it, and this function is the object it builds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.eos import components as _components
from azoth.process.reference import capacity

#: The expander's and compressor's tip-speed caps, m/s.
MAX_EXPANDER_TIP_SPEED = 450.0
MAX_COMPRESSOR_TIP_SPEED = 350.0
#: The casing design margins.
DESIGN_PRESSURE_MARGIN = 1.10
DESIGN_TEMPERATURE_MARGIN = 30.0
#: Steel's density, kg/m3.
STEEL_DENSITY = 7850.0
#: The radial inflow turbine's velocity ratio.
VELOCITY_RATIO = 0.7
#: The compressor's stage head cap, kJ/kg.
MAX_HEAD_PER_STAGE = 30.0
#: The shaft's allowable shear stress, Pa.
ALLOWABLE_SHEAR_STRESS = 50.0e6
#: The casing's base allowable stress, MPa.
BASE_ALLOWABLE_STRESS_MPA = 137.9


@dataclass(frozen=True, slots=True)
class TurboExpanderCompressorMechanicalDesign:
    """``calcDesign``'s stored result, in SI."""

    expander_wheel_diameter: float
    expander_tip_speed: float
    expander_type: str
    expander_stages: int
    expander_design_inlet_pressure: float
    expander_design_outlet_pressure: float
    expander_casing_pressure: float
    expander_casing_temperature: float
    expander_power: float
    expander_enthalpy_drop: float
    compressor_impeller_diameter: float
    compressor_tip_speed: float
    compressor_stages: int
    compressor_polytropic_head: float
    compressor_head_per_stage: float
    compressor_casing_type: str
    compressor_design_suction_pressure: float
    compressor_design_discharge_pressure: float
    compressor_casing_pressure: float
    compressor_casing_temperature: float
    compressor_power: float
    shaft_diameter: float
    shaft_weight: float
    operating_speed: float
    max_continuous_speed: float
    trip_speed: float
    first_critical_speed: float
    second_critical_speed: float
    bearing_span: float
    bearing_type: str
    expander_seal_type: str
    compressor_seal_type: str
    gear_ratio: float
    design_torque: float
    design_shaft_stress: float
    design_surge_flow: float
    design_stonewall_flow: float
    expander_nozzle_throat_area: float
    expander_wheel_weight: float
    expander_casing_weight: float
    compressor_impeller_weight: float
    compressor_casing_weight: float
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
    expander_components: list[str],
    expander_n: float,
    expander_z: list[float],
    expander_t: float,
    expander_p: float,
    expander_outlet_t: float,
    expander_outlet_p: float,
    compressor_components: list[str],
    compressor_n: float,
    compressor_z: list[float],
    compressor_t: float,
    compressor_p: float,
    compressor_outlet_t: float,
    compressor_outlet_p: float,
    expander_power: float,
    compressor_power: float,
    compressor_head: float,
    operating_speed: float,
    gear_ratio: float = 1.0,
) -> TurboExpanderCompressorMechanicalDesign:
    """``calcDesign`` on the coupled machine's run state, in SI."""
    (expander_mixture, _) = _components.mixture_of(expander_components, eos="pr")
    expander_molar_mass = capacity._molar_mass(expander_mixture, list(expander_z))
    expander_mass_flow = expander_n * expander_molar_mass * 3600.0
    expander_volume = (
        capacity.system_volumetric_flow(
            expander_components, expander_t, expander_p, expander_z, expander_n, from_feed=True
        )
        * 3600.0
    )
    compressor_volume = (
        capacity.system_volumetric_flow(
            compressor_components,
            compressor_t,
            compressor_p,
            compressor_z,
            compressor_n,
            from_feed=True,
        )
        * 3600.0
    )
    expander_power_kw = abs(expander_power) / 1000.0
    compressor_power_kw = abs(compressor_power) / 1000.0

    expander_inlet_pressure = expander_p / 1.0e5
    expander_outlet_pressure = expander_outlet_p / 1.0e5
    expander_inlet_temperature_c = expander_t - 273.15
    expander_outlet_temperature_c = expander_outlet_t - 273.15
    compressor_inlet_pressure = compressor_p / 1.0e5
    compressor_outlet_pressure = compressor_outlet_p / 1.0e5
    compressor_inlet_temperature_c = compressor_t - 273.15
    compressor_outlet_temperature_c = compressor_outlet_t - 273.15

    expander_casing_pressure = expander_inlet_pressure * DESIGN_PRESSURE_MARGIN
    expander_casing_temperature = (
        max(expander_inlet_temperature_c, expander_outlet_temperature_c) + DESIGN_TEMPERATURE_MARGIN
    )
    compressor_casing_pressure = compressor_outlet_pressure * DESIGN_PRESSURE_MARGIN
    compressor_casing_temperature = (
        max(compressor_inlet_temperature_c, compressor_outlet_temperature_c)
        + DESIGN_TEMPERATURE_MARGIN
    )

    expander_enthalpy_drop = (
        expander_power_kw * 3600.0 / expander_mass_flow if expander_mass_flow > 0.0 else 0.0
    )

    expander = _size_expander(expander_volume, expander_enthalpy_drop, operating_speed)
    compressor = _size_compressor(
        compressor_volume, compressor_head, operating_speed, compressor_casing_pressure
    )
    shaft = _size_shaft(
        expander_power_kw + compressor_power_kw,
        operating_speed,
        expander.wheel_diameter,
        compressor.impeller_diameter,
        compressor.stages,
    )
    bearing, expander_seal, compressor_seal = _select_bearings_and_seals(
        operating_speed, expander_casing_pressure, compressor_casing_pressure
    )
    expander_casing_wall = _casing_wall(
        expander.wheel_diameter, expander_casing_pressure, expander_casing_temperature
    )
    compressor_casing_wall = _casing_wall(
        compressor.impeller_diameter, compressor_casing_pressure, compressor_casing_temperature
    )
    nozzle_area = _nozzle_throat_area(expander_volume, expander_enthalpy_drop)
    weights = _weights(
        expander.wheel_diameter,
        compressor.stages,
        compressor.impeller_diameter,
        shaft.diameter,
        shaft.bearing_span,
        expander_casing_wall,
        compressor_casing_wall,
        expander_power_kw,
        compressor_power_kw,
        compressor_seal,
    )
    module_length, module_width, module_height = _module_dimensions(
        expander.wheel_diameter,
        expander_casing_wall,
        compressor.stages,
        compressor.impeller_diameter,
        compressor_casing_wall,
    )

    return TurboExpanderCompressorMechanicalDesign(
        expander_wheel_diameter=expander.wheel_diameter / 1000.0,
        expander_tip_speed=expander.tip_speed,
        expander_type=expander.kind,
        expander_stages=expander.stages,
        expander_design_inlet_pressure=expander_inlet_pressure * 1.0e5,
        expander_design_outlet_pressure=expander_outlet_pressure * 1.0e5,
        expander_casing_pressure=expander_casing_pressure * 1.0e5,
        expander_casing_temperature=expander_casing_temperature + 273.15,
        expander_power=expander_power_kw * 1000.0,
        expander_enthalpy_drop=expander_enthalpy_drop,
        compressor_impeller_diameter=compressor.impeller_diameter / 1000.0,
        compressor_tip_speed=compressor.tip_speed,
        compressor_stages=compressor.stages,
        compressor_polytropic_head=compressor_head,
        compressor_head_per_stage=compressor.head_per_stage,
        compressor_casing_type=compressor.casing,
        compressor_design_suction_pressure=compressor_inlet_pressure * 1.0e5,
        compressor_design_discharge_pressure=compressor_outlet_pressure * 1.0e5,
        compressor_casing_pressure=compressor_casing_pressure * 1.0e5,
        compressor_casing_temperature=compressor_casing_temperature + 273.15,
        compressor_power=compressor_power_kw * 1000.0,
        shaft_diameter=shaft.diameter / 1000.0,
        shaft_weight=weights.shaft,
        operating_speed=operating_speed,
        max_continuous_speed=shaft.max_continuous_speed,
        trip_speed=shaft.trip_speed,
        first_critical_speed=shaft.first_critical_speed,
        second_critical_speed=shaft.second_critical_speed,
        bearing_span=shaft.bearing_span / 1000.0,
        bearing_type=bearing,
        expander_seal_type=expander_seal,
        compressor_seal_type=compressor_seal,
        gear_ratio=gear_ratio,
        design_torque=shaft.torque,
        design_shaft_stress=shaft.stress,
        design_surge_flow=compressor_volume * 0.65,
        design_stonewall_flow=compressor_volume * 1.25,
        expander_nozzle_throat_area=nozzle_area,
        expander_wheel_weight=weights.expander_wheel,
        expander_casing_weight=weights.expander_casing,
        compressor_impeller_weight=weights.compressor_impeller,
        compressor_casing_weight=weights.compressor_casing,
        weight_vessel=weights.vessel,
        weight_internals=weights.internals,
        weight_piping=weights.piping,
        weight_electro_instrument=weights.electrical,
        weight_structural_steel=weights.structural,
        weight_total=weights.total,
        module_length=module_length,
        module_width=module_width,
        module_height=module_height,
    )


@dataclass(frozen=True, slots=True)
class _Expander:
    wheel_diameter: float
    tip_speed: float
    kind: str
    stages: int


def _size_expander(
    volume_flow_m3_hr: float, enthalpy_drop_kj_kg: float, operating_speed: float
) -> _Expander:
    """``sizeExpanderWheel``, in millimetres."""
    if enthalpy_drop_kj_kg <= 0.0 or volume_flow_m3_hr <= 0.0:
        return _Expander(300.0, 300.0, "RADIAL_INFLOW", 1)
    if volume_flow_m3_hr > 100_000.0:
        kind = "AXIAL"
        stages = max(1, math.ceil(enthalpy_drop_kj_kg / 50.0))
    else:
        kind = "RADIAL_INFLOW"
        stages = 1
    spouting = math.sqrt(2.0 * enthalpy_drop_kj_kg * 1000.0)
    tip_speed = min(VELOCITY_RATIO * spouting, MAX_EXPANDER_TIP_SPEED)
    volume_m3_s = volume_flow_m3_hr / 3600.0
    enthalpy_per_stage = enthalpy_drop_kj_kg / stages
    calc_speed = 70.0 * enthalpy_per_stage**0.75 / math.sqrt(volume_m3_s * 60.0)
    calc_speed = min(max(calc_speed, 3000.0), 100_000.0)
    speed = operating_speed if operating_speed > 0.0 else calc_speed
    wheel_diameter = tip_speed * 60.0 / (math.pi * speed) * 1000.0
    wheel_diameter = min(max(wheel_diameter, 100.0), 1500.0)
    tip_speed = math.pi * (wheel_diameter / 1000.0) * speed / 60.0
    return _Expander(wheel_diameter, tip_speed, kind, stages)


@dataclass(frozen=True, slots=True)
class _Compressor:
    impeller_diameter: float
    tip_speed: float
    stages: int
    head_per_stage: float
    casing: str


def _size_compressor(
    volume_flow_m3_hr: float,
    polytropic_head: float,
    operating_speed: float,
    casing_design_pressure: float,
) -> _Compressor:
    """``sizeCompressorImpeller``, in millimetres."""
    if casing_design_pressure > 100.0:
        casing = "BARREL"
    elif casing_design_pressure > 40.0:
        casing = "HORIZONTALLY_SPLIT"
    else:
        casing = "VERTICALLY_SPLIT"
    if polytropic_head <= 0.0 or volume_flow_m3_hr <= 0.0:
        return _Compressor(300.0, 0.0, 1, 0.0, casing)
    stages = max(1, math.ceil(polytropic_head / MAX_HEAD_PER_STAGE))
    head_per_stage = polytropic_head / stages
    tip_speed = min(math.sqrt(head_per_stage * 1000.0 / 0.50), MAX_COMPRESSOR_TIP_SPEED)
    speed = operating_speed if operating_speed > 0.0 else 10_000.0
    impeller_diameter = tip_speed * 60.0 / (math.pi * speed) * 1000.0
    impeller_diameter = min(max(impeller_diameter, 100.0), 1500.0)
    tip_speed = math.pi * (impeller_diameter / 1000.0) * speed / 60.0
    return _Compressor(impeller_diameter, tip_speed, stages, head_per_stage, casing)


@dataclass(frozen=True, slots=True)
class _Shaft:
    diameter: float
    bearing_span: float
    torque: float
    stress: float
    max_continuous_speed: float
    trip_speed: float
    first_critical_speed: float
    second_critical_speed: float


def _size_shaft(
    total_power_kw: float,
    speed_rpm: float,
    expander_wheel_diameter: float,
    compressor_impeller_diameter: float,
    compressor_stages: int,
) -> _Shaft:
    """``sizeShaft`` and ``calculateRotorDynamics``, in millimetres, rpm and Nm."""
    torque = (
        total_power_kw * 1000.0 * 60.0 / (2.0 * math.pi * speed_rpm) if speed_rpm > 0.0 else 0.0
    )
    if speed_rpm <= 0.0 or total_power_kw <= 0.0:
        diameter = 50.0
    else:
        d = (
            float((16.0 * torque / (math.pi * ALLOWABLE_SHEAR_STRESS)) ** (1.0 / 3.0))
            * 1000.0
            * 1.5
        )
        diameter = max(math.ceil(d / 5.0) * 5.0, 30.0)
    d_m = diameter / 1000.0
    polar = math.pi * d_m**4 / 32.0
    stress = (torque * d_m / 2.0) / polar / 1.0e6 if polar > 0.0 else 0.0
    expander_end = expander_wheel_diameter * 1.3
    compressor_end = (
        compressor_stages * compressor_impeller_diameter * 0.8 + compressor_impeller_diameter
    )
    bearing_span = expander_end + compressor_end + 200.0
    max_continuous_speed = speed_rpm * 1.05
    trip_speed = max_continuous_speed * 1.05
    l_span = bearing_span / 1000.0
    mass_per_length = STEEL_DENSITY * math.pi * d_m * d_m / 4.0
    if l_span > 0.0 and mass_per_length > 0.0:
        stiffness = 48.0 * 2.1e11 * math.pi * d_m**4 / 64.0
        first_critical_speed = (
            (946.0 / (2.0 * math.pi)) * math.sqrt(stiffness / (mass_per_length * l_span**3)) * 60.0
        )
    else:
        first_critical_speed = 0.0
    return _Shaft(
        diameter,
        bearing_span,
        torque,
        stress,
        max_continuous_speed,
        trip_speed,
        first_critical_speed,
        first_critical_speed * 4.0,
    )


def _select_bearings_and_seals(
    speed_rpm: float, expander_casing_pressure: float, compressor_casing_pressure: float
) -> tuple[str, str, str]:
    """``selectBearingsAndSeals``, by speed and the two casing pressures."""
    if speed_rpm > 50_000.0:
        bearing = "Active Magnetic"
    elif speed_rpm > 20_000.0:
        bearing = "Tilting Pad"
    else:
        bearing = "Sleeve"
    expander_seal = "Dry Gas" if expander_casing_pressure > 50.0 else "Labyrinth"
    compressor_seal = "Dry Gas" if compressor_casing_pressure > 20.0 else "Labyrinth"
    return bearing, expander_seal, compressor_seal


def _derated_stress(temperature_c: float) -> float:
    """``getTemperatureDeratedStress``, MPa."""
    if temperature_c > 200.0:
        return BASE_ALLOWABLE_STRESS_MPA * 0.85
    if temperature_c < -50.0:
        return BASE_ALLOWABLE_STRESS_MPA * 0.90
    return BASE_ALLOWABLE_STRESS_MPA


def _casing_wall(diameter_mm: float, pressure_bara: float, temperature_c: float) -> float:
    """``calculateExpanderCasing``/``calculateCompressorCasing``, in millimetres."""
    inner = diameter_mm * 1.3
    pressure_mpa = pressure_bara * 0.1
    allowable = _derated_stress(temperature_c)
    joint = 0.85
    thickness = (pressure_mpa * inner) / (2.0 * allowable * joint - 0.6 * pressure_mpa) + 3.0
    return max(thickness, 10.0)


def _nozzle_throat_area(volume_flow_m3_hr: float, enthalpy_drop_kj_kg: float) -> float:
    """``estimateExpanderNozzleThroatArea``, m2."""
    if volume_flow_m3_hr <= 0.0 or enthalpy_drop_kj_kg <= 0.0:
        return 0.01
    velocity = 0.95 * math.sqrt(2.0 * enthalpy_drop_kj_kg * 1000.0)
    return (volume_flow_m3_hr / 3600.0) / max(velocity, 1.0) * 1.2


@dataclass(frozen=True, slots=True)
class _Weights:
    expander_wheel: float
    expander_casing: float
    compressor_impeller: float
    compressor_casing: float
    shaft: float
    vessel: float
    internals: float
    piping: float
    electrical: float
    structural: float
    total: float


def _weights(
    expander_wheel_diameter: float,
    compressor_stages: int,
    compressor_impeller_diameter: float,
    shaft_diameter: float,
    bearing_span: float,
    expander_casing_wall: float,
    compressor_casing_wall: float,
    expander_power_kw: float,
    compressor_power_kw: float,
    compressor_seal: str,
) -> _Weights:
    """``calculateWeights``, in kilograms."""
    wheel_volume = (
        math.pi * (expander_wheel_diameter / 2000.0) ** 2 * (expander_wheel_diameter / 3000.0)
    )
    expander_wheel = wheel_volume * STEEL_DENSITY * 0.7
    impeller_volume = (
        compressor_stages
        * math.pi
        * (compressor_impeller_diameter / 2000.0) ** 2
        * (compressor_impeller_diameter / 3000.0)
    )
    compressor_impeller = impeller_volume * STEEL_DENSITY * 0.7
    shaft_length = bearing_span * 1.5
    shaft_volume = math.pi * (shaft_diameter / 2000.0) ** 2 * (shaft_length / 1000.0)
    shaft = shaft_volume * STEEL_DENSITY

    exp_out = expander_wheel_diameter * 1.3 + 2.0 * expander_casing_wall
    exp_in = expander_wheel_diameter * 1.3
    exp_length = expander_wheel_diameter * 2.0
    exp_volume = (
        math.pi / 4.0 * ((exp_out / 1000.0) ** 2 - (exp_in / 1000.0) ** 2) * (exp_length / 1000.0)
    )
    expander_casing = exp_volume * STEEL_DENSITY * 1.5

    comp_out = compressor_impeller_diameter * 1.3 + 2.0 * compressor_casing_wall
    comp_in = compressor_impeller_diameter * 1.3
    comp_length = compressor_stages * compressor_impeller_diameter * 1.2
    comp_volume = (
        math.pi
        / 4.0
        * ((comp_out / 1000.0) ** 2 - (comp_in / 1000.0) ** 2)
        * (comp_length / 1000.0)
    )
    compressor_casing = comp_volume * STEEL_DENSITY * 1.5

    internals = (expander_wheel + compressor_impeller + shaft) * 1.3
    lube_system = 200.0 + (expander_power_kw + compressor_power_kw) * 0.15
    seal_system = 500.0 if compressor_seal == "Dry Gas" else 100.0
    baseplate = (expander_casing + compressor_casing + internals) * 0.2
    piping_weight = (expander_casing + compressor_casing) * 0.15
    electrical = 50.0
    structural = (expander_casing + compressor_casing + internals) * 0.1
    vessel = expander_casing + compressor_casing
    piping = piping_weight + lube_system + seal_system
    structural_steel = structural + baseplate
    return _Weights(
        expander_wheel=expander_wheel,
        expander_casing=expander_casing,
        compressor_impeller=compressor_impeller,
        compressor_casing=compressor_casing,
        shaft=shaft,
        vessel=vessel,
        internals=internals,
        piping=piping,
        electrical=electrical,
        structural=structural_steel,
        total=vessel + internals + piping + electrical + structural_steel,
    )


def _module_dimensions(
    expander_wheel_diameter: float,
    expander_casing_wall: float,
    compressor_stages: int,
    compressor_impeller_diameter: float,
    compressor_casing_wall: float,
) -> tuple[float, float, float]:
    """``calculateModuleDimensions``, in metres."""
    expander_len = expander_wheel_diameter / 1000.0 * 2.5
    compressor_len = compressor_stages * compressor_impeller_diameter / 1000.0 * 1.2 + 1.0
    module_length = expander_len + compressor_len + 3.0
    max_casing_od = max(
        expander_wheel_diameter * 1.3 + 2.0 * expander_casing_wall,
        compressor_impeller_diameter * 1.3 + 2.0 * compressor_casing_wall,
    )
    module_width = max(2.5, max_casing_od / 1000.0 * 2.0 + 1.5)
    module_height = max(2.0, max_casing_od / 1000.0 + 1.5)
    return module_length, module_width, module_height
