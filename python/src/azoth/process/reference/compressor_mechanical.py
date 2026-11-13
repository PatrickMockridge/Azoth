"""``CompressorMechanicalDesign.calcDesign`` - the API 617 sizing surface.

The twin of ``crates/azoth-process/src/rotating_design/compressor.rs``.

The stored sizing, weights and module are transcribed in the order ``calcDesign`` runs
them. ``CompressorCasingDesignCalculator``'s audit - API 617/ASME/NACE checks and their
prose - is a separate report and is not published here.

**The impeller screen is not clipped.** The diameter and the flow coefficient are computed
unclipped and the band violations are reported as issues; the head-based candidate stands
when no stage count satisfies them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.process.reference import capacity

#: The stage head cap, kJ/kg.
MAX_HEAD_PER_STAGE = 30.0
#: The work coefficient.
WORK_COEFFICIENT = 0.50
#: The flow-coefficient band.
FLOW_COEFFICIENT_MIN = 0.01
FLOW_COEFFICIENT_MAX = 0.15
#: The impeller-diameter band, mm.
MIN_IMPELLER_DIAMETER = 100.0
MAX_IMPELLER_DIAMETER = 1500.0
#: The driver margins: below 150 kW, below 750 kW, above.
DRIVER_MARGIN_SMALL = 1.25
DRIVER_MARGIN_MEDIUM = 1.15
DRIVER_MARGIN_LARGE = 1.10
#: The casing design margins.
DESIGN_PRESSURE_MARGIN = 1.10
DESIGN_TEMPERATURE_MARGIN = 30.0
#: Steel's density, kg/m3.
STEEL_DENSITY = 7850.0
#: The shaft's allowable shear stress, Pa.
SHAFT_ALLOWABLE_SHEAR_PA = 50.0e6
#: The gravity the class converts a chart head to metres with.
HEAD_TO_METRES_GRAVITY = 9.81

#: ``CompressorMechanicalDesign.CasingType``, as the spec spells it.
CASING_TYPES = ("vertically_split", "horizontally_split", "barrel")


@dataclass(frozen=True, slots=True)
class CompressorMechanicalDesign:
    """``calcDesign``'s stored result, in SI."""

    number_of_stages: int
    impeller_diameter: float
    shaft_diameter: float
    tip_speed: float
    flow_coefficient: float
    impeller_sizing_speed: float
    impeller_sizing_flow: float
    impeller_sizing_head: float
    impeller_sizing_ok: bool
    driver_power: float
    driver_margin: float
    design_pressure: float
    design_temperature: float
    max_continuous_speed: float
    trip_speed: float
    first_critical_speed: float
    casing_type: str
    head_per_stage: float
    bearing_span: float
    casing_weight: float
    bundle_weight: float
    rotor_weight: float
    outer_diameter: float
    inner_diameter: float
    wall_thickness: float
    tantan_length: float
    weight_vessel_shell: float
    weight_internals: float
    weight_nozzle: float
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
    discharge_pressure: float,
    outlet_t: float,
    power: float,
    polytropic_head: float,
    speed: float,
    max_stages_per_casing: int = 10,
) -> CompressorMechanicalDesign:
    """``calcDesign`` on the compressor's run state, in SI.

    ``polytropic_head`` is the class's ``getPolytropicFluidHead``, kJ/kg. The mass flow is
    derived from the inlet.
    """
    volume_flow_m3_hr = (
        capacity.system_volumetric_flow(
            components, inlet_t, inlet_p, inlet_z, inlet_n, from_feed=True
        )
        * 3600.0
    )
    discharge_pressure_bara = discharge_pressure / 1.0e5
    discharge_temperature_c = outlet_t - 273.15

    design_pressure_bara = discharge_pressure_bara * DESIGN_PRESSURE_MARGIN
    design_temperature_c = discharge_temperature_c + DESIGN_TEMPERATURE_MARGIN
    casing_type = _select_casing_type(design_pressure_bara)

    head_based_stages, _ = _number_of_stages(polytropic_head)
    impeller = _impeller_sizing(
        volume_flow_m3_hr,
        polytropic_head,
        speed,
        head_based_stages,
        max_stages_per_casing,
    )
    if not math.isfinite(impeller.diameter):
        raise ValueError("the impeller screen could not size the machine")

    power_kw = power / 1000.0
    shaft_diameter = _shaft_diameter(power_kw, speed)
    driver_margin = _driver_margin(power_kw)
    driver_power_kw = power_kw * driver_margin
    rotor = _rotor_dynamics(impeller.number_of_stages, impeller.diameter, shaft_diameter, speed)
    weights = _weights(
        impeller.number_of_stages,
        impeller.diameter,
        shaft_diameter,
        rotor.bearing_span,
        casing_type,
        design_pressure_bara,
        driver_power_kw,
    )
    outer_diameter = impeller.diameter * 1.5 / 1000.0
    inner_diameter = impeller.diameter / 1000.0
    module_length, module_width, module_height = _module_dimensions(
        rotor.bearing_span, outer_diameter, driver_power_kw
    )

    return CompressorMechanicalDesign(
        number_of_stages=impeller.number_of_stages,
        impeller_diameter=impeller.diameter / 1000.0,
        shaft_diameter=shaft_diameter / 1000.0,
        tip_speed=impeller.tip_speed,
        flow_coefficient=impeller.flow_coefficient,
        impeller_sizing_speed=speed,
        impeller_sizing_flow=volume_flow_m3_hr,
        impeller_sizing_head=polytropic_head,
        impeller_sizing_ok=not impeller.issues,
        driver_power=driver_power_kw * 1000.0,
        driver_margin=driver_margin,
        design_pressure=design_pressure_bara * 1.0e5,
        design_temperature=design_temperature_c + 273.15,
        max_continuous_speed=rotor.max_continuous_speed,
        trip_speed=rotor.trip_speed,
        first_critical_speed=rotor.first_critical_speed,
        casing_type=casing_type,
        head_per_stage=impeller.head_per_stage,
        bearing_span=rotor.bearing_span / 1000.0,
        casing_weight=weights.casing,
        bundle_weight=weights.bundle,
        rotor_weight=weights.rotor,
        outer_diameter=outer_diameter,
        inner_diameter=inner_diameter,
        wall_thickness=(outer_diameter - inner_diameter) / 2.0,
        tantan_length=rotor.bearing_span * 1.3 / 1000.0,
        weight_vessel_shell=weights.casing,
        weight_internals=weights.bundle,
        weight_nozzle=0.0,
        weight_piping=weights.piping,
        weight_electro_instrument=weights.electrical,
        weight_structural_steel=weights.structural,
        weight_total=weights.total,
        module_length=module_length,
        module_width=module_width,
        module_height=module_height,
    )


def head_from_metres(head_m: float) -> float:
    """The class's fluid head, kJ/kg, from its metre head."""
    return head_m * HEAD_TO_METRES_GRAVITY / 1000.0


def head_from_power(power_w: float, mass_flow_kg_s: float, efficiency: float) -> float:
    """The class's fluid head in the no-chart polytropic branch, kJ/kg."""
    return power_w / mass_flow_kg_s / 1000.0 * efficiency


def _select_casing_type(pressure_bara: float) -> str:
    """``selectCasingType``: barrel above 100 bara, horizontal above 40, vertical below."""
    if pressure_bara > 100.0:
        return "barrel"
    if pressure_bara > 40.0:
        return "horizontally_split"
    return "vertically_split"


def _number_of_stages(total_head: float) -> tuple[int, float]:
    """``calculateNumberOfStages``, in kJ/kg."""
    if not math.isfinite(total_head) or total_head <= 0.0:
        return 1, 0.0
    stages = max(1, math.ceil(total_head / MAX_HEAD_PER_STAGE))
    return stages, total_head / stages


@dataclass(frozen=True, slots=True)
class _Impeller:
    number_of_stages: int
    head_per_stage: float
    diameter: float
    tip_speed: float
    flow_coefficient: float
    issues: list[str]


def _impeller_sizing(
    volume_flow_m3_hr: float,
    total_head: float,
    speed_rpm: float,
    head_based_stages: int,
    max_stages: int,
) -> _Impeller:
    """``calculateImpellerSizing`` and its ``setImpellerCandidate``, in millimetres."""
    if (
        not math.isfinite(speed_rpm)
        or speed_rpm <= 0.0
        or not math.isfinite(volume_flow_m3_hr)
        or volume_flow_m3_hr <= 0.0
        or not math.isfinite(total_head)
        or total_head <= 0.0
    ):
        return _Impeller(0, math.nan, math.nan, math.nan, math.nan, ["not positive"])
    stages = max(1, head_based_stages)
    while 0 < stages <= max(1, max_stages):
        candidate = _impeller_candidate(volume_flow_m3_hr, total_head, speed_rpm, stages)
        if not candidate.issues:
            return candidate
        if (
            candidate.diameter < MIN_IMPELLER_DIAMETER
            or candidate.flow_coefficient > FLOW_COEFFICIENT_MAX
        ):
            break
        stages += 1
    return _impeller_candidate(volume_flow_m3_hr, total_head, speed_rpm, max(1, head_based_stages))


def _impeller_candidate(
    volume_flow_m3_hr: float, total_head: float, speed_rpm: float, stages: int
) -> _Impeller:
    """``setImpellerCandidate``: the coupled geometry and the issue list, unclipped."""
    head_per_stage = total_head / stages
    required_tip_speed = math.sqrt(head_per_stage * 1000.0 / WORK_COEFFICIENT)
    diameter = required_tip_speed * 60.0 / (math.pi * speed_rpm) * 1000.0
    tip_speed = math.pi * (diameter / 1000.0) * speed_rpm / 60.0
    flow_coefficient = (volume_flow_m3_hr / 3600.0) / ((diameter / 1000.0) ** 2 * tip_speed)
    issues = []
    if not MIN_IMPELLER_DIAMETER <= diameter <= MAX_IMPELLER_DIAMETER:
        issues.append("the impeller diameter is outside 100-1500 mm")
    if not FLOW_COEFFICIENT_MIN <= flow_coefficient <= FLOW_COEFFICIENT_MAX:
        issues.append("the inlet flow coefficient is outside 0.01-0.15")
    return _Impeller(stages, head_per_stage, diameter, tip_speed, flow_coefficient, issues)


def _shaft_diameter(power_kw: float, speed_rpm: float) -> float:
    """``calculateShaftDiameter``, in millimetres."""
    if speed_rpm <= 0.0 or power_kw <= 0.0:
        return 50.0
    torque_nm = power_kw * 1000.0 * 60.0 / (2.0 * math.pi * speed_rpm)
    diameter = (
        float((16.0 * torque_nm / (math.pi * SHAFT_ALLOWABLE_SHEAR_PA)) ** (1.0 / 3.0))
        * 1000.0
        * 1.5
    )
    return min(max(diameter, 30.0), 300.0)


def _driver_margin(shaft_power_kw: float) -> float:
    """``calculateDriverSizing``'s margin."""
    if shaft_power_kw < 150.0:
        return DRIVER_MARGIN_SMALL
    if shaft_power_kw < 750.0:
        return DRIVER_MARGIN_MEDIUM
    return DRIVER_MARGIN_LARGE


@dataclass(frozen=True, slots=True)
class _Rotor:
    max_continuous_speed: float
    trip_speed: float
    bearing_span: float
    first_critical_speed: float


def _rotor_dynamics(
    stages: int, impeller_diameter: float, shaft_diameter: float, operating_speed: float
) -> _Rotor:
    """``calculateRotorDynamics``, in rpm and millimetres."""
    max_continuous_speed = operating_speed * 1.05
    trip_speed = max_continuous_speed * 1.05
    bearing_span = stages * (impeller_diameter * 0.8) + impeller_diameter
    stiffness = 48.0 * 2.1e11 * math.pi * (shaft_diameter / 2000.0) ** 4 / 64.0
    mass_per_length = STEEL_DENSITY * math.pi * (shaft_diameter / 2000.0) ** 2
    first_critical_speed = (
        (946.0 / (2.0 * math.pi))
        * math.sqrt(stiffness / (mass_per_length * (bearing_span / 1000.0) ** 3))
        * 60.0
        if bearing_span > 0.0
        else math.nan
    )
    return _Rotor(max_continuous_speed, trip_speed, bearing_span, first_critical_speed)


@dataclass(frozen=True, slots=True)
class _Weights:
    rotor: float
    casing: float
    bundle: float
    piping: float
    electrical: float
    structural: float
    total: float


def _weights(
    stages: int,
    impeller_diameter: float,
    shaft_diameter: float,
    bearing_span: float,
    casing_type: str,
    design_pressure_bara: float,
    driver_power_kw: float,
) -> _Weights:
    """``calculateWeights``, in kilograms."""
    impeller_weight = stages * 0.5 * (impeller_diameter / 100.0) ** 2.5
    shaft_weight = bearing_span / 1000.0 * STEEL_DENSITY * math.pi * (shaft_diameter / 2000.0) ** 2
    rotor = impeller_weight + shaft_weight

    casing_thickness = max(design_pressure_bara * impeller_diameter / (2.0 * 150.0), 10.0)
    casing_length = bearing_span * 1.3
    casing_outer_diameter = impeller_diameter * 1.5
    casing = (
        math.pi
        * (casing_outer_diameter / 1000.0)
        * (casing_length / 1000.0)
        * (casing_thickness / 1000.0)
        * STEEL_DENSITY
    )
    casing *= 1.2
    if casing_type == "barrel":
        casing *= 1.3

    internals = stages * 20.0 * (impeller_diameter / 300.0) ** 2
    bundle = rotor + internals
    seal_system = 100.0 * (shaft_diameter / 100.0)
    lube_system = 200.0 + driver_power_kw * 0.1
    base = casing * 0.3
    empty_vessel = casing + base
    piping = empty_vessel * 0.2
    electrical = driver_power_kw * 0.5
    structural = empty_vessel * 0.15
    total = empty_vessel + bundle + seal_system + lube_system + piping + electrical + structural
    return _Weights(rotor, casing, bundle, piping, electrical, structural, total)


def _module_dimensions(
    bearing_span: float, outer_diameter: float, driver_power_kw: float
) -> tuple[float, float, float]:
    """``calculateModuleDimensions``, in metres."""
    compressor_length = bearing_span * 1.5 / 1000.0
    driver_length = 0.5 + driver_power_kw * 0.001
    coupling_space = 0.5
    auxiliary_space = 2.0
    module_length = max(compressor_length + driver_length + coupling_space + auxiliary_space, 4.0)
    module_width = max(outer_diameter + 3.0, 3.0)
    module_height = max(outer_diameter + 2.0, 3.0)
    return module_length, module_width, module_height
