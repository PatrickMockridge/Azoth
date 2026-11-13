"""``PumpMechanicalDesign.calcDesign`` - the API 610 sizing surface.

The twin of ``crates/azoth-process/src/rotating_design/pump.rs``.

The class's own gravity is ``9.81``, which is not the chart kernel's ``9.80665`` and not
the API 610 calculator's; that is why the design's head is a distinct quantity from the
model's. The API 610 assessment's driver selection is here because its selected rating
overwrites the design's driver power, and every weight and module dimension is built on
the intermediate that precedes it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from azoth.process.reference import capacity

#: The class's own gravity, ``PumpMechanicalDesign.GRAVITY``.
GRAVITY = 9.81
#: The API 610 calculator's gravity, which the assessment's pressure basis uses.
API610_GRAVITY = 9.80665
#: The driver margins, from absorbed power: below 22 kW, at or below 55 kW, above.
DRIVER_MARGIN_SMALL = 1.25
DRIVER_MARGIN_MEDIUM = 1.15
DRIVER_MARGIN_LARGE = 1.10
#: The design temperature's margin, degC.
DESIGN_TEMPERATURE_MARGIN = 30.0
#: The class's default pump efficiency.
DEFAULT_PUMP_EFFICIENCY = 0.75
#: Steel's density, kg/m3.
STEEL_DENSITY = 7850.0
#: The casing's allowable stress, MPa.
ALLOWABLE_STRESS_MPA = 137.9
#: The shaft's allowable shear stress, Pa.
SHAFT_ALLOWABLE_SHEAR_PA = 40.0e6

#: ``PumpMechanicalDesign.SealType``, as the spec spells it.
SEAL_TYPES = (
    "packed",
    "single_mechanical",
    "double_mechanical",
    "magnetic_drive",
    "canned_motor",
)

#: ``PumpApi610DesignCalculator.Api610PumpType``, as the spec spells it.
API610_PUMP_TYPES = (
    "unspecified",
    "oh1",
    "oh2",
    "oh3",
    "oh4",
    "oh5",
    "oh6",
    "bb1",
    "bb2",
    "bb3",
    "bb4",
    "bb5",
    "vs1",
    "vs2",
    "vs3",
    "vs4",
    "vs5",
    "vs6",
    "vs7",
)


@dataclass(frozen=True, slots=True)
class PumpMechanicalDesign:
    """``calcDesign``'s stored result, in SI."""

    pump_type: str
    api610_pump_type: str
    seal_type: str
    rated_flow: float
    rated_head: float
    hydraulic_power: float
    pump_efficiency: float
    rated_speed: float
    specific_speed: float
    suction_specific_speed: float
    number_of_stages: int
    bep_flow: float | None
    bep_head: float | None
    bep_source: str
    npshr_source: str
    shutoff_source: str
    shutoff_head: float | None
    design_pressure: float
    design_temperature: float
    impeller_diameter: float
    impeller_width: float
    shaft_diameter: float
    casing_wall_thickness: float
    max_allowable_working_pressure: float
    driver_power: float
    driver_margin: float
    suction_nozzle_size: float
    discharge_nozzle_size: float
    minimum_flow: float
    maximum_flow: float
    casing_weight: float
    motor_weight: float
    baseplate_weight: float
    inner_diameter: float
    outer_diameter: float
    tantan_length: float
    wall_thickness: float
    weight_piping: float
    weight_electro_instrument: float
    weight_structural_steel: float
    weight_total: float
    module_length: float
    module_width: float
    module_height: float
    npsh_margin: float
    npsh_margin_ratio: float
    assessment_status: str
    assessment_basis: str
    operating_flow_ratio: float
    operating_region: str
    governing_shutoff_head: float
    required_casing_pressure: float
    preliminary_hydrostatic_test_pressure: float
    head_rise_to_shutoff: float
    required_driver_power: float
    pump_type_source: str


def design(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_t: float,
    inlet_p: float,
    discharge_pressure: float,
    outlet_t: float,
    power: float,
    efficiency: float,
    npsh_required: float,
    npsh_available: float,
    speed: float,
    chart: Any | None,
    api610_pump_type: str | None,
    seal_type: str | None,
    maximum_suction_pressure: float | None,
    furnished_casing_mawp: float | None,
    shutoff_head: float | None,
    npsh_margin_factor: float = 1.20,
    minimum_npsh_margin: float = 1.0,
    por_low: float = 0.70,
    por_high: float = 1.20,
    aor_low: float = 0.65,
    aor_high: float = 1.25,
    rated_low: float = 0.80,
    rated_high: float = 1.10,
) -> PumpMechanicalDesign:
    """``calcDesign`` on the pump's run state, in SI.

    ``discharge_pressure`` and ``outlet_t`` are the run's outlet, ``power`` its shaft
    power, and the NPSH pair the model's own values.

    Example:
        >>> design(
        ...     ["water"], 2775.9, [1.0], 293.15, 1e5, 60e5, 293.5, 132.1, 0.75,
        ...     8.0, 10.1, 3000.0, None, None, None, None, None,
        ... ).pump_type
        'between_bearings'
    """
    density = capacity.system_mass_density(components, inlet_t, inlet_p, inlet_z)
    volume_flow_m3_hr = (
        capacity.system_volumetric_flow(components, inlet_t, inlet_p, inlet_z, inlet_n) * 3600.0
    )
    suction_pressure_bara = inlet_p / 1.0e5
    discharge_pressure_bara = discharge_pressure / 1.0e5
    suction_temperature_c = inlet_t - 273.15
    discharge_temperature_c = outlet_t - 273.15

    differential_pressure = discharge_pressure_bara - suction_pressure_bara
    pump_head = (
        differential_pressure * 1.0e5 / (density * GRAVITY) if density > 0.0 else math.nan
    )
    hydraulic_power_kw = (
        density * GRAVITY * volume_flow_m3_hr / 3600.0 * pump_head / 1000.0
        if density > 0.0 and volume_flow_m3_hr > 0.0 and pump_head > 0.0
        else math.nan
    )

    pump_efficiency = _normalize_efficiency(efficiency)
    if not (0.0 < pump_efficiency <= 1.0):
        pump_efficiency = DEFAULT_PUMP_EFFICIENCY
    absorbed_power_kw = power / 1000.0
    if not absorbed_power_kw > 0.0 and hydraulic_power_kw > 0.0 and pump_efficiency > 0.0:
        absorbed_power_kw = hydraulic_power_kw / pump_efficiency

    bep_flow, bep_head, bep_source = _vendor_bep(chart, speed)
    shutoff, shutoff_source = _vendor_shutoff(chart, speed, shutoff_head)
    npshr_source = (
        "vendor_curve"
        if chart is not None and chart.has_npsh_curve()
        else "screening_estimate"
        if npsh_required > 0.0
        else "not_available"
    )

    specific_flow = bep_flow if bep_flow is not None and bep_flow > 0.0 else volume_flow_m3_hr
    specific_head = bep_head if bep_head is not None and bep_head > 0.0 else pump_head
    specific_speed, suction_specific_speed = _specific_speeds(
        specific_flow, specific_head, speed, npsh_required, api610_pump_type
    )

    pump_type, api610, number_of_stages = _select_pump_type(api610_pump_type, pump_head)

    maximum_suction_bara = (
        maximum_suction_pressure / 1.0e5
        if maximum_suction_pressure is not None
        else suction_pressure_bara
    )
    governing_shutoff = shutoff if shutoff is not None and shutoff > 0.0 else pump_head * 1.10
    design_pressure_bara = max(
        maximum_suction_bara + density * GRAVITY * governing_shutoff / 1.0e5,
        discharge_pressure_bara,
    )
    design_temperature_c = (
        max(suction_temperature_c, discharge_temperature_c) + DESIGN_TEMPERATURE_MARGIN
    )

    impeller_diameter, impeller_width = _impeller_sizing(
        volume_flow_m3_hr, pump_head, speed, number_of_stages
    )
    shaft_diameter = _shaft_sizing(absorbed_power_kw, speed)
    casing = _casing_design(impeller_diameter, design_pressure_bara, design_temperature_c)
    driver_margin = _driver_margin(absorbed_power_kw)
    suction_nozzle, discharge_nozzle = _nozzle_sizes(volume_flow_m3_hr)
    intermediate_driver = _calculate_driver_sizing(absorbed_power_kw)
    weights = _weights(
        impeller_diameter,
        impeller_width,
        shaft_diameter,
        casing.wall_thickness,
        pump_type,
        number_of_stages,
        intermediate_driver,
    )
    module_length, module_width, module_height = _module_dimensions(
        impeller_diameter, intermediate_driver
    )
    driver_power_kw = _api610_driver_rating(absorbed_power_kw)
    npsh_margin = npsh_available - npsh_required
    npsh_margin_ratio = npsh_available / npsh_required if npsh_required > 0.0 else math.nan
    basis = bep_flow if bep_flow is not None and bep_flow > 0.0 else volume_flow_m3_hr

    # **The per-check list is a presentation surface and is not published.** Its statuses
    # aggregate into ``assessment_status``; the scalars below are its numeric consequences.
    bep_flow_value = bep_flow if bep_flow is not None else math.nan
    operating_flow_ratio = (
        volume_flow_m3_hr / bep_flow_value
        if volume_flow_m3_hr > 0.0 and bep_flow_value > 0.0
        else math.nan
    )
    if math.isnan(operating_flow_ratio):
        operating_region = "NOT_EVALUATED"
        operating_status = _NOT_EVALUATED
    elif rated_low <= operating_flow_ratio <= rated_high:
        operating_region = "POR"
        operating_status = _PASS
    elif por_low <= operating_flow_ratio <= por_high:
        operating_region = "POR"
        operating_status = _WARNING
    elif aor_low <= operating_flow_ratio <= aor_high:
        operating_region = "AOR"
        operating_status = _WARNING
    else:
        operating_region = "OUTSIDE_AOR"
        operating_status = _FAIL

    if math.isfinite(npsh_available) and npsh_required > 0.0:
        npsh_status = (
            _PASS
            if npsh_margin_ratio >= npsh_margin_factor and npsh_margin >= minimum_npsh_margin
            else _FAIL
        )
    else:
        npsh_status = _NOT_EVALUATED

    required_casing_pressure_bara = (
        maximum_suction_bara + density * API610_GRAVITY * governing_shutoff / 1.0e5
    )
    if maximum_suction_bara > 0.0 and pump_head > 0.0 and density > 0.0:
        if furnished_casing_mawp is None:
            pressure_status = _NOT_EVALUATED
        elif furnished_casing_mawp / 1.0e5 >= required_casing_pressure_bara:
            pressure_status = _PASS
        else:
            pressure_status = _FAIL
    else:
        pressure_status = _NOT_EVALUATED
    hydrostatic_basis_bara = (
        furnished_casing_mawp / 1.0e5
        if furnished_casing_mawp is not None and furnished_casing_mawp > 0.0
        else required_casing_pressure_bara
    )
    preliminary_hydrostatic_test_pressure_bara = hydrostatic_basis_bara * 1.50

    if shutoff is not None and shutoff > 0.0 and pump_head > 0.0 and shutoff_source in (
        "vendor_curve",
        "purchaser_input",
    ):
        head_rise_to_shutoff = shutoff / pump_head - 1.0
        head_rise_status = _PASS if head_rise_to_shutoff >= 0.10 else _FAIL
    else:
        head_rise_to_shutoff = math.nan
        head_rise_status = _NOT_EVALUATED

    required_driver_power_kw = absorbed_power_kw * driver_margin
    if absorbed_power_kw > 0.0:
        driver_status = _PASS if required_driver_power_kw <= 10000.0 else _FAIL
    else:
        driver_status = _NOT_EVALUATED

    if api610_pump_type is None:
        pump_type_status = _WARNING
    elif api610_pump_type == "unspecified":
        pump_type_status = _NOT_EVALUATED
    else:
        pump_type_status = _PASS
    statuses = [
        operating_status,
        npsh_status,
        pressure_status,
        _NOT_EVALUATED,
        head_rise_status,
        driver_status,
        _NOT_EVALUATED,
        _NOT_EVALUATED,
        _NOT_EVALUATED,
        _NOT_EVALUATED,
        _NOT_EVALUATED,
        pump_type_status,
    ]
    if seal_type in ("magnetic_drive", "canned_motor"):
        statuses.append(_WARNING)
    status = _assessment_status(statuses)
    pump_type_source = "purchaser_input" if api610_pump_type is not None else "screening_estimate"

    return PumpMechanicalDesign(
        pump_type=pump_type,
        api610_pump_type=api610,
        seal_type="single_mechanical" if seal_type is None else seal_type,
        rated_flow=volume_flow_m3_hr,
        rated_head=pump_head,
        hydraulic_power=hydraulic_power_kw * 1000.0,
        pump_efficiency=pump_efficiency,
        rated_speed=speed,
        specific_speed=specific_speed,
        suction_specific_speed=suction_specific_speed,
        number_of_stages=number_of_stages,
        bep_flow=bep_flow,
        bep_head=bep_head,
        bep_source=bep_source,
        npshr_source=npshr_source,
        shutoff_source=shutoff_source,
        shutoff_head=shutoff,
        design_pressure=design_pressure_bara * 1.0e5,
        design_temperature=design_temperature_c + 273.15,
        impeller_diameter=impeller_diameter / 1000.0,
        impeller_width=impeller_width / 1000.0,
        shaft_diameter=shaft_diameter / 1000.0,
        casing_wall_thickness=casing.wall_thickness / 1000.0,
        max_allowable_working_pressure=casing.mawp_bara * 1.0e5,
        driver_power=driver_power_kw * 1000.0,
        driver_margin=driver_margin,
        suction_nozzle_size=suction_nozzle * 0.0254,
        discharge_nozzle_size=discharge_nozzle * 0.0254,
        minimum_flow=basis * 0.65,
        maximum_flow=basis * 1.15,
        casing_weight=weights.casing,
        motor_weight=weights.motor,
        baseplate_weight=weights.baseplate,
        inner_diameter=casing.inner_diameter / 1000.0,
        outer_diameter=casing.outer_diameter / 1000.0,
        tantan_length=impeller_width * number_of_stages / 1000.0,
        wall_thickness=casing.wall_thickness / 1000.0,
        weight_piping=weights.piping,
        weight_electro_instrument=weights.electro_instrument,
        weight_structural_steel=weights.structural_steel,
        weight_total=weights.total,
        module_length=module_length,
        module_width=module_width,
        module_height=module_height,
        npsh_margin=npsh_margin,
        npsh_margin_ratio=npsh_margin_ratio,
        assessment_status=status,
        assessment_basis="API 610 screening; vendor and purchaser verification required",
        operating_flow_ratio=operating_flow_ratio,
        operating_region=operating_region,
        governing_shutoff_head=governing_shutoff,
        required_casing_pressure=required_casing_pressure_bara * 1.0e5,
        preliminary_hydrostatic_test_pressure=preliminary_hydrostatic_test_pressure_bara * 1.0e5,
        head_rise_to_shutoff=head_rise_to_shutoff,
        required_driver_power=required_driver_power_kw * 1000.0,
        pump_type_source=pump_type_source,
    )


def _normalize_efficiency(efficiency: float) -> float:
    """The class's ``normalizeEfficiency``: a percentage becomes a fraction."""
    if 1.0 < efficiency <= 100.0:
        return efficiency / 100.0
    return efficiency


def _vendor_bep(chart: Any | None, speed: float) -> tuple[float | None, float | None, str]:
    """``populateVendorCurveData``'s BEP pair and provenance."""
    if chart is None or not speed > 0.0:
        return None, None, "not_available"
    flow = chart.best_efficiency_flow_rate(speed)
    if not flow > 0.0:
        return None, None, "not_available"
    head = chart.head(flow, speed)
    if not head > 0.0:
        return None, None, "not_available"
    return flow, head, "vendor_curve"


def _vendor_shutoff(
    chart: Any | None, speed: float, purchaser: float | None
) -> tuple[float | None, str]:
    """``populateVendorCurveData``'s shutoff head and provenance."""
    if purchaser is not None:
        return purchaser, "purchaser_input"
    if chart is None or not speed > 0.0:
        return None, "not_available"
    head = chart.head(0.0, speed)
    return (head, "vendor_curve") if head > 0.0 else (None, "not_available")


def _specific_speeds(
    flow_m3_hr: float,
    head_m: float,
    speed_rpm: float,
    npsh_required_m: float,
    api610: str | None,
) -> tuple[float, float]:
    """``calculateSpecificSpeeds``, in the class's two conventions."""
    flow_m3_s = flow_m3_hr / 3600.0
    specific = (
        speed_rpm * math.sqrt(flow_m3_s) / head_m**0.75
        if flow_m3_s > 0.0 and head_m > 0.0 and speed_rpm > 0.0
        else math.nan
    )
    if npsh_required_m > 0.0:
        flow_gpm = flow_m3_hr * 4.402867
        if api610 == "bb1":
            flow_gpm /= 2.0
        npsh_feet = npsh_required_m * 3.28084
        suction = speed_rpm * math.sqrt(flow_gpm) / npsh_feet**0.75
    else:
        suction = math.nan
    return specific, suction


def _select_pump_type(api610: str | None, head_m: float) -> tuple[str, str, int]:
    """``selectPumpType``: the purchaser's API 610 code, or the head-driven choice."""
    if api610 is not None:
        if api610.startswith("bb"):
            stages = 2 if api610 in ("bb3", "bb4", "bb5") else 1
            return "between_bearings", api610, stages
        if api610.startswith("vs"):
            return "vertically_suspended", api610, 1
        return "overhung", api610, 1
    if head_m > 250.0:
        stages = max(2, math.ceil(head_m / 200.0))
        return "between_bearings", "bb3" if stages > 2 else "bb2", stages
    return "overhung", "oh2", 1


def _impeller_sizing(
    flow_m3_hr: float, head_m: float, speed_rpm: float, stages: int
) -> tuple[float, float]:
    """``calculateImpellerSizing``, in millimetres."""
    if speed_rpm <= 0.0 or flow_m3_hr <= 0.0 or head_m <= 0.0:
        return 200.0, 30.0
    head_per_stage = head_m / max(1, stages)
    tip_speed = math.sqrt(2.0 * GRAVITY * head_per_stage)
    diameter = min(max(60.0 * tip_speed / (math.pi * speed_rpm) * 1000.0, 50.0), 1000.0)
    flow_m3_s = flow_m3_hr / 3600.0
    d2 = diameter / 1000.0
    cm2 = 0.1 * tip_speed
    width = (
        min(max(flow_m3_s / (math.pi * d2 * cm2) * 1000.0, 10.0), 200.0)
        if cm2 > 0.0 and d2 > 0.0
        else 10.0
    )
    return diameter, width


def _shaft_sizing(power_kw: float, speed_rpm: float) -> float:
    """``calculateShaftSizing``, in millimetres."""
    if speed_rpm <= 0.0:
        return 40.0
    torque_nm = power_kw * 60.0 / (2.0 * math.pi * speed_rpm) * 1000.0
    diameter_m = (16.0 * torque_nm / (math.pi * SHAFT_ALLOWABLE_SHEAR_PA)) ** (1.0 / 3.0)
    diameter_mm = diameter_m * 1000.0 * 1.2
    return max(math.ceil(diameter_mm / 5.0) * 5.0, 20.0)


@dataclass(frozen=True, slots=True)
class _Casing:
    inner_diameter: float
    outer_diameter: float
    wall_thickness: float
    mawp_bara: float


def _casing_design(
    impeller_diameter: float, design_pressure_bara: float, design_temperature_c: float
) -> _Casing:
    """``calculateCasingDesign``, in millimetres and bar."""
    inner_diameter = impeller_diameter * 1.15
    pressure_mpa = design_pressure_bara * 0.1
    allowable_stress = ALLOWABLE_STRESS_MPA
    if design_temperature_c > 200.0:
        allowable_stress *= 0.85
    joint_efficiency = 1.0
    wall_thickness = max(
        pressure_mpa * inner_diameter / (2.0 * allowable_stress * joint_efficiency - 0.6 * pressure_mpa)
        + 3.0,
        6.0,
    )
    mawp_bara = (
        2.0
        * allowable_stress
        * joint_efficiency
        * (wall_thickness - 3.0)
        / (inner_diameter + 0.6 * (wall_thickness - 3.0))
        * 10.0
    )
    return _Casing(
        inner_diameter,
        inner_diameter + 2.0 * wall_thickness,
        wall_thickness,
        mawp_bara,
    )


def _driver_margin(absorbed_power_kw: float) -> float:
    """The driver margin's three bands."""
    if absorbed_power_kw < 22.0:
        return DRIVER_MARGIN_SMALL
    if absorbed_power_kw <= 55.0:
        return DRIVER_MARGIN_MEDIUM
    return DRIVER_MARGIN_LARGE


#: ``calculateDriverSizing``'s own 25-rating IEC list.
IEC_RATINGS = (
    0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0, 45.0, 55.0,
    75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 400.0,
)

#: The API 610 assessment's longer candidate list.
API610_RATINGS = (
    0.75, 1.1, 1.5, 2.2, 3.0, 4.0, 5.5, 7.5, 11.0, 15.0, 18.5, 22.0, 30.0, 37.0, 45.0, 55.0,
    75.0, 90.0, 110.0, 132.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0,
)


def _select_driver_rating(required_kw: float, candidates: tuple[float, ...]) -> float | None:
    """``selectDriverRating``: the smallest candidate covering the required power."""
    if not (math.isfinite(required_kw) and required_kw > 0.0):
        return None
    covering = [candidate for candidate in candidates if candidate >= required_kw]
    return min(covering) if covering else None


def _calculate_driver_sizing(absorbed_power_kw: float) -> float:
    """``calculateDriverSizing``'s intermediate, which the weights and module read."""
    required = absorbed_power_kw * _driver_margin(absorbed_power_kw)
    selected = _select_driver_rating(required, IEC_RATINGS)
    return required if selected is None else selected


def _api610_driver_rating(absorbed_power_kw: float) -> float:
    """The API 610 assessment's selected rating, from its longer list."""
    required = absorbed_power_kw * _driver_margin(absorbed_power_kw)
    selected = _select_driver_rating(required, API610_RATINGS)
    return required if selected is None else selected


def _nozzle_sizes(flow_m3_hr: float) -> tuple[float, float]:
    """``calculateNozzleSizes``, in the class's inches."""
    standard = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0)
    flow_m3_s = flow_m3_hr / 3600.0
    suction_diameter_mm = 2.0 * math.sqrt(flow_m3_s / 2.5 / math.pi) * 1000.0
    discharge_diameter_mm = 2.0 * math.sqrt(flow_m3_s / 5.0 / math.pi) * 1000.0
    suction = 4.0
    for size in standard:
        if size * 25.4 >= suction_diameter_mm:
            suction = size
            break
    discharge = 3.0
    for size in standard:
        if size * 25.4 >= discharge_diameter_mm:
            discharge = size
            break
    if discharge >= suction and suction in standard:
        index = standard.index(suction)
        if index > 0:
            discharge = standard[index - 1]
    return suction, discharge


@dataclass(frozen=True, slots=True)
class _Weights:
    casing: float
    motor: float
    baseplate: float
    piping: float
    electro_instrument: float
    structural_steel: float
    total: float


def _weights(
    impeller_diameter: float,
    impeller_width: float,
    shaft_diameter: float,
    casing_wall_thickness: float,
    pump_type: str,
    number_of_stages: int,
    driver_power_kw: float,
) -> _Weights:
    """``calculateWeights``, in kilograms."""
    impeller_volume = (
        math.pi * (impeller_diameter / 2000.0) ** 2 * (impeller_width / 1000.0) * 0.3
    )
    impeller_weight = impeller_volume * STEEL_DENSITY

    shaft_length = 500.0 if pump_type == "overhung" else 1000.0
    shaft_volume = (
        math.pi * (shaft_diameter / 2000.0) ** 2 * (shaft_length / 1000.0)
    )
    shaft_weight = shaft_volume * STEEL_DENSITY

    casing_od = (impeller_diameter * 1.15 + 2.0 * casing_wall_thickness) / 1000.0
    casing_id = impeller_diameter * 1.15 / 1000.0
    casing_length = (impeller_width + 100.0) / 1000.0 * number_of_stages
    casing_volume = math.pi / 4.0 * (casing_od**2 - casing_id**2) * casing_length
    casing_weight = casing_volume * STEEL_DENSITY * 2.0

    if driver_power_kw <= 0.0:
        motor_weight = 50.0
    elif driver_power_kw <= 10.0:
        motor_weight = 20.0 + driver_power_kw * 3.0
    elif driver_power_kw <= 100.0:
        motor_weight = 50.0 + driver_power_kw * 2.5
    else:
        motor_weight = 300.0 + driver_power_kw * 1.5
    baseplate_weight = (casing_weight + motor_weight) * 0.6
    seal_weight = 5.0 + shaft_diameter * 0.5
    coupling_weight = 10.0 + driver_power_kw * 0.1

    empty_weight = (
        casing_weight
        + impeller_weight
        + shaft_weight
        + motor_weight
        + baseplate_weight
        + seal_weight
        + coupling_weight
    )
    piping = empty_weight * 0.15
    electro_instrument = empty_weight * 0.10
    structural_steel = empty_weight * 0.05
    return _Weights(
        casing=casing_weight,
        motor=motor_weight,
        baseplate=baseplate_weight,
        piping=piping,
        electro_instrument=electro_instrument,
        structural_steel=structural_steel,
        total=empty_weight + piping + electro_instrument + structural_steel,
    )


def _module_dimensions(impeller_diameter: float, driver_power_kw: float) -> tuple[float, float, float]:
    """``calculateModuleDimensions``, in metres."""
    module_width = max(1.0, impeller_diameter / 500.0 + 0.5)
    pump_length = impeller_diameter / 1000.0 * 2.0
    motor_length = 0.3 + math.sqrt(driver_power_kw) * 0.05
    module_length = pump_length + motor_length + 1.0
    module_height = max(0.5, impeller_diameter / 1000.0 + 0.3)
    return module_length, module_width, module_height


_PASS = "pass"
_WARNING = "warning"
_FAIL = "fail"
_NOT_EVALUATED = "not_evaluated"


def _assessment_status(statuses: list[str]) -> str:
    """``updateAssessmentStatus``: a failure anywhere defeats a warning or a gap."""
    evaluated = False
    warning_or_gap = False
    for status in statuses:
        if status == _FAIL:
            return "FAIL"
        if status == _PASS:
            evaluated = True
        else:
            warning_or_gap = True
    if not evaluated:
        return "NOT_EVALUATED"
    return "PASS_WITH_WARNINGS" if warning_or_gap else "PASS"
