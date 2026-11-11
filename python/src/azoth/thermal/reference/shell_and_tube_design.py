"""``thermal.shell_and_tube_design`` - the detailed TEMA/ASME VIII/NACE design loop.

Spec: ``specs/calcs/thermal/shell_and_tube_design.toml``

The Python twin of ``crates/azoth-thermal/src/shell_and_tube_design.rs``, and both mirror
``ShellAndTubeDesignCalculator``. The pipeline is single-pass: TEMA configuration, material
lookup, tube bundle, shell dimensions, baffles, UHX-13, UG-37, UG-27, UG-99, the optional NACE
assessment, the optional thermal and vibration block, weights and costs.

**Two upstream defects are reproduced**: the NACE test is ``Yes`` against a ``true``/``false``
column, and the price test looks for ``STAINLESS`` in a specification string that only ever
holds ``ASME-II-A``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import AppliedStandardRow, NaceIssueRow
from azoth.core.result import ShellAndTubeDesignResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.thermal.reference.shell_and_tube_rating import _evaluate_si
from azoth.thermal.reference.shell_and_tube_vibration import screen
from azoth.thermal.reference.tema_standard import (
    CLASS_LIMITS,
    TUBE_OUTSIDE_DIAMETER,
    configuration_cost_factor,
    max_unsupported_span,
    min_baffle_spacing,
)
from azoth.thermal.reference.tube_materials import lookup

CALC_ID = "thermal.shell_and_tube_design"

TUBE_WALL_MM = 2.108
CORROSION_SHELL_MM = 3.175
CORROSION_TUBE_MM = 0.0
PITCH_RATIO = 1.25
BAFFLE_CUT = 0.25
SHELL_WALL_JOINT_EFFICIENCY = 0.85
STEEL_DENSITY = 7850.0
WATER_DENSITY = 1000.0
SOUR_THRESHOLD_BAR = 0.003
HARDNESS_HV = 248.0
CARBON_PRICE = 1.5
STAINLESS_PRICE = 5.0

SHELL_SIZES_MM = (
    168.3,
    219.1,
    273.0,
    323.9,
    355.6,
    406.4,
    457.2,
    508.0,
    558.8,
    609.6,
    660.4,
    711.2,
    762.0,
    812.8,
    863.6,
    914.4,
    965.2,
    1016.0,
    1066.8,
    1117.6,
    1168.4,
    1219.2,
    1270.0,
    1320.8,
    1371.6,
    1422.4,
    1473.2,
    1524.0,
    1574.8,
    1625.6,
    1676.4,
    1727.2,
    1778.0,
    1828.8,
    1879.6,
    1930.4,
    1981.2,
    2032.0,
)

PLATE_THICKNESSES_MM = (
    3.175,
    4.763,
    6.35,
    7.938,
    9.525,
    12.7,
    15.875,
    19.05,
    22.225,
    25.4,
    31.75,
    38.1,
    44.45,
    50.8,
)

FRONT_HEADS = frozenset("ABCDN")
SHELL_LETTERS = frozenset("EFGHJKX")
REAR_HEADS = frozenset("LMNPSTUW")


@dataclass(frozen=True, slots=True)
class Material:
    """One material's properties, as the class loads them."""

    name: str
    allowable_mpa: float
    smys_mpa: float
    conductivity: float
    grade: str


def _load_material(grade: str, is_shell: bool) -> Material:
    row = lookup(grade)
    if row is not None:
        return Material(
            name=f"{row.specification} {grade}",
            allowable_mpa=row.allowable_stress_mpa,
            smys_mpa=row.smys_mpa,
            conductivity=row.thermal_conductivity_wpm_k,
            grade=grade,
        )
    return Material(
        name="Carbon Steel SA-516-70" if is_shell else "Carbon Steel SA-179",
        allowable_mpa=138.0 if is_shell else 117.0,
        smys_mpa=260.0 if is_shell else 228.0,
        conductivity=52.0,
        grade=grade,
    )


def _parse_designation(designation: str) -> tuple[str, str, str]:
    if len(designation) != 3:
        raise InvalidInputError("tema_designation", f"{designation} is not three letters")
    front, shell, rear = designation[0], designation[1], designation[2]
    if front not in FRONT_HEADS:
        raise InvalidInputError("tema_designation", f"{front} is not a front head")
    if shell not in SHELL_LETTERS:
        raise InvalidInputError("tema_designation", f"{shell} is not a shell type")
    if rear not in REAR_HEADS:
        raise InvalidInputError("tema_designation", f"{rear} is not a rear head")
    return front, shell, rear


def _rear_clearance_mm(rear: str) -> float:
    if rear == "T":
        return 90.0
    if rear == "S":
        return 35.0
    if rear == "U":
        return 15.0
    return 10.0


def _round_shell_size(diameter: float) -> float:
    for size in SHELL_SIZES_MM:
        if size >= diameter:
            return size
    return SHELL_SIZES_MM[-1]


def _round_plate(thickness: float) -> float:
    for size in PLATE_THICKNESSES_MM:
        if size >= thickness:
            return size
    return PLATE_THICKNESSES_MM[-1]


def _allowable_stress(
    shell: Material, tube: Material, material: str, temperature_c: float
) -> float:
    if material == shell.name:
        base = shell.allowable_mpa
    elif material == tube.name:
        base = tube.allowable_mpa
    else:
        base = 137.9
    if temperature_c > 400.0:
        base *= 0.55
    elif temperature_c > 300.0:
        base *= 0.70
    elif temperature_c > 200.0:
        base *= 0.85
    return base


def _whole(value: float, field: str) -> int:
    if not math.isfinite(value) or value <= 0.0 or value != math.floor(value):
        raise InvalidInputError(field, f"{value} is not a positive whole count")
    return int(value)


def _group_present(values: list[float | None], side: str) -> bool:
    present = sum(value is not None for value in values)
    if present == 0:
        return False
    if present != len(values):
        raise InvalidInputError(side, f"the {side} thermal properties are all-or-none")
    return True


def _estimate_tubesheet_thickness(
    shell_inside_diameter: float,
    design_temperature_c: float,
    shell_pressure: float,
    tube_pressure: float,
    tema_class: str,
    shell_material: Material,
    tube_material: Material,
) -> float:
    pressure_mpa = max(shell_pressure, tube_pressure) * 0.1
    allowable = _allowable_stress(
        shell_material, tube_material, shell_material.name, design_temperature_c
    )
    thickness = (
        shell_inside_diameter * math.sqrt(0.785 * pressure_mpa / allowable) / 0.7
        + CORROSION_SHELL_MM
    )
    minimum = {"C": 15.875, "B": 12.7, "R": 19.05}[tema_class]
    return max(thickness, minimum)


def shell_and_tube_design(
    tema_designation: str,
    tema_class: str,
    required_area: Q,
    shell_side_pressure: Q,
    tube_side_pressure: Q,
    design_temperature: Q,
    shell_material_grade: str,
    tube_material_grade: str,
    tube_size: str,
    tube_length: Q,
    tube_passes: float,
    shell_joint_efficiency: float,
    nozzle_diameter: Q,
    nozzle_wall_thickness: Q,
    sour_service_assessment: bool,
    h2s_partial_pressure: Q | None = None,
    fouling_tube: float | None = None,
    fouling_shell: float | None = None,
    shell_side_method: str | None = None,
    tube_side_density: Q | None = None,
    tube_side_viscosity: Q | None = None,
    tube_side_cp: Q | None = None,
    tube_side_conductivity: Q | None = None,
    tube_side_mass_flow: Q | None = None,
    shell_side_density: Q | None = None,
    shell_side_viscosity: Q | None = None,
    shell_side_cp: Q | None = None,
    shell_side_conductivity: Q | None = None,
    shell_side_mass_flow: Q | None = None,
) -> ShellAndTubeDesignResult:
    """The detailed design.

    Raises:
        InvalidInputError: for a designation that is not three declared letters, a
            non-integral pass count, a partial thermal side, or a state the ASME formulas
            cannot take.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    raw: dict[str, Any] = {
        "required_area": input_to_si(spec, "required_area", required_area),
        "shell_side_pressure": input_to_si(spec, "shell_side_pressure", shell_side_pressure) / 1e5,
        "tube_side_pressure": input_to_si(spec, "tube_side_pressure", tube_side_pressure) / 1e5,
        "design_temperature": input_to_si(spec, "design_temperature", design_temperature),
        "tube_length": input_to_si(spec, "tube_length", tube_length) * 1000.0,
        "tube_passes": tube_passes,
        "shell_joint_efficiency": shell_joint_efficiency,
        "nozzle_diameter": input_to_si(spec, "nozzle_diameter", nozzle_diameter) * 1000.0,
        "nozzle_wall_thickness": input_to_si(spec, "nozzle_wall_thickness", nozzle_wall_thickness)
        * 1000.0,
        "h2s_partial_pressure": None
        if h2s_partial_pressure is None
        else input_to_si(spec, "h2s_partial_pressure", h2s_partial_pressure) / 1e5,
        "fouling_tube": fouling_tube,
        "fouling_shell": fouling_shell,
    }
    apply_checks(checks.on_input, raw.get, warnings)

    front, shell_letter, rear = _parse_designation(tema_designation)
    passes = _whole(tube_passes, "tube_passes")
    if raw["required_area"] <= 0.0:
        raise InvalidInputError("required_area", "the required area must be positive")
    design_temperature_c = raw["design_temperature"] - 273.15

    shell_material = _load_material(shell_material_grade, True)
    tube_material = _load_material(tube_material_grade, False)

    tube_od = TUBE_OUTSIDE_DIAMETER[tube_size] * 1000.0
    tube_id = tube_od - 2.0 * TUBE_WALL_MM
    tube_pitch = tube_od * PITCH_RATIO
    area_per_tube = math.pi * tube_od / 1000.0 * raw["tube_length"] / 1000.0
    tube_count = math.ceil(raw["required_area"] / area_per_tube)
    if tube_count % passes != 0:
        tube_count = (tube_count // passes + 1) * passes
    actual_area = tube_count * area_per_tube

    area_per_tube_layout = tube_pitch * tube_pitch * math.sqrt(3.0) / 2.0
    pass_lane_factor = 1.0 - 0.05 * (passes - 1.0)
    bundle_area = tube_count * area_per_tube_layout / pass_lane_factor
    bundle_diameter = 2.0 * math.sqrt(bundle_area / math.pi)
    shell_inside_diameter = _round_shell_size(
        bundle_diameter + 2.0 * tube_od + _rear_clearance_mm(rear)
    )
    design_pressure_mpa = max(raw["shell_side_pressure"], 1.0) * 1.1 * 0.1
    shell_allowable = _allowable_stress(
        shell_material, tube_material, shell_material.name, design_temperature_c
    )
    shell_radius = shell_inside_diameter / 2.0
    shell_wall = (
        design_pressure_mpa
        * shell_radius
        / (shell_allowable * SHELL_WALL_JOINT_EFFICIENCY - 0.6 * design_pressure_mpa)
        + CORROSION_SHELL_MM
    )
    min_wall_mm = CLASS_LIMITS[tema_class][1] * 1000.0
    shell_wall = _round_plate(max(shell_wall, min_wall_mm))
    shell_outside_diameter = shell_inside_diameter + 2.0 * shell_wall

    minimum_spacing = min_baffle_spacing(shell_inside_diameter / 1000.0, tema_class) * 1000.0
    maximum_spacing = shell_inside_diameter
    maximum_span = max_unsupported_span(tube_od / 1000.0, tube_material.name) * 1000.0
    target_spacing = shell_inside_diameter * 0.4
    baffle_spacing = max(minimum_spacing, min(target_spacing, maximum_spacing))
    baffle_spacing = min(baffle_spacing, maximum_span)
    effective_length = raw["tube_length"] - 2.0 * baffle_spacing
    baffle_count = max(1, math.floor(effective_length / baffle_spacing))

    pressure_mpa = max(raw["shell_side_pressure"], raw["tube_side_pressure"]) * 0.1
    ligament = max(1.0 - tube_od / (tube_od * PITCH_RATIO), 0.25)
    f_factor = 0.785 / ligament
    tubesheet = (
        shell_inside_diameter * math.sqrt(f_factor * pressure_mpa / shell_allowable)
        + CORROSION_SHELL_MM
    )
    tubesheet = max(tubesheet, {"C": 15.875, "B": 12.7, "R": 19.05}[tema_class])
    standards = [AppliedStandardRow(name="ASME VIII Div.1 UHX-13 (Tubesheet)")]

    joint = raw["shell_joint_efficiency"]
    t_r = (pressure_mpa * shell_inside_diameter / 2.0) / (
        shell_allowable * joint - 0.6 * pressure_mpa
    )
    required_reinforcement = raw["nozzle_diameter"] * t_r
    t_shell_net = shell_wall - CORROSION_SHELL_MM
    a1 = max(t_shell_net - t_r, 0.0) * raw["nozzle_diameter"] * 2.0
    t_rn = (pressure_mpa * raw["nozzle_diameter"] / 2.0) / (
        shell_allowable * joint - 0.6 * pressure_mpa
    )
    a2 = max(raw["nozzle_wall_thickness"] - t_rn, 0.0) * 2.0 * 2.5 * t_shell_net
    available_area = a1 + a2
    nozzle_adequate = available_area >= required_reinforcement
    standards.append(AppliedStandardRow(name="ASME VIII Div.1 UG-37 (Nozzle Reinforcement)"))

    mawp_shell_mpa = (shell_allowable * joint * t_shell_net) / (shell_radius + 0.6 * t_shell_net)
    tube_allowable = _allowable_stress(
        shell_material, tube_material, tube_material.name, design_temperature_c
    )
    t_tube = TUBE_WALL_MM - CORROSION_TUBE_MM
    tube_radius = (tube_od - 2.0 * TUBE_WALL_MM) / 2.0
    mawp_tube_mpa = (tube_allowable * t_tube) / (tube_radius + 0.6 * t_tube)
    mawp_shell = mawp_shell_mpa * 10.0
    mawp_tube = mawp_tube_mpa * 10.0
    standards.append(AppliedStandardRow(name="ASME VIII Div.1 UG-27 (MAWP)"))

    stress_ratio = 1.0
    if design_temperature_c > 100.0:
        stress_ratio = _allowable_stress(
            shell_material, tube_material, shell_material.name, 20.0
        ) / _allowable_stress(
            shell_material, tube_material, shell_material.name, design_temperature_c
        )
    hydro_shell = 1.3 * mawp_shell * stress_ratio
    hydro_tube = 1.3 * mawp_tube * stress_ratio
    standards.append(AppliedStandardRow(name="ASME VIII Div.1 UG-99 (Hydrostatic Test)"))

    nace_issues: list[NaceIssueRow] = []
    sour_required = False
    h2s_output: float | None = None
    shell_compliant = False
    tube_compliant = False
    if sour_service_assessment:
        h2s = 0.0 if raw["h2s_partial_pressure"] is None else raw["h2s_partial_pressure"]
        h2s_output = h2s
        sour_required = h2s > SOUR_THRESHOLD_BAR
        shell_row = lookup(shell_material_grade)
        tube_row = lookup(tube_material_grade)
        shell_compliant = shell_row is not None and shell_row.is_nace_compliant()
        tube_compliant = tube_row is not None and tube_row.is_nace_compliant()
        if not sour_required:
            standards.append(
                AppliedStandardRow(
                    name="NACE MR0175/ISO 15156 (Not Required - H2S pp below threshold)"
                )
            )
        else:
            if not shell_compliant:
                nace_issues.append(
                    NaceIssueRow(
                        message=(
                            f"Shell material {shell_material_grade} is not NACE MR0175 "
                            "compliant for sour service"
                        )
                    )
                )
            if not tube_compliant:
                nace_issues.append(
                    NaceIssueRow(
                        message=(
                            f"Tube material {tube_material_grade} is not NACE MR0175 "
                            "compliant for sour service"
                        )
                    )
                )
            if HARDNESS_HV > 248.0:
                nace_issues.append(
                    NaceIssueRow(
                        message=f"Shell hardness {HARDNESS_HV:.0f} HV exceeds 248 HV limit"
                    )
                )
            if HARDNESS_HV > 248.0:
                nace_issues.append(
                    NaceIssueRow(message=f"Tube hardness {HARDNESS_HV:.0f} HV exceeds 248 HV limit")
                )
            standards.append(
                AppliedStandardRow(name="NACE MR0175/ISO 15156 (Sour Service Assessment)")
            )
            if nace_issues:
                standards.append(
                    AppliedStandardRow(name="NORSOK M-001 Rev.6 (Material Hardness Limits)")
                )

    tube_values = [
        tube_side_density,
        tube_side_viscosity,
        tube_side_cp,
        tube_side_conductivity,
        tube_side_mass_flow,
    ]
    shell_values = [
        shell_side_density,
        shell_side_viscosity,
        shell_side_cp,
        shell_side_conductivity,
        shell_side_mass_flow,
    ]
    tube_present = _group_present([None if v is None else 1.0 for v in tube_values], "tube_side")
    shell_present = _group_present([None if v is None else 1.0 for v in shell_values], "shell_side")
    if tube_present != shell_present:
        raise InvalidInputError("thermal", "the thermal block needs both sides stated")

    thermal: dict[str, float] | None = None
    vibration: dict[str, float | bool] | None = None
    if tube_present and shell_present:
        fouling_t = 1.76e-4 if fouling_tube is None else fouling_tube
        fouling_s = 1.76e-4 if fouling_shell is None else fouling_shell
        method = "kern" if shell_side_method is None else shell_side_method
        rating_raw: dict[str, Any] = {
            "tube_outside_diameter": tube_od / 1000.0,
            "tube_inside_diameter": tube_id / 1000.0,
            "tube_length": raw["tube_length"] / 1000.0,
            "tube_count": float(tube_count),
            "tube_passes": float(passes),
            "tube_pitch": tube_pitch / 1000.0,
            "triangular_pitch": True,
            "shell_inside_diameter": shell_inside_diameter / 1000.0,
            "baffle_spacing": baffle_spacing / 1000.0,
            "baffle_count": float(baffle_count),
            "baffle_cut": BAFFLE_CUT,
            "tube_wall_conductivity": tube_material.conductivity,
            "tube_density": input_to_si(spec, "tube_side_density", tube_side_density),
            "tube_viscosity": input_to_si(spec, "tube_side_viscosity", tube_side_viscosity),
            "tube_cp": input_to_si(spec, "tube_side_cp", tube_side_cp),
            "tube_conductivity": input_to_si(
                spec, "tube_side_conductivity", tube_side_conductivity
            ),
            "tube_mass_flow": input_to_si(spec, "tube_side_mass_flow", tube_side_mass_flow),
            "shell_density": input_to_si(spec, "shell_side_density", shell_side_density),
            "shell_viscosity": input_to_si(spec, "shell_side_viscosity", shell_side_viscosity),
            "shell_cp": input_to_si(spec, "shell_side_cp", shell_side_cp),
            "shell_conductivity": input_to_si(
                spec, "shell_side_conductivity", shell_side_conductivity
            ),
            "shell_mass_flow": input_to_si(spec, "shell_side_mass_flow", shell_side_mass_flow),
            "fouling_tube": fouling_t,
            "fouling_shell": fouling_s,
            "shell_viscosity_wall": None,
            "sealing_pairs": None,
            "tube_to_baffle_clearance": None,
            "shell_to_baffle_clearance": None,
            "bypass_area": None,
            "tube_fouling_layer_thickness": None,
            "tube_fouling_layer_conductivity": None,
        }
        thermal = _evaluate_si(rating_raw, method, None, True)
        standards.append(AppliedStandardRow(name="Gnielinski correlation (Tube-side HTC)"))
        standards.append(
            AppliedStandardRow(
                name="Bell-Delaware method (Shell-side HTC)"
                if method == "bell_delaware"
                else "Kern method (Shell-side HTC)"
            )
        )

        area = (
            baffle_spacing
            / 1000.0
            * (shell_inside_diameter / 1000.0)
            * (1.0 - (tube_od / 1000.0) / (tube_pitch / 1000.0))
        )
        shell_density = rating_raw["shell_density"]
        velocity = rating_raw["shell_mass_flow"] / (shell_density * area) if area > 0.0 else 0.0
        sonic = math.sqrt(1.4 * 8314.0 * 300.0 / 29.0)
        damping = 0.03 if shell_density > 100.0 else 0.01
        vibration = screen(
            tube_od / 1000.0,
            tube_id / 1000.0,
            baffle_spacing / 1000.0,
            tube_pitch / 1000.0,
            200e9,
            STEEL_DENSITY,
            velocity,
            shell_density,
            rating_raw["tube_density"],
            shell_inside_diameter / 1000.0,
            sonic,
            damping,
            True,
        )
        standards.append(
            AppliedStandardRow(
                name="TEMA RCB-4.6 (Vibration Screening - PASS)"
                if vibration["passed"]
                else "TEMA RCB-4.6 (Vibration Screening - FAIL)"
            )
        )

    shell_volume = (
        math.pi
        / 4.0
        * ((shell_outside_diameter / 1000.0) ** 2 - (shell_inside_diameter / 1000.0) ** 2)
        * raw["tube_length"]
        / 1000.0
    )
    shell_weight = shell_volume * STEEL_DENSITY
    tube_volume = (
        math.pi
        / 4.0
        * ((tube_od / 1000.0) ** 2 - (tube_id / 1000.0) ** 2)
        * raw["tube_length"]
        / 1000.0
    )
    tube_weight = tube_count * tube_volume * STEEL_DENSITY
    estimated_tubesheet = _estimate_tubesheet_thickness(
        shell_inside_diameter,
        design_temperature_c,
        raw["shell_side_pressure"],
        raw["tube_side_pressure"],
        tema_class,
        shell_material,
        tube_material,
    )
    tubesheet_volume = (
        math.pi / 4.0 * (shell_outside_diameter / 1000.0) ** 2 * estimated_tubesheet / 1000.0
    )
    tubesheet_count = 1.0 if rear == "U" else 2.0
    tubesheet_weight = tubesheet_count * tubesheet_volume * STEEL_DENSITY
    head_thickness = shell_wall * 0.5
    head_volume = (
        2.0
        * math.pi
        / 3.0
        * (
            ((shell_outside_diameter / 2.0 + head_thickness) / 1000.0) ** 3
            - ((shell_outside_diameter / 2.0) / 1000.0) ** 3
        )
    )
    head_weight = 2.0 * head_volume * STEEL_DENSITY
    baffle_area = math.pi / 4.0 * (shell_inside_diameter / 1000.0) ** 2 * (1.0 - BAFFLE_CUT)
    baffle_volume = baffle_area * 6.35 / 1000.0
    baffle_weight = baffle_count * baffle_volume * STEEL_DENSITY
    total_dry = shell_weight + tube_weight + tubesheet_weight + head_weight + baffle_weight
    shell_side_volume = (
        math.pi / 4.0 * (shell_inside_diameter / 1000.0) ** 2 * raw["tube_length"] / 1000.0
        - tube_count * math.pi / 4.0 * (tube_od / 1000.0) ** 2 * raw["tube_length"] / 1000.0
    )
    tube_side_volume = (
        tube_count * math.pi / 4.0 * (tube_id / 1000.0) ** 2 * raw["tube_length"] / 1000.0
    )
    operating_weight = total_dry + (shell_side_volume + tube_side_volume) * WATER_DENSITY

    shell_price = STAINLESS_PRICE if "STAINLESS" in shell_material.name.upper() else CARBON_PRICE
    tube_price = STAINLESS_PRICE if "STAINLESS" in tube_material.name.upper() else CARBON_PRICE
    material_cost = (
        shell_weight * shell_price
        + tube_weight * tube_price
        + (tubesheet_weight + head_weight + baffle_weight) * shell_price
    )
    fabrication_factor = (
        2.5 * configuration_cost_factor(front, shell_letter, rear) * CLASS_LIMITS[tema_class][0]
    )
    fabrication_cost = material_cost * fabrication_factor

    def optional_thermal(name: str) -> float | None:
        return None if thermal is None else thermal[name]

    return ShellAndTubeDesignResult(
        tema_designation=tema_designation,
        tema_class=tema_class,
        shell_inside_diameter=from_si(shell_inside_diameter / 1000.0, "mm"),
        shell_outside_diameter=from_si(shell_outside_diameter / 1000.0, "mm"),
        shell_wall_thickness=from_si(shell_wall / 1000.0, "mm"),
        shell_material=shell_material.name,
        shell_material_grade=shell_material_grade,
        shell_allowable_stress=from_si(shell_allowable * 1e6, "MPa"),
        shell_smys=from_si(shell_material.smys_mpa * 1e6, "MPa"),
        shell_joint_efficiency=joint,
        tube_count=float(tube_count),
        tube_outer_diameter=from_si(tube_od / 1000.0, "mm"),
        tube_wall_thickness=from_si(TUBE_WALL_MM / 1000.0, "mm"),
        tube_length=from_si(raw["tube_length"] / 1000.0, "mm"),
        tube_passes=float(passes),
        tube_pitch=from_si(tube_pitch / 1000.0, "mm"),
        tube_pitch_pattern="Triangular 30°",
        tube_material=tube_material.name,
        tube_material_grade=tube_material_grade,
        tube_allowable_stress=from_si(tube_allowable * 1e6, "MPa"),
        tube_thermal_conductivity=from_si(tube_material.conductivity, "W/(m*K)"),
        baffle_type="Single Segmental",
        baffle_count=float(baffle_count),
        baffle_spacing=from_si(baffle_spacing / 1000.0, "mm"),
        baffle_cut=BAFFLE_CUT,
        required_area=from_si(raw["required_area"], "m**2"),
        actual_area=from_si(actual_area, "m**2"),
        area_margin=(actual_area - raw["required_area"]) / raw["required_area"],
        tubesheet_thickness_uhx=from_si(tubesheet / 1000.0, "mm"),
        mawp_shell_side=from_si(mawp_shell * 1e5, "bar"),
        mawp_tube_side=from_si(mawp_tube * 1e5, "bar"),
        hydro_test_pressure_shell=from_si(hydro_shell * 1e5, "bar"),
        hydro_test_pressure_tube=from_si(hydro_tube * 1e5, "bar"),
        nozzle_reinforcement_area=from_si(available_area / 1e6, "m**2"),
        nozzle_reinforcement_adequate=nozzle_adequate,
        sour_service_required=sour_required,
        h2s_partial_pressure=None if h2s_output is None else from_si(h2s_output * 1e5, "bar"),
        shell_nace_compliant=shell_compliant,
        tube_nace_compliant=tube_compliant,
        shell_max_hardness=HARDNESS_HV,
        tube_max_hardness=HARDNESS_HV,
        nace_issues=tuple(nace_issues),
        applied_standards=tuple(standards),
        tube_side_heat_transfer_coefficient=None
        if thermal is None
        else from_si(thermal["tube_htc"], "W/(m**2*K)"),
        shell_side_heat_transfer_coefficient=None
        if thermal is None
        else from_si(thermal["shell_htc"], "W/(m**2*K)"),
        overall_heat_transfer_coefficient=None
        if thermal is None
        else from_si(thermal["overall_u"], "W/(m**2*K)"),
        tube_side_pressure_drop=None if thermal is None else from_si(thermal["tube_dp"], "Pa"),
        shell_side_pressure_drop=None if thermal is None else from_si(thermal["shell_dp"], "Pa"),
        tube_side_velocity=None if thermal is None else from_si(thermal["tube_velocity"], "m/s"),
        shell_side_velocity=None if thermal is None else from_si(thermal["shell_velocity"], "m/s"),
        tube_side_reynolds=optional_thermal("tube_re"),
        shell_side_reynolds=optional_thermal("shell_re"),
        effective_tube_inside_diameter=None
        if thermal is None
        else from_si(thermal["effective_id"], "m"),
        outside_heat_transfer_area=None
        if thermal is None
        else from_si(thermal["outside_area"], "m**2"),
        tube_fouling_layer_resistance_outside=optional_thermal("deposit_resistance"),
        jc=optional_thermal("jc"),
        jl=optional_thermal("jl"),
        jb=optional_thermal("jb"),
        js=optional_thermal("js"),
        jr=optional_thermal("jr"),
        passed=None if vibration is None else bool(vibration["passed"]),
        natural_frequency=None if vibration is None else from_si(float(vibration["fn"]), "1/s"),
        vortex_shedding_frequency=None
        if vibration is None
        else from_si(float(vibration["fvs"]), "1/s"),
        vortex_shedding_ratio=None if vibration is None else float(vibration["vortex_ratio"]),
        vortex_shedding_critical=None if vibration is None else bool(vibration["vortex_critical"]),
        critical_velocity=None if vibration is None else from_si(float(vibration["vcrit"]), "m/s"),
        velocity_ratio=None if vibration is None else float(vibration["velocity_ratio"]),
        fluid_elastic_critical=None
        if vibration is None
        else bool(vibration["fluid_elastic_critical"]),
        acoustic_frequency=None if vibration is None else from_si(float(vibration["fac"]), "1/s"),
        acoustic_ratio=None if vibration is None else float(vibration["acoustic_ratio"]),
        acoustic_critical=None if vibration is None else bool(vibration["acoustic_critical"]),
        shell_weight=from_si(shell_weight, "kg"),
        tube_weight=from_si(tube_weight, "kg"),
        tubesheet_weight=from_si(tubesheet_weight, "kg"),
        head_weight=from_si(head_weight, "kg"),
        baffle_weight=from_si(baffle_weight, "kg"),
        total_dry_weight=from_si(total_dry, "kg"),
        operating_weight=from_si(operating_weight, "kg"),
        material_cost=material_cost,
        fabrication_cost=fabrication_cost,
        total_cost=material_cost + fabrication_cost,
        warnings=tuple(warnings),
    )


__all__ = ["CALC_ID", "shell_and_tube_design"]
