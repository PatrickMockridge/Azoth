"""``thermal.tema_standard`` - TEMA configuration, class limits and layout rules.

Spec: ``specs/calcs/thermal/tema_standard.toml``

The Python twin of ``crates/azoth-thermal/src/tema_standard.rs``, and both mirror
``TEMAStandard``: the configuration cost and head predicates, the construction-class limits,
the tube-size, pitch and baffle tables, and the six numeric rules the design calculator reads.

**Every length is SI metres**; the class's millimetre constants are converted here. The tube
count is an int truncation, and a material name is matched by upper-cased substring.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TemaStandardResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.tema_standard"

SHELL_PRESSURE_DROP = {
    "E": 1.0,
    "F": 0.8,
    "G": 0.6,
    "H": 0.5,
    "J": 0.65,
    "K": 0.7,
    "X": 0.3,
}

CLASS_LIMITS = {
    "R": (1.0, 0.0127, 0.00476),
    "C": (0.8, 0.0095, 0.00318),
    "B": (0.6, 0.00635, 0.00211),
}

TUBE_OUTSIDE_DIAMETER = {
    "tube_3_8_inch": 0.009525,
    "tube_1_2_inch": 0.0127,
    "tube_5_8_inch": 0.015875,
    "tube_3_4_inch": 0.01905,
    "tube_1_inch": 0.0254,
}

PITCH_LIMITS = {
    "triangular_30": (30.0, 1.25, 1.1),
    "triangular_60": (60.0, 1.25, 1.15),
    "square_90": (90.0, 1.25, 1.0),
    "square_45": (45.0, 1.25, 1.05),
}

BAFFLE_FACTORS = {
    "single_segmental": (0.75, 1.0),
    "double_segmental": (0.5, 0.6),
    "triple_segmental": (0.35, 0.4),
    "no_tubes_in_window": (0.6, 0.5),
    "disc_and_doughnut": (0.5, 0.55),
    "rod_baffles": (0.2, 0.3),
}

FLOATING_REAR_HEADS = ("P", "S", "T", "W")


def configuration_cost_factor(front: str, shell: str, rear: str) -> float:
    """The class's step additions from one, in its order."""
    factor = 1.0
    if front == "D":
        factor += 0.4
    elif front == "A":
        factor += 0.15
    elif front in ("C", "N"):
        factor += 0.1
    if shell == "K":
        factor += 0.3
    elif shell == "F":
        factor += 0.2
    if rear == "T":
        factor += 0.25
    elif rear == "S":
        factor += 0.2
    elif rear in ("P", "W"):
        factor += 0.15
    return factor


def rear_is_removable(rear: str) -> bool:
    """The class's removable predicate: a floating rear head, or the U-tube."""
    return rear in FLOATING_REAR_HEADS or rear == "U"


def recommend_configuration(
    needs_cleanable: bool,
    large_temperature_difference: bool,
    high_pressure: bool,
    hazardous: bool,
) -> str:
    """The class's four-branch recommendation."""
    if hazardous:
        return "AEW" if large_temperature_difference else "BEM"
    if high_pressure:
        return "NEU" if large_temperature_difference else "NEN"
    if large_temperature_difference:
        return "AES" if needs_cleanable else "AEU"
    return "AES" if needs_cleanable else "BEM"


def min_tube_pitch(diameter: float, pattern: str) -> float:
    """``OD * ratio``."""
    return diameter * PITCH_LIMITS[pattern][1]


def estimate_tube_count(
    shell_inside_diameter: float,
    tube_diameter: float,
    tube_pitch: float,
    pattern: str,
    tube_passes: int,
) -> float:
    """The truncated bundle-area estimate."""
    effective_diameter = shell_inside_diameter - 2.0 * tube_diameter
    pass_lane_factor = 1.0 - 0.05 * (tube_passes - 1)
    effective_area = math.pi * math.pow(effective_diameter / 2.0, 2) * pass_lane_factor
    if pattern in ("triangular_30", "triangular_60"):
        area_per_tube = tube_pitch * tube_pitch * math.sqrt(3.0) / 2.0
    else:
        area_per_tube = tube_pitch * tube_pitch
    return float(int(effective_area / area_per_tube))


def material_factor(material: str | None) -> float:
    """The upper-cased substring match on the material name."""
    if material is None:
        return 1.0
    upper = material.upper()
    if "STAINLESS" in upper or "SS" in upper:
        return 0.95
    if "COPPER" in upper or "CU" in upper:
        return 0.75
    if "TITANIUM" in upper or "TI" in upper:
        return 0.9
    return 1.0


def max_unsupported_span(tube_diameter: float, material: str | None) -> float:
    """The class's diameter band times the material factor."""
    if tube_diameter <= 0.0127:
        base = 0.660
    elif tube_diameter <= 0.01905:
        base = 0.889
    elif tube_diameter <= 0.0254:
        base = 1.067
    else:
        base = 1.321
    return base * material_factor(material)


def min_baffle_spacing(shell_inside_diameter: float, tema_class: str) -> float:
    """The greater of a fifth of the shell and the class's floor."""
    floor = 0.0381 if tema_class == "B" else 0.0508
    return max(shell_inside_diameter / 5.0, floor)


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
    """The composed TEMA calculation.

    Raises:
        InvalidInputError: for a non-integral pass count or an undeclared enum value.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    values = {
        "shell_inside_diameter": input_to_si(spec, "shell_inside_diameter", shell_inside_diameter),
        "tube_pitch": input_to_si(spec, "tube_pitch", tube_pitch),
        "tube_passes": tube_passes,
    }
    apply_checks(checks.on_input, values.get, warnings)

    if not math.isfinite(tube_passes) or tube_passes != math.floor(tube_passes) or tube_passes < 1:
        raise InvalidInputError("tube_passes", f"{tube_passes} is not a whole pass count")
    passes = int(tube_passes)

    try:
        class_cost, min_wall, min_corrosion = CLASS_LIMITS[tema_class]
        diameter = TUBE_OUTSIDE_DIAMETER[tube_size]
        angle, ratio, pitch_heat = PITCH_LIMITS[pitch_pattern]
        baffle_heat, baffle_drop = BAFFLE_FACTORS[baffle_type]
        if front_head not in ("A", "B", "C", "N", "D"):
            raise KeyError(front_head)
        if shell not in SHELL_PRESSURE_DROP:
            raise KeyError(shell)
        if rear_head not in ("L", "M", "N", "P", "S", "T", "U", "W"):
            raise KeyError(rear_head)
    except KeyError as exc:
        raise InvalidInputError("enum", f"{exc.args[0]} is not a declared value") from None

    removable = rear_is_removable(rear_head)
    return TemaStandardResult(
        configuration_cost_factor=configuration_cost_factor(front_head, shell, rear_head),
        shell_pressure_drop_factor=SHELL_PRESSURE_DROP[shell],
        bundle_removable=removable,
        good_thermal_expansion=removable,
        class_cost_factor=class_cost,
        min_tube_wall_thickness=from_si(min_wall, "m"),
        min_corrosion_allowance=from_si(min_corrosion, "m"),
        tube_outer_diameter=from_si(diameter, "m"),
        tube_outer_diameter_inch=from_si(diameter, "in"),
        layout_angle=angle,
        min_pitch_ratio=ratio,
        pitch_heat_transfer_factor=pitch_heat,
        baffle_heat_transfer_factor=baffle_heat,
        baffle_pressure_drop_factor=baffle_drop,
        min_tube_pitch=from_si(min_tube_pitch(diameter, pitch_pattern), "m"),
        estimated_tube_count=estimate_tube_count(
            values["shell_inside_diameter"], diameter, values["tube_pitch"], pitch_pattern, passes
        ),
        max_unsupported_span=from_si(max_unsupported_span(diameter, tube_material), "m"),
        min_baffle_spacing=from_si(
            min_baffle_spacing(values["shell_inside_diameter"], tema_class), "m"
        ),
        max_baffle_spacing=from_si(values["shell_inside_diameter"], "m"),
        recommended_designation=recommend_configuration(
            needs_cleanable, large_temperature_difference, high_pressure, hazardous
        ),
        warnings=tuple(warnings),
    )
