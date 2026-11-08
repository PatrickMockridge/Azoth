"""``thermal.tube_insert_model`` - twisted tape, wire matrix and coiled wire enhancement.

Spec: ``specs/calcs/thermal/tube_insert_model.toml``

The Python twin of ``crates/azoth-thermal/src/tube_insert_model.rs``, and both mirror
``TubeInsertModel``: the enhancement and friction ratios of the three inserts, the performance
evaluation criterion, and the enhancement applied to a plain-tube state.

**The tape thickness cancels in the blockage factor** (the class writes ``t/(t/0.05)``), and
**the coiled wire's pitch cancels against itself** in the ``p/D``; both are carried and
unread. The twist ratio is floored at two, and the laminar branches read neither Reynolds
nor Prandtl number.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TubeInsertModelResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.tube_insert_model"
TWIST_RATIO = 5.0
TAPE_THICKNESS = 0.001
TAPE_WIDTH_FRACTION = 0.95
WIRE_DIAMETER = 0.001
WIRE_COIL_PITCH = 0.010
MATRIX_DENSITY = 0.5
ROUGHNESS_RATIO = 0.02
HELIX_ANGLE = 45.0


def tape_blockage(thickness: float, width: float) -> float:
    """The class's blockage factor, whose ``t/(t/0.05)`` cancels the thickness."""
    factor = math.pi / (math.pi - 4.0 * thickness * width / (thickness / 0.05))
    return min(max(factor, 1.0), 1.5)


def tape_swirl(twist_ratio: float) -> float:
    """``sqrt(1 + (pi/(2y))**2)``."""
    return math.sqrt(1.0 + math.pow(math.pi / (2.0 * twist_ratio), 2))


def twisted_tape_heat(reynolds: float, twist_ratio: float, thickness: float, width: float) -> float:
    """The twisted tape's heat-transfer enhancement ratio."""
    y = max(twist_ratio, 2.0)
    blockage = tape_blockage(thickness, width)
    if reynolds < 2300.0:
        return math.pow(1.0 + 0.769 / y, 1.0 / 3.0) * blockage
    swirl = tape_swirl(y)
    enhancement = math.pow(swirl, 0.8) * math.pow(blockage, 0.8)
    return max(1.0, enhancement * (1.0 + 2.0 / y) / swirl)


def twisted_tape_friction(
    reynolds: float, twist_ratio: float, thickness: float, width: float
) -> float:
    """The twisted tape's friction-penalty ratio."""
    y = max(twist_ratio, 2.0)
    blockage = tape_blockage(thickness, width)
    if reynolds < 2300.0:
        return (1.0 + 2.75 / y) * math.pow(blockage, 1.29)
    swirl = tape_swirl(y)
    return max(1.0, (1.0 + 2.75 / y) * math.pow(swirl, 1.8) / swirl)


def wire_matrix_heat(reynolds: float, density: float) -> float:
    """The wire matrix's heat-transfer enhancement ratio."""
    if reynolds < 500.0:
        return 1.0 + 5.0 * density
    if reynolds < 2300.0:
        return 1.0 + 4.0 * density * math.pow(reynolds / 2300.0, 0.3)
    return max(1.0, 1.0 + 2.5 * density / math.pow(reynolds / 10000.0, 0.2))


def wire_matrix_friction(reynolds: float, density: float) -> float:
    """The wire matrix's friction-penalty ratio."""
    base = 1.0 + 8.0 * density
    if reynolds < 2300.0:
        return base * 1.5
    return max(1.0, base / math.pow(reynolds / 10000.0, 0.15))


def coiled_wire_heat(
    reynolds: float,
    prandtl: float,
    roughness: float,
    pitch: float,
    helix_angle: float,
) -> float:
    """The coiled wire's heat-transfer enhancement ratio."""
    if reynolds < 2300.0:
        return 1.0 + 0.5 * roughness * 20.0
    p_over_d = pitch / (pitch / roughness)
    a = (
        2.64
        * math.pow(reynolds, 0.036)
        * math.pow(roughness, 0.212)
        * math.pow(max(p_over_d, 0.1), -0.21)
        * math.pow(helix_angle / 90.0, 0.29)
        * math.pow(prandtl, -0.024)
    )
    return math.pow(1.0 + math.pow(a, 7), 1.0 / 7.0)


def coiled_wire_friction(reynolds: float, roughness: float, helix_angle: float) -> float:
    """The coiled wire's friction-penalty ratio."""
    if reynolds < 2300.0:
        return 1.0 + 10.0 * roughness
    return max(
        1.0,
        1.0
        + 29.1
        * math.pow(roughness, 0.67)
        * math.pow(reynolds, -0.06)
        * math.pow(helix_angle / 90.0, 0.49),
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
    """The composed insert calculation.

    Raises:
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    plain_htc = input_to_si(
        spec, "plain_heat_transfer_coefficient", plain_heat_transfer_coefficient
    )
    plain_drop = input_to_si(spec, "plain_pressure_drop", plain_pressure_drop)
    values = {
        "reynolds_number": reynolds_number,
        "prandtl_number": prandtl_number,
        "twist_ratio": twist_ratio,
        "tape_thickness": None
        if tape_thickness is None
        else input_to_si(spec, "tape_thickness", tape_thickness),
        "tape_width_fraction": tape_width_fraction,
        "wire_diameter": None
        if wire_diameter is None
        else input_to_si(spec, "wire_diameter", wire_diameter),
        "wire_coil_pitch": None
        if wire_coil_pitch is None
        else input_to_si(spec, "wire_coil_pitch", wire_coil_pitch),
        "matrix_density": matrix_density,
        "roughness_ratio": roughness_ratio,
        "helix_angle": helix_angle,
        "plain_heat_transfer_coefficient": plain_htc,
        "plain_pressure_drop": plain_drop,
    }
    apply_checks(checks.on_input, values.get, warnings)

    y = TWIST_RATIO if values["twist_ratio"] is None else values["twist_ratio"]
    thickness = TAPE_THICKNESS if values["tape_thickness"] is None else values["tape_thickness"]
    width = (
        TAPE_WIDTH_FRACTION
        if values["tape_width_fraction"] is None
        else values["tape_width_fraction"]
    )
    density = MATRIX_DENSITY if values["matrix_density"] is None else values["matrix_density"]
    roughness = ROUGHNESS_RATIO if values["roughness_ratio"] is None else values["roughness_ratio"]
    pitch = WIRE_COIL_PITCH if values["wire_coil_pitch"] is None else values["wire_coil_pitch"]
    helix = HELIX_ANGLE if values["helix_angle"] is None else values["helix_angle"]

    if insert == "none":
        heat_ratio = 1.0
        friction_ratio = 1.0
    elif insert == "twisted_tape":
        heat_ratio = twisted_tape_heat(reynolds_number, y, thickness, width)
        friction_ratio = twisted_tape_friction(reynolds_number, y, thickness, width)
    elif insert == "wire_matrix":
        heat_ratio = wire_matrix_heat(reynolds_number, density)
        friction_ratio = wire_matrix_friction(reynolds_number, density)
    elif insert == "coiled_wire":
        heat_ratio = coiled_wire_heat(reynolds_number, prandtl_number, roughness, pitch, helix)
        friction_ratio = coiled_wire_friction(reynolds_number, roughness, helix)
    else:
        raise InvalidInputError(
            "insert", f"{insert} is not one of none, twisted_tape, wire_matrix or coiled_wire"
        )

    return TubeInsertModelResult(
        heat_transfer_enhancement_ratio=heat_ratio,
        pressure_drop_penalty_ratio=friction_ratio,
        performance_evaluation_criteria=heat_ratio / math.pow(friction_ratio, 1.0 / 3.0),
        enhanced_heat_transfer_coefficient=from_si(plain_htc * heat_ratio, "W/(m**2*K)"),
        enhanced_pressure_drop=from_si(plain_drop * friction_ratio, "Pa"),
        warnings=tuple(warnings),
    )
