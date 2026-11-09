"""``thermal.shell_and_tube_vibration`` - the tube-vibration screen.

Spec: ``specs/calcs/thermal/shell_and_tube_vibration.toml``

The Python twin of ``crates/azoth-thermal/src/shell_and_tube_vibration.rs``, and both
mirror ``VibrationAnalysis.performScreening``. **The screen fixes the natural frequency's
end condition at ``pinned``**, as the class does; the ``fixed`` and ``clamped-pinned``
branches are the helper's and are measured in the capture as evidence.

The Strouhal number is Fitzhugh's pitch correction, not the nominal 0.22: on the captured
state the ratio 1.2499 gives 0.3178 and a 33.37 Hz shedding frequency the class's own test
calls "~23 Hz". The three critical bands are strict ``0.8 < ratio < 1.2``.
"""

from __future__ import annotations

from typing import Any

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ShellAndTubeVibrationResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "thermal.shell_and_tube_vibration"

STROUHAL_DEFAULT = 0.22
CONNORS_TRIANGULAR = 3.3
CONNORS_SQUARE = 3.4
FLUID_ELASTIC_SAFETY_FACTOR = 0.5
CN_PINNED = 9.87


def natural_frequency(
    tube_od: float,
    tube_id: float,
    unsupported_span: float,
    material_e: float,
    tube_density: float,
    tube_fluid_density: float,
    shell_fluid_density: float,
) -> float:
    """The first natural frequency of a pinned span, Hz."""
    if tube_od <= 0.0 or tube_id >= tube_od or unsupported_span <= 0.0 or material_e <= 0.0:
        return 0.0
    second_moment = math.pi / 64.0 * (tube_od**4 - tube_id**4)
    tube_area = math.pi / 4.0 * (tube_od * tube_od - tube_id * tube_id)
    core_area = math.pi / 4.0 * tube_id * tube_id
    metal_mass = tube_density * tube_area
    internal_mass = tube_fluid_density * core_area
    hydrodynamic_mass = shell_fluid_density * math.pi / 4.0 * tube_od * tube_od
    effective_mass = metal_mass + internal_mass + hydrodynamic_mass
    if effective_mass <= 0.0:
        return 0.0
    return (
        CN_PINNED
        / (2.0 * math.pi)
        * math.sqrt(material_e * second_moment / (effective_mass * unsupported_span**4))
    )


def strouhal_number(pitch_ratio: float) -> float:
    """Fitzhugh's pitch correction inside ``1 < r < 3``, the class's default otherwise."""
    if 1.0 < pitch_ratio < 3.0:
        return STROUHAL_DEFAULT * (1.0 + 0.2 / (pitch_ratio - 0.8))
    return STROUHAL_DEFAULT


def vortex_shedding_frequency(
    crossflow_velocity: float, tube_od: float, tube_pitch: float
) -> float:
    """The Strouhal shedding frequency, Hz."""
    if tube_od <= 0.0 or crossflow_velocity <= 0.0:
        return 0.0
    return strouhal_number(tube_pitch / tube_od) * crossflow_velocity / tube_od


def critical_velocity_connors(
    natural_frequency_hz: float,
    tube_od: float,
    damping_ratio: float,
    effective_mass_per_length: float,
    shell_fluid_density: float,
    triangular: bool,
) -> float:
    """Connors' critical velocity; ``math.inf`` is the class's "not computable" sentinel."""
    if natural_frequency_hz <= 0.0 or tube_od <= 0.0 or shell_fluid_density <= 0.0:
        return math.inf
    k = CONNORS_TRIANGULAR if triangular else CONNORS_SQUARE
    mass_ratio = effective_mass_per_length / (shell_fluid_density * tube_od * tube_od)
    return (
        k * natural_frequency_hz * tube_od * math.sqrt(2.0 * math.pi * damping_ratio * mass_ratio)
    )


def acoustic_frequency(shell_id: float, sonic_velocity: float, mode: int) -> float:
    """The acoustic resonance frequency of the shell cavity, Hz."""
    if shell_id <= 0.0 or sonic_velocity <= 0.0 or mode < 1:
        return 0.0
    return mode * sonic_velocity / (2.0 * shell_id)


def _effective_mass(
    tube_od: float,
    tube_id: float,
    tube_density: float,
    tube_fluid_density: float,
    shell_fluid_density: float,
) -> float:
    tube_area = math.pi / 4.0 * (tube_od * tube_od - tube_id * tube_id)
    core_area = math.pi / 4.0 * tube_id * tube_id
    return (
        tube_density * tube_area
        + tube_fluid_density * core_area
        + shell_fluid_density * math.pi / 4.0 * tube_od * tube_od
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
    """The tube-vibration screen.

    Raises:
        InvalidInputError: for an inside diameter not smaller than the outside one, a
            negative damping ratio or velocity.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    values = {
        "tube_outside_diameter": input_to_si(spec, "tube_outside_diameter", tube_outside_diameter),
        "tube_inside_diameter": input_to_si(spec, "tube_inside_diameter", tube_inside_diameter),
        "unsupported_span": input_to_si(spec, "unsupported_span", unsupported_span),
        "tube_pitch": input_to_si(spec, "tube_pitch", tube_pitch),
        "tube_material_e": input_to_si(spec, "tube_material_e", tube_material_e),
        "tube_density": input_to_si(spec, "tube_density", tube_density),
        "crossflow_velocity": input_to_si(spec, "crossflow_velocity", crossflow_velocity),
        "shell_fluid_density": input_to_si(spec, "shell_fluid_density", shell_fluid_density),
        "tube_fluid_density": input_to_si(spec, "tube_fluid_density", tube_fluid_density),
        "shell_inside_diameter": input_to_si(spec, "shell_inside_diameter", shell_inside_diameter),
        "sonic_velocity": input_to_si(spec, "sonic_velocity", sonic_velocity),
        "damping_ratio": damping_ratio,
    }
    apply_checks(checks.on_input, values.get, warnings)

    if values["tube_inside_diameter"] >= values["tube_outside_diameter"]:
        raise InvalidInputError(
            "tube_inside_diameter",
            "the screen requires an inside diameter smaller than the outside diameter",
        )
    if damping_ratio < 0.0:
        raise InvalidInputError(
            "damping_ratio",
            "Connors' criterion takes the square root of the damping, so it cannot be negative",
        )
    if values["crossflow_velocity"] < 0.0:
        raise InvalidInputError("crossflow_velocity", "a velocity cannot be negative")

    screened = screen(
        values["tube_outside_diameter"],
        values["tube_inside_diameter"],
        values["unsupported_span"],
        values["tube_pitch"],
        values["tube_material_e"],
        values["tube_density"],
        values["crossflow_velocity"],
        values["shell_fluid_density"],
        values["tube_fluid_density"],
        values["shell_inside_diameter"],
        values["sonic_velocity"],
        damping_ratio,
        triangular_pitch,
    )
    return ShellAndTubeVibrationResult(
        natural_frequency=from_si(screened["fn"], "1/s"),
        vortex_shedding_frequency=from_si(screened["fvs"], "1/s"),
        vortex_shedding_ratio=screened["vortex_ratio"],
        vortex_shedding_critical=screened["vortex_critical"],
        critical_velocity=from_si(screened["vcrit"], "m/s"),
        velocity_ratio=screened["velocity_ratio"],
        fluid_elastic_critical=screened["fluid_elastic_critical"],
        acoustic_frequency=from_si(screened["fac"], "1/s"),
        acoustic_ratio=screened["acoustic_ratio"],
        acoustic_critical=screened["acoustic_critical"],
        passed=screened["passed"],
        warnings=tuple(warnings),
    )


def screen(
    tube_outside_diameter: float,
    tube_inside_diameter: float,
    unsupported_span: float,
    tube_pitch: float,
    tube_material_e: float,
    tube_density: float,
    crossflow_velocity: float,
    shell_fluid_density: float,
    tube_fluid_density: float,
    shell_inside_diameter: float,
    sonic_velocity: float,
    damping_ratio: float,
    triangular_pitch: bool,
) -> dict[str, Any]:
    """`performScreening` without the input bounds, in SI magnitudes."""
    fn_value = natural_frequency(
        tube_outside_diameter,
        tube_inside_diameter,
        unsupported_span,
        tube_material_e,
        tube_density,
        tube_fluid_density,
        shell_fluid_density,
    )
    fvs = vortex_shedding_frequency(
        crossflow_velocity,
        tube_outside_diameter,
        tube_pitch,
    )
    vortex_ratio = fvs / fn_value if fn_value > 0.0 else 0.0
    vortex_critical = 0.8 < vortex_ratio < 1.2

    mass = _effective_mass(
        tube_outside_diameter,
        tube_inside_diameter,
        tube_density,
        tube_fluid_density,
        shell_fluid_density,
    )
    vcrit = critical_velocity_connors(
        fn_value,
        tube_outside_diameter,
        damping_ratio,
        mass,
        shell_fluid_density,
        triangular_pitch,
    )
    velocity_ratio = (
        crossflow_velocity / (vcrit * FLUID_ELASTIC_SAFETY_FACTOR)
        if 0.0 < vcrit < math.inf
        else 0.0
    )
    fluid_elastic_critical = velocity_ratio > 1.0

    fac = acoustic_frequency(shell_inside_diameter, sonic_velocity, 1)
    acoustic_ratio = fvs / fac if fac > 0.0 else 0.0
    acoustic_critical = 0.8 < acoustic_ratio < 1.2

    return {
        "fn": fn_value,
        "fvs": fvs,
        "vortex_ratio": vortex_ratio,
        "vortex_critical": vortex_critical,
        "vcrit": vcrit,
        "velocity_ratio": velocity_ratio,
        "fluid_elastic_critical": fluid_elastic_critical,
        "fac": fac,
        "acoustic_ratio": acoustic_ratio,
        "acoustic_critical": acoustic_critical,
        "passed": not vortex_critical and not fluid_elastic_critical and not acoustic_critical,
    }
