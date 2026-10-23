"""``eos.pipeline_cooldown`` - NeqSim's ``PipelineCooldownCalculator``.

The lumped thermal model the class marches with explicit Euler: the fluid, steel and
insulation are one thermal mass losing heat through the layers (or a supplied U) to the
ambient temperature.

Spec: ``specs/models/eos/pipeline_cooldown.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PipelineCooldownResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.pipeline_cooldown"

STEEL_CONDUCTIVITY = 50.0
INSULATION_DENSITY = 600.0
INSULATION_SPECIFIC_HEAT = 1700.0
DEFAULT_U_VALUE = 3.0
DEFAULT_EXTERNAL_HTC = 300.0
MINIMUM_INSULATION = 0.001


def pipeline_cooldown(
    internal_diameter: Q,
    wall_thickness: Q,
    insulation_thickness: Q,
    coating_thickness: Q,
    steel_density: Q,
    steel_specific_heat: Q,
    insulation_conductivity: Q,
    fluid_density: Q,
    fluid_specific_heat: Q,
    initial_fluid_temperature: Q,
    ambient_temperature: Q,
    time_step: Q,
    total_time: Q,
    overall_u_value: Q | None = None,
    use_layer_calculation: bool | None = None,
    external_htc: Q | None = None,
    target_temperature: Q | None = None,
) -> PipelineCooldownResult:
    """The explicit-Euler pipeline cooldown march.

    Args:
        internal_diameter: the pipe's internal diameter.
        wall_thickness: the steel wall thickness.
        insulation_thickness: the insulation thickness, zero for a bare pipe.
        coating_thickness: the external coating thickness.
        steel_density: the steel's density.
        steel_specific_heat: the steel's specific heat capacity.
        insulation_conductivity: the insulation's thermal conductivity.
        fluid_density: the fluid's density at the operating point.
        fluid_specific_heat: the fluid's specific heat capacity.
        initial_fluid_temperature: the fluid temperature at shutdown.
        ambient_temperature: the seawater or soil temperature.
        time_step: the explicit-Euler step.
        total_time: the simulation horizon.
        overall_u_value: the supplied U value; absent takes 3.0.
        use_layer_calculation: whether U is computed from the layers; absent is false.
        external_htc: the external convection coefficient; absent takes 300.
        target_temperature: a temperature whose crossing time is wanted.

    Returns:
        The time and temperature profiles, the time constant, both U values and the
        target crossing.

    Raises:
        OutOfRangeError: if an input is outside the spec's hard range.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pipeline_cooldown(
        ...     q(0.254, "m"), q(0.0127, "m"), q(0.050, "m"), q(0.006, "m"),
        ...     q(7850.0, "kg/m**3"), q(500.0, "J/(kg*K)"), q(0.17, "W/(m*K)"),
        ...     q(750.0, "kg/m**3"), q(2200.0, "J/(kg*K)"),
        ...     q(353.15, "K"), q(277.15, "K"), q(3600.0, "s"), q(36000.0, "s"),
        ...     q(3.0, "W/(m**2*K)"))
        >>> len(r.time)
        11
    """
    from azoth._models_gen import model as _model_for

    spec: dict[str, Any] = _model_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values: dict[str, float] = {
        name: input_to_si(spec, name, value)
        for name, value in {
            "internal_diameter": internal_diameter,
            "wall_thickness": wall_thickness,
            "insulation_thickness": insulation_thickness,
            "coating_thickness": coating_thickness,
            "steel_density": steel_density,
            "steel_specific_heat": steel_specific_heat,
            "insulation_conductivity": insulation_conductivity,
            "fluid_density": fluid_density,
            "fluid_specific_heat": fluid_specific_heat,
            "initial_fluid_temperature": initial_fluid_temperature,
            "ambient_temperature": ambient_temperature,
            "time_step": time_step,
            "total_time": total_time,
        }.items()
    }
    optional = {
        "overall_u_value": overall_u_value,
        "external_htc": external_htc,
        "target_temperature": target_temperature,
    }
    for name, value in optional.items():
        if value is not None:
            values[name] = input_to_si(spec, name, value)
    apply_checks(checks.on_input, values.get, warnings)

    dt = values["time_step"]
    step_count = math.ceil(values["total_time"] / dt)

    radius_inner = values["internal_diameter"] / 2.0
    radius_outer = radius_inner + values["wall_thickness"]
    radius_insulation = radius_outer + values["insulation_thickness"]
    radius_external = radius_insulation + values["coating_thickness"]

    steel_resistance = math.log(radius_outer / radius_inner) / (2.0 * math.pi * STEEL_CONDUCTIVITY)
    insulation_resistance = (
        math.log(radius_insulation / radius_outer)
        / (2.0 * math.pi * values["insulation_conductivity"])
        if values["insulation_thickness"] > MINIMUM_INSULATION
        else 0.0
    )
    external_resistance = 1.0 / (
        2.0 * math.pi * radius_external * values.get("external_htc", DEFAULT_EXTERNAL_HTC)
    )
    calculated_u = 1.0 / (
        2.0
        * math.pi
        * radius_external
        * (steel_resistance + insulation_resistance + external_resistance)
    )

    supplied_u = values.get("overall_u_value", DEFAULT_U_VALUE)
    use_layers = bool(use_layer_calculation)
    effective_u = calculated_u if use_layers else supplied_u
    conductance = effective_u * 2.0 * math.pi * radius_external

    fluid_area = math.pi * radius_inner * radius_inner
    steel_area = math.pi * (radius_outer * radius_outer - radius_inner * radius_inner)
    insulation_area = math.pi * (
        radius_insulation * radius_insulation - radius_outer * radius_outer
    )
    heat_capacity = (
        values["fluid_density"] * fluid_area * values["fluid_specific_heat"]
        + values["steel_density"] * steel_area * values["steel_specific_heat"]
        + INSULATION_DENSITY * insulation_area * INSULATION_SPECIFIC_HEAT
    )

    time = [0.0]
    fluid_temperature = [values["initial_fluid_temperature"]]
    for index in range(step_count):
        previous = fluid_temperature[index]
        heat_loss = conductance * (previous - values["ambient_temperature"])
        fluid_temperature.append(previous - heat_loss / heat_capacity * dt)
        time.append((index + 1) * dt)

    target = values.get("target_temperature")
    time_to_reach_target: float | None = None
    if target is not None:
        for index, temperature in enumerate(fluid_temperature):
            if temperature <= target:
                if index == 0:
                    time_to_reach_target = 0.0
                else:
                    t1, t2 = time[index - 1], time[index]
                    temperature_1 = fluid_temperature[index - 1]
                    temperature_2 = fluid_temperature[index]
                    time_to_reach_target = t1 + (t2 - t1) * (temperature_1 - target) / (
                        temperature_1 - temperature_2
                    )
                break

    return PipelineCooldownResult(
        time=tuple(from_si(value, "s") for value in time),
        fluid_temperature=tuple(from_si(value, "K") for value in fluid_temperature),
        time_constant=from_si(heat_capacity / conductance, "s"),
        calculated_u_value=from_si(calculated_u, "W/(m**2*K)"),
        effective_u_value=from_si(effective_u, "W/(m**2*K)"),
        time_to_reach_target=(
            None if time_to_reach_target is None else from_si(time_to_reach_target, "s")
        ),
        warnings=tuple(warnings),
    )
