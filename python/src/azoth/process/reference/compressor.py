"""``process.compressor`` - the isentropic route or the machine's chart.

Spec: ``specs/models/process/compressor.toml``

The Python twin of ``crates/azoth-process/src/models/compressor.rs`` over
``crates/azoth-process/src/kernels/compressor.rs``, written to mirror them.

# The step is isentropic, and the inlet's entropy is derived

NeqSim reads ``getEntropy()`` off the inlet system, sets the outlet pressure, runs a
``PSflash`` on that entropy and takes the enthalpy it lands on as the reversible outlet. The
record has no entropy field - ``s = f(T, P, z)``, the same flash that gives ``h`` - so a port
that agrees with the captured efficiency-of-one rows has computed it rather than carried it.

# The efficiency divides the isentropic step

``h_out = h_in + (h_isentropic - h_in) / eta``, so an efficiency below one costs more
work than the reversible step. An expansion's difference is negative, so the same rule with a
multiplication is `process.expander`'s - dividing there would beat the reversible machine.

# A curve selects the fixed-speed chart branch

``CompressorChart``'s reduced head and polytropic efficiency, then the target's
``chartPressureForHead``: the historical kappa ratio brackets the pressure and the selected
polytropic reference - Schultz by default - is inverted to the chart head. **The chart's head
crosses in metres** and the branch converts it with ``9.81`` - the class's own
``polytropicHeadMeter`` factor.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, replace
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CompressorResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._heat_capacity import phase_capacity, volume_derivatives
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.eos.reference.sound_speed import sound_speed as _sound_speed
from azoth.process.reference import _unported, capacity, compressor_chart
from azoth.process.reference.compressor_mechanical import (
    design as design_compressor,
)
from azoth.process.reference.compressor_mechanical import (
    head_from_metres,
    head_from_power,
)
from azoth.process.reference.surge_curve import SurgeCurve

#: ``Compressor.R``, J/(mol*K).
R = 8.314
#: ``Compressor``'s own metre-to-specific-energy factor, not standard gravity.
METRE_TO_JOULE_PER_KG = 9.81
#: ``Compressor.MIN_EFFICIENCY``.
MIN_EFFICIENCY = 1.0e-6
#: ``Compressor.speed``'s initialiser, rpm.
DEFAULT_SPEED = 3000.0
#: ``Compressor.run``'s speed-solve loop: ``while (iteration < 200)``.
SPEED_SOLVE_ITERATIONS = 200
#: The numerical derivative's speed increment, ``deltaSpeed``, rpm.
SPEED_SOLVE_DELTA_RPM = 100.0
#: The derivative floor, ``1e-6`` bar/rpm, in pascals per rpm.
SPEED_SOLVE_DERIVATIVE_FLOOR_PA = 1.0e-4
#: The fallback slope, ``0.01`` bar/rpm, in pascals per rpm.
SPEED_SOLVE_DERIVATIVE_FALLBACK_PA = 1000.0
#: ``polytropicFluidHead = 0.0001`` kJ/kg, in joules per kilogram.
MIN_POLYTROPIC_HEAD_J_PER_KG = 0.1


def compressor(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    isentropic_efficiency: float,
    use_polytropic_calc: bool | None = None,
    polytropic_efficiency: float | None = None,
    polytropic_method: str | None = None,
    compression_ratio: float | None = None,
    max_outlet_pressure: Q | None = None,
    power: Q | None = None,
    deposit: list[dict[str, Q]] | None = None,
    wetted_area: Q | None = None,
    passage_half_height: Q | None = None,
    outlet_pressure: Q | None = None,
    curve: list[dict[str, Q]] | None = None,
    speed: float | None = None,
    solve_speed: bool | None = None,
    min_speed: float | None = None,
    max_speed: float | None = None,
    limit_speed: bool | None = None,
    chart_type: str | None = None,
    chart_conditions: list[float] | None = None,
    impeller_diameter: Q | None = None,
    surge_curve: list[dict[str, Q]] | None = None,
    surge_control_margin: float | None = None,
    anti_surge_active: bool | None = None,
    surge_control_factor: float | None = None,
    mechanical_design: bool | None = None,
    max_stages_per_casing: float | None = None,
) -> CompressorResult:
    """Raise a stream's pressure along an isentrope or the machine's chart.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        isentropic_efficiency: the machine's isentropic efficiency, in ``(0, 1]``.
        outlet_pressure: the pressure the isentropic branch takes the stream to.
        curve: the chart's measured points; an empty list is the isentropic branch.
        speed: the machine's speed, rpm; absent takes ``Compressor.speed``'s 3000. In the
            speed-solve mode it is the Newton's starting point.
        solve_speed: solve the speed that reaches ``outlet_pressure`` through the chart.
        min_speed: the speed-solve mode's floor, rpm; absent takes the class's 0.
        max_speed: the speed-solve mode's ceiling, rpm; absent takes the class's 30000.
        limit_speed: peg the solve at a bound and report the chart's own pressure there.
        chart_type: which map class; absent is ``fan_law``, the others refused by name.
        surge_curve: the surge boundary's ``(flow, head)`` rows; absent is no boundary.
        surge_control_margin: the control line as a fraction above the surge flow.
        anti_surge_active: whether the run-time recycle block runs.
        surge_control_factor: the recycle target's multiple of the surge flow.

    Returns:
        The outlet's record, the chart head and efficiency, and the shaft power.

    Raises:
        InvalidInputError: where the shapes disagree, the efficiency is outside ``(0, 1]``, a
            compressor with no curve has no outlet pressure, or the chart type is refused.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = compressor(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     isentropic_efficiency=0.75,
        ...     outlet_pressure=q(60.0, "bar"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 3)
        380.212
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    if len(inlet_z) != len(components):
        raise InvalidInputError(
            "z", f"a feed for {len(components)} components has {len(inlet_z)} entries"
        )

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    p_out = (
        None if outlet_pressure is None else input_to_si(spec, "outlet_pressure", outlet_pressure)
    )
    if compression_ratio is not None:
        ratio_in = float(compression_ratio)
        if not (math.isfinite(ratio_in) and ratio_in > 0.0):
            raise InvalidInputError(
                "compression_ratio", f"a compression ratio is positive, and {ratio_in} is not"
            )
        p_out = p * ratio_in
        if max_outlet_pressure is not None:
            p_out = min(p_out, input_to_si(spec, "max_outlet_pressure", max_outlet_pressure))
    power_si = None if power is None else input_to_si(spec, "power", power)
    deposit_state = _deposit_state(
        spec,
        deposit,
        None if wetted_area is None else input_to_si(spec, "wetted_area", wetted_area),
        None
        if passage_half_height is None
        else input_to_si(spec, "passage_half_height", passage_half_height),
    )

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "isentropic_efficiency": isentropic_efficiency,
            "speed": speed,
            "surge_control_margin": surge_control_margin,
        }.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the machine's isentropic efficiency must be in (0, 1], and "
            f"{isentropic_efficiency} is not",
        )

    chart_type = "fan_law" if chart_type is None else chart_type
    kind = {
        "fan_law": compressor_chart.FAN_LAW,
        "interpolate": compressor_chart.MAP_LOOKUP,
        "interpolate_and_extrapolate": compressor_chart.MAP_LOOKUP_EXTRAPOLATE,
        "khader_2015": compressor_chart.KHADER_2015,
    }.get(chart_type)
    if kind is None:
        raise InvalidInputError(
            "chart_type", f"'{chart_type}' is not one of the four declared chart types"
        )

    rows: list[dict[str, float]] = (
        [] if curve is None else input_records_to_si(spec, "curve", curve)
    )
    points = [
        compressor_chart.CurvePoint(
            speed=float(row["speed"]),
            flow=float(row["flow"]) * 3600.0,
            head=float(row["head"]),
            efficiency=float(row["efficiency"]),
        )
        for row in rows
    ]

    surge_rows: list[dict[str, float]] = (
        [] if surge_curve is None else input_records_to_si(spec, "surge_curve", surge_curve)
    )
    settings = _SurgeSettings(
        active=bool(anti_surge_active),
        control_factor=1.05 if surge_control_factor is None else float(surge_control_factor),
        control_margin=0.0 if surge_control_margin is None else float(surge_control_margin),
    )

    polytropic_exponent: float | None = None
    shaft_power: float = 0.0
    if not points:
        if p_out is None:
            raise InvalidInputError(
                "outlet_pressure",
                "a compressor without a curve reaches only the pressure it is given, and this "
                "one has neither",
            )
        exponent = None
        if power_si is not None:
            stated = (
                (1.0 if polytropic_efficiency is None else float(polytropic_efficiency))
                * (1.0 if deposit_state is None else deposit_state.efficiency_multiplier)
                if use_polytropic_calc
                else None
            )
            h_out, outlet_t, efficiency = _power_states(
                components, t, p, inlet_z, p_out, power_si, n, stated
            )
            shaft_power = power_si
        elif use_polytropic_calc:
            method = "default" if polytropic_method is None else polytropic_method
            if method == "detailed":
                raise _unported.refuse("polytropic_method=detailed")
            if method == "schultz":
                raise _unported.refuse("polytropic_method=schultz")
            if method != "default":
                raise InvalidInputError(
                    "polytropic_method",
                    f"`{method}` is not one of default, detailed or schultz",
                )
            eta = 1.0 if polytropic_efficiency is None else float(polytropic_efficiency)
            if deposit_state is not None:
                eta *= deposit_state.efficiency_multiplier
            if not 0.0 < eta <= 1.0:
                raise InvalidInputError(
                    "polytropic_efficiency",
                    f"the compressor's polytropic efficiency is in (0, 1], and {eta} is not",
                )
            h_in, h_out, outlet_t, exponent = _polytropic_states(
                components, t, p, inlet_z, p_out, eta
            )
            efficiency = eta
            polytropic_exponent = exponent
            shaft_power = n * (h_out - h_in)
        else:
            states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)
            efficiency = 0.0
            shaft_power = n * states.dh_actual
            outlet_t = states.outlet_t
            h_out = states.h_out
        head = 0.0
        out_p = p_out
        state = _ChartState(
            head=0.0,
            efficiency=0.0,
            power=shaft_power,
            out_p=p_out,
            h_out=h_out,
            outlet_t=outlet_t,
            read=None,
            recycle_fraction=None,
            speed=DEFAULT_SPEED if speed is None else float(speed),
            speed_limited=False,
        )
    else:
        if kind == compressor_chart.KHADER_2015:
            if chart_conditions is None or len(chart_conditions) != 4:
                raise InvalidInputError(
                    "chart_conditions",
                    "the Khader chart builds its reference fluid at the chart conditions, and "
                    "this machine states none",
                )
            reference = _reference_sound_speed(chart_conditions)
            feed_speed = _sound_speed(
                components, quantity(t, "K"), quantity(p, "Pa"), list(inlet_z)
            ).system_sound_speed
            if feed_speed is None:
                raise InvalidInputError(
                    "inlet_t", "the feed's sound speed is not a number, and the map needs one"
                )
            process = feed_speed.to("m/s").magnitude
            diameter = (
                0.3
                if impeller_diameter is None
                else input_to_si(spec, "impeller_diameter", impeller_diameter)
            )
            chart = compressor_chart.CompressorChart.from_khader(
                points, reference, process, diameter
            )
        else:
            chart = compressor_chart.CompressorChart(points, kind)
        if surge_rows:
            chart.set_surge_curve(
                SurgeCurve.from_points(
                    [float(row["flow"]) * 3600.0 for row in surge_rows],
                    [float(row["head"]) for row in surge_rows],
                )
            )
        if solve_speed:
            if p_out is None:
                raise InvalidInputError(
                    "outlet_pressure",
                    "the speed-solve mode drives the chart to a stated pressure, and this "
                    "machine states none",
                )
            state = _speed_solve(
                components,
                n,
                inlet_z,
                t,
                p,
                p_out,
                chart,
                DEFAULT_SPEED if speed is None else float(speed),
                0.0 if min_speed is None else float(min_speed),
                30000.0 if max_speed is None else float(max_speed),
                bool(limit_speed),
                settings,
            )
        else:
            state = _chart_states(
                components,
                n,
                inlet_z,
                t,
                p,
                chart,
                DEFAULT_SPEED if speed is None else float(speed),
                settings,
            )
        head, efficiency, shaft_power = state.head, state.efficiency, state.power
        out_p, h_out, outlet_t = state.out_p, state.h_out, state.outlet_t

    design = None
    if mechanical_design:
        if points:
            fluid_head = head_from_metres(head)
        elif use_polytropic_calc:
            density = capacity.system_mass_density(components, t, p, inlet_z, from_feed=True)
            volume = capacity.system_volumetric_flow(components, t, p, inlet_z, n, from_feed=True)
            fluid_head = head_from_power(shaft_power, density * volume, efficiency)
        else:
            fluid_head = 0.0
        if not (math.isfinite(fluid_head) and fluid_head > 0.0):
            raise InvalidInputError(
                "mechanical_design",
                "the class's polytropic fluid head is zero on this branch, and `calcDesign` "
                "returns before it sizes anything",
            )
        design = design_compressor(
            components,
            n,
            list(inlet_z),
            t,
            p,
            out_p,
            outlet_t.to("K").magnitude,
            shaft_power,
            fluid_head,
            state.speed,
            10 if max_stages_per_casing is None else int(max_stages_per_casing),
        )
    d = design

    return CompressorResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(out_p, "Pa"),
        outlet_t=outlet_t,
        outlet_h=from_si(h_out, "J/mol"),
        polytropic_head=from_si(head, "m"),
        polytropic_efficiency=efficiency,
        polytropic_exponent=polytropic_exponent,
        total_deposit_mass=(
            None if deposit_state is None else from_si(deposit_state.total_mass, "kg")
        ),
        deposit_volume=(None if deposit_state is None else from_si(deposit_state.volume, "m**3")),
        film_thickness=(
            None if deposit_state is None else from_si(deposit_state.film_thickness, "m")
        ),
        area_blockage_fraction=None if deposit_state is None else deposit_state.blockage,
        efficiency_multiplier=(
            None if deposit_state is None else deposit_state.efficiency_multiplier
        ),
        head_multiplier=None if deposit_state is None else deposit_state.head_multiplier,
        degradation_factor=(None if deposit_state is None else deposit_state.efficiency_multiplier),
        fouling_factor=(None if deposit_state is None else 1.0 - deposit_state.head_multiplier),
        power=from_si(shaft_power, "W"),
        shaft_speed=state.speed,
        surge_flow_rate=_flow_or_none(state.read, "surge_flow_rate"),
        surge_flow_rate_margin=_flow_or_none(state.read, "surge_flow_rate_margin"),
        surge_flow_rate_std=_flow_or_none(state.read, "surge_flow_rate_std"),
        distance_to_surge=_plain_or_none(state.read, "distance_to_surge"),
        control_line_flow=_flow_or_none(state.read, "control_line_flow"),
        distance_to_control_line=_plain_or_none(state.read, "distance_to_control_line"),
        required_recycle_fraction=_plain_or_none(state.read, "required_recycle_fraction"),
        anti_surge_recycle_power=_power_or_none(state.read, "anti_surge_recycle_power"),
        recycle_fraction=state.recycle_fraction,
        design_number_of_stages=None if d is None else float(d.number_of_stages),
        design_impeller_diameter=None if d is None else from_si(d.impeller_diameter, "m"),
        design_shaft_diameter=None if d is None else from_si(d.shaft_diameter, "m"),
        design_tip_speed=None if d is None else from_si(d.tip_speed, "m/s"),
        design_flow_coefficient=None if d is None else d.flow_coefficient,
        design_impeller_sizing_speed=None if d is None else d.impeller_sizing_speed,
        design_impeller_sizing_flow=None
        if d is None
        else from_si(d.impeller_sizing_flow / 3600.0, "m**3/h"),
        design_impeller_sizing_head=None if d is None else d.impeller_sizing_head,
        design_impeller_sizing_ok=None if d is None else d.impeller_sizing_ok,
        design_driver_power=None if d is None else from_si(d.driver_power, "W"),
        design_driver_margin=None if d is None else d.driver_margin,
        design_pressure=None if d is None else from_si(d.design_pressure, "Pa"),
        design_temperature=None if d is None else from_si(d.design_temperature, "K"),
        design_max_continuous_speed=None if d is None else d.max_continuous_speed,
        design_trip_speed=None if d is None else d.trip_speed,
        design_first_critical_speed=None if d is None else d.first_critical_speed,
        design_casing_type=None if d is None else d.casing_type,
        design_head_per_stage=None if d is None else d.head_per_stage,
        design_bearing_span=None if d is None else from_si(d.bearing_span, "m"),
        design_casing_weight=None if d is None else from_si(d.casing_weight, "kg"),
        design_bundle_weight=None if d is None else from_si(d.bundle_weight, "kg"),
        design_rotor_weight=None if d is None else from_si(d.rotor_weight, "kg"),
        design_outer_diameter=None if d is None else from_si(d.outer_diameter, "m"),
        design_inner_diameter=None if d is None else from_si(d.inner_diameter, "m"),
        design_wall_thickness=None if d is None else from_si(d.wall_thickness, "m"),
        design_tantan_length=None if d is None else from_si(d.tantan_length, "m"),
        design_weight_vessel_shell=None if d is None else from_si(d.weight_vessel_shell, "kg"),
        design_weight_internals=None if d is None else from_si(d.weight_internals, "kg"),
        design_weight_nozzle=None if d is None else from_si(d.weight_nozzle, "kg"),
        design_weight_piping=None if d is None else from_si(d.weight_piping, "kg"),
        design_weight_electro_instrument=None
        if d is None
        else from_si(d.weight_electro_instrument, "kg"),
        design_weight_structural_steel=None
        if d is None
        else from_si(d.weight_structural_steel, "kg"),
        design_weight_total=None if d is None else from_si(d.weight_total, "kg"),
        design_module_length=None if d is None else from_si(d.module_length, "m"),
        design_module_width=None if d is None else from_si(d.module_width, "m"),
        design_module_height=None if d is None else from_si(d.module_height, "m"),
        warnings=tuple(warnings),
    )


MAX_SURGE_PASSES = 20


@dataclass(frozen=True, slots=True)
class _SurgeSettings:
    """``AntiSurge``'s state that the chart branch reads."""

    active: bool
    control_factor: float
    control_margin: float


@dataclass(frozen=True, slots=True)
class _SurgeRead:
    """One row of the anti-surge read-outs."""

    surge_flow_rate: float
    surge_flow_rate_margin: float
    surge_flow_rate_std: float
    distance_to_surge: float
    control_line_flow: float
    distance_to_control_line: float | None
    required_recycle_fraction: float
    anti_surge_recycle_power: float


@dataclass(frozen=True, slots=True)
class _ChartState:
    """The chart branch's answer, its surge read-outs and the recycle fraction."""

    head: float
    efficiency: float
    power: float
    out_p: float
    h_out: float
    outlet_t: Q
    read: _SurgeRead | None
    recycle_fraction: float | None
    speed: float
    speed_limited: bool


def _flow_or_none(read: _SurgeRead | None, name: str) -> Q | None:
    """One optional flow read-out; the twin holds it in m3/hr, and the boundary crosses SI."""
    if read is None:
        return None
    value = getattr(read, name)
    return None if value is None else from_si(float(value) / 3600.0, "m**3/h")


def _power_or_none(read: _SurgeRead | None, name: str) -> Q | None:
    """One optional power read-out, from SI."""
    if read is None:
        return None
    value = getattr(read, name)
    return None if value is None else from_si(value, "W")


def _plain_or_none(read: _SurgeRead | None, name: str) -> float | None:
    """One optional dimensionless read-out."""
    if read is None:
        return None
    value = getattr(read, name)
    return None if value is None else float(value)


def _reference_sound_speed(conditions: list[float]) -> float:
    """``createDefaultFluid``'s phase sound speed: methane/ethane/propane, rescaled to MW."""
    names = ["methane", "ethane", "propane"]
    mixture, _ideal_gas = _components.mixture_of(names, eos="pr")
    masses = [
        (component.molar_mass.to("kg/mol").magnitude * 1000.0)
        if component.molar_mass is not None
        else 0.0
        for component in mixture.components
    ]
    x1, x2, x3 = 0.90, 0.05, 0.05
    if conditions[3] > 0.0:
        a = (conditions[3] - masses[2]) / (
            x1 * (masses[0] - masses[2]) + x2 * (masses[1] - masses[2])
        )
        first = a * x1
        second = a * x2
        z = [first, second, 1.0 - first - second]
    else:
        z = [x1, x2, x3]
    state = _sound_speed(
        names,
        quantity(conditions[0] + 273.15, "K"),
        quantity(conditions[1] * 1.0e5, "Pa"),
        z,
    )
    return float(state.sound_speed.to("m/s").magnitude)


def _feed_state(
    mixture: Any, ideal_gas: Any, t: float, p: float, z: list[float]
) -> tuple[float, float, float]:
    """``getZ()``, ``getCP()/n`` and the molar volume at the feed's own composition."""
    flash = pt_flash_solve(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))
    reduced = reduced_parameters(mixture, t, p)
    kij = mixture.effective_kij(reduced.t_kelvin)

    def single(composition: list[float], liquid: bool) -> tuple[float, float, float]:
        state = phase_state(reduced, kij, list(composition), liquid=liquid)
        root = float(state.z)
        cp = float(
            molar_enthalpy_entropy(
                mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), list(composition), root
            )
            .cp.to_base_units()
            .magnitude
        )
        volume = root * R * t / p - mixture.volume_shift(list(composition))
        return root, cp, volume

    if flash.phase == "two_phase":
        beta = flash.vapour_fraction if flash.vapour_fraction is not None else 1.0
        vapour = single(list(flash.y), False)
        liquid = single(list(flash.x), True)
        return (
            beta * vapour[0] + (1.0 - beta) * liquid[0],
            beta * vapour[1] + (1.0 - beta) * liquid[1],
            beta * vapour[2] + (1.0 - beta) * liquid[2],
        )
    return single(list(z), flash.phase != "all_vapour")


def _feed_volume_derivatives(
    mixture: Any, ideal_gas: Any, t: float, p: float, z: list[float]
) -> tuple[float, float]:
    """``getdVdTpn`` and ``getdVdPtn`` at the feed's own root, m3/(mol K) and m3/Pa."""
    flash = pt_flash_solve(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))
    reduced = reduced_parameters(mixture, t, p)
    kij = mixture.effective_kij(reduced.t_kelvin)

    def single(composition: list[float], liquid: bool) -> tuple[float, float]:
        state = phase_state(reduced, kij, list(composition), liquid=liquid)
        capacity = phase_capacity(
            mixture,
            ideal_gas,
            reduced,
            kij,
            t,
            p,
            list(composition),
            float(state.z),
        )
        return volume_derivatives(capacity)

    if flash.phase == "two_phase":
        beta = flash.vapour_fraction if flash.vapour_fraction is not None else 1.0
        vapour = single(list(flash.y), False)
        liquid = single(list(flash.x), True)
        return (
            beta * vapour[0] + (1.0 - beta) * liquid[0],
            beta * vapour[1] + (1.0 - beta) * liquid[1],
        )
    return single(list(z), flash.phase != "all_vapour")


@dataclass(frozen=True, slots=True)
class _SchultzHead:
    """``createSchultzHeadFunction``'s four state numbers, evaluated at an efficiency."""

    schultz_x: float
    schultz_y: float
    isentropic_exponent: float
    head_factor: float
    pressure_ratio: float

    def head_j_per_kg(self, efficiency: float) -> float:
        inverse = 1.0 / efficiency
        exponent = (1.0 + self.schultz_x) / (
            (inverse + self.schultz_x) / self.isentropic_exponent - self.schultz_y * (inverse - 1.0)
        )
        term = exponent / (exponent - 1.0)
        return float(term * self.head_factor * (self.pressure_ratio ** (1.0 / term) - 1.0) * 1000.0)


def _schultz_head(
    mixture: Any,
    ideal_gas: Any,
    t: float,
    p: float,
    inlet_z: list[float],
    n: float,
    molar_mass: float,
    inlet_density: float,
    molar_volume: float,
    d_v_d_t: float,
    d_v_d_p: float,
    outlet_p: float,
) -> _SchultzHead:
    """Target ``createSchultzHeadFunction`` at ``outlet_p``."""
    if outlet_p <= p:
        raise InvalidInputError(
            "outlet_pressure",
            f"the Schultz reference compresses, and it is given {outlet_p} Pa against {p} Pa",
        )
    entropy = float(entropy_at(mixture, ideal_gas, t, p, list(inlet_z))[0])
    isentropic = ps_flash_solve(
        mixture,
        ideal_gas,
        quantity(outlet_p, "Pa"),
        from_si(entropy, "J/(mol*K)"),
        list(inlet_z),
    )
    isentropic_t = float(isentropic.T.to("K").magnitude)
    density_volume = _feed_state(mixture, ideal_gas, isentropic_t, outlet_p, list(inlet_z))[2]
    isentropic_density = molar_mass / density_volume
    h_isentropic = enthalpy_at(
        mixture,
        ideal_gas,
        isentropic_t,
        outlet_p,
        list(inlet_z),
    )[0]
    h_in = enthalpy_at(mixture, ideal_gas, t, p, list(inlet_z))[0]
    pressure_ratio = outlet_p / p
    density_ratio = isentropic_density / inlet_density
    if pressure_ratio <= 1.0 or density_ratio <= 1.0:
        raise InvalidInputError(
            "outlet_pressure",
            f"the Schultz reference has pressure ratio {pressure_ratio} and density ratio "
            f"{density_ratio}, which the isentropic exponent cannot divide",
        )
    isentropic_exponent = math.log(pressure_ratio) / math.log(density_ratio)
    correction_factor = (
        (h_isentropic - h_in)
        * n
        / (
            1.0e5
            * (outlet_p / isentropic_density - p / inlet_density)
            * isentropic_exponent
            / (isentropic_exponent - 1.0)
        )
    )
    head_factor = correction_factor * 1.0e5 * p / inlet_density / (n * molar_mass) / 1000.0
    return _SchultzHead(
        schultz_x=t / molar_volume * d_v_d_t - 1.0,
        schultz_y=-p / molar_volume * d_v_d_p,
        isentropic_exponent=isentropic_exponent,
        head_factor=head_factor,
        pressure_ratio=pressure_ratio,
    )


def _chart_pressure_for_head(
    components: list[str],
    inlet_z: list[float],
    t: float,
    p: float,
    head_j: float,
    efficiency: float,
    molar_mass: float,
    n: float,
) -> float:
    """``chartPressureForHead``: the discharge pressure whose Schultz head is the chart's.

    The target's chart branch brackets the historical kappa-ratio pressure and inverts the
    selected polytropic reference's head; its default method is Schultz. The reference's
    residual is ``getPolytropicHead("kJ/kg") - head``, which is the function below.
    """
    if head_j == 0.0:
        return p
    if not 0.0 < efficiency <= 1.0:
        raise InvalidInputError(
            "curve", f"the chart head's efficiency is in (0, 1], and {efficiency} is not"
        )
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    z, cp, molar_volume = _feed_state(mixture, ideal_gas, t, p, list(inlet_z))
    inlet_density = molar_mass / molar_volume
    d_v_d_t, d_v_d_p = _feed_volume_derivatives(mixture, ideal_gas, t, p, list(inlet_z))
    kappa = cp / (cp - R)
    ideal_pressure = p * _polytropic_ratio(head_j, z, t, molar_mass, kappa, efficiency)

    def head_at(outlet_p: float) -> float:
        return _schultz_head(
            mixture,
            ideal_gas,
            t,
            p,
            list(inlet_z),
            n,
            molar_mass,
            inlet_density,
            molar_volume,
            d_v_d_t,
            d_v_d_p,
            outlet_p,
        ).head_j_per_kg(efficiency)

    lower = 0.0
    lower_residual = -head_j
    upper = math.log(ideal_pressure / p) if ideal_pressure > p else math.log(2.0)
    upper_residual = head_at(p * math.exp(upper)) - head_j
    for _ in range(12):
        if upper_residual >= 0.0:
            break
        lower = upper
        lower_residual = upper_residual
        upper *= 2.0
        upper_residual = head_at(p * math.exp(upper)) - head_j
    if upper_residual < 0.0:
        raise InvalidInputError(
            "outlet_pressure",
            f"the chart head {head_j} J/kg cannot be bracketed above {p} Pa",
        )
    tolerance = max(1.0e-12, 1.0e-8 * head_j)
    for iteration in range(72):
        trial = upper - upper_residual * (upper - lower) / (upper_residual - lower_residual)
        if (
            iteration >= 12
            or not math.isfinite(trial)
            or trial <= lower
            or trial >= upper
            or iteration % 5 == 4
        ):
            trial = 0.5 * (lower + upper)
        residual = head_at(p * math.exp(trial)) - head_j
        if abs(residual) <= tolerance:
            return p * math.exp(trial)
        if residual < 0.0:
            lower = trial
            lower_residual = residual
        else:
            upper = trial
            upper_residual = residual
    raise InvalidInputError("outlet_pressure", "the chart head-pressure solve did not converge")


def _chart_states(
    components: list[str],
    n: float,
    inlet_z: list[float],
    t: float,
    p: float,
    chart: compressor_chart.CompressorChart,
    speed: float,
    settings: _SurgeSettings,
) -> _ChartState:
    """The chart branch: head, efficiency, power, outlet pressure, enthalpy and temperature."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))
    # **The feed's own root, at the feed's composition** - `getZ()` and `getCP()/n` on the one
    # phase the class holds, and the molar volume behind `getFlowRate("m3/sec")`.
    z, _cp, molar_volume = _feed_state(mixture, ideal_gas, t, p, list(inlet_z))
    fresh_flow = n * molar_volume * 3600.0
    if not n > 0.0 or not fresh_flow > 0.0 or not math.isfinite(n):
        raise InvalidInputError(
            "inlet_n",
            f"the chart branch needs a positive flow to scale, and {n} mol/s gives "
            f"{fresh_flow} m3/hr",
        )
    # ``getZvolcorr`` = ``P v/(R T)`` at the Peneloux-shifted molar volume.
    zvolcorr = z - p * mixture.volume_shift(list(inlet_z)) / (R * t)

    boosted = n
    passes = 0
    while True:
        flow_m3hr = fresh_flow * (boosted / n)
        efficiency = min(1.0, max(MIN_EFFICIENCY, chart.polytropic_efficiency(flow_m3hr, speed)))
        head = chart.polytropic_head(flow_m3hr, speed)
        head_j = head * METRE_TO_JOULE_PER_KG

        if settings.active and chart.is_surge(head, flow_m3hr):
            target = chart.surge_flow(head)
            if target is not None and math.isfinite(target) and target > 0.0:
                boosted = n * settings.control_factor * target / fresh_flow
                passes += 1
                if passes < MAX_SURGE_PASSES:
                    continue
        break

    out_p = _chart_pressure_for_head(
        components,
        inlet_z,
        t,
        p,
        head_j,
        efficiency,
        molar_mass,
        boosted,
    )
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    dh = head_j * molar_mass / efficiency
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, quantity(out_p, "Pa"), quantity(h_out, "J/mol"), inlet_z
    ).T
    power = boosted * dh
    recycle_fraction = boosted / n - 1.0

    read = _surge_read(
        chart,
        settings,
        p,
        t,
        zvolcorr,
        fresh_flow,
        power / (n * molar_mass) * efficiency / METRE_TO_JOULE_PER_KG,
        power,
        boosted,
        n,
    )
    return _ChartState(
        head=head,
        efficiency=efficiency,
        power=power,
        out_p=out_p,
        h_out=h_out,
        outlet_t=outlet_t,
        read=read,
        recycle_fraction=recycle_fraction if settings.active else None,
        speed=speed,
        speed_limited=False,
    )


def _surge_read(
    chart: compressor_chart.CompressorChart,
    settings: _SurgeSettings,
    p: float,
    t: float,
    zvolcorr: float,
    fresh_flow: float,
    head_for_curve: float,
    power: float,
    boosted: float,
    n: float,
) -> _SurgeRead | None:
    """The anti-surge read-outs at ``head_for_curve`` and the fresh flow."""
    surge_flow = chart.surge_flow(head_for_curve)
    if surge_flow is None or not surge_flow > 0.0:
        return None
    control = surge_flow * (1.0 + settings.control_margin)
    required = max(0.0, (control - fresh_flow) / control) if control > 0.0 else 0.0
    return _SurgeRead(
        surge_flow_rate=surge_flow,
        surge_flow_rate_margin=fresh_flow - surge_flow,
        surge_flow_rate_std=surge_flow * (p / 1.0e5) / 1.01325 * 288.15 / t / zvolcorr,
        distance_to_surge=fresh_flow / surge_flow - 1.0,
        control_line_flow=control,
        distance_to_control_line=(fresh_flow / control - 1.0) if control > 0.0 else None,
        required_recycle_fraction=required,
        # ``getAntiSurgeRecyclePower(currentSurgeFraction)``: the added power.
        anti_surge_recycle_power=power * max(boosted / n - 1.0, 0.0),
    )


def _polytropic_ratio(
    head_j: float,
    z: float,
    t: float,
    molar_mass: float,
    kappa: float,
    efficiency: float,
) -> float:
    """The class's polytropic pressure ratio: ``((H + W)/W)^(n/(n-1))``."""
    exponent = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency)
    flow_work = exponent / (exponent - 1.0) * z * R * t / molar_mass
    if not math.isfinite(flow_work) or not flow_work > 0.0 or not math.isfinite(exponent):
        raise InvalidInputError(
            "curve",
            f"the chart's head and the inlet's state give no polytropic ratio: n={exponent}, "
            f"flow work={flow_work}",
        )
    return float(((head_j + flow_work) / flow_work) ** (exponent / (exponent - 1.0)))


def _speed_solve(
    components: list[str],
    n: float,
    inlet_z: list[float],
    t: float,
    p: float,
    target: float,
    chart: compressor_chart.CompressorChart,
    initial_speed: float,
    min_speed: float,
    max_speed: float,
    limit_speed: bool,
    settings: _SurgeSettings,
) -> _ChartState:
    """``Compressor.run``'s ``solveSpeed`` branch: the speed that reaches a stated pressure."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))
    z, _cp, molar_volume = _feed_state(mixture, ideal_gas, t, p, list(inlet_z))
    fresh_flow = n * molar_volume * 3600.0
    if not n > 0.0 or not fresh_flow > 0.0 or not math.isfinite(n):
        raise InvalidInputError(
            "inlet_n",
            f"the chart branch needs a positive flow to scale, and {n} mol/s gives "
            f"{fresh_flow} m3/hr",
        )
    zvolcorr = z - p * mixture.volume_shift(list(inlet_z)) / (R * t)
    inlet_density = molar_mass / molar_volume
    d_v_d_t, d_v_d_p = _feed_volume_derivatives(mixture, ideal_gas, t, p, list(inlet_z))
    schultz = _schultz_head(
        mixture,
        ideal_gas,
        t,
        p,
        list(inlet_z),
        n,
        molar_mass,
        inlet_density,
        molar_volume,
        d_v_d_t,
        d_v_d_p,
        target,
    )

    if not initial_speed > 0.0 or not math.isfinite(initial_speed):
        raise InvalidInputError(
            "speed",
            f"the speed-solve mode starts from the machine's speed, and {initial_speed} rpm "
            "does not",
        )
    chart_min = chart.min_curve_speed()
    chart_max = chart.max_curve_speed()
    if not math.isfinite(chart_min) or not math.isfinite(chart_max) or chart_max <= chart_min:
        raise InvalidInputError(
            "curve",
            f"the speed-solve mode needs a measured speed range, and the chart's is "
            f"{chart_min} to {chart_max} rpm",
        )
    speed = initial_speed
    if speed < chart_min * 0.5 or speed > chart_max * 2.0:
        speed = 0.5 * (chart_min + chart_max)
    delta = SPEED_SOLVE_DELTA_RPM
    efficiency = 0.0
    head = 0.0
    dh = 0.0
    converged: float | None = None
    limited: float | None = None

    for iteration in range(1, SPEED_SOLVE_ITERATIONS):
        percent = chart.polytropic_efficiency(fresh_flow, speed) * 100.0
        if percent <= 0.0:
            stored, formula = 0.01, 0.01
        elif percent > 100.0:
            stored, formula = 100.0, 1.0
        else:
            stored = formula = percent / 100.0
        head = chart.polytropic_head(fresh_flow, speed)
        head_j = head * METRE_TO_JOULE_PER_KG
        if head_j <= 0.0:
            head_j = MIN_POLYTROPIC_HEAD_J_PER_KG
        current_residual = head_j - schultz.head_j_per_kg(formula)
        residual_tolerance = 1.0e-8 * max(1.0, head_j)
        # The derivative's read is never clamped, as the class's second read is not.
        percent_delta = chart.polytropic_efficiency(fresh_flow, speed + delta) * 100.0
        head_delta = chart.polytropic_head(fresh_flow, speed + delta) * METRE_TO_JOULE_PER_KG
        delta_residual = head_delta - schultz.head_j_per_kg(percent_delta / 100.0)
        derivative = (delta_residual - current_residual) / delta
        if not math.isfinite(derivative):
            raise InvalidInputError(
                "curve",
                f"the chart's derivative at {speed} rpm is not a number: efficiency "
                f"{percent_delta}%, head {head_delta} J/kg",
            )
        if abs(derivative) < SPEED_SOLVE_DERIVATIVE_FLOOR_PA:
            derivative = SPEED_SOLVE_DERIVATIVE_FALLBACK_PA
        relaxation = min(0.8, iteration / (iteration + 3.0))
        cap = max(500.0, 0.2 * speed)
        update = -current_residual / derivative
        if abs(update) > cap:
            update = math.copysign(cap, update)
        speed += relaxation * update
        efficiency = stored
        dh = head_j * molar_mass / stored
        if speed < 0.0:
            speed = min_speed if min_speed > 1.0 else chart_min
        if iteration % 10 == 0 and delta > 10.0:
            delta /= 2.0
        if (speed < min_speed or speed > max_speed) and limit_speed:
            limited = max_speed if speed > max_speed else min_speed
            break
        if abs(current_residual) <= residual_tolerance:
            converged = speed
            break

    if limited is not None:
        state = _chart_states(components, n, list(inlet_z), t, p, chart, limited, settings)
        return replace(state, speed_limited=True)
    if converged is not None:
        speed = converged

    # The class's post-convergence surge check: one pass, no re-solve of the speed.
    boosted = n
    if settings.active and chart.is_surge(head, fresh_flow):
        target_flow = chart.surge_flow(head)
        if target_flow is not None and math.isfinite(target_flow) and target_flow > 0.0:
            boosted = n * settings.control_factor * target_flow / fresh_flow
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, quantity(target, "Pa"), quantity(h_out, "J/mol"), inlet_z
    ).T
    power = boosted * dh
    read = _surge_read(chart, settings, p, t, zvolcorr, fresh_flow, head, power, boosted, n)
    return _ChartState(
        head=head,
        efficiency=efficiency,
        power=power,
        out_p=target,
        h_out=h_out,
        outlet_t=outlet_t,
        read=read,
        recycle_fraction=(boosted / n - 1.0) if settings.active else None,
        speed=speed,
        speed_limited=False,
    )


@dataclass(frozen=True, slots=True)
class CompressorStates:
    """The machine's intermediates, in the order it computes them, in SI."""

    #: The inlet's molar enthalpy and entropy, J/mol and J/(mol*K).
    h_in: float
    s_in: float
    #: The reversible outlet: the temperature at the outlet pressure on the inlet's entropy,
    #: and its molar enthalpy.
    t_isentropic: Q
    h_isentropic: float
    #: The isentropic step, and the step the efficiency actually asks for.
    dh_isentropic: float
    dh_actual: float
    #: The outlet: its molar enthalpy, its temperature, and its molar entropy - which rises
    #: across an irreversible step and is unchanged across a reversible one.
    h_out: float
    outlet_t: Q
    s_out: float


def _polytropic_states(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    efficiency: float,
) -> tuple[float, float, Q, float]:
    """``Compressor.run``'s no-chart polytropic default branch, in SI."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    molar_mass = capacity._molar_mass(mixture, inlet_z)
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)
    # The class reads `getDensity()`, the cubic's own ``M/(Z R T/P)`` - the capture's
    # ``inlet_density_kg_m3`` (14.62717) is that number, not the translated one.
    shift = mixture.volume_shift(inlet_z)
    _, _, volume_in = _feed_state(mixture, ideal_gas, t, p, inlet_z)
    density_in = molar_mass / (volume_in + shift)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )
    _, _, volume_out = _feed_state(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )
    density_out = molar_mass / (volume_out + shift)

    ratio = p_out / p
    exponent = math.log(ratio) / math.log(density_out / density_in)
    if (
        not math.isfinite(exponent)
        or exponent <= 0.0
        or abs(exponent - 1.0) < sys.float_info.epsilon
    ):
        raise InvalidInputError(
            "outlet_pressure",
            f"the isenthalpic volume exponent is {exponent} from densities {density_out} and "
            f"{density_in} kg/m3, which the polytropic path cannot divide by",
        )
    term = exponent / (exponent - 1.0) * efficiency
    term2 = p_out / density_out - p / density_in
    if term2 == 0.0 or not math.isfinite(term2):
        raise InvalidInputError(
            "outlet_pressure", f"the polytropic path's pressure-over-density term is {term2}"
        )
    correction = (h_isentropic - h_in) / (term * term2)
    dh = term * correction * (p / density_in) * (ratio ** (1.0 / term) - 1.0)
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    ).T
    # **The class reports a different exponent from the one it divided by**: after `run` it
    # recomputes it from the final outlet density (Compressor.java:1784).
    _, _, volume_actual = _feed_state(
        mixture, ideal_gas, outlet_t.to("K").magnitude, p_out, inlet_z
    )
    density_actual = molar_mass / (volume_actual + shift)
    density_ratio = density_actual / density_in
    exponent_polytropic = (
        math.log(ratio) / math.log(density_ratio)
        if ratio > 1.0 and density_ratio > 1.0
        else exponent
    )
    return h_in, h_out, outlet_t, exponent_polytropic


_DEPOSIT_DENSITIES = {
    "sulfur_s8": 2070.0,
    "salt_nacl": 2160.0,
    "scale_caco3": 2710.0,
    "scale_baso4": 4500.0,
    "scale_caso4": 2960.0,
    "iron_sulfide": 4840.0,
    "iron_oxide": 5240.0,
    "hydrocarbon_wax": 900.0,
    "particulate": 2500.0,
    "generic": 2000.0,
}


@dataclass(frozen=True, slots=True)
class _DepositState:
    """`CompressorDeposit`'s screening multipliers, in SI."""

    total_mass: float
    volume: float
    film_thickness: float
    blockage: float
    efficiency_multiplier: float
    head_multiplier: float


def _deposit_state(
    spec: dict[str, Any],
    deposit: list[dict[str, Q]] | None,
    wetted_area: float | None,
    passage_half_height: float | None,
) -> _DepositState | None:
    if not deposit:
        return None
    rows: list[dict[str, Any]] = input_records_to_si(spec, "deposit", deposit)
    masses: dict[str, float] = {}
    for row in rows:
        mechanism = str(row["mechanism"])
        mass = float(row["mass"])
        if mass != mass or mass <= 0.0:
            continue
        masses[mechanism] = masses.get(mechanism, 0.0) + mass
    volume = 0.0
    for mechanism, mass in masses.items():
        density = _DEPOSIT_DENSITIES.get(mechanism)
        if density is None:
            raise InvalidInputError(
                "deposit", f"`{mechanism}` is not one of the ten deposit mechanisms"
            )
        volume += mass / density
    film = volume / wetted_area if wetted_area and wetted_area > 0.0 else 0.0
    blockage = (
        min(film / passage_half_height, 0.95)
        if passage_half_height and passage_half_height > 0.0
        else 0.0
    )
    return _DepositState(
        total_mass=sum(masses.values()),
        volume=volume,
        film_thickness=film,
        blockage=blockage,
        efficiency_multiplier=max(1.0 - 1.6 * blockage, 0.02),
        head_multiplier=max((1.0 - blockage) ** 2, 0.05),
    )


def _power_states(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    power: float,
    n: float,
    polytropic_efficiency: float | None,
) -> tuple[float, Q, float]:
    """``setPower``'s fixed-shaft-work branch: ``h_out = h_in + power / n``.

    The step is per mole, so a non-positive flow has no step and is refused on `inlet_n`.
    """
    if not (math.isfinite(power) and power > 0.0):
        raise InvalidInputError("power", f"a fixed shaft power is positive, and {power} W is not")
    if p_out <= p:
        raise InvalidInputError(
            "outlet_pressure",
            f"the power branch compresses, and it is given {p_out} Pa against the inlet's {p} Pa",
        )
    if not math.isfinite(n) or not n > 0.0:
        raise InvalidInputError(
            "inlet_n",
            f"the power branch divides the shaft work by the flow, and {n} mol/s is not a "
            "positive flow",
        )
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    enthalpy_change = power / n
    h_out = h_in + enthalpy_change
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    ).T
    if polytropic_efficiency is not None:
        efficiency = polytropic_efficiency
    else:
        entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)
        isentropic = ps_flash_solve(
            mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
        )
        h_isentropic, _ = enthalpy_at(
            mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
        )
        efficiency = (h_isentropic - h_in) / enthalpy_change
    return h_out, outlet_t, efficiency


def _states(
    components: list[str],
    t: float,
    p: float,
    inlet_z: list[float],
    p_out: float,
    isentropic_efficiency: float,
) -> CompressorStates:
    """The machine's intermediates, in SI, in the order it computes them."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )

    dh_isentropic = h_isentropic - h_in
    dh_actual = dh_isentropic / isentropic_efficiency
    h_out = h_in + dh_actual
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
    )

    s_out, _ = entropy_at(mixture, ideal_gas, outlet_t.T.to("K").magnitude, p_out, inlet_z)

    return CompressorStates(
        h_in=h_in,
        s_in=entropy,
        t_isentropic=isentropic.T,
        h_isentropic=h_isentropic,
        dh_isentropic=dh_isentropic,
        dh_actual=dh_actual,
        h_out=h_out,
        outlet_t=outlet_t.T,
        s_out=s_out,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.compressor")
