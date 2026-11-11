"""``process.heater`` - the heater's kernel.

Spec: ``specs/models/process/heater.toml``

The Python twin of ``crates/azoth-process/src/models/heater.rs`` over
``crates/azoth-process/src/kernels/heater.rs``, written to mirror them.

# The duty is the state's, and it is in W against an enthalpy in J/mol

``Heater.run`` adds ``getDuty()`` (W) to ``system.getEnthalpy()`` (a total, J) directly, and
the two are the same dimension only while a process fluid's total moles equal its mol/s rate.
Here the conversion is written out - ``duty / inlet_n`` is the shift in J/mol - and the duty
*reported* is ``inlet_n * (h_out - h_in)``, which is ``run``'s own recomputation after the
flash rather than whatever the caller set.

# A temperature and a duty together are refused

``setOutletTemperature`` sets the temperature flag and clears the duty flag; ``setDuty`` does
the reverse. The class's answer to both is therefore *the order they were called in*, which
these inputs cannot express, so the pair is refused rather than resolved.

# With neither, the drop is isothermal

``run``'s else branch is ``T_in + dT`` with ``dT`` defaulting to zero, and ``dT`` is not a
palette parameter - so a heater that specifies nothing still moves the pressure.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import HeaterInstrument, HeaterMechanicalType
from azoth.core.result import HeaterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve

#: `HeatExchangerInstrumentDesign.HeatExchangerType`'s spellings.
INSTRUMENT_TYPES = ("electric_heater", "air_cooler", "shell_and_tube")
#: `ElectricalMotor.efficiencyClass`'s spellings.
IE_CLASSES = ("ie1", "ie2", "ie3", "ie4")


def heater(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
    max_outlet_temperature: Q | None = None,
    min_outlet_temperature: Q | None = None,
    instrument_design: bool | None = None,
    include_safety_instruments: bool | None = None,
    default_sil_level: float | None = None,
    heat_exchanger_type: str | None = None,
    electrical_design: bool | None = None,
    motor_sizing_margin: float | None = None,
    motor_poles: float | None = None,
    efficiency_class: str | None = None,
    frequency_hz: float | None = None,
    phases: float | None = None,
    number_of_fans: float | None = None,
    instrumentation_power: Q | None = None,
    cooling_water_pump_power: Q | None = None,
    cable_length: Q | None = None,
    use_vfd: bool | None = None,
    utility_specification: bool | None = None,
    utility_supply_temperature: Q | None = None,
    utility_return_temperature: Q | None = None,
    utility_approach_temperature: Q | None = None,
    utility_heat_capacity_rate: Q | None = None,
    utility_overall_heat_transfer_coefficient: Q | None = None,
    mechanical_design: bool | None = None,
    mechanical_type: str | None = None,
    mechanical_selection: str | None = None,
) -> HeaterResult:
    """Heat or cool a stream to a stated temperature, or by a stated duty.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_temperature: the temperature to reach, or ``None`` for the enthalpy route.
        duty: the heat to add, or ``None`` to state a temperature instead. Negative removes
            heat.
        pressure_drop: the pressure loss across the heater, or ``None`` for none.
        max_outlet_temperature: ``setMaxOutletTemperature``; absent is the class's unset
            sentinel, and supplying it initializes the capacity map.
        min_outlet_temperature: ``setMinOutletTemperature``; absent is 0, and only a positive
            limit is a constraint.

    Returns:
        The outlet's record, the duty it moved and the outlet-temperature capacity map.

    Raises:
        InvalidInputError: if a temperature and a duty are both given, if the pressure drop
            leaves a non-positive pressure, or if a duty is asked of a stream carrying no
            flow.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = heater(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     duty=q(5000.0, "W"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 3)
        420.458
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    t_out = (
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature)
    )
    q = None if duty is None else input_to_si(spec, "duty", duty)
    drop = None if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop)
    max_t = (
        None
        if max_outlet_temperature is None
        else input_to_si(spec, "max_outlet_temperature", max_outlet_temperature)
    )
    min_t = (
        None
        if min_outlet_temperature is None
        else input_to_si(spec, "min_outlet_temperature", min_outlet_temperature)
    )

    apply_checks(
        checks.on_input,
        {"inlet_t": t, "outlet_temperature": t_out}.get,
        warnings,
    )

    if t_out is not None and q is not None:
        raise InvalidInputError(
            "duty",
            "a heater holds a temperature specification or a duty specification and not "
            "both: `Heater.setOutletTemperature` clears the duty flag and `setDuty` clears "
            "the temperature flag, so which one the class runs is the order they were set "
            "in and not a function of these inputs",
        )

    instrument_requested = bool(instrument_design)
    for name, stated in (
        ("include_safety_instruments", include_safety_instruments is not None),
        ("default_sil_level", default_sil_level is not None),
    ):
        if stated and not instrument_requested:
            raise InvalidInputError(
                name, "an instrument-design input is stated without `instrument_design`"
            )
    sil = 2.0 if default_sil_level is None else float(default_sil_level)
    if not math.isfinite(sil) or sil < 0.0 or not sil.is_integer():
        raise InvalidInputError(
            "default_sil_level",
            f"the SIL level is {sil}, and the class's setter takes an integer",
        )
    kind = "electric_heater" if heat_exchanger_type is None else heat_exchanger_type
    if kind not in INSTRUMENT_TYPES:
        raise InvalidInputError(
            "heat_exchanger_type",
            f"{kind} is not one of electric_heater, air_cooler or shell_and_tube",
        )

    electrical_requested = bool(electrical_design)
    if heat_exchanger_type is not None and not instrument_requested and not electrical_requested:
        raise InvalidInputError(
            "heat_exchanger_type",
            "the type is shared by the two designs, and neither flag is set",
        )
    for name, stated in (
        ("motor_sizing_margin", motor_sizing_margin is not None),
        ("motor_poles", motor_poles is not None),
        ("efficiency_class", efficiency_class is not None),
        ("frequency_hz", frequency_hz is not None),
        ("phases", phases is not None),
        ("number_of_fans", number_of_fans is not None),
        ("instrumentation_power", instrumentation_power is not None),
        (
            "cooling_water_pump_power",
            cooling_water_pump_power is not None,
        ),
        ("cable_length", cable_length is not None),
        ("use_vfd", use_vfd is not None),
    ):
        if stated and not electrical_requested:
            raise InvalidInputError(
                name, "an electrical-design input is stated without `electrical_design`"
            )
    instrumentation_si = (
        None
        if instrumentation_power is None
        else input_to_si(spec, "instrumentation_power", instrumentation_power)
    )
    pump_si = (
        None
        if cooling_water_pump_power is None
        else input_to_si(spec, "cooling_water_pump_power", cooling_water_pump_power)
    )
    cable_length_si = (
        None if cable_length is None else input_to_si(spec, "cable_length", cable_length)
    )
    for name, value in (
        ("motor_sizing_margin", motor_sizing_margin),
        ("frequency_hz", frequency_hz),
        ("instrumentation_power", instrumentation_si),
        ("cooling_water_pump_power", pump_si),
        ("cable_length", cable_length_si),
    ):
        if value is not None and not math.isfinite(float(value)):
            raise InvalidInputError(name, f"{value} is not a finite number")
    efficiency = "ie3" if efficiency_class is None else efficiency_class
    if efficiency not in IE_CLASSES:
        raise InvalidInputError(
            "efficiency_class",
            f"{efficiency} is not one of ie1, ie2, ie3 or ie4",
        )
    poles = 4.0 if motor_poles is None else float(motor_poles)
    phase_count = 3.0 if phases is None else float(phases)
    fans = 2.0 if number_of_fans is None else float(number_of_fans)
    for name, value in (
        ("motor_poles", poles),
        ("phases", phase_count),
        ("number_of_fans", fans),
    ):
        if value < 1.0 or not value.is_integer():
            raise InvalidInputError(name, f"{value} is not an integer of at least one")
    instruments = (
        _instruments(kind, include_safety_instruments is not False, int(sil))
        if instrument_requested
        else ()
    )

    utility_requested = bool(utility_specification)
    utility_inputs = (
        ("utility_supply_temperature", utility_supply_temperature),
        ("utility_return_temperature", utility_return_temperature),
        ("utility_approach_temperature", utility_approach_temperature),
        ("utility_heat_capacity_rate", utility_heat_capacity_rate),
        (
            "utility_overall_heat_transfer_coefficient",
            utility_overall_heat_transfer_coefficient,
        ),
    )
    for name, given in utility_inputs:
        if given is not None and not utility_requested:
            raise InvalidInputError(
                name,
                "a utility-specification input is stated without `utility_specification`",
            )
    utility_supply_si = (
        None
        if utility_supply_temperature is None
        else input_to_si(spec, "utility_supply_temperature", utility_supply_temperature)
    )
    utility_return_si = (
        None
        if utility_return_temperature is None
        else input_to_si(spec, "utility_return_temperature", utility_return_temperature)
    )
    utility_approach_si = (
        None
        if utility_approach_temperature is None
        else input_to_si(spec, "utility_approach_temperature", utility_approach_temperature)
    )
    utility_capacity_si = (
        None
        if utility_heat_capacity_rate is None
        else input_to_si(spec, "utility_heat_capacity_rate", utility_heat_capacity_rate)
    )
    utility_u_si = (
        None
        if utility_overall_heat_transfer_coefficient is None
        else input_to_si(
            spec,
            "utility_overall_heat_transfer_coefficient",
            utility_overall_heat_transfer_coefficient,
        )
    )
    for name, value in (
        ("utility_supply_temperature", utility_supply_si),
        ("utility_return_temperature", utility_return_si),
        ("utility_approach_temperature", utility_approach_si),
        ("utility_heat_capacity_rate", utility_capacity_si),
        ("utility_overall_heat_transfer_coefficient", utility_u_si),
    ):
        if value is not None and not math.isfinite(float(value)):
            raise InvalidInputError(name, f"{value} is not a finite number")

    mechanical_requested = bool(mechanical_design)
    for name, stated in (
        ("mechanical_type", mechanical_type is not None),
        ("mechanical_selection", mechanical_selection is not None),
    ):
        if stated and not mechanical_requested:
            raise InvalidInputError(
                name,
                "a mechanical-design input is stated without `mechanical_design`",
            )
    manual = None
    if mechanical_type is not None:
        if mechanical_type not in MECHANICAL_NAMES:
            raise InvalidInputError(
                "mechanical_type",
                f"{mechanical_type} is not one of the five mechanical configurations",
            )
        manual = mechanical_type
    criterion = "min_area" if mechanical_selection is None else mechanical_selection
    if criterion not in MECHANICAL_CRITERIA:
        raise InvalidInputError(
            "mechanical_selection",
            f"{criterion} is not one of min_area, min_weight or min_pressure_drop",
        )

    states = _states(components, n, inlet_z, p, t, t_out, q, drop)
    capacity = _capacity(states.outlet_t.to("K").magnitude, max_t, min_t)
    electrical = (
        None
        if not electrical_requested
        else _electrical(
            kind,
            states.duty,
            1.10 if motor_sizing_margin is None else float(motor_sizing_margin),
            int(poles),
            efficiency,
            50.0 if frequency_hz is None else float(frequency_hz),
            int(phase_count),
            int(fans),
            1500.0 if instrumentation_si is None else float(instrumentation_si),
            0.0 if pump_si is None else float(pump_si),
            50.0 if cable_length_si is None else float(cable_length_si),
            bool(use_vfd),
        )
    )
    analog_input = sum(1 for row in instruments if row.io_type == "AI")
    analog_output = sum(1 for row in instruments if row.io_type == "AO")
    digital_input = sum(1 for row in instruments if row.io_type == "DI")
    digital_output = sum(1 for row in instruments if row.io_type == "DO")
    utility_state = {
        "supply": utility_supply_si,
        "return": utility_return_si,
        "approach": utility_approach_si,
        "rate": utility_capacity_si,
        "u": utility_u_si,
    }
    mechanical = (
        _mechanical(
            states.duty,
            t,
            states.outlet_t.to("K").magnitude,
            utility_state,
            manual,
            criterion,
        )
        if mechanical_requested
        else None
    )
    counts = (
        None
        if not instrument_requested
        else (
            float(len(instruments)),
            float(analog_input),
            float(analog_output),
            float(digital_input),
            float(digital_output),
            float(analog_input + analog_output + digital_input + digital_output),
            float(sum(1 for row in instruments if row.safety_related)),
        )
    )

    vfd = None if electrical is None else electrical.vfd
    return HeaterResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=states.outlet_t,
        outlet_h=from_si(states.h_out, "J/mol"),
        outlet_duty=from_si(states.duty, "W"),
        max_outlet_temperature_utilization=capacity.max_outlet_temperature_utilization,
        min_outlet_temperature_utilization=capacity.min_outlet_temperature_utilization,
        max_utilization=capacity.max_utilization,
        capacity_exceeded=capacity.capacity_exceeded,
        hard_limit_exceeded=capacity.hard_limit_exceeded,
        bottleneck_constraint=capacity.bottleneck_constraint,
        instrument_count=None if counts is None else counts[0],
        analog_input_count=None if counts is None else counts[1],
        analog_output_count=None if counts is None else counts[2],
        digital_input_count=None if counts is None else counts[3],
        digital_output_count=None if counts is None else counts[4],
        total_io_count=None if counts is None else counts[5],
        safety_instrument_count=None if counts is None else counts[6],
        instruments=instruments,
        shaft_power=(None if electrical is None else from_si(electrical.shaft_power_w, "W")),
        electrical_input_power=(
            None if electrical is None else from_si(electrical.electrical_input_w, "W")
        ),
        apparent_power=None if electrical is None else electrical.apparent_power_va,
        reactive_power=None if electrical is None else electrical.reactive_power_var,
        power_factor=None if electrical is None else electrical.power_factor,
        rated_voltage=None if electrical is None else electrical.rated_voltage_v,
        frequency_hz=None if electrical is None else electrical.frequency_hz,
        full_load_current=(None if electrical is None else electrical.full_load_current_a),
        starting_current=(None if electrical is None else electrical.starting_current_a),
        total_electrical_losses=(
            None if electrical is None else from_si(electrical.total_electrical_losses_w, "W")
        ),
        motor_rated_power=(
            None if electrical is None else from_si(electrical.motor_rated_power_w, "W")
        ),
        motor_efficiency=None if electrical is None else electrical.motor_efficiency,
        motor_power_factor=(None if electrical is None else electrical.motor_power_factor),
        motor_rated_speed=(None if electrical is None else electrical.motor_rated_speed_rpm),
        motor_locked_rotor_multiplier=(
            None if electrical is None else electrical.motor_locked_rotor_multiplier
        ),
        total_auxiliary_power=(
            None if electrical is None else from_si(electrical.total_auxiliary_w, "W")
        ),
        total_connected_power=(
            None if electrical is None else from_si(electrical.total_connected_w, "W")
        ),
        cable_length=(None if electrical is None else from_si(electrical.cable_length_m, "m")),
        cable_cross_section=(None if electrical is None else electrical.cable_cross_section_mm2),
        cable_ampacity=None if electrical is None else electrical.cable_ampacity_a,
        cable_voltage_drop=(None if electrical is None else electrical.cable_voltage_drop_percent),
        cable_short_circuit_withstand=(
            None if electrical is None else electrical.cable_short_circuit_withstand_ka
        ),
        switchgear_starter=None if electrical is None else electrical.switchgear_starter,
        switchgear_type=None if electrical is None else electrical.switchgear_type,
        switchgear_circuit_breaker=(
            None if electrical is None else electrical.switchgear_circuit_breaker_a
        ),
        switchgear_fuse=None if electrical is None else electrical.switchgear_fuse_a,
        switchgear_rated_current=(
            None if electrical is None else electrical.switchgear_rated_current_a
        ),
        switchgear_short_circuit=(
            None if electrical is None else electrical.switchgear_short_circuit_ka
        ),
        switchgear_weight=(
            None if electrical is None else from_si(electrical.switchgear_weight_kg, "kg")
        ),
        vfd_rated_power=(None if vfd is None else from_si(vfd.rated_power_kw * 1000.0, "W")),
        vfd_rated_current=None if vfd is None else vfd.rated_current_a,
        vfd_input_voltage=None if vfd is None else vfd.input_voltage_v,
        vfd_topology=None if vfd is None else vfd.topology,
        vfd_pulse_configuration=None if vfd is None else vfd.pulse_configuration,
        vfd_efficiency=None if vfd is None else vfd.efficiency_percent / 100.0,
        vfd_input_power_factor=None if vfd is None else vfd.input_power_factor,
        vfd_thd_current=None if vfd is None else vfd.thd_current_percent,
        vfd_has_active_rectifier=None if vfd is None else vfd.has_active_rectifier,
        vfd_requires_input_filter=None if vfd is None else vfd.requires_input_filter,
        vfd_heat_dissipation=(
            None if vfd is None else from_si(vfd.heat_dissipation_kw * 1000.0, "W")
        ),
        vfd_cooling_method=None if vfd is None else vfd.cooling_method,
        vfd_weight=None if vfd is None else from_si(vfd.weight_kg, "kg"),
        vfd_enclosure_rating=None if vfd is None else vfd.enclosure_rating,
        has_supply_temperature=(None if not utility_requested else utility_supply_si is not None),
        has_return_temperature=(None if not utility_requested else utility_return_si is not None),
        has_approach_temperature=(
            None if not utility_requested else utility_approach_si is not None
        ),
        has_heat_capacity_rate=(
            None
            if not utility_requested
            else utility_capacity_si is not None and utility_capacity_si > 0.0
        ),
        has_overall_heat_transfer_coefficient=(
            None if not utility_requested else utility_u_si is not None and utility_u_si > 0.0
        ),
        mechanical_selected_type=(
            None if mechanical is None else mechanical.selected.mechanical_type
        ),
        mechanical_selected_area=(
            None if mechanical is None else mechanical.selected.required_area
        ),
        log_mean_temperature_difference=(
            None if mechanical is None else from_si(mechanical.lmtd, "K")
        ),
        mechanical_approach_temperature=(
            None if mechanical is None else from_si(mechanical.approach, "K")
        ),
        used_overall_heat_transfer_coefficient=(
            None if mechanical is None else from_si(mechanical.used_u, "W/(m**2*K)")
        ),
        calculated_ua=(None if mechanical is None else from_si(mechanical.ua, "W/K")),
        mechanical_types=() if mechanical is None else mechanical.types,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class Vfd:
    """`VariableFrequencyDrive.sizeVFD`'s table and ratings."""

    rated_power_kw: float
    rated_current_a: float
    input_voltage_v: float
    topology: str
    pulse_configuration: str
    efficiency_percent: float
    input_power_factor: float
    thd_current_percent: float
    has_active_rectifier: bool
    requires_input_filter: bool
    heat_dissipation_kw: float
    cooling_method: str
    weight_kg: float
    enclosure_rating: str


@dataclass(frozen=True, slots=True)
class Electrical:
    """`HeatExchangerElectricalDesign.calcDesign`'s numbers, in SI."""

    shaft_power_w: float
    electrical_input_w: float
    apparent_power_va: float
    reactive_power_var: float
    power_factor: float
    rated_voltage_v: float
    frequency_hz: float
    full_load_current_a: float
    starting_current_a: float
    total_electrical_losses_w: float
    motor_rated_power_w: float
    motor_efficiency: float
    motor_power_factor: float
    motor_rated_speed_rpm: float
    motor_locked_rotor_multiplier: float
    total_auxiliary_w: float
    total_connected_w: float
    cable_length_m: float
    cable_cross_section_mm2: float
    cable_ampacity_a: float
    cable_voltage_drop_percent: float
    cable_short_circuit_withstand_ka: float
    switchgear_starter: str
    switchgear_type: str
    switchgear_circuit_breaker_a: float
    switchgear_fuse_a: float
    switchgear_rated_current_a: float
    switchgear_short_circuit_ka: float
    switchgear_weight_kg: float
    vfd: Vfd | None


#: `ElectricalMotor.IEC_STANDARD_POWERS_KW`.
IEC_POWERS_KW = (
    0.37,
    0.55,
    0.75,
    1.1,
    1.5,
    2.2,
    3.0,
    4.0,
    5.5,
    7.5,
    11.0,
    15.0,
    18.5,
    22.0,
    30.0,
    37.0,
    45.0,
    55.0,
    75.0,
    90.0,
    110.0,
    132.0,
    160.0,
    200.0,
    250.0,
    315.0,
    355.0,
    400.0,
    450.0,
    500.0,
    560.0,
    630.0,
    710.0,
    800.0,
    900.0,
    1000.0,
    1120.0,
    1250.0,
    1400.0,
    1600.0,
    1800.0,
    2000.0,
    2240.0,
    2500.0,
    2800.0,
    3150.0,
    3550.0,
    4000.0,
    4500.0,
    5000.0,
    5600.0,
    6300.0,
    7100.0,
    8000.0,
    9000.0,
    10000.0,
)

#: The efficiency-class offsets, in points.
IE_OFFSETS = {"ie1": -3.0, "ie2": -1.5, "ie3": 0.0, "ie4": 1.0}

#: `ElectricalCable.STANDARD_CROSS_SECTIONS`, mm2.
CABLE_CROSS_SECTIONS_MM2 = (
    1.5,
    2.5,
    4.0,
    6.0,
    10.0,
    16.0,
    25.0,
    35.0,
    50.0,
    70.0,
    95.0,
    120.0,
    150.0,
    185.0,
    240.0,
    300.0,
    400.0,
    500.0,
    630.0,
)
#: `ElectricalCable.BASE_AMPACITY_XLPE_TRAY`, A.
CABLE_BASE_AMPACITY_A = (
    18.0,
    25.0,
    33.0,
    43.0,
    60.0,
    80.0,
    106.0,
    131.0,
    158.0,
    200.0,
    241.0,
    278.0,
    315.0,
    360.0,
    420.0,
    480.0,
    560.0,
    640.0,
    740.0,
)
#: `Switchgear.STANDARD_RATINGS_A`.
SWITCHGEAR_RATINGS_A = (
    100.0,
    160.0,
    250.0,
    400.0,
    630.0,
    800.0,
    1000.0,
    1250.0,
    1600.0,
    2000.0,
    2500.0,
    3150.0,
    4000.0,
    5000.0,
)
#: `Switchgear.selectStandardFuseRating`'s ratings.
FUSE_RATINGS_A = (
    16.0,
    20.0,
    25.0,
    32.0,
    40.0,
    50.0,
    63.0,
    80.0,
    100.0,
    125.0,
    160.0,
    200.0,
    250.0,
    315.0,
    400.0,
    500.0,
    630.0,
    800.0,
    1000.0,
    1250.0,
)


def _next_cable_size(min_mm2: float) -> float:
    return next(
        (standard for standard in CABLE_CROSS_SECTIONS_MM2 if standard >= min_mm2),
        CABLE_CROSS_SECTIONS_MM2[-1],
    )


def _voltage_drop_percent(
    current_a: float, voltage_v: float, length_m: float, section: float
) -> float:
    if section <= 0.0 or voltage_v <= 0.0 or length_m <= 0.0:
        return 0.0
    resistance = 0.0175 / section
    reactance = 0.00008
    cos_phi = 0.85
    sin_phi = math.sqrt(1.0 - cos_phi * cos_phi)
    drop = math.sqrt(3.0) * current_a * length_m * (resistance * cos_phi + reactance * sin_phi)
    return drop / voltage_v * 100.0


def _size_cable(
    current_a: float, voltage_v: float, length_m: float
) -> tuple[float, float, float, float]:
    """`ElectricalCable.sizeCable` on the heater path's Tray/40 C/copper defaults."""
    ambient = min(1.0, max(0.5, math.sqrt((90.0 - 40.0) / (90.0 - 30.0))))
    derating = ambient * 0.85
    required = current_a / derating if derating > 0.0 else current_a
    section = CABLE_CROSS_SECTIONS_MM2[0]
    base_ampacity = CABLE_BASE_AMPACITY_A[0]
    for index, standard in enumerate(CABLE_CROSS_SECTIONS_MM2):
        if CABLE_BASE_AMPACITY_A[index] >= required:
            section = standard
            base_ampacity = CABLE_BASE_AMPACITY_A[index]
            break
        if index == len(CABLE_CROSS_SECTIONS_MM2) - 1:
            section = standard
            base_ampacity = CABLE_BASE_AMPACITY_A[index]
    ampacity = base_ampacity * derating
    drop = _voltage_drop_percent(current_a, voltage_v, length_m, section)
    while drop > 5.0 and section < CABLE_CROSS_SECTIONS_MM2[-1]:
        section = _next_cable_size(section + 1.0)
        drop = _voltage_drop_percent(current_a, voltage_v, length_m, section)
    withstand = 143.0 * section / 1000.0
    return section, ampacity, drop, withstand


def _size_vfd(rated_power_kw: float, motor_rated_current_a: float, input_voltage_v: float) -> Vfd:
    """`VariableFrequencyDrive.sizeVFD`, on the motor's rating and rated current."""
    rated_current_a = motor_rated_current_a * 1.10
    if input_voltage_v > 3300.0:
        if rated_power_kw > 2000.0:
            topology, pulse, active, thd = ("Multi-level", "AFE", True, 3.0)
        else:
            topology, pulse, active, thd = ("3-level", "12-pulse", False, 10.0)
    elif rated_power_kw > 250.0:
        topology, pulse, active, thd = ("2-level", "AFE", True, 5.0)
    else:
        topology, pulse, active, thd = ("2-level", "6-pulse", False, 35.0)
    # The 3-level arm is tested before the Multi-level one, and a Multi-level drive always has
    # an active rectifier, so the class's ``"Multi-level"`` efficiency arm is unreachable.
    if active:
        efficiency_percent, input_power_factor = (96.5, 0.98)
    elif topology == "3-level":
        efficiency_percent, input_power_factor = (97.0, 0.95)
    elif topology == "Multi-level":
        efficiency_percent, input_power_factor = (97.5, 0.99)
    else:
        efficiency_percent, input_power_factor = (97.5, 0.95)
    heat_dissipation_kw = rated_power_kw * (1.0 - efficiency_percent / 100.0)
    cooling_method = "Water" if rated_power_kw > 500.0 or input_voltage_v > 3300.0 else "Air"
    if input_voltage_v > 3300.0:
        weight_kg = rated_power_kw * 1.5 + 500.0
    elif input_voltage_v > 690.0:
        weight_kg = rated_power_kw + 200.0
    else:
        weight_kg = rated_power_kw * 0.6 + 50.0
    if input_voltage_v > 1000.0:
        enclosure_rating = "IP42"
    elif cooling_method == "Water":
        enclosure_rating = "IP54"
    else:
        enclosure_rating = "IP21"
    return Vfd(
        rated_power_kw=rated_power_kw,
        rated_current_a=rated_current_a,
        input_voltage_v=input_voltage_v,
        topology=topology,
        pulse_configuration=pulse,
        efficiency_percent=efficiency_percent,
        input_power_factor=input_power_factor,
        thd_current_percent=thd,
        has_active_rectifier=active,
        requires_input_filter=not active and pulse == "6-pulse" and rated_power_kw > 100.0,
        heat_dissipation_kw=heat_dissipation_kw,
        cooling_method=cooling_method,
        weight_kg=weight_kg,
        enclosure_rating=enclosure_rating,
    )


def _size_switchgear(
    current_a: float, rated_power_kw: float, voltage_v: float, use_vfd: bool
) -> tuple[str, str, float, float, float, float, float]:
    """`Switchgear.sizeSwitchgear`."""
    if use_vfd:
        starter = "VFD"
    elif rated_power_kw > 200.0:
        starter = "Soft Starter"
    elif rated_power_kw > 11.0:
        starter = "Star-Delta"
    else:
        starter = "DOL"
    switchgear_type = "MV Switchgear" if voltage_v > 1000.0 else "MCC"
    required = current_a * 1.25
    breaker = next(
        (rating for rating in SWITCHGEAR_RATINGS_A if rating >= required),
        math.ceil(required / 100.0) * 100.0,
    )
    target = current_a * 1.6
    fuse = next(
        (rating for rating in FUSE_RATINGS_A if rating >= target),
        FUSE_RATINGS_A[-1],
    )
    if voltage_v > 6000.0:
        short_circuit = 40.0
    elif voltage_v > 1000.0:
        short_circuit = 31.5
    else:
        short_circuit = 25.0
    weight = breaker * 0.5 + 200.0 if voltage_v > 1000.0 else breaker * 0.15 + 50.0
    return starter, switchgear_type, breaker, fuse, breaker, short_circuit, weight


def _electrical(
    kind: str,
    process_duty_w: float,
    margin: float,
    poles: int,
    efficiency_class: str,
    frequency: float,
    phases: int,
    fans: int,
    instrumentation_w: float,
    pump_w: float,
    cable_length_m: float,
    use_vfd: bool,
) -> Electrical:
    """`HeatExchangerElectricalDesign.calcDesign`, in SI."""
    duty_kw = abs(process_duty_w) / 1000.0
    if kind == "electric_heater":
        shaft_kw = duty_kw
    elif kind == "air_cooler":
        shaft_kw = max(2.0 * fans, duty_kw * 0.01)
    else:
        shaft_kw = 0.0
    if kind == "electric_heater":
        voltage = 6600.0 if shaft_kw > 200.0 else (690.0 if shaft_kw > 75.0 else 400.0)
    elif kind == "air_cooler":
        voltage = 690.0 if shaft_kw > 200.0 else 400.0
    else:
        voltage = 400.0

    if shaft_kw > 0.0:
        required = shaft_kw * margin
        rated = next(
            (size for size in IEC_POWERS_KW if size >= required * (1.0 - 1.0e-9)),
            None,
        )
        if rated is None:
            rated = math.ceil(required / 100.0) * 100.0
        # The motor is sized before the design's frequency is copied to it.
        synchronous = 120.0 * 50.0 / max(1, poles)
        if rated <= 7.5:
            slip = 0.05
        elif rated <= 75.0:
            slip = 0.03
        elif rated <= 375.0:
            slip = 0.02
        else:
            slip = 0.015
        if rated <= 1.1:
            base = 82.5
        elif rated <= 4.0:
            base = 87.7
        elif rated <= 11.0:
            base = 91.0
        elif rated <= 37.0:
            base = 93.3
        elif rated <= 110.0:
            base = 95.0
        elif rated <= 375.0:
            base = 96.0
        else:
            base = 96.5
        pole_offset = {2: -0.5, 6: -1.0, 8: -1.5}.get(poles, 0.0)
        efficiency = min(99.0, base + IE_OFFSETS[efficiency_class] + pole_offset)
        if rated <= 4.0:
            factor = 0.78
        elif rated <= 15.0:
            factor = 0.83
        elif rated <= 55.0:
            factor = 0.86
        elif rated <= 200.0:
            factor = 0.88
        elif rated <= 1000.0:
            factor = 0.90
        else:
            factor = 0.91
        if rated <= 7.5:
            locked = 6.0
        elif rated <= 75.0:
            locked = 6.5
        elif rated <= 500.0:
            locked = 7.0
        else:
            locked = 6.5
        speed = synchronous * (1.0 - slip)
    else:
        rated = 0.0
        efficiency = 95.0
        factor = 0.86
        speed = 0.0
        locked = 7.0

    auxiliary = instrumentation_w + (pump_w if kind == "shell_and_tube" else 0.0)
    if kind == "shell_and_tube":
        input_w = auxiliary
        base_factor = 0.85
        apparent = input_w / base_factor if base_factor > 0.0 else 0.0
        full_load = _current(input_w / 1000.0, voltage, base_factor, phases)
        return Electrical(
            shaft_power_w=0.0,
            electrical_input_w=input_w,
            apparent_power_va=apparent,
            reactive_power_var=apparent * math.sqrt(1.0 - base_factor * base_factor),
            power_factor=base_factor,
            rated_voltage_v=voltage,
            frequency_hz=frequency,
            full_load_current_a=full_load,
            starting_current_a=full_load * locked,
            total_electrical_losses_w=input_w,
            motor_rated_power_w=0.0,
            motor_efficiency=efficiency / 100.0,
            motor_power_factor=factor,
            motor_rated_speed_rpm=speed,
            motor_locked_rotor_multiplier=locked,
            total_auxiliary_w=auxiliary,
            total_connected_w=input_w + auxiliary,
            cable_length_m=cable_length_m,
            cable_cross_section_mm2=0.0,
            cable_ampacity_a=0.0,
            cable_voltage_drop_percent=0.0,
            cable_short_circuit_withstand_ka=0.0,
            switchgear_starter="DOL",
            switchgear_type="MCC",
            switchgear_circuit_breaker_a=0.0,
            switchgear_fuse_a=0.0,
            switchgear_rated_current_a=0.0,
            switchgear_short_circuit_ka=25.0,
            switchgear_weight_kg=0.0,
            vfd=None,
        )
    efficiency_fraction = efficiency / 100.0 if efficiency > 0.0 else 0.0
    motor_input_kw = shaft_kw / efficiency_fraction if efficiency_fraction > 0.0 else shaft_kw
    # The motor's rated current is computed inside ``sizeMotor``, before ``calcDesign`` copies
    # the design's tap to the motor, so it is always at the motor's own 400 V default.
    motor_current_denominator = math.sqrt(3.0) * 400.0 * efficiency_fraction * factor
    motor_rated_current = (
        rated * 1000.0 / motor_current_denominator if motor_current_denominator > 0.0 else 0.0
    )
    vfd = _size_vfd(rated, motor_rated_current, voltage) if use_vfd else None
    if vfd is None:
        input_kw = motor_input_kw
        effective_factor = factor
    else:
        input_kw = motor_input_kw / (vfd.efficiency_percent / 100.0)
        effective_factor = vfd.input_power_factor
    apparent = input_kw / effective_factor if effective_factor > 0.0 else 0.0
    full_load = _current(input_kw, voltage, effective_factor, phases)
    cable_section, cable_ampacity, cable_drop, cable_withstand = _size_cable(
        full_load, voltage, cable_length_m
    )
    starter, switchgear_type, breaker, fuse, rated_current, short_circuit, weight = (
        _size_switchgear(full_load, rated, voltage, use_vfd)
    )
    return Electrical(
        shaft_power_w=shaft_kw * 1000.0,
        electrical_input_w=input_kw * 1000.0,
        apparent_power_va=apparent * 1000.0,
        reactive_power_var=apparent * 1000.0 * math.sqrt(1.0 - effective_factor * effective_factor),
        power_factor=effective_factor,
        rated_voltage_v=voltage,
        frequency_hz=frequency,
        full_load_current_a=full_load,
        starting_current_a=full_load if use_vfd else full_load * locked,
        total_electrical_losses_w=(input_kw - shaft_kw) * 1000.0,
        motor_rated_power_w=rated * 1000.0,
        motor_efficiency=efficiency_fraction,
        motor_power_factor=factor,
        motor_rated_speed_rpm=speed,
        motor_locked_rotor_multiplier=locked,
        total_auxiliary_w=auxiliary,
        total_connected_w=input_kw * 1000.0 + auxiliary,
        cable_length_m=cable_length_m,
        cable_cross_section_mm2=cable_section,
        cable_ampacity_a=cable_ampacity,
        cable_voltage_drop_percent=cable_drop,
        cable_short_circuit_withstand_ka=cable_withstand,
        switchgear_starter=starter,
        switchgear_type=switchgear_type,
        switchgear_circuit_breaker_a=breaker,
        switchgear_fuse_a=fuse,
        switchgear_rated_current_a=rated_current,
        switchgear_short_circuit_ka=short_circuit,
        switchgear_weight_kg=weight,
        vfd=vfd,
    )


def _current(input_kw: float, voltage: float, factor: float, phases: int) -> float:
    if voltage <= 0.0 or factor <= 0.0:
        return 0.0
    if phases == 3:
        return input_kw * 1000.0 / (math.sqrt(3.0) * voltage * factor)
    return input_kw * 1000.0 / (voltage * factor)


@dataclass(frozen=True, slots=True)
class HeaterStates:
    """The heater's intermediates, in the order it computes them, in SI."""

    #: The inlet's molar enthalpy, J/mol.
    h_in: float
    #: The outlet pressure: the inlet's less the drop.
    p_out: float
    #: The outlet molar enthalpy, J/mol.
    h_out: float
    #: The outlet temperature.
    outlet_t: Q
    #: The duty moved, W.
    duty: float


def _states(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: float,
    inlet_t: float,
    outlet_temperature: float | None,
    duty: float | None,
    pressure_drop: float | None,
) -> HeaterStates:
    """The heater's intermediates, in SI, in the order it computes them.

    **One arithmetic, two consumers**: :func:`heater` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    h_in, _ = enthalpy_at(mixture, ideal_gas, inlet_t, inlet_p, inlet_z)

    p_out = inlet_p - (0.0 if pressure_drop is None else pressure_drop)
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {0.0 if pressure_drop is None else pressure_drop} Pa takes "
            f"a {inlet_p} Pa inlet to {p_out} Pa, which is not a pressure",
        )

    if outlet_temperature is not None:
        # The state at the stated temperature. `Heater.run` sets the temperature and runs a
        # `TPflash`, which is what `enthalpy_at` computes the enthalpy of.
        h_out, _ = enthalpy_at(mixture, ideal_gas, outlet_temperature, p_out, inlet_z)
        outlet_t = from_si(outlet_temperature, "K")
    elif duty is not None:
        if inlet_n == 0.0:
            raise InvalidInputError(
                "duty",
                "a duty is moved per mole of flow, and this stream carries none",
            )
        h_out = h_in + duty / inlet_n
        outlet_t = ph_flash_solve(
            mixture, ideal_gas, from_si(p_out, "Pa"), from_si(h_out, "J/mol"), inlet_z
        ).T
    else:
        # `run`'s else branch: `T_in + dT`, with `dT` zero. **An isothermal drop and not a
        # throttling**, which is the distinction an earlier version of this file got wrong
        # and the cross-implementation test caught: the outlet holds the inlet's
        # temperature, and its enthalpy is the *state's* at the outlet pressure, which for a
        # real fluid is not the inlet's. Measured on the captured fluid, two bar of drop at
        # 320 K moves it `49.6` J/mol - and NeqSim reports that movement as the duty.
        h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, p_out, inlet_z)
        outlet_t = from_si(inlet_t, "K")

    return HeaterStates(
        h_in=h_in,
        p_out=p_out,
        h_out=h_out,
        outlet_t=outlet_t,
        duty=inlet_n * (h_out - h_in),
    )


@dataclass(frozen=True, slots=True)
class HeaterCapacity:
    """The outlet-temperature constraints and their summaries, in SI."""

    max_outlet_temperature_utilization: float | None
    min_outlet_temperature_utilization: float | None
    max_utilization: float
    capacity_exceeded: bool
    hard_limit_exceeded: bool
    bottleneck_constraint: str | None


def _capacity(
    outlet_temperature: float,
    max_outlet_temperature: float | None,
    min_outlet_temperature: float | None,
) -> HeaterCapacity:
    """`Heater`'s outlet-temperature capacity map, evaluated on the outlet temperature.

    **The maximum's hard-limit check is a no-op in the class.** `isHardLimitExceeded` reads
    `maxValue`, while `setMaxOutletTemperature` sets `designValue`; only the minimum can trip
    it. A maximum of `f64::MAX` is the unset sentinel, and a minimum of zero or less is not a
    constraint at all.
    """
    max_limit = (
        None
        if max_outlet_temperature is None or max_outlet_temperature == 1.7976931348623157e308
        else max_outlet_temperature
    )
    min_limit = (
        None
        if min_outlet_temperature is None or not min_outlet_temperature > 0.0
        else min_outlet_temperature
    )

    max_utilization = None
    if max_limit is not None:
        max_utilization = 0.0 if max_limit <= 0.0 else min(outlet_temperature / max_limit, 9.99)
    min_utilization = None
    if min_limit is not None:
        min_utilization = (
            9.99 if outlet_temperature <= 0.0 else min(min_limit / outlet_temperature, 9.99)
        )

    # The map is insertion-ordered: the maximum first, the minimum second.
    max_utilization_value = 0.0
    bottleneck_constraint = None
    if max_utilization is not None and max_utilization > max_utilization_value:
        max_utilization_value = max_utilization
        bottleneck_constraint = "maxOutletTemperature"
    if min_utilization is not None and min_utilization > max_utilization_value:
        max_utilization_value = min_utilization
        bottleneck_constraint = "minOutletTemperature"

    return HeaterCapacity(
        max_outlet_temperature_utilization=max_utilization,
        min_outlet_temperature_utilization=min_utilization,
        max_utilization=max_utilization_value,
        capacity_exceeded=(max_utilization is not None and max_utilization > 1.0)
        or (min_utilization is not None and min_utilization > 1.0),
        hard_limit_exceeded=(min_limit is not None and outlet_temperature < min_limit),
        bottleneck_constraint=bottleneck_constraint,
    )


#: The five mechanical configurations and their (typical U, allowable approach).
MECHANICAL_TYPES = (
    ("shell_and_tube", 580.0, 5.0),
    ("plate_and_frame", 1150.0, 2.0),
    ("air_cooler", 90.0, 12.0),
    ("double_pipe", 350.0, 8.0),
    ("plate_fin", 600.0, 1.5),
)
#: The five configurations' spellings.
MECHANICAL_NAMES = tuple(kind for kind, _, _ in MECHANICAL_TYPES)
#: `HeatExchangerMechanicalDesign.SelectionCriterion`'s spellings.
MECHANICAL_CRITERIA = ("min_area", "min_weight", "min_pressure_drop")


def _lmtd(delta_t1: float, delta_t2: float) -> float:
    safe1 = max(delta_t1, 1.0e-6)
    safe2 = max(delta_t2, 1.0e-6)
    if abs(safe1 - safe2) < 1.0e-6:
        return (safe1 + safe2) / 2.0
    return abs((safe1 - safe2) / math.log(safe1 / safe2))


def _geometry(name: str, area: float, ua: float, approach: float) -> HeaterMechanicalType:
    """One `HeatExchangerType.createSizingResult`, in SI."""
    length = width = height = inner = outer = thickness = 0.0
    tubes = passes = fin_area = drop = weight = 0.0
    if name == "shell_and_tube":
        length = max(6.0, min(area / (math.pi * 0.019 * 40.0), 12.0))
        tubes = max(math.ceil(area / (math.pi * 0.019 * length)), 2)
        passes = 4.0 if tubes > 100 else 2.0
        bundle = math.sqrt(tubes) * (0.019 * 1.25)
        inner = max(bundle, 0.3)
        outer = inner + 2.0 * 0.012
        thickness = 0.012
        drop = 0.05 * passes * length
        weight = math.pi * inner * length * 0.012 * 7850.0
        width = height = outer + 0.8
    elif name == "plate_and_frame":
        passes = max(math.ceil(area / 2.2), 5)
        tubes = passes * 2.0 + 1.0
        length = 0.6 + tubes * 0.003
        inner = math.sqrt(area / math.pi)
        outer = inner + 0.016
        thickness = 0.008
        drop = 0.02 * passes
        weight = tubes * 12.0
        width = 1.0
        height = 1.8
    elif name == "air_cooler":
        plan = max(area / 2.5, 5.0)
        width = math.sqrt(plan)
        length = plan / width
        inner = math.sqrt(4.0 * plan / math.pi)
        outer = inner + 0.012
        thickness = 0.006
        fin_area = area * 2.5
        drop = 0.01 * area
        weight = plan * 45.0
        height = 3.2
    elif name == "double_pipe":
        length = max(area / (math.pi * 0.05), 2.5)
        inner = 0.05
        outer = 0.11
        thickness = 0.005
        tubes = 2.0
        passes = 2.0
        drop = 0.08 * length
        weight = math.pi * outer * length * 0.005 * 7850.0
        width = height = outer + 0.6
    else:  # plate_fin
        volume = area / 1000.0
        wh = (volume / 6.0) ** (1.0 / 3.0)
        length = 6.0 * wh
        inner = wh
        outer = wh + 2.0 * 0.0016
        thickness = 0.0016
        fin_area = area
        drop = 0.03 * length
        weight = volume * 0.30 * 2730.0
        width = wh + 0.8
        height = wh + 0.5
    return HeaterMechanicalType(
        mechanical_type=name,
        required_area=from_si(area, "m**2"),
        required_ua=from_si(ua, "W/K"),
        overall_heat_transfer_coefficient=from_si(
            next(typical for kind, typical, _ in MECHANICAL_TYPES if kind == name),
            "W/(m**2*K)",
        ),
        approach_temperature=from_si(approach, "K"),
        estimated_length=from_si(length, "m"),
        inner_diameter=from_si(inner, "m"),
        outer_diameter=from_si(outer, "m"),
        wall_thickness=from_si(thickness, "m"),
        tube_count=float(tubes),
        tube_passes=float(passes),
        fin_surface_area=from_si(fin_area, "m**2"),
        estimated_pressure_drop=drop,
        estimated_weight=from_si(weight, "kg"),
        module_length=from_si(
            length
            + (1.5 if name == "shell_and_tube" else 0.0)
            + (
                0.0
                if name not in ("double_pipe", "plate_fin")
                else (0.5 if name == "double_pipe" else 2.0)
            ),
            "m",
        ),
        module_width=from_si(width, "m"),
        module_height=from_si(height, "m"),
    )


@dataclass(frozen=True, slots=True)
class _Mechanical:
    lmtd: float
    approach: float
    used_u: float
    ua: float
    types: tuple[HeaterMechanicalType, ...]
    selected: HeaterMechanicalType


def _mechanical(
    duty_w: float,
    process_in: float,
    process_out: float,
    utility: dict[str, float | None],
    manual: str | None,
    criterion: str,
) -> _Mechanical:
    """`HeatExchangerMechanicalDesign.calcDesign`'s heater branch."""
    duty = abs(duty_w)
    used_u = 500.0
    if utility["u"] is not None and utility["u"] > 0.0:
        used_u = utility["u"]
    heating = process_out >= process_in
    supply = utility["supply"]
    ret = utility["return"]
    if supply is None and utility["approach"] is not None:
        supply = process_out + utility["approach"] if heating else process_out - utility["approach"]
    if (
        ret is None
        and utility["rate"] is not None
        and utility["rate"] > 0.0
        and duty > 0.0
        and supply is not None
    ):
        delta = duty / utility["rate"]
        ret = supply - delta if heating else supply + delta
    if ret is None and supply is not None:
        ret = supply
    if supply is None or ret is None:
        lmtd = None
        approach = None
        delta = abs(process_out - process_in)
        if delta > 0.0:
            lmtd = delta
            approach = delta
    else:
        temperatures = (
            (supply, ret, process_in, process_out)
            if heating
            else (process_in, process_out, supply, ret)
        )
        delta_t1 = abs(temperatures[0] - temperatures[3])
        delta_t2 = abs(temperatures[1] - temperatures[2])
        lmtd = _lmtd(delta_t1, delta_t2)
        approach = min(delta_t1, delta_t2)
    ua = None
    if duty > 0.0 and lmtd is not None and lmtd > 0.0:
        ua = duty / lmtd
    if approach is None or approach < 0.0:
        approach = 0.0
    if lmtd is None or lmtd <= 0.0:
        if approach > 0.0:
            lmtd = approach
        elif duty > 0.0 and used_u > 0.0:
            lmtd = duty / used_u
        else:
            lmtd = 1.0
    if (ua is None or not ua > 0.0) and duty > 0.0 and lmtd > 0.0:
        ua = duty / lmtd
    if ua is None or ua <= 0.0:
        ua = used_u
    rows = []
    for name, typical, allowable in MECHANICAL_TYPES:
        type_approach = max(approach, allowable) if approach > 0.0 else allowable
        ua_for_type = ua
        if not ua_for_type > 0.0 and duty > 0.0:
            ua_for_type = duty / max(type_approach, 1.0)
        if ua_for_type is None or ua_for_type <= 0.0:
            ua_for_type = typical
        area = max(ua_for_type / max(used_u, typical), 0.1)
        rows.append(_geometry(name, area, ua_for_type, type_approach))
    if manual is not None:
        selected = next(row for row in rows if row.mechanical_type == manual)
    else:

        def metric(row: HeaterMechanicalType) -> float:
            if criterion == "min_weight":
                return row.estimated_weight.to("kg").magnitude
            if criterion == "min_pressure_drop":
                return row.estimated_pressure_drop
            return row.required_area.to("m**2").magnitude

        selected = min(rows, key=metric)
    return _Mechanical(
        lmtd=lmtd,
        approach=approach,
        used_u=used_u,
        ua=ua,
        types=tuple(rows),
        selected=selected,
    )


def _instruments(kind: str, include_safety: bool, sil: int) -> tuple[HeaterInstrument, ...]:
    """`HeatExchangerInstrumentDesign.calcDesign`'s list, in the class's order."""
    rows: list[tuple[str, str, float, float, str | None, str, int, bool]] = [
        ("TT", "Process Inlet Temperature", -50.0, 400.0, "degC", "AI", 0, False),
        ("TT", "Process Outlet Temperature", -50.0, 400.0, "degC", "AI", 0, False),
        ("PT", "Process Inlet Pressure", 0.0, 200.0, "bara", "AI", 0, False),
        ("PT", "Process Outlet Pressure", 0.0, 200.0, "bara", "AI", 0, False),
        ("PDT", "Process Side dP", 0.0, 10.0, "bar", "AI", 0, False),
    ]
    if include_safety:
        rows.append(("TSHH", "Overtemperature Trip", 0.0, 0.0, None, "DI", sil, True))
    if kind == "electric_heater":
        rows.append(("TT", "Heating Element Temperature", 0.0, 600.0, "degC", "AI", 0, False))
        if include_safety:
            rows.append(("TSHH", "Element Overtemperature Trip", 0.0, 0.0, None, "DI", sil, True))
    elif kind == "shell_and_tube":
        rows.extend(
            [
                ("TT", "Utility Inlet Temperature", -50.0, 400.0, "degC", "AI", 0, False),
                ("TT", "Utility Outlet Temperature", -50.0, 400.0, "degC", "AI", 0, False),
                ("PT", "Utility Inlet Pressure", 0.0, 50.0, "bara", "AI", 0, False),
            ]
        )
    else:
        rows.extend(
            [
                ("TT", "Ambient Air Temperature", -40.0, 60.0, "degC", "AI", 0, False),
                ("ST", "Fan Speed", 0.0, 1800.0, "rpm", "AI", 0, False),
                ("VT", "Fan Vibration", 0.0, 100.0, "mm/s", "AI", 0, False),
            ]
        )
    return tuple(
        HeaterInstrument(
            tag=f"{symbol}-{1000 + index + 1}",
            isa_symbol=symbol,
            service=service,
            io_type=io_type,
            range_unit=unit,
            range_min=low,
            range_max=high,
            sil_rating=float(sil_rating),
            safety_related=safety,
        )
        for index, (
            symbol,
            service,
            low,
            high,
            unit,
            io_type,
            sil_rating,
            safety,
        ) in enumerate(rows)
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.heater")
