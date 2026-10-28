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

``CompressorChart``'s reduced head and polytropic efficiency, the ideal-gas kappa, the
polytropic pressure ratio, and the head's shaft work. **The chart's head crosses in metres**
and the branch converts it with ``9.81`` - the class's own ``polytropicHeadMeter`` factor - and
**the pressure is solved**, so a stated outlet pressure is the isentropic branch's.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CompressorResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.process.reference import _unported, capacity, compressor_chart
from azoth.process.reference.surge_curve import SurgeCurve

#: ``Compressor.R``, J/(mol*K).
R = 8.314
#: ``Compressor``'s own metre-to-specific-energy factor, not standard gravity.
METRE_TO_JOULE_PER_KG = 9.81
#: ``Compressor.MIN_EFFICIENCY``.
MIN_EFFICIENCY = 1.0e-6
#: ``Compressor.speed``'s initialiser, rpm.
DEFAULT_SPEED = 3000.0


def compressor(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    isentropic_efficiency: float,
    outlet_pressure: Q | None = None,
    curve: list[dict[str, Q]] | None = None,
    speed: float | None = None,
    chart_type: str | None = None,
    surge_curve: list[dict[str, Q]] | None = None,
    surge_control_margin: float | None = None,
    anti_surge_active: bool | None = None,
    surge_control_factor: float | None = None,
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
        speed: the machine's speed, rpm; absent takes ``Compressor.speed``'s 3000.
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

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    p_out = (
        None if outlet_pressure is None else input_to_si(spec, "outlet_pressure", outlet_pressure)
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
    }.get(chart_type)
    if kind is None:
        if chart_type == "khader_2015":
            _unported.refuse("chart_type=khader_2015")
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

    if not points:
        if p_out is None:
            raise InvalidInputError(
                "outlet_pressure",
                "a compressor without a curve reaches only the pressure it is given, and this "
                "one has neither",
            )
        states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)
        head = 0.0
        efficiency = 0.0
        power = n * states.dh_actual
        outlet_t = states.outlet_t
        h_out = states.h_out
        out_p = p_out
        state = _ChartState(
            head=0.0,
            efficiency=0.0,
            power=power,
            out_p=p_out,
            h_out=h_out,
            outlet_t=outlet_t,
            read=None,
            recycle_fraction=None,
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
        head, efficiency, power = state.head, state.efficiency, state.power
        out_p, h_out, outlet_t = state.out_p, state.h_out, state.outlet_t

    return CompressorResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(out_p, "Pa"),
        outlet_t=outlet_t,
        outlet_h=from_si(h_out, "J/mol"),
        polytropic_head=from_si(head, "m"),
        polytropic_efficiency=efficiency,
        power=from_si(power, "W"),
        surge_flow_rate=_flow_or_none(state.read, "surge_flow_rate"),
        surge_flow_rate_margin=_flow_or_none(state.read, "surge_flow_rate_margin"),
        surge_flow_rate_std=_flow_or_none(state.read, "surge_flow_rate_std"),
        distance_to_surge=_plain_or_none(state.read, "distance_to_surge"),
        control_line_flow=_flow_or_none(state.read, "control_line_flow"),
        distance_to_control_line=_plain_or_none(state.read, "distance_to_control_line"),
        required_recycle_fraction=_plain_or_none(state.read, "required_recycle_fraction"),
        anti_surge_recycle_power=_power_or_none(state.read, "anti_surge_recycle_power"),
        recycle_fraction=state.recycle_fraction,
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
    fresh_flow = capacity.system_volumetric_flow(components, t, p, inlet_z, n) * 3600.0
    if not n > 0.0 or not fresh_flow > 0.0 or not math.isfinite(n):
        raise InvalidInputError(
            "inlet_n",
            f"the chart branch needs a positive flow to scale, and {n} mol/s gives "
            f"{fresh_flow} m3/hr",
        )
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))
    flash = pt_flash_solve(mixture, quantity(t, "K"), quantity(p, "Pa"), list(inlet_z))
    root = flash.z_vapour if flash.phase == "all_vapour" else flash.z_liquid
    volume = pr_molar_volume(root, quantity(t, "K"), quantity(p, "Pa")).v.to("m**3/mol").magnitude
    z = p * volume / (R * t)
    state = molar_enthalpy_entropy(
        mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), list(inlet_z), root
    )
    cp = float(state.cp.to_base_units().magnitude)
    kappa = cp / (cp - R)
    # ``getZvolcorr`` = ``P v/(R T)`` at the Peneloux-shifted molar volume.
    zvolcorr = z - p * mixture.volume_shift(list(inlet_z)) / (R * t)

    boosted = n
    passes = 0
    while True:
        flow_m3hr = fresh_flow * (boosted / n)
        efficiency = min(1.0, max(MIN_EFFICIENCY, chart.polytropic_efficiency(flow_m3hr, speed)))
        head = chart.polytropic_head(flow_m3hr, speed)
        head_j = head * METRE_TO_JOULE_PER_KG
        exponent = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency)
        flow_work = exponent / (exponent - 1.0) * z * R * t / molar_mass
        if not flow_work > 0.0 or not math.isfinite(flow_work):
            raise InvalidInputError(
                "curve",
                f"the chart's head and the inlet's state give no polytropic ratio: n={exponent}, "
                f"flow work={flow_work}",
            )
        pressure_ratio = ((head_j + flow_work) / flow_work) ** (exponent / (exponent - 1.0))

        if settings.active and chart.is_surge(head, flow_m3hr):
            target = chart.surge_flow(head)
            if target is not None and math.isfinite(target) and target > 0.0:
                boosted = n * settings.control_factor * target / fresh_flow
                passes += 1
                if passes < MAX_SURGE_PASSES:
                    continue
        break

    out_p = p * pressure_ratio
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    dh = head_j * molar_mass / efficiency
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, quantity(out_p, "Pa"), quantity(h_out, "J/mol"), inlet_z
    ).T
    power = boosted * dh
    recycle_fraction = boosted / n - 1.0

    surge_flow = chart.surge_flow(power / (n * molar_mass) * efficiency / METRE_TO_JOULE_PER_KG)
    read: _SurgeRead | None = None
    if surge_flow is not None and surge_flow > 0.0:
        control = surge_flow * (1.0 + settings.control_margin)
        required = max(0.0, (control - fresh_flow) / control) if control > 0.0 else 0.0
        read = _SurgeRead(
            surge_flow_rate=surge_flow,
            surge_flow_rate_margin=fresh_flow - surge_flow,
            surge_flow_rate_std=surge_flow * (p / 1.0e5) / 1.01325 * 288.15 / t / zvolcorr,
            distance_to_surge=fresh_flow / surge_flow - 1.0,
            control_line_flow=control,
            distance_to_control_line=(fresh_flow / control - 1.0) if control > 0.0 else None,
            required_recycle_fraction=required,
            # ``getAntiSurgeRecyclePower(currentSurgeFraction)``: the added power.
            anti_surge_recycle_power=power * max(recycle_fraction, 0.0),
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
