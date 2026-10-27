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
        {"inlet_t": t, "isentropic_efficiency": isentropic_efficiency, "speed": speed}.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the machine's isentropic efficiency must be in (0, 1], and "
            f"{isentropic_efficiency} is not",
        )

    chart_type = "fan_law" if chart_type is None else chart_type
    if chart_type == "interpolate":
        _unported.refuse("chart_type=interpolate")
    if chart_type == "interpolate_and_extrapolate":
        _unported.refuse("chart_type=interpolate_and_extrapolate")
    if chart_type == "khader_2015":
        _unported.refuse("chart_type=khader_2015")
    if chart_type != "fan_law":
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
    else:
        chart = compressor_chart.CompressorChart(points)
        head, efficiency, power, out_p, h_out, outlet_t = _chart_states(
            components,
            n,
            inlet_z,
            t,
            p,
            chart,
            DEFAULT_SPEED if speed is None else float(speed),
        )

    return CompressorResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(out_p, "Pa"),
        outlet_t=outlet_t,
        outlet_h=from_si(h_out, "J/mol"),
        polytropic_head=from_si(head, "m"),
        polytropic_efficiency=efficiency,
        power=from_si(power, "W"),
        warnings=tuple(warnings),
    )


def _chart_states(
    components: list[str],
    n: float,
    inlet_z: list[float],
    t: float,
    p: float,
    chart: compressor_chart.CompressorChart,
    speed: float,
) -> tuple[float, float, float, float, float, Q]:
    """The chart branch: head, efficiency, power, outlet pressure, enthalpy and temperature."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flow_m3hr = capacity.system_volumetric_flow(components, t, p, inlet_z, n) * 3600.0
    molar_mass = capacity._molar_mass(mixture, list(inlet_z))
    flash = pt_flash_solve(mixture, quantity(t, "K"), quantity(p, "Pa"), list(inlet_z))
    if flash.phase == "all_vapour":
        root = flash.z_vapour
    elif flash.phase in ("all_liquid", "trivial"):
        root = flash.z_liquid
    else:
        root = flash.z_liquid
    volume = pr_molar_volume(root, quantity(t, "K"), quantity(p, "Pa")).v.to("m**3/mol").magnitude
    z = p * volume / (R * t)
    state = molar_enthalpy_entropy(
        mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), list(inlet_z), root
    )
    cp = float(state.cp.to_base_units().magnitude)
    kappa = cp / (cp - R)

    efficiency = min(1.0, max(MIN_EFFICIENCY, chart.polytropic_efficiency(flow_m3hr, speed)))
    head = chart.polytropic_head(flow_m3hr, speed)
    head_j = head * METRE_TO_JOULE_PER_KG
    exponent = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency)
    flow_work = exponent / (exponent - 1.0) * z * R * t / molar_mass
    if not flow_work > 0.0:
        raise InvalidInputError(
            "curve",
            f"the chart's head and the inlet's state give no polytropic ratio: n={exponent}, "
            f"flow work={flow_work}",
        )
    pressure_ratio = ((head_j + flow_work) / flow_work) ** (exponent / (exponent - 1.0))
    out_p = p * pressure_ratio
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    dh = head_j * molar_mass / efficiency
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, quantity(out_p, "Pa"), quantity(h_out, "J/mol"), inlet_z
    ).T
    return head, efficiency, n * dh, out_p, h_out, outlet_t


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
