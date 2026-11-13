"""``process.pump`` - the pump's kernel.

Spec: ``specs/models/process/pump.toml``

The Python twin of ``crates/azoth-process/src/models/pump.rs`` over
``crates/azoth-process/src/kernels/pump.rs``, written to mirror them.

# The work is an isentropic head, not ``v (P_out - P_in)``

NeqSim's ``Pump.run`` clones the fluid to the outlet pressure and runs a ``PSflash`` at the
inlet entropy, because ``calculateAsCompressor`` defaults to ``true`` and selects that
branch. So the work is ``(H(P_out, s_in) - H(P_in)) / eta`` - the enthalpy a real fluid
needs along an isentrope - and the incompressible form is its linearisation. Measured on the
captured fluid the two differ by ``1.2e-3``, which is why this is the isentropic route.

# The chart is a second branch

A non-empty ``curve`` selects ``PumpChart``'s reduced quadratic: the head becomes
``dP = H rho g``, the shaft power is ``Q dP / eta``, and the outlet is the state at the raised
pressure carrying the inlet enthalpy plus the shaft work per mole. The curve overrides
``outlet_pressure``, as the class's ``isUsePumpChart`` flag does.

# The stream record

An inlet carries ``n, z, P, T`` and an outlet carries ``n, z, P, T, h``. ``h`` is a state
function of ``(T, P, z)`` - the same flash decides it - so accepting one as well would let a
case hand over a state that does not exist.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PumpResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.bubble_pressure import bubble_pressure

# **The modules' functions, imported by name.** `azoth.eos.reference` re-exports `ph_flash`
# and `ps_flash` as *functions*, so reaching the modules through the package gives the
# function where a module is meant. The `from ... import` form resolves the module path
# itself and is unambiguous.
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve
from azoth.eos.reference.pure_saturation import pure_saturation
from azoth.process.reference import capacity
from azoth.process.reference.pump_mechanical import (
    API610_PUMP_TYPES,
    SEAL_TYPES,
    design as design_pump,
)
from azoth.process.reference.pump_chart import (
    GRAVITY,
    CurvePoint,
    PumpChart,
)

#: `Pump.speed`'s initialiser, rpm.
DEFAULT_SPEED = 1000.0
#: The class's water-like default when the fluid carries no oil or aqueous phase, cSt.
DEFAULT_VISCOSITY_CST = 1.0
#: `Pump.npshMargin`'s initialiser.
DEFAULT_NPSH_MARGIN = 1.3
#: `Pump.run`'s threshold for the fully corrected branch, cSt.
VISCOSITY_CORRECTION_THRESHOLD_CST = 4.0


def pump(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    isentropic_efficiency: float,
    outlet_pressure: Q | None = None,
    curve: list[dict[str, Q]] | None = None,
    speed: float | None = None,
    reference_density: Q | None = None,
    viscosity: Q | None = None,
    use_viscosity_correction: bool | None = None,
    check_npsh: bool | None = None,
    npsh_margin: float | None = None,
    mechanical_design: bool | None = None,
    api610_pump_type: str | None = None,
    seal_type: str | None = None,
    maximum_suction_pressure: Q | None = None,
    shutoff_head: Q | None = None,
) -> PumpResult:
    """Raise a stream's pressure, adding the pump's work as enthalpy.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the pump raises the stream to, in the isentropic branch.
            A non-empty ``curve`` overrides it.
        isentropic_efficiency: the pump's isentropic efficiency, in ``(0, 1]``.
        curve: the manufacturer's measured points, one mapping each; an empty list selects the
            isentropic branch.
        speed: the pump's speed, rpm. Absent takes ``Pump.speed``'s 1000.
        reference_density: the density the curve was measured at; absent is no correction.
        viscosity: the inlet's kinematic viscosity; absent takes the class's 1 cSt.
        use_viscosity_correction: whether the Hydraulic Institute correction is enabled.
        check_npsh: whether cavitation is assessed.
        npsh_margin: the factor NPSH available must exceed; absent takes 1.3.

    Returns:
        The outlet's record, the chart's head and efficiency, the two NPSH values and the
        shaft power.

    Raises:
        InvalidInputError: where the shapes disagree, the efficiency is outside ``(0, 1]``, a
            pump with no curve has no outlet pressure, or the curve is not a curve.
        OutOfRangeError: for a temperature or flow the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pump(
        ...     ["n-butane"],
        ...     1.0,
        ...     q(1.0, "mol/s"),
        ...     q(5.0, "bar"),
        ...     q(250.0, "K"),
        ...     isentropic_efficiency=0.75,
        ...     outlet_pressure=q(20.0, "bar"),
        ...     curve=[],
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 4)
        250.7621
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
    reference_density_si = (
        None
        if reference_density is None
        else input_to_si(spec, "reference_density", reference_density)
    )
    viscosity_si = None if viscosity is None else input_to_si(spec, "viscosity", viscosity)
    max_suction_si = (
        None
        if maximum_suction_pressure is None
        else input_to_si(spec, "maximum_suction_pressure", maximum_suction_pressure)
    )
    shutoff_si = (
        None if shutoff_head is None else input_to_si(spec, "shutoff_head", shutoff_head)
    )

    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "isentropic_efficiency": isentropic_efficiency,
            "inlet_n": n,
            "speed": speed,
        }.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the pump's isentropic efficiency must be in (0, 1], and {isentropic_efficiency} "
            f"is not",
        )

    speed_rpm = DEFAULT_SPEED if speed is None else float(speed)
    viscosity_cst = DEFAULT_VISCOSITY_CST if viscosity_si is None else float(viscosity_si) * 1.0e6
    check = bool(check_npsh)
    margin = DEFAULT_NPSH_MARGIN if npsh_margin is None else float(npsh_margin)

    design_requested = bool(mechanical_design)
    if not design_requested:
        for name, stated in (
            ("api610_pump_type", api610_pump_type is not None),
            ("seal_type", seal_type is not None),
            ("maximum_suction_pressure", max_suction_si is not None),
            ("shutoff_head", shutoff_si is not None),
        ):
            if stated:
                raise InvalidInputError(
                    name, f"`{name}` is a `mechanical_design` input, and the design is off"
                )
    if api610_pump_type is not None and api610_pump_type not in API610_PUMP_TYPES:
        raise InvalidInputError(
            "api610_pump_type", f"`{api610_pump_type}` is not an API 610 construction type"
        )
    if seal_type is not None and seal_type not in SEAL_TYPES:
        raise InvalidInputError("seal_type", f"`{seal_type}` is not a seal type")

    rows: list[dict[str, float]] = (
        [] if curve is None else input_records_to_si(spec, "curve", curve)
    )
    points = [
        CurvePoint(
            speed=float(row["speed"]),
            flow=float(row["flow"]) * 3600.0,
            head=float(row["head"]),
            efficiency=float(row["efficiency"]),
            npsh=None if row.get("npsh_required") is None else float(row["npsh_required"]),
        )
        for row in rows
    ]
    chart = (
        None
        if not points
        else PumpChart(
            points,
            reference_density_si,
            bool(use_viscosity_correction),
        )
    )

    if chart is None:
        if p_out is None:
            raise InvalidInputError(
                "outlet_pressure",
                "a pump without a curve reaches only the pressure it is given, and this one has "
                "neither",
            )
        states = _states(components, t, p, inlet_z, p_out, isentropic_efficiency)
        out_p = p_out
        h_out = states.h_out
        outlet_t = states.outlet_t
        # **The isentropic branch reports the head its work corresponds to**, and the stated
        # efficiency - the class's chart getters have nothing to read without a curve, but
        # the work per mole over ``M g`` is the head the machine developed.
        mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
        molar_mass = capacity._molar_mass(mixture, list(inlet_z))
        head = states.dh_actual / (molar_mass * GRAVITY)
        efficiency = isentropic_efficiency
        power = n * states.dh_actual
        flow_m3hr = capacity.system_volumetric_flow(components, t, p, inlet_z, n) * 3600.0
    else:
        density = capacity.system_mass_density(components, t, p, inlet_z)
        volumetric = capacity.system_volumetric_flow(components, t, p, inlet_z, n)
        flow_m3hr = volumetric * 3600.0
        corrected = (
            chart.uses_viscosity_correction() and viscosity_cst > VISCOSITY_CORRECTION_THRESHOLD_CST
        )
        if corrected:
            head = chart.fully_corrected_head(flow_m3hr, speed_rpm, density, viscosity_cst)
            efficiency = chart.corrected_efficiency(flow_m3hr, speed_rpm, viscosity_cst)
        else:
            head = chart.corrected_head(flow_m3hr, speed_rpm, density)
            efficiency = chart.efficiency(flow_m3hr, speed_rpm)
        if efficiency <= 0.0:
            raise InvalidInputError(
                "curve",
                f"the chart's efficiency at {flow_m3hr} m3/hr and {speed_rpm} rpm is "
                f"{efficiency}, and the shaft power divides by it",
            )

        mixture, ideal_gas = _components.mixture_of(components, eos="pr")
        h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
        delta_p = head * density * GRAVITY
        shaft_power = volumetric * delta_p / efficiency
        dh = shaft_power / n
        h_out = h_in + dh
        out_p = p + delta_p
        power = shaft_power
        outlet_t = ph_flash_solve(
            mixture, ideal_gas, quantity(out_p, "Pa"), quantity(h_out, "J/mol"), inlet_z
        ).T

    density_for_npsh = capacity.system_mass_density(components, t, p, inlet_z)
    npsh_required = _npsh_required(chart, flow_m3hr, speed_rpm)
    npsh_available = _npsh_available(components, t, p, inlet_z, density_for_npsh)
    cavitating = _cavitating(check, npsh_available, npsh_required, margin)

    design = (
        None
        if not design_requested
        else design_pump(
            components,
            n,
            list(inlet_z),
            t,
            p,
            out_p,
            outlet_t.to("K").magnitude,
            power,
            efficiency,
            npsh_required,
            npsh_available,
            speed_rpm,
            chart,
            api610_pump_type,
            seal_type,
            max_suction_si,
            shutoff_si,
        )
    )
    d = design

    return PumpResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(out_p, "Pa"),
        outlet_t=outlet_t,
        outlet_h=from_si(h_out, "J/mol"),
        head=from_si(head, "m"),
        efficiency=efficiency,
        npsh_required=from_si(npsh_required, "m"),
        npsh_available=from_si(npsh_available, "m"),
        cavitating=cavitating,
        power=from_si(power, "W"),
        design_pump_type=None if d is None else d.pump_type,
        design_api610_pump_type=None if d is None else d.api610_pump_type,
        design_seal_type=None if d is None else d.seal_type,
        design_rated_flow=None if d is None else from_si(d.rated_flow / 3600.0, "m**3/h"),
        design_rated_head=None if d is None else from_si(d.rated_head, "m"),
        design_hydraulic_power=None if d is None else from_si(d.hydraulic_power, "W"),
        design_rated_speed=None if d is None else d.rated_speed,
        design_specific_speed=None if d is None else d.specific_speed,
        design_suction_specific_speed=None if d is None else d.suction_specific_speed,
        design_number_of_stages=None if d is None else float(d.number_of_stages),
        design_bep_flow=(
            None if d is None or d.bep_flow is None else from_si(d.bep_flow / 3600.0, "m**3/h")
        ),
        design_bep_head=(
            None if d is None or d.bep_head is None else from_si(d.bep_head, "m")
        ),
        design_bep_source=None if d is None else d.bep_source,
        design_npshr_source=None if d is None else d.npshr_source,
        design_shutoff_source=None if d is None else d.shutoff_source,
        design_shutoff_head=(
            None if d is None or d.shutoff_head is None else from_si(d.shutoff_head, "m")
        ),
        design_pressure=None if d is None else from_si(d.design_pressure, "Pa"),
        design_temperature=None if d is None else from_si(d.design_temperature, "K"),
        design_impeller_diameter=(
            None if d is None else from_si(d.impeller_diameter, "m")
        ),
        design_impeller_width=None if d is None else from_si(d.impeller_width, "m"),
        design_shaft_diameter=None if d is None else from_si(d.shaft_diameter, "m"),
        design_casing_wall_thickness=(
            None if d is None else from_si(d.casing_wall_thickness, "m")
        ),
        design_max_allowable_working_pressure=(
            None if d is None else from_si(d.max_allowable_working_pressure, "Pa")
        ),
        design_driver_power=None if d is None else from_si(d.driver_power, "W"),
        design_driver_margin=None if d is None else d.driver_margin,
        design_suction_nozzle_size=(
            None if d is None else from_si(d.suction_nozzle_size, "m")
        ),
        design_discharge_nozzle_size=(
            None if d is None else from_si(d.discharge_nozzle_size, "m")
        ),
        design_minimum_flow=None if d is None else from_si(d.minimum_flow / 3600.0, "m**3/h"),
        design_maximum_flow=None if d is None else from_si(d.maximum_flow / 3600.0, "m**3/h"),
        design_casing_weight=None if d is None else from_si(d.casing_weight, "kg"),
        design_motor_weight=None if d is None else from_si(d.motor_weight, "kg"),
        design_baseplate_weight=(
            None if d is None else from_si(d.baseplate_weight, "kg")
        ),
        design_inner_diameter=None if d is None else from_si(d.inner_diameter, "m"),
        design_outer_diameter=None if d is None else from_si(d.outer_diameter, "m"),
        design_tantan_length=None if d is None else from_si(d.tantan_length, "m"),
        design_wall_thickness=None if d is None else from_si(d.wall_thickness, "m"),
        design_weight_piping=None if d is None else from_si(d.weight_piping, "kg"),
        design_weight_electro_instrument=(
            None if d is None else from_si(d.weight_electro_instrument, "kg")
        ),
        design_weight_structural_steel=(
            None if d is None else from_si(d.weight_structural_steel, "kg")
        ),
        design_weight_total=None if d is None else from_si(d.weight_total, "kg"),
        design_module_length=None if d is None else from_si(d.module_length, "m"),
        design_module_width=None if d is None else from_si(d.module_width, "m"),
        design_module_height=None if d is None else from_si(d.module_height, "m"),
        design_npsh_margin=None if d is None else from_si(d.npsh_margin, "m"),
        design_npsh_margin_ratio=None if d is None else d.npsh_margin_ratio,
        warnings=tuple(warnings),
    )


def _npsh_required(chart: PumpChart | None, flow_m3hr: float, speed: float) -> float:
    """``getNPSHRequired``: the curve where the chart has one, else the flow-band estimate."""
    if chart is not None:
        required = chart.npsh_required(flow_m3hr, speed)
        if required > 0.0:
            return required
    if flow_m3hr < 10.0:
        return 2.0
    if flow_m3hr < 100.0:
        return 5.0
    return 8.0


def _npsh_available(
    components: list[str], t: float, p: float, z: list[float], density: float
) -> float:
    """``(P_suction - P_vapour)/(rho g)``, m; ``NaN`` where the bubble point does not solve."""
    mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
    try:
        if len(mixture.components) == 1:
            component = mixture.components[0]
            vapor = pure_saturation(
                component.Tc, component.Pc, component.omega, quantity(t, "K")
            ).p_sat
            pressure = vapor.to("Pa").magnitude
        else:
            pressure = (
                bubble_pressure(mixture, quantity(t, "K"), list(z)).pressure.to("Pa").magnitude
            )
    except Exception:
        return math.nan
    return (p - pressure) / (density * GRAVITY)


def _cavitating(
    check_npsh: bool, npsh_available: float, npsh_required: float, npsh_margin: float
) -> bool:
    """``isCavitating``: ``NPSHa < margin * NPSHr``, with an unknown status cavitating."""
    if not check_npsh:
        return False
    if not math.isfinite(npsh_available) or not math.isfinite(npsh_required):
        return True
    if npsh_required <= 0.0:
        return True
    return npsh_available < npsh_margin * npsh_required


@dataclass(frozen=True, slots=True)
class PumpStates:
    """The pump's intermediates, in the order it computes them, in SI.

    A record rather than a dict so that a layer the harness reads and a field the result
    carries are the same typed thing, and neither can acquire a name the other lacks.
    """

    #: The inlet's molar enthalpy and entropy, J/mol and J/(mol*K).
    h_in: float
    s_in: float
    #: The isentropic outlet: the temperature at the raised pressure on the inlet's
    #: entropy, and its molar enthalpy.
    t_isentropic: Q
    h_isentropic: float
    #: The isentropic head, and the head the efficiency actually asks for.
    dh_isentropic: float
    dh_actual: float
    #: The outlet: its molar enthalpy, its temperature, and its molar entropy - which is
    #: the second-law statement of the step, since it rises iff the step was irreversible.
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
) -> PumpStates:
    """The pump's intermediates, in the order it computes them, in SI.

    **One arithmetic, two consumers.** :func:`pump` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a
    layer the harness compares is a layer the model actually took - which is the rule
    `azoth.eos.layers` states and the reason this is factored rather than restated.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    # The inlet's own state, which is where the entropy and the enthalpy come from.
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

    # The outlet's own entropy, which is the second-law statement of the head: it rises
    # across an irreversible step and is unchanged across a reversible one. `Pump` exposes
    # it as `getEntropyProduction` and nothing else about the machine's interior, which is
    # why the harness localises a wrong head through it.
    s_out, _ = entropy_at(mixture, ideal_gas, outlet_t.T.to("K").magnitude, p_out, inlet_z)

    return PumpStates(
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

    return model("process.pump")
