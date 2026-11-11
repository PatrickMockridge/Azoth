"""``process.turbo_expander_compressor`` - an expander and a compressor on one shaft.

Spec: ``specs/models/process/turbo_expander_compressor.toml``

The Python twin of ``crates/azoth-process/src/models/turbo_expander_compressor.rs`` over
``crates/azoth-process/src/kernels/turbo_expander.rs``, written to mirror them.

# The speed is the unknown

The expander's isentropic head comes from one ``PSflash`` to the stated outlet pressure; its
efficiency is the velocity-ratio parabola scaled by the design efficiency, or the 2-D map read
at the *raw* ``U/C`` when one is supplied - the class's own asymmetry, because the parabola
takes the ratio normalised by the design value. The compressor's absorbed power is its design
polytropic head scaled by ``(N/N_design)**2`` and a Q/N head factor, over a Q/N efficiency
factor times the design value. A damped Newton on ``N`` drives the expander's power to the
compressor's plus the speed-scaled bearing loss.

# The replay is two machines

The expander outlet is the isentropic expansion at the solved efficiency, which the
outlet-temperature inverse may push above one. The compressor outlet is the polytropic
pressure ratio the absorbed power implies - the class's ``setPower`` + ``setCalcPressureOut``
route.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TurboExpanderResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.process.reference import capacity, expander_chart

#: ``ThermodynamicConstantsInterface.R``, J/(mol*K).
R = 8.314
#: ``Compressor``'s own metre-to-specific-energy factor.
METRE_TO_JOULE_PER_KG = 9.81
#: The class's no-fit parabola.
LEGACY_PARABOLA = -3.56
#: ``run``'s hard-coded speed bounds, rpm.
MIN_SPEED = 1000.0
MAX_SPEED = 9000.0
#: ``run``'s Newton caps, iterations.
MAX_ITERATIONS = 50
MIN_ITERATIONS = 3
#: ``TurboExpanderCompressor``'s defaults.
DEFAULT_IMPELLER_DIAMETER = 0.424
DEFAULT_MAXIMUM_IGV_AREA = 1.637e-2
DEFAULT_IGV_AREA_INCREASE_FACTOR = 1.14
DEFAULT_BEARING_LOSS = 10.0
DEFAULT_EXPANDER_DESIGN_EFFICIENCY = 1.0
DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY = 0.81
DEFAULT_COMPRESSOR_DESIGN_HEAD_M = 20.47 * 1000.0 / METRE_TO_JOULE_PER_KG
DEFAULT_IGV_OPENING = 1.0


class _Pass:
    """One pass's quantities, in SI."""

    __slots__ = (
        "bearing",
        "compressor_efficiency",
        "efficiency",
        "head",
        "power_compressor",
        "power_expander",
        "qn_ratio",
        "total",
        "uc_norm",
    )

    def __init__(
        self,
        efficiency: float,
        total: float,
        power_expander: float,
        qn_ratio: float,
        head: float,
        compressor_efficiency: float,
        power_compressor: float,
        bearing: float,
        uc_norm: float,
    ) -> None:
        self.efficiency = efficiency
        self.total = total
        self.power_expander = power_expander
        self.qn_ratio = qn_ratio
        self.head = head
        self.compressor_efficiency = compressor_efficiency
        self.power_compressor = power_compressor
        self.bearing = bearing
        self.uc_norm = uc_norm

    def residual(self) -> float:
        """The shaft's imbalance, W."""
        return self.power_expander - (self.power_compressor + self.bearing)


def turbo_expander_compressor(
    expander_components: list[str],
    compressor_components: list[str],
    expander_n: Q,
    expander_z: list[float],
    expander_p: Q,
    expander_t: Q,
    compressor_n: Q,
    compressor_z: list[float],
    compressor_p: Q,
    compressor_t: Q,
    expander_out_pressure: Q,
    impeller_diameter: Q,
    design_speed: float,
    design_uc: float,
    design_qn: float,
    igv_opening: float | None = None,
    igv_control_mode: bool | None = None,
    maximum_igv_area: Q | None = None,
    igv_area_increase_factor: float | None = None,
    bearing_loss: Q | None = None,
    expander_design_isentropic_efficiency: float | None = None,
    compressor_design_polytropic_efficiency: float | None = None,
    compressor_design_polytropic_head: Q | None = None,
    expander_out_temperature: Q | None = None,
    uc_curve: list[dict[str, float]] | None = None,
    qn_efficiency_curve: list[dict[str, float]] | None = None,
    qn_head_curve: list[dict[str, float]] | None = None,
    chart: list[dict[str, float]] | None = None,
    igv_penalty_curve: list[dict[str, float]] | None = None,
) -> TurboExpanderResult:
    """Balance an expander and a compressor on one shaft.

    Args:
        expander_components: the expander feed's substances, by name.
        compressor_components: the compressor feed's substances, by name.
        expander_n: the expander feed's molar flow.
        expander_z: the expander feed's composition.
        expander_p: the expander feed's pressure.
        expander_t: the expander feed's temperature.
        compressor_n: the compressor feed's molar flow.
        compressor_z: the compressor feed's composition.
        compressor_p: the compressor feed's pressure.
        compressor_t: the compressor feed's temperature.
        expander_out_pressure: the expander's stated outlet pressure.
        impeller_diameter: the impeller's outer diameter.
        design_speed: the design speed, rpm.
        design_uc: the design velocity ratio.
        design_qn: the design flow coefficient, m3/rev.
        igv_opening: the stated inlet-guide-vane opening, 0..1.
        igv_control_mode: whether the stated opening is held.
        maximum_igv_area: the maximum IGV area.
        igv_area_increase_factor: the expanded-area factor.
        bearing_loss: the bearing loss at the design speed.
        expander_design_isentropic_efficiency: the expander's design efficiency.
        compressor_design_polytropic_efficiency: the compressor's design efficiency.
        compressor_design_polytropic_head: the compressor's design polytropic head.
        expander_out_temperature: a target expander outlet temperature, selecting the inverse.
        uc_curve: the legacy velocity-ratio efficiency curve.
        qn_efficiency_curve: the compressor's Q/N efficiency correction.
        qn_head_curve: the compressor's Q/N head correction.
        chart: the expander's 2-D map; absent is the legacy parabola.
        igv_penalty_curve: the IGV efficiency penalty.

    Returns:
        Both outlet records, the solved speed, both powers and the machine's ratios.

    Raises:
        InvalidInputError: for a non-positive design or geometry value, an efficiency outside
            ``(0, 1]``, a non-positive isentropic drop or a map that is not a map; and whatever
            the flash refuses.

    See :func:`azoth.process.turbo_expander_compressor`.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    _check_feed(expander_components, expander_z)
    _check_feed(compressor_components, compressor_z)

    e_n = input_to_si(spec, "expander_n", expander_n)
    e_p = input_to_si(spec, "expander_p", expander_p)
    e_t = input_to_si(spec, "expander_t", expander_t)
    c_n = input_to_si(spec, "compressor_n", compressor_n)
    c_p = input_to_si(spec, "compressor_p", compressor_p)
    c_t = input_to_si(spec, "compressor_t", compressor_t)
    p_out = input_to_si(spec, "expander_out_pressure", expander_out_pressure)
    diameter = input_to_si(spec, "impeller_diameter", impeller_diameter)

    apply_checks(
        checks.on_input,
        {
            "design_speed": design_speed,
            "design_uc": design_uc,
            "design_qn": design_qn,
            "impeller_diameter": diameter,
            "expander_design_isentropic_efficiency": (expander_design_isentropic_efficiency),
            "compressor_design_polytropic_efficiency": (compressor_design_polytropic_efficiency),
        }.get,
        warnings,
    )

    design_eff = (
        DEFAULT_EXPANDER_DESIGN_EFFICIENCY
        if expander_design_isentropic_efficiency is None
        else float(expander_design_isentropic_efficiency)
    )
    compressor_eff = (
        DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY
        if compressor_design_polytropic_efficiency is None
        else float(compressor_design_polytropic_efficiency)
    )
    _validate(
        diameter,
        design_speed,
        design_uc,
        design_qn,
        DEFAULT_MAXIMUM_IGV_AREA
        if maximum_igv_area is None
        else input_to_si(spec, "maximum_igv_area", maximum_igv_area),
        DEFAULT_IGV_AREA_INCREASE_FACTOR
        if igv_area_increase_factor is None
        else float(igv_area_increase_factor),
        DEFAULT_BEARING_LOSS
        if bearing_loss is None
        else input_to_si(spec, "bearing_loss", bearing_loss),
        DEFAULT_COMPRESSOR_DESIGN_HEAD_M
        if compressor_design_polytropic_head is None
        else input_to_si(
            spec, "compressor_design_polytropic_head", compressor_design_polytropic_head
        ),
        design_eff,
        compressor_eff,
    )

    opening = DEFAULT_IGV_OPENING if igv_opening is None else float(igv_opening)
    hold = False if igv_control_mode is None else bool(igv_control_mode)
    area = (
        DEFAULT_MAXIMUM_IGV_AREA
        if maximum_igv_area is None
        else input_to_si(spec, "maximum_igv_area", maximum_igv_area)
    )
    increase = (
        DEFAULT_IGV_AREA_INCREASE_FACTOR
        if igv_area_increase_factor is None
        else float(igv_area_increase_factor)
    )
    bearing_loss_w = (
        DEFAULT_BEARING_LOSS
        if bearing_loss is None
        else input_to_si(spec, "bearing_loss", bearing_loss)
    )
    design_head_m = (
        DEFAULT_COMPRESSOR_DESIGN_HEAD_M
        if compressor_design_polytropic_head is None
        else input_to_si(
            spec, "compressor_design_polytropic_head", compressor_design_polytropic_head
        )
    )
    target = (
        None
        if expander_out_temperature is None
        else input_to_si(spec, "expander_out_temperature", expander_out_temperature)
    )

    expander_mixture, expander_ideal_gas = _components.mixture_of(expander_components, eos="pr")
    compressor_mixture, compressor_ideal_gas = _components.mixture_of(
        compressor_components, eos="pr"
    )
    expander_mass = capacity._molar_mass(expander_mixture, list(expander_z))
    compressor_mass = capacity._molar_mass(compressor_mixture, list(compressor_z))
    expander_mass_flow = e_n * expander_mass
    compressor_mass_flow = c_n * compressor_mass
    expander_state = _feed_state(expander_components, e_t, e_p, list(expander_z))
    compressor_state = _feed_state(compressor_components, c_t, c_p, list(compressor_z))
    compressor_flow = c_n * compressor_state.molar_volume
    expander_density = expander_mass / expander_state.molar_volume

    h_in, _ = enthalpy_at(expander_mixture, expander_ideal_gas, e_t, e_p, list(expander_z))
    entropy, _ = entropy_at(expander_mixture, expander_ideal_gas, e_t, e_p, list(expander_z))
    isentropic = ps_flash_solve(
        expander_mixture,
        expander_ideal_gas,
        from_si(p_out, "Pa"),
        from_si(entropy, "J/(mol*K)"),
        list(expander_z),
    )
    h_isentropic, _ = enthalpy_at(
        expander_mixture,
        expander_ideal_gas,
        isentropic.T.to("K").magnitude,
        p_out,
        list(expander_z),
    )
    isentropic_drop = (h_in - h_isentropic) / expander_mass
    if not isentropic_drop > 0.0:
        raise InvalidInputError(
            "expander_out_pressure",
            f"the expander's isentropic drop is {isentropic_drop} J/kg at {p_out} Pa, and the "
            f"machine needs a positive one",
        )
    spouting = math.sqrt(2.0 * isentropic_drop)

    required: float | None = None
    if target is not None:
        h_target, _ = enthalpy_at(
            expander_mixture, expander_ideal_gas, target, p_out, list(expander_z)
        )
        actual = (h_in - h_target) / expander_mass
        required = actual / isentropic_drop

    uc_rows = _rows(spec, "uc_curve", uc_curve)
    qn_efficiency_rows = _rows(spec, "qn_efficiency_curve", qn_efficiency_curve)
    qn_head_rows = _rows(spec, "qn_head_curve", qn_head_curve)
    penalty_rows = _rows(spec, "igv_penalty_curve", igv_penalty_curve)
    chart_rows = _rows(spec, "chart", chart)

    parabola = _parabola_coefficient([(row["uc"], row["efficiency"]) for row in uc_rows])
    qn_efficiency = sorted((row["qn"], row["factor"]) for row in qn_efficiency_rows)
    qn_head = sorted((row["qn"], row["factor"]) for row in qn_head_rows)
    penalty = sorted((row["igv"], row["factor"]) for row in penalty_rows)
    map_chart = (
        None
        if not chart_rows
        else expander_chart.ExpanderChart(
            [
                expander_chart.CurvePoint(
                    igv=row["igv"], uc=row["uc"], efficiency=row["efficiency"], head=0.0
                )
                for row in chart_rows
            ],
            1.0,
        )
    )

    def evaluate(speed: float) -> _Pass:
        tip = math.pi * diameter * speed / 60.0
        uc_raw = tip / spouting
        uc_norm = uc_raw / design_uc
        if map_chart is None:
            efficiency = design_eff * _efficiency_from_uc(parabola, uc_norm)
        else:
            efficiency = map_chart.efficiency(uc_raw, opening)
        efficiency *= _igv_penalty(penalty, opening)
        if efficiency < 0.0:
            efficiency = 0.0
        total = max(0.0, efficiency / design_eff) if design_eff > 1.0e-9 else 0.0
        efficiency = efficiency if required is None else required
        power_expander = expander_mass_flow * isentropic_drop * efficiency

        qn_ratio = (compressor_flow * 60.0 / speed) / design_qn
        efficiency_factor = max(1.0e-6, _spline(qn_efficiency, qn_ratio))
        head_factor = max(1.0e-6, _spline(qn_head, qn_ratio))
        head = (
            design_head_m
            * METRE_TO_JOULE_PER_KG
            * (speed / design_speed)
            * (speed / design_speed)
            * head_factor
        )
        compressor_efficiency = efficiency_factor * compressor_eff
        power_compressor = compressor_mass_flow * head / compressor_efficiency
        bearing = bearing_loss_w * (speed / design_speed) * (speed / design_speed)
        return _Pass(
            efficiency=efficiency,
            total=total,
            power_expander=power_expander,
            qn_ratio=qn_ratio,
            head=head,
            compressor_efficiency=compressor_efficiency,
            power_compressor=power_compressor,
            bearing=bearing,
            uc_norm=uc_norm,
        )

    speed = design_speed
    delta = 10.0
    iteration = 0
    last = evaluate(speed)
    while True:
        first = evaluate(speed)
        second = evaluate(speed + delta)
        residual = first.residual()
        derivative = (second.residual() - residual) / delta
        if abs(derivative) < 1.0e-8:
            delta += 10.0
            speed += 10.0
        else:
            damping = (1.0 + iteration) / (iteration + 5.0)
            speed -= damping * residual / derivative
        last = first
        if speed > MAX_SPEED:
            speed = MAX_SPEED
            break
        if speed < MIN_SPEED:
            speed = MIN_SPEED
            break
        iteration += 1
        if not (
            (abs(residual) * 100.0 > 1.0e-3 and iteration < MAX_ITERATIONS)
            or iteration < MIN_ITERATIONS
        ):
            break

    design_efficiency = (
        required / last.total if required is not None and last.total > 1.0e-9 else design_eff
    )

    if hold:
        report_opening = min(1.0, max(0.0, opening))
    else:
        report_opening = _evaluate_igv(
            isentropic_drop,
            expander_mass_flow,
            expander_density,
            area,
            increase,
        )

    if not last.efficiency > 0.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the replay needs a positive efficiency, and {last.efficiency} is not",
        )
    expander_outlet_h = h_in + (h_isentropic - h_in) * last.efficiency
    expander_outlet_t = ph_flash_solve(
        expander_mixture,
        expander_ideal_gas,
        from_si(p_out, "Pa"),
        from_si(expander_outlet_h, "J/mol"),
        list(expander_z),
    ).T

    (
        compressor_outlet_p_value,
        compressor_outlet_h,
        compressor_outlet_t,
    ) = _compressor_replay(
        compressor_mixture,
        compressor_ideal_gas,
        list(compressor_z),
        c_t,
        c_p,
        c_n,
        compressor_mass,
        compressor_state,
        last.power_compressor,
        last.compressor_efficiency,
    )

    return TurboExpanderResult(
        expander_outlet_n=from_si(e_n, "mol/s"),
        expander_outlet_z=tuple(expander_z),
        expander_outlet_p=from_si(p_out, "Pa"),
        expander_outlet_t=expander_outlet_t,
        expander_outlet_h=from_si(expander_outlet_h, "J/mol"),
        compressor_outlet_n=from_si(c_n, "mol/s"),
        compressor_outlet_z=tuple(compressor_z),
        compressor_outlet_p=from_si(compressor_outlet_p_value, "Pa"),
        compressor_outlet_t=compressor_outlet_t,
        compressor_outlet_h=from_si(compressor_outlet_h, "J/mol"),
        speed=speed,
        power_expander=from_si(last.power_expander, "W"),
        power_compressor=from_si(last.power_compressor, "W"),
        expander_isentropic_efficiency=last.efficiency,
        expander_design_isentropic_efficiency=design_efficiency,
        compressor_polytropic_efficiency=last.compressor_efficiency,
        compressor_polytropic_head=from_si(last.head / METRE_TO_JOULE_PER_KG, "m"),
        uc_ratio_expander=last.uc_norm,
        qn_ratio_compressor=last.qn_ratio,
        igv_opening=report_opening,
        warnings=tuple(warnings),
    )


def _validate(
    diameter: float,
    design_speed: float,
    design_uc: float,
    design_qn: float,
    maximum_area: float,
    increase_factor: float,
    bearing_loss: float,
    design_head_m: float,
    design_eff: float,
    compressor_eff: float,
) -> None:
    """The class's own refusal of a design number the machine cannot run with."""
    for name, value, positive in (
        ("impeller_diameter", diameter, True),
        ("design_speed", design_speed, True),
        ("design_uc", design_uc, True),
        ("design_qn", design_qn, True),
        ("maximum_igv_area", maximum_area, True),
        ("igv_area_increase_factor", increase_factor, True),
        ("bearing_loss", bearing_loss, False),
        ("compressor_design_head_m", design_head_m, True),
    ):
        if (
            not math.isfinite(value)
            or (positive and value <= 0.0)
            or (not positive and value < 0.0)
        ):
            raise InvalidInputError(name, f"the machine needs a valid {name}, and {value} is not")
    for name, value in (
        ("expander_design_isentropic_efficiency", design_eff),
        ("compressor_design_polytropic_efficiency", compressor_eff),
    ):
        if not 0.0 < value <= 1.0:
            raise InvalidInputError(name, f"{name} must be in (0, 1], and {value} is not")


def _parabola_coefficient(curve: list[tuple[float, float]]) -> float:
    """``setUCcurve``'s constrained parabola: ``a (uc-1)**2 + 1``, the vertex at one."""
    if len(curve) < 2:
        return 0.0
    numerator = 0.0
    denominator = 0.0
    for uc, efficiency in curve:
        dx = uc - 1.0
        dx2 = dx * dx
        numerator += (efficiency - 1.0) * dx2
        denominator += dx2 * dx2
    return numerator / denominator if denominator != 0.0 else 0.0


def _efficiency_from_uc(coefficient: float, uc: float) -> float:
    """``getEfficiencyFromUC``."""
    if coefficient != 0.0:
        return coefficient * (uc - 1.0) * (uc - 1.0) + 1.0
    return LEGACY_PARABOLA * (uc - 1.0) * (uc - 1.0) + 1.0


def _spline(curve: list[tuple[float, float]], qn: float) -> float:
    """``getEfficiencyFromQN``/``getHeadFromQN``: cubic Hermite, linear extrapolation."""
    if len(curve) < 2:
        return 1.0
    n = len(curve)
    if qn <= curve[0][0]:
        slope = (curve[1][1] - curve[0][1]) / (curve[1][0] - curve[0][0])
        return curve[0][1] + slope * (qn - curve[0][0])
    if qn >= curve[n - 1][0]:
        slope = (curve[n - 1][1] - curve[n - 2][1]) / (curve[n - 1][0] - curve[n - 2][0])
        return curve[n - 1][1] + slope * (qn - curve[n - 1][0])
    index = 0
    for i in range(n - 1):
        if curve[i][0] <= qn <= curve[i + 1][0]:
            index = i
            break
    x0, y0 = curve[index]
    x1, y1 = curve[index + 1]
    h = x1 - x0
    t = (qn - x0) / h
    m0 = (
        (y1 - y0) / (x1 - x0)
        if index == 0
        else (y1 - curve[index - 1][1]) / (x1 - curve[index - 1][0])
    )
    m1 = (
        (y1 - y0) / (x1 - x0)
        if index == n - 2
        else (curve[index + 2][1] - y0) / (curve[index + 2][0] - x0)
    )
    h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t)
    h10 = t * (1.0 - t) * (1.0 - t)
    h01 = t * t * (3.0 - 2.0 * t)
    h11 = t * t * (t - 1.0)
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1


def _igv_penalty(curve: list[tuple[float, float]], igv: float) -> float:
    """``getIgvEfficiencyPenalty``: edge-clamped linear, unity without a curve."""
    if len(curve) < 2:
        return 1.0
    n = len(curve)
    if igv <= curve[0][0]:
        return curve[0][1]
    if igv >= curve[n - 1][0]:
        return curve[n - 1][1]
    for i in range(n - 1):
        if curve[i][0] <= igv <= curve[i + 1][0]:
            t = (igv - curve[i][0]) / (curve[i + 1][0] - curve[i][0])
            return curve[i][1] + t * (curve[i + 1][1] - curve[i][1])
    return curve[n - 1][1]


def _evaluate_igv(
    stage_drop: float,
    mass_flow: float,
    density: float,
    maximum_area: float,
    increase_factor: float,
) -> float:
    """``evaluateIGV``: the opening the stage drop and flow imply."""
    if not mass_flow > 0.0 or not stage_drop > 0.0 or not density > 0.0:
        return 0.0
    half = 0.5 * stage_drop
    if half <= 0.0:
        half = max(stage_drop * 0.5, 1.0)
    velocity = math.sqrt(2.0 * max(half, 1.0e-6))
    required = mass_flow / (density * velocity)
    expanded = maximum_area * increase_factor
    available = expanded if required > maximum_area and expanded > maximum_area else maximum_area
    if available > 0.0:
        return min(1.0, max(0.0, required / available))
    return 0.0


def _compressor_replay(
    mixture: Any,
    ideal_gas: Any,
    z_feed: list[float],
    t: float,
    p: float,
    n: float,
    molar_mass: float,
    state: _FeedState,
    power: float,
    efficiency: float,
) -> tuple[float, float, Q]:
    """The polytropic pressure ratio a stated power implies, and the outlet's state."""
    if not 0.0 < efficiency <= 1.0:
        raise InvalidInputError(
            "compressor_polytropic_efficiency",
            f"the replay efficiency is {efficiency}, and it must be in (0, 1]",
        )
    z = state.z
    cp = state.cp
    kappa = cp / (cp - R)
    exponent = 1.0 / (1.0 - (kappa - 1.0) / kappa / efficiency)
    dh = power / n
    head = dh * efficiency / molar_mass
    flow_work = exponent / (exponent - 1.0) * z * R * t / molar_mass
    if not flow_work > 0.0 or not math.isfinite(flow_work) or not math.isfinite(exponent):
        raise InvalidInputError(
            "power_compressor",
            f"the replay has no polytropic ratio: n={exponent}, flow work={flow_work}",
        )
    pressure_ratio = ((head + flow_work) / flow_work) ** (exponent / (exponent - 1.0))
    out_p = p * pressure_ratio
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, list(z_feed))
    h_out = h_in + dh
    outlet_t = ph_flash_solve(
        mixture, ideal_gas, from_si(out_p, "Pa"), from_si(h_out, "J/mol"), list(z_feed)
    ).T
    return out_p, h_out, outlet_t


def _rows(
    spec: dict[str, object], name: str, rows: list[dict[str, float]] | None
) -> list[dict[str, float]]:
    """One record list, converted at the boundary; absent is empty."""
    if not rows:
        return []
    converted: list[dict[str, float]] = input_records_to_si(spec, name, rows)
    return converted


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.turbo_expander_compressor")


class _FeedState:
    """The feed's own state, as the class reads it: `getZ()`, `getCP()/n` and the molar volume."""

    __slots__ = ("cp", "molar_volume", "z")

    def __init__(self, z: float, cp: float, molar_volume: float) -> None:
        self.z = z
        self.cp = cp
        self.molar_volume = molar_volume


def _feed_state(components: list[str], t: float, p: float, z: list[float]) -> _FeedState:
    """Read the feed's state at its own composition.

    **A single phase is evaluated at the feed's composition, not the incipient phase's.**
    ``pt_flash`` on a superheated stream reports the *incipient* phase's root, and with a large
    negative vapour fraction that composition sits far from the feed's; the class holds one
    phase whose composition is the feed's. A genuinely two-phase feed is the beta-weighted mean
    of its phases, which is what ``SystemThermo.getZ()`` and ``getCP()`` sum over.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flash = pt_flash_solve(mixture, quantity(t, "K"), quantity(p, "Pa"), list(z))

    def at(composition: list[float], root: float) -> tuple[float, float]:
        cp = float(
            molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                quantity(t, "K"),
                quantity(p, "Pa"),
                list(composition),
                root,
            )
            .cp.to_base_units()
            .magnitude
        )
        volume = root * R * t / p - mixture.volume_shift(list(composition))
        return cp, volume

    if flash.phase == "two_phase":
        beta = flash.vapour_fraction if flash.vapour_fraction is not None else 1.0
        cp_vapour, volume_vapour = at(list(flash.y), flash.z_vapour)
        cp_liquid, volume_liquid = at(list(flash.x), flash.z_liquid)
        state = _FeedState(
            z=beta * flash.z_vapour + (1.0 - beta) * flash.z_liquid,
            cp=beta * cp_vapour + (1.0 - beta) * cp_liquid,
            molar_volume=beta * volume_vapour + (1.0 - beta) * volume_liquid,
        )
    else:
        reduced = reduced_parameters(mixture, t, p)
        kij = mixture.effective_kij(reduced.t_kelvin)
        liquid = flash.phase != "all_vapour"
        root = float(phase_state(reduced, kij, list(z), liquid=liquid).z)
        cp, volume = at(list(z), root)
        state = _FeedState(z=root, cp=cp, molar_volume=volume)
    if not state.molar_volume > 0.0:
        raise InvalidInputError(
            "stream",
            f"the feed's molar volume is {state.molar_volume} m3/mol, and a state needs a "
            f"positive one",
        )
    return state


def _check_feed(components: list[str], z: list[float]) -> None:
    """Refuse a composition that is not one entry per component, as `Stream::from_pt` does."""
    if len(z) != len(components):
        raise InvalidInputError(
            "z", f"a feed for {len(components)} components has {len(z)} entries"
        )
