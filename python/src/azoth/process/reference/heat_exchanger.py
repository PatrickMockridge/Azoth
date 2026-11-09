"""``process.heat_exchanger`` - the exchanger's kernel.

Spec: ``specs/models/process/heat_exchanger.toml``

The Python twin of ``crates/azoth-process/src/models/heat_exchanger.rs`` over
``crates/azoth-process/src/kernels/heat_exchanger.rs``, written to mirror them.

# The rating, not a duty

``HeatExchanger.run``'s default branch specifies ``UAvalue`` and runs an effectiveness-NTU
rating: each side's capacity is estimated by flashing it at the **other side's inlet
temperature** and dividing the total enthalpy change by that temperature change, then
``NTU = UA / C_min`` and ``calcThermalEffectivenes(NTU, C_min / C_max)`` scales that swing.
``duty`` is an *output* of this class - ``energyInput`` appears nowhere in it - so a kernel
that moved a caller's duty was a different unit operation under the same name.

``C`` is a difference of *totals* over a difference of temperatures, which is W/K; that is
what makes the NTU dimensionless. The absolute value matters: the hot side's seeded swing
is negative because it cools, and taking it at face value would make ``C_min`` and the NTU
negative, which is an exponential in the wrong direction.

At an effectiveness of one ``run`` leaves the capacity-limiting side at its seeded state -
"this side at the other's inlet temperature". That is the same *enthalpy* the arithmetic
gives, so the branch is an optimisation rather than a second answer, and this reaches it
through the flash either way.

# The other mode pins an outlet

``hot_outlet_temperature`` or ``cold_outlet_temperature`` pins that side at its inlet
pressure and energy-balances the other against it. **NeqSim's own ``runSpecifiedStream``
does not reach the state it names**: it clones the stream's already-flashed fluid, sets the
temperature and flashes again, and the enthalpy it lands on is one NeqSim's own ``PHflash``
places at another temperature. This reaches the pinned temperature directly, which is why
the case is pinned to the capture's ``reference_pin_*`` rows rather than to its
``out_temperature_*`` ones.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.numerics import divide
from azoth.core.range import apply_checks, checks_for
from azoth.core.record_outputs_gen import HeaterMechanicalType, HeatExchangerMechanicalType
from azoth.core.result import HeatExchangerResult, Phase
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import phase_label
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.phase_transport import phase_transport
from azoth.eos.reference.pr_mass_density import pr_mass_density
from azoth.eos.reference.pr_molar_volume import pr_molar_volume
from azoth.eos.reference.pt_flash import pt_flash as pt_flash_solve
from azoth.process.reference.heater import (
    MECHANICAL_CRITERIA,
    MECHANICAL_NAMES,
    MECHANICAL_TYPES,
    _geometry,
    _lmtd,
)
from azoth.thermal.reference.shell_and_tube_rating import (
    shell_and_tube_rating as _rating_calculator,
)

#: The arrangements ``calcThermalEffectivenes`` names, in azoth's own spelling.
FLOW_ARRANGEMENTS = ("counterflow", "parallelflow", "shell_and_tube")


def heat_exchanger(
    hot_components: list[str],
    hot_in_n: Q,
    hot_in_z: list[float],
    hot_in_p: Q,
    hot_in_t: Q,
    cold_components: list[str],
    cold_in_n: Q,
    cold_in_z: list[float],
    cold_in_p: Q,
    cold_in_t: Q,
    ua: Q | None = None,
    flow_arrangement: str = "counterflow",
    hot_outlet_temperature: Q | None = None,
    cold_outlet_temperature: Q | None = None,
    use_delta_t: bool | None = None,
    delta_t: Q | None = None,
    guess_out_temperature: Q | None = None,
    mechanical_design: bool | None = None,
    mechanical_type: str | None = None,
    mechanical_selection: str | None = None,
    use_rating_pressure_drop: bool | None = None,
    rating_area: Q | None = None,
    rating_tube_side: str | None = None,
    rating_tube_outside_diameter: Q | None = None,
    rating_tube_inside_diameter: Q | None = None,
    rating_tube_length: Q | None = None,
    rating_tube_count: float | None = None,
    rating_tube_passes: float | None = None,
    rating_tube_pitch: Q | None = None,
    rating_triangular_pitch: bool | None = None,
    rating_shell_inside_diameter: Q | None = None,
    rating_baffle_spacing: Q | None = None,
    rating_baffle_count: float | None = None,
    rating_baffle_cut: float | None = None,
    rating_tube_wall_conductivity: Q | None = None,
    rating_fouling_tube: float | None = None,
    rating_fouling_shell: float | None = None,
    rating_tube_fouling_layer_thickness: Q | None = None,
    rating_tube_fouling_layer_conductivity: Q | None = None,
    rating_shell_side_method: str | None = None,
    rating_shell_viscosity_wall: Q | None = None,
    rating_has_sealing: bool | None = None,
    rating_sealing_pairs: float | None = None,
    rating_tube_to_baffle_clearance: Q | None = None,
    rating_shell_to_baffle_clearance: Q | None = None,
    rating_bypass_area: Q | None = None,
    time_step: Q | None = None,
    dynamic_model_enabled: bool | None = None,
    transient_tube_side: str | None = None,
    wall_mass: Q | None = None,
    wall_cp: Q | None = None,
    wall_temperature: Q | None = None,
    shell_holdup_volume: Q | None = None,
    tube_holdup_volume: Q | None = None,
    shell_side_htc: Q | None = None,
    tube_side_htc: Q | None = None,
    heat_transfer_area: Q | None = None,
    shell_fluid_temperature: Q | None = None,
    tube_fluid_temperature: Q | None = None,
) -> HeatExchangerResult:
    """Exchange heat between two streams.

    Args:
        hot_components: the hot side's substances, by name.
        hot_in_n: the hot inlet's molar flow.
        hot_in_z: the hot inlet's composition.
        hot_in_p: the hot inlet's pressure.
        hot_in_t: the hot inlet's temperature.
        cold_components: the cold side's substances, by name.
        cold_in_n: the cold inlet's molar flow.
        cold_in_z: the cold inlet's composition.
        cold_in_p: the cold inlet's pressure.
        cold_in_t: the cold inlet's temperature.
        ua: the overall conductance, for the rating. Unread when an outlet is pinned.
        flow_arrangement: ``counterflow``, ``parallelflow`` or ``shell_and_tube``.
        hot_outlet_temperature: pin the hot outlet at this temperature.
        cold_outlet_temperature: pin the cold outlet at this temperature.
        mechanical_design: run the mechanical design's thermal estimate.
        mechanical_type: the manual type selection.
        mechanical_selection: the automatic ranking's metric.

    Returns:
        Both outlets' records.

    Raises:
        InvalidInputError: for the mode and arrangement refusals the spec lists.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = heat_exchanger(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.8, 0.2],
        ...     q(20.0, "bar"),
        ...     q(400.0, "K"),
        ...     ["n-butane", "n-pentane"],
        ...     q(1.0, "mol/s"),
        ...     [0.5, 0.5],
        ...     q(5.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(100.0, "W/K"),
        ...     "counterflow",
        ... )
        >>> round(r.hot_out_t.to("K").magnitude, 4)
        313.7571
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    ua_si = None if ua is None else input_to_si(spec, "ua", ua)
    hot_t_si = input_to_si(spec, "hot_in_t", hot_in_t)
    cold_t_si = input_to_si(spec, "cold_in_t", cold_in_t)
    hot_out_t_si = (
        None
        if hot_outlet_temperature is None
        else input_to_si(spec, "hot_outlet_temperature", hot_outlet_temperature)
    )
    cold_out_t_si = (
        None
        if cold_outlet_temperature is None
        else input_to_si(spec, "cold_outlet_temperature", cold_outlet_temperature)
    )
    delta_t_si = 1.0 if delta_t is None else input_to_si(spec, "delta_t", delta_t)
    # The class's `guessOutTemperature` seeds a pass the answering branch overwrites.
    _ = (
        None
        if guess_out_temperature is None
        else input_to_si(spec, "guess_out_temperature", guess_out_temperature)
    )
    time_step_si = None if time_step is None else input_to_si(spec, "time_step", time_step)
    wall_mass_si = None if wall_mass is None else input_to_si(spec, "wall_mass", wall_mass)
    wall_cp_si = None if wall_cp is None else input_to_si(spec, "wall_cp", wall_cp)
    wall_temperature_si = (
        None
        if wall_temperature is None
        else input_to_si(spec, "wall_temperature", wall_temperature)
    )
    shell_holdup_si = (
        None
        if shell_holdup_volume is None
        else input_to_si(spec, "shell_holdup_volume", shell_holdup_volume)
    )
    tube_holdup_si = (
        None
        if tube_holdup_volume is None
        else input_to_si(spec, "tube_holdup_volume", tube_holdup_volume)
    )
    shell_side_htc_si = (
        None if shell_side_htc is None else input_to_si(spec, "shell_side_htc", shell_side_htc)
    )
    tube_side_htc_si = (
        None if tube_side_htc is None else input_to_si(spec, "tube_side_htc", tube_side_htc)
    )
    area_si = (
        None
        if heat_transfer_area is None
        else input_to_si(spec, "heat_transfer_area", heat_transfer_area)
    )
    shell_fluid_temperature_si = (
        None
        if shell_fluid_temperature is None
        else input_to_si(spec, "shell_fluid_temperature", shell_fluid_temperature)
    )
    tube_fluid_temperature_si = (
        None
        if tube_fluid_temperature is None
        else input_to_si(spec, "tube_fluid_temperature", tube_fluid_temperature)
    )
    # Every declared bound's quantity is fed here: a closure that supplies some and not
    # the rest emits a `RANGE_CHECK_SKIPPED` warning for a check that could have run.
    apply_checks(
        checks.on_input,
        {
            "ua": ua_si,
            "hot_in_t": hot_t_si,
            "cold_in_t": cold_t_si,
            "rating_area": _optional_si(spec, "rating_area", rating_area),
            "rating_tube_outside_diameter": _optional_si(
                spec, "rating_tube_outside_diameter", rating_tube_outside_diameter
            ),
            "rating_tube_inside_diameter": _optional_si(
                spec, "rating_tube_inside_diameter", rating_tube_inside_diameter
            ),
            "rating_tube_length": _optional_si(spec, "rating_tube_length", rating_tube_length),
            "rating_tube_count": rating_tube_count,
            "rating_tube_passes": rating_tube_passes,
            "rating_tube_pitch": _optional_si(spec, "rating_tube_pitch", rating_tube_pitch),
            "rating_shell_inside_diameter": _optional_si(
                spec, "rating_shell_inside_diameter", rating_shell_inside_diameter
            ),
            "rating_baffle_spacing": _optional_si(
                spec, "rating_baffle_spacing", rating_baffle_spacing
            ),
            "rating_baffle_count": rating_baffle_count,
            "rating_baffle_cut": rating_baffle_cut,
            "rating_tube_wall_conductivity": _optional_si(
                spec, "rating_tube_wall_conductivity", rating_tube_wall_conductivity
            ),
            "rating_fouling_tube": rating_fouling_tube,
            "rating_fouling_shell": rating_fouling_shell,
            "rating_tube_fouling_layer_thickness": _optional_si(
                spec, "rating_tube_fouling_layer_thickness", rating_tube_fouling_layer_thickness
            ),
            "rating_tube_fouling_layer_conductivity": _optional_si(
                spec,
                "rating_tube_fouling_layer_conductivity",
                rating_tube_fouling_layer_conductivity,
            ),
            "rating_sealing_pairs": rating_sealing_pairs,
            "rating_tube_to_baffle_clearance": _optional_si(
                spec, "rating_tube_to_baffle_clearance", rating_tube_to_baffle_clearance
            ),
            "rating_shell_to_baffle_clearance": _optional_si(
                spec, "rating_shell_to_baffle_clearance", rating_shell_to_baffle_clearance
            ),
            "rating_bypass_area": _optional_si(spec, "rating_bypass_area", rating_bypass_area),
            "time_step": time_step_si,
            "wall_mass": wall_mass_si,
            "wall_cp": wall_cp_si,
            "shell_holdup_volume": shell_holdup_si,
            "tube_holdup_volume": tube_holdup_si,
            "shell_side_htc": shell_side_htc_si,
            "tube_side_htc": tube_side_htc_si,
            "heat_transfer_area": area_si,
        }.get,
        warnings,
    )
    if use_delta_t and (hot_out_t_si is not None or cold_out_t_si is not None):
        warnings.append(
            Warning(
                WarningCode.FALLBACK_TO_INPUT,
                "a pinned outlet wins in `HeatExchanger.run`: `runDeltaT` delegates to "
                "`runSpecifiedStream`, so `use_delta_t` is ignored on this state",
                field="use_delta_t",
            )
        )

    mechanical_requested = bool(mechanical_design)
    for name, stated in (
        ("mechanical_type", mechanical_type is not None),
        ("mechanical_selection", mechanical_selection is not None),
    ):
        if stated and not mechanical_requested:
            raise InvalidInputError(
                name, "a mechanical-design input is stated without `mechanical_design`"
            )
    if mechanical_type is not None and mechanical_type not in MECHANICAL_NAMES:
        raise InvalidInputError(
            "mechanical_type",
            f"{mechanical_type} is not one of the five mechanical configurations",
        )
    criterion = "min_area" if mechanical_selection is None else mechanical_selection
    if criterion not in MECHANICAL_CRITERIA:
        raise InvalidInputError(
            "mechanical_selection",
            f"{mechanical_selection} is not one of min_area, min_weight or min_pressure_drop",
        )

    hot_n = input_to_si(spec, "hot_in_n", hot_in_n)
    cold_n = input_to_si(spec, "cold_in_n", cold_in_n)
    # A zero flow is not a state, and every duty-per-flow divides by one.
    for field, flow in (("hot_in_n", hot_n), ("cold_in_n", cold_n)):
        if flow <= 0.0:
            raise OutOfRangeError(field, flow, "an exchanger side needs a positive molar flow")

    # **The transient gate.** `time_step`'s presence is the caller's statement that this is
    # one `runTransient` step; the class's `dynamicModelEnabled`, `wallMass > 0` and
    # `heatTransferArea > 0` then choose between the enabled branch and the steady fallback.
    transient_requested = time_step_si is not None
    if not transient_requested:
        for name, stated in (
            ("dynamic_model_enabled", dynamic_model_enabled is not None),
            ("transient_tube_side", transient_tube_side is not None),
            ("wall_mass", wall_mass_si is not None),
            ("wall_cp", wall_cp_si is not None),
            ("wall_temperature", wall_temperature_si is not None),
            ("shell_holdup_volume", shell_holdup_si is not None),
            ("tube_holdup_volume", tube_holdup_si is not None),
            ("shell_side_htc", shell_side_htc_si is not None),
            ("tube_side_htc", tube_side_htc_si is not None),
            ("heat_transfer_area", area_si is not None),
            ("shell_fluid_temperature", shell_fluid_temperature_si is not None),
            ("tube_fluid_temperature", tube_fluid_temperature_si is not None),
        ):
            if stated:
                raise InvalidInputError(name, "a transient input is stated without `time_step`")
    if transient_tube_side not in (None, "hot", "cold"):
        raise InvalidInputError(
            "transient_tube_side",
            f"{transient_tube_side} is not one of hot or cold",
        )
    if transient_requested and transient_tube_side is None:
        raise InvalidInputError(
            "transient_tube_side",
            "a stated `time_step` needs the side that is the tube",
        )
    for name, value in (
        ("time_step", time_step_si),
        ("wall_mass", wall_mass_si),
        ("wall_cp", wall_cp_si),
        ("shell_holdup_volume", shell_holdup_si),
        ("tube_holdup_volume", tube_holdup_si),
        ("shell_side_htc", shell_side_htc_si),
        ("tube_side_htc", tube_side_htc_si),
        ("heat_transfer_area", area_si),
    ):
        if value is not None and not math.isfinite(value):
            raise InvalidInputError(name, f"a finite value is required, and {value} is not")
    for name, value in (
        ("wall_temperature", wall_temperature_si),
        ("shell_fluid_temperature", shell_fluid_temperature_si),
        ("tube_fluid_temperature", tube_fluid_temperature_si),
    ):
        if value is not None and math.isinf(value):
            raise InvalidInputError(
                name,
                f"a carried temperature is finite or NaN, and {value} is neither",
            )
    wall_mass_kg = 0.0 if wall_mass_si is None else max(0.0, wall_mass_si)
    holdup_shell = 0.0 if shell_holdup_si is None else max(0.0, shell_holdup_si)
    holdup_tube = 0.0 if tube_holdup_si is None else max(0.0, tube_holdup_si)
    area_m2 = 0.0 if area_si is None else max(0.0, area_si)
    active_transient = (
        transient_requested and bool(dynamic_model_enabled) and wall_mass_kg > 0.0 and area_m2 > 0.0
    )

    diagnostics: RatingDiagnostics | None = None
    transient_state: TransientState | None = None
    hot_p_si = input_to_si(spec, "hot_in_p", hot_in_p)
    cold_p_si = input_to_si(spec, "cold_in_p", cold_in_p)
    if active_transient:
        for name, stated in (
            ("hot_outlet_temperature", hot_out_t_si is not None),
            ("cold_outlet_temperature", cold_out_t_si is not None),
            ("use_delta_t", bool(use_delta_t)),
            ("use_rating_pressure_drop", bool(use_rating_pressure_drop)),
        ):
            if stated:
                raise InvalidInputError(
                    name, "the dynamic branch has no steady outlet to pin or rate"
                )
        time_step_value = 0.0 if time_step_si is None else time_step_si
        states = _transient_states(
            hot_components,
            hot_n,
            hot_in_z,
            hot_p_si,
            hot_t_si,
            cold_components,
            cold_n,
            cold_in_z,
            cold_p_si,
            cold_t_si,
            time_step_value,
            wall_mass_kg,
            500.0 if wall_cp_si is None else wall_cp_si,
            wall_temperature_si,
            holdup_shell,
            holdup_tube,
            500.0 if shell_side_htc_si is None else shell_side_htc_si,
            1000.0 if tube_side_htc_si is None else tube_side_htc_si,
            area_m2,
            shell_fluid_temperature_si,
            tube_fluid_temperature_si,
            tube_is_hot=transient_tube_side == "hot",
        )
        transient_state = states.transient
    elif use_rating_pressure_drop:
        if hot_out_t_si is not None or cold_out_t_si is not None or use_delta_t:
            raise InvalidInputError(
                "use_rating_pressure_drop",
                "the thermal-hydraulic rating requires a UA specification without a pinned "
                "outlet or delta-T",
            )
        if rating_area is None:
            raise InvalidInputError(
                "rating_area", "the thermal-hydraulic rating needs a positive outside area"
            )
        if rating_tube_side not in ("hot", "cold"):
            raise InvalidInputError(
                "rating_tube_side",
                "the thermal-hydraulic rating needs the side that is the tube",
            )
        geometry = {
            "tube_outside_diameter": _optional_si(
                spec, "rating_tube_outside_diameter", rating_tube_outside_diameter
            ),
            "tube_inside_diameter": _optional_si(
                spec, "rating_tube_inside_diameter", rating_tube_inside_diameter
            ),
            "tube_length": _optional_si(spec, "rating_tube_length", rating_tube_length),
            "tube_count": rating_tube_count,
            "tube_passes": rating_tube_passes,
            "tube_pitch": _optional_si(spec, "rating_tube_pitch", rating_tube_pitch),
            "triangular_pitch": rating_triangular_pitch,
            "shell_inside_diameter": _optional_si(
                spec, "rating_shell_inside_diameter", rating_shell_inside_diameter
            ),
            "baffle_spacing": _optional_si(spec, "rating_baffle_spacing", rating_baffle_spacing),
            "baffle_count": rating_baffle_count,
            "baffle_cut": rating_baffle_cut,
            "tube_wall_conductivity": _optional_si(
                spec, "rating_tube_wall_conductivity", rating_tube_wall_conductivity
            ),
            "fouling_tube": rating_fouling_tube,
            "fouling_shell": rating_fouling_shell,
            "tube_fouling_layer_thickness": _optional_si(
                spec, "rating_tube_fouling_layer_thickness", rating_tube_fouling_layer_thickness
            ),
            "tube_fouling_layer_conductivity": _optional_si(
                spec,
                "rating_tube_fouling_layer_conductivity",
                rating_tube_fouling_layer_conductivity,
            ),
            "shell_side_method": rating_shell_side_method,
            "shell_viscosity_wall": _optional_si(
                spec, "rating_shell_viscosity_wall", rating_shell_viscosity_wall
            ),
            "has_sealing": rating_has_sealing,
            "sealing_pairs": rating_sealing_pairs,
            "tube_to_baffle_clearance": _optional_si(
                spec, "rating_tube_to_baffle_clearance", rating_tube_to_baffle_clearance
            ),
            "shell_to_baffle_clearance": _optional_si(
                spec, "rating_shell_to_baffle_clearance", rating_shell_to_baffle_clearance
            ),
            "bypass_area": _optional_si(spec, "rating_bypass_area", rating_bypass_area),
        }
        states, diagnostics = _rating_states(
            hot_components,
            hot_n,
            hot_in_z,
            input_to_si(spec, "hot_in_p", hot_in_p),
            hot_t_si,
            cold_components,
            cold_n,
            cold_in_z,
            input_to_si(spec, "cold_in_p", cold_in_p),
            cold_t_si,
            rating_tube_side == "hot",
            input_to_si(spec, "rating_area", rating_area),
            geometry,
            flow_arrangement,
        )
    else:
        states = _states(
            hot_components,
            hot_n,
            hot_in_z,
            input_to_si(spec, "hot_in_p", hot_in_p),
            hot_t_si,
            cold_components,
            cold_n,
            cold_in_z,
            input_to_si(spec, "cold_in_p", cold_in_p),
            cold_t_si,
            ua_si,
            flow_arrangement,
            hot_out_t_si,
            cold_out_t_si,
            bool(use_delta_t),
            delta_t_si,
        )
    if transient_requested and not active_transient:
        # **The class's fields survive `run`**, so a steady fallback echoes what the caller
        # stated rather than inventing a state the step did not advance.
        transient_state = TransientState(
            wall_temperature=wall_temperature_si,
            shell_fluid_temperature=shell_fluid_temperature_si,
            tube_fluid_temperature=tube_fluid_temperature_si,
        )

    mechanical = None
    if mechanical_requested:
        pinned = hot_out_t_si is not None or cold_out_t_si is not None
        # **The class's `getDuty()` basis, not this result's `duty`**: the pinned branch
        # writes neither, so `determineDuty` sees zero there.
        duty_basis = 0.0 if pinned else abs(states.duty)
        ua_basis = (
            states.ua_value
            if states.ua_value is not None
            else (
                diagnostics.rating_ua
                if diagnostics is not None
                else (ua_si if ua_si is not None else 500.0)
            )
        )
        effectiveness = 0.0 if states.rating is None else states.rating.effectiveness
        mechanical = _two_stream_mechanical(
            hot_t_si,
            states.hot_out_t.to("K").magnitude,
            cold_t_si,
            states.cold_out_t.to("K").magnitude,
            duty_basis,
            ua_basis,
            effectiveness,
            mechanical_type,
            criterion,
        )

    return HeatExchangerResult(
        hot_out_n=from_si(hot_n, "mol/s"),
        hot_out_z=tuple(hot_in_z),
        hot_out_p=from_si(states.hot_out_p, "Pa"),
        hot_out_t=states.hot_out_t,
        hot_out_h=from_si(states.hot_out_h, "J/mol"),
        cold_out_n=from_si(cold_n, "mol/s"),
        cold_out_z=tuple(cold_in_z),
        cold_out_p=from_si(states.cold_out_p, "Pa"),
        cold_out_t=states.cold_out_t,
        cold_out_h=from_si(states.cold_out_h, "J/mol"),
        # **The hot side's release, and not `states.duty`.** That one is the *kept* side's own
        # change - NeqSim's `getDuty`, whose sign depends on which side `run` happened to keep -
        # so the same exchanger could report either sign for one physical direction of heat. This
        # is `n_hot * (h_hot_in - h_hot_out)`, positive whenever the hot side cools, which is the
        # frame a reader comparing two runs needs and the one the Rust half publishes.
        duty=from_si((states.hot_in_h - states.hot_out_h) * hot_n, "W"),
        ntu=None if states.rating is None else states.rating.ntu,
        effectiveness=None if states.rating is None else states.rating.effectiveness,
        c_min=None if states.rating is None else from_si(states.rating.c_min, "W/K"),
        c_max=None if states.rating is None else from_si(states.rating.c_max, "W/K"),
        capacity_ratio=None if states.rating is None else states.rating.capacity_ratio,
        ua_value=None if states.ua_value is None else from_si(states.ua_value, "W/K"),
        delta_t_side=states.delta_t_side,
        wall_temperature=(
            None
            if transient_state is None
            or transient_state.wall_temperature is None
            or math.isnan(transient_state.wall_temperature)
            else from_si(transient_state.wall_temperature, "K")
        ),
        shell_fluid_temperature=(
            None
            if transient_state is None
            or transient_state.shell_fluid_temperature is None
            or math.isnan(transient_state.shell_fluid_temperature)
            else from_si(transient_state.shell_fluid_temperature, "K")
        ),
        tube_fluid_temperature=(
            None
            if transient_state is None
            or transient_state.tube_fluid_temperature is None
            or math.isnan(transient_state.tube_fluid_temperature)
            else from_si(transient_state.tube_fluid_temperature, "K")
        ),
        mechanical_selected_type=(
            None if mechanical is None else mechanical.selected.mechanical_type
        ),
        mechanical_selected_area=None if mechanical is None else mechanical.selected.required_area,
        log_mean_temperature_difference=(
            None if mechanical is None else from_si(mechanical.lmtd, "K")
        ),
        mechanical_approach_temperature=(
            None if mechanical is None else from_si(mechanical.approach, "K")
        ),
        used_overall_heat_transfer_coefficient=(
            None if mechanical is None else from_si(mechanical.used_u, "W/(m**2*K)")
        ),
        calculated_ua=None if mechanical is None else from_si(mechanical.ua, "W/K"),
        rating_u=None if diagnostics is None else from_si(diagnostics.rating_u, "W/(m**2*K)"),
        rating_ua=None if diagnostics is None else from_si(diagnostics.rating_ua, "W/K"),
        rating_area=None if diagnostics is None else from_si(diagnostics.rating_area, "m**2"),
        tube_side_pressure_drop=(
            None if diagnostics is None else from_si(diagnostics.tube_side_pressure_drop, "Pa")
        ),
        shell_side_pressure_drop=(
            None if diagnostics is None else from_si(diagnostics.shell_side_pressure_drop, "Pa")
        ),
        tube_side_heat_transfer_coefficient=(
            None
            if diagnostics is None
            else from_si(diagnostics.tube_side_heat_transfer_coefficient, "W/(m**2*K)")
        ),
        shell_side_heat_transfer_coefficient=(
            None
            if diagnostics is None
            else from_si(diagnostics.shell_side_heat_transfer_coefficient, "W/(m**2*K)")
        ),
        tube_side_reynolds=None if diagnostics is None else diagnostics.tube_side_reynolds,
        shell_side_reynolds=None if diagnostics is None else diagnostics.shell_side_reynolds,
        tube_side_velocity=(
            None if diagnostics is None else from_si(diagnostics.tube_side_velocity, "m/s")
        ),
        shell_side_velocity=(
            None if diagnostics is None else from_si(diagnostics.shell_side_velocity, "m/s")
        ),
        jc=None if diagnostics is None else diagnostics.jc,
        jl=None if diagnostics is None else diagnostics.jl,
        jb=None if diagnostics is None else diagnostics.jb,
        js=None if diagnostics is None else diagnostics.js,
        jr=None if diagnostics is None else diagnostics.jr,
        mechanical_types=() if mechanical is None else mechanical.types,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class Rating:
    """The rating's own numbers, in the order ``run`` computes them, in SI.

    ``C`` is a total over a temperature difference, so W/K; that is what makes ``NTU``
    dimensionless. ``NTU`` is a package-private field of ``HeatExchanger`` with no getter, so
    a capture cannot carry it - but ``duty`` and ``effectiveness`` are two numbers it does
    carry, and those two pin this record: the kept side is the larger swing, so
    ``C_max = duty / (effectiveness * span)``, and the effectiveness relation then leaves
    ``C_min`` as the only unknown. A capacity this got wrong cannot hide behind the other.
    """

    hot_capacity: float
    cold_capacity: float
    c_min: float
    c_max: float
    capacity_ratio: float
    ntu: float
    effectiveness: float


@dataclass(frozen=True, slots=True)
class TransientState:
    """The three temperatures one dynamic step hands to the next, in K.

    ``None`` means the class's field is NaN, which is its first-step sentinel and not a
    temperature; a fallback echoes what the caller stated.
    """

    wall_temperature: float | None
    shell_fluid_temperature: float | None
    tube_fluid_temperature: float | None


@dataclass(frozen=True, slots=True)
class ExchangerStates:
    """Both sides' states on the way through, and the rating's interior, in SI."""

    hot_in_h: float
    cold_in_h: float
    hot_out_h: float
    cold_out_h: float
    hot_out_t: Q
    cold_out_t: Q
    #: Each outlet's pressure, Pa. The ordinary and pinned modes keep the inlet's; the
    #: thermal-hydraulic rating reduces it by the correlated drop.
    hot_out_p: float
    cold_out_p: float
    #: The kept side's own total enthalpy change, W - negative when that side cools. It is
    #: what ``run`` publishes as ``getDuty``, and the pinned branch reports the pinned
    #: side's change the same way.
    duty: float
    #: ``None`` when an outlet was pinned instead: that branch does not rate anything.
    rating: Rating | None
    #: ``runDeltaT``'s computed conductance, W/K; ``None`` outside that mode.
    ua_value: float | None
    #: Which side ``runDeltaT`` set; ``None`` outside that mode.
    delta_t_side: str | None
    #: The carried state after a ``runTransient`` step; ``None`` on every steady branch.
    transient: TransientState | None = None


def _states(
    hot_components: list[str],
    hot_in_n: float,
    hot_in_z: list[float],
    hot_in_p: float,
    hot_in_t: float,
    cold_components: list[str],
    cold_in_n: float,
    cold_in_z: list[float],
    cold_in_p: float,
    cold_in_t: float,
    ua: float | None,
    flow_arrangement: str,
    hot_outlet_temperature: float | None,
    cold_outlet_temperature: float | None,
    use_delta_t: bool,
    delta_t: float,
) -> ExchangerStates:
    """The exchanger's intermediates, in the order it computes them, in SI.

    **One arithmetic, two consumers.** :func:`heat_exchanger` builds the result from this and
    :mod:`azoth.process.layers` reports the same numbers against a NeqSim capture, so a layer
    the harness compares is a layer the model actually took. The refusals are here rather
    than in the caller for the same reason the Rust kernel holds them: they are about which
    branch runs, and a caller that had to know would be a second copy of that decision.
    """
    if hot_outlet_temperature is not None and cold_outlet_temperature is not None:
        raise InvalidInputError(
            "hot_outlet_temperature",
            "an exchanger pins one outlet at most: with both stated the other side has "
            "nothing left to solve for",
        )
    pinned = hot_outlet_temperature is not None or cold_outlet_temperature is not None
    if not pinned and not use_delta_t and ua is None:
        raise InvalidInputError(
            "ua",
            "the rating needs an overall conductance. NeqSim's `UAvalue` defaults to "
            "500 W/K, and this refuses to guess it: a default that only changes the "
            "answer is not one a caller can see",
        )
    if flow_arrangement not in FLOW_ARRANGEMENTS:
        raise InvalidInputError(
            "flow_arrangement",
            f"{flow_arrangement} is not one of {list(FLOW_ARRANGEMENTS)}. NeqSim's `run` "
            f"falls through to the counterflow relation for an arrangement it does not "
            f"know, which would make a misspelling a plausible answer",
        )
    span = abs(hot_in_t - cold_in_t)
    if not pinned and span <= 2.220446049250313e-16 * max(abs(hot_in_t), abs(cold_in_t)):
        raise InvalidInputError(
            "hot_in_t",
            f"both inlets are at {hot_in_t} K, so there is no driving force for the "
            f"rating to size against",
        )

    hot_mixture, hot_gas = _components.mixture_of(hot_components, eos="pr")
    cold_mixture, cold_gas = _components.mixture_of(cold_components, eos="pr")
    hot_h, _ = enthalpy_at(hot_mixture, hot_gas, hot_in_t, hot_in_p, hot_in_z)
    cold_h, _ = enthalpy_at(cold_mixture, cold_gas, cold_in_t, cold_in_p, cold_in_z)

    rating: Rating | None = None
    ua_value: float | None = None
    delta_t_side: str | None = None
    if pinned:
        hot_pinned = hot_outlet_temperature is not None
        pinned_t = hot_outlet_temperature if hot_pinned else cold_outlet_temperature
        if pinned_t is None:
            raise InvalidInputError(
                "mode", "a pinned outlet was requested and neither outlet temperature was given"
            )
        if hot_pinned:
            specified, _ = enthalpy_at(hot_mixture, hot_gas, pinned_t, hot_in_p, hot_in_z)
            duty = hot_in_n * (specified - hot_h)
            hot_out_h, cold_out_h = specified, cold_h - divide(duty, cold_in_n)
        else:
            specified, _ = enthalpy_at(cold_mixture, cold_gas, pinned_t, cold_in_p, cold_in_z)
            duty = cold_in_n * (specified - cold_h)
            cold_out_h, hot_out_h = specified, hot_h - divide(duty, hot_in_n)
    elif use_delta_t:
        if delta_t <= 0.0:
            raise InvalidInputError(
                "delta_t",
                f"`delta_t` is a positive temperature step, and {delta_t} K is not",
            )
        sign = 1.0 if hot_in_t > cold_in_t else -1.0
        hot_seeded, _ = enthalpy_at(
            hot_mixture, hot_gas, cold_in_t - sign * delta_t, hot_in_p, hot_in_z
        )
        cold_seeded, _ = enthalpy_at(
            cold_mixture, cold_gas, hot_in_t + sign * delta_t, cold_in_p, cold_in_z
        )
        hot_swing = hot_in_n * (hot_seeded - hot_h)
        cold_swing = cold_in_n * (cold_seeded - cold_h)
        hot_capacity = divide(abs(hot_swing), abs(cold_in_t - sign * delta_t - hot_in_t))
        cold_capacity = divide(abs(cold_swing), abs(hot_in_t + sign * delta_t - cold_in_t))
        _ = (hot_capacity, cold_capacity)
        if abs(cold_swing) < abs(hot_swing):
            d_h = hot_swing
            hot_out_h, cold_out_h = hot_seeded, cold_h - divide(d_h, cold_in_n)
            set_t, set_in_t = cold_in_t - sign * delta_t, hot_in_t
            delta_t_side = "hot"
        else:
            d_h = cold_swing
            cold_out_h, hot_out_h = cold_seeded, hot_h - divide(d_h, hot_in_n)
            set_t, set_in_t = hot_in_t + sign * delta_t, cold_in_t
            delta_t_side = "cold"
        ua_value = divide(d_h, set_t - set_in_t)
        duty = d_h
    else:
        if ua is None:
            raise InvalidInputError("ua", "the rating mode needs an overall conductance")
        # Each side flashed at the *other's* inlet temperature, which is `run`'s own
        # estimate of its capacity and is what makes the number of passes not matter.
        hot_seeded, _ = enthalpy_at(hot_mixture, hot_gas, cold_in_t, hot_in_p, hot_in_z)
        cold_seeded, _ = enthalpy_at(cold_mixture, cold_gas, hot_in_t, cold_in_p, cold_in_z)
        hot_swing = hot_in_n * (hot_seeded - hot_h)
        cold_swing = cold_in_n * (cold_seeded - cold_h)
        hot_capacity = divide(abs(hot_swing), span)
        cold_capacity = divide(abs(cold_swing), span)
        c_min = min(hot_capacity, cold_capacity)
        c_max = max(hot_capacity, cold_capacity)
        capacity_ratio = divide(c_min, c_max)
        ntu = divide(ua, c_min)
        effectiveness = _effectiveness(ntu, capacity_ratio, flow_arrangement)
        rating = Rating(
            hot_capacity=hot_capacity,
            cold_capacity=cold_capacity,
            c_min=c_min,
            c_max=c_max,
            capacity_ratio=capacity_ratio,
            ntu=ntu,
            effectiveness=effectiveness,
        )
        # The side whose seeded swing is the larger is the one `run` keeps, and the other
        # is energy-balanced against it.
        if abs(cold_swing) > abs(hot_swing):
            duty = effectiveness * hot_swing
            hot_out_h, cold_out_h = (
                hot_h + divide(duty, hot_in_n),
                cold_h - divide(duty, cold_in_n),
            )
        else:
            duty = effectiveness * cold_swing
            cold_out_h, hot_out_h = (
                cold_h + divide(duty, cold_in_n),
                hot_h - divide(duty, hot_in_n),
            )

    return ExchangerStates(
        hot_in_h=hot_h,
        cold_in_h=cold_h,
        hot_out_h=hot_out_h,
        cold_out_h=cold_out_h,
        hot_out_t=_temperature(hot_mixture, hot_gas, hot_in_p, hot_out_h, hot_in_z),
        cold_out_t=_temperature(cold_mixture, cold_gas, cold_in_p, cold_out_h, cold_in_z),
        hot_out_p=hot_in_p,
        cold_out_p=cold_in_p,
        duty=duty,
        rating=rating,
        ua_value=ua_value,
        delta_t_side=delta_t_side,
    )


def _transient_phase(
    components: list[str],
    n: float,
    z: list[float],
    pressure: float,
    temperature: float,
) -> tuple[float, float, float]:
    """A single-phase state's density, mass Cp and mass flow, without transport.

    **There is no transport here on purpose**: the wall and CSTR ODEs read only the inlet's
    density and mass heat capacity, and `_phase_properties`' transport step refuses a pure
    component the way NeqSim's effective-diffusivity denominator does.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flash = pt_flash_solve(mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), z)
    if flash.phase == Phase.TWO_PHASE:
        raise InvalidInputError(
            "stream", "the exchanger's transient step requires single-phase streams"
        )
    root = flash.z_vapour if flash.phase == Phase.ALL_VAPOUR else flash.z_liquid
    composition = list(z)
    molar_mass = sum(
        (component.molar_mass.to("kg/mol").magnitude if component.molar_mass else 0.0) * fraction
        for component, fraction in zip(mixture.components, composition, strict=True)
    )
    volume = (
        pr_molar_volume(root, quantity(temperature, "K"), quantity(pressure, "Pa"))
        .v.to("m**3/mol")
        .magnitude
    )
    density = (
        pr_mass_density(
            quantity(molar_mass, "kg/mol"),
            quantity(volume - mixture.volume_shift(composition), "m**3/mol"),
        )
        .rho.to("kg/m**3")
        .magnitude
    )
    cp_mass = (
        molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            composition,
            root,
        )
        .cp.to("J/(mol*K)")
        .magnitude
        / molar_mass
    )
    return density, cp_mass, n * molar_mass


def _transient_states(
    hot_components: list[str],
    hot_in_n: float,
    hot_in_z: list[float],
    hot_in_p: float,
    hot_in_t: float,
    cold_components: list[str],
    cold_in_n: float,
    cold_in_z: list[float],
    cold_in_p: float,
    cold_in_t: float,
    time_step: float,
    wall_mass: float,
    wall_cp: float,
    wall_temperature: float | None,
    shell_holdup: float,
    tube_holdup: float,
    shell_side_htc: float,
    tube_side_htc: float,
    area: float,
    shell_fluid_temperature: float | None,
    tube_fluid_temperature: float | None,
    *,
    tube_is_hot: bool,
) -> ExchangerStates:
    """One `HeatExchanger.runTransient` step of the enabled branch.

    The class's stream 0 is the shell and stream 1 the tube, the reverse of the rating; the
    outlet shift applies the heat **rate** as though it were a joule while the temperatures
    integrate by `time_step`.
    """
    hot_mixture, hot_gas = _components.mixture_of(hot_components, eos="pr")
    cold_mixture, cold_gas = _components.mixture_of(cold_components, eos="pr")
    hot_h, _ = enthalpy_at(hot_mixture, hot_gas, hot_in_t, hot_in_p, hot_in_z)
    cold_h, _ = enthalpy_at(cold_mixture, cold_gas, cold_in_t, cold_in_p, cold_in_z)
    hot = (hot_components, hot_in_n, hot_in_z, hot_in_p, hot_in_t, hot_mixture, hot_gas, hot_h)
    cold = (
        cold_components,
        cold_in_n,
        cold_in_z,
        cold_in_p,
        cold_in_t,
        cold_mixture,
        cold_gas,
        cold_h,
    )
    shell, tube = (cold, hot) if tube_is_hot else (hot, cold)
    shell_names, shell_n, shell_z, shell_p, shell_t, _shell_mix, _shell_gas, shell_h = shell
    tube_names, tube_n, tube_z, tube_p, tube_t, _tube_mix, _tube_gas, tube_h = tube

    # The refusal is independent of the holdup: a two-phase side has no one density or Cp.
    _transient_phase(shell_names, shell_n, shell_z, shell_p, shell_t)
    _transient_phase(tube_names, tube_n, tube_z, tube_p, tube_t)

    wall = (
        0.5 * (shell_t + tube_t)
        if wall_temperature is None or math.isnan(wall_temperature)
        else wall_temperature
    )

    if shell_holdup > 0.0:
        shell_state = (
            shell_t
            if shell_fluid_temperature is None or math.isnan(shell_fluid_temperature)
            else shell_fluid_temperature
        )
        shell_density, shell_cp, shell_mass_flow = _transient_phase(
            shell_names, shell_n, shell_z, shell_p, shell_t
        )
        fluid_mass = shell_density * shell_holdup
        if fluid_mass > 0.0 and shell_cp > 0.0:
            mixing = shell_mass_flow * shell_cp * (shell_t - shell_state)
            to_wall = shell_side_htc * area * (shell_state - wall)
            shell_state += (mixing - to_wall) / (fluid_mass * shell_cp) * time_step
        shell_for_wall = shell_state
        shell_out_state = shell_state
    else:
        shell_for_wall = shell_t
        shell_out_state = None

    if tube_holdup > 0.0:
        tube_state = (
            tube_t
            if tube_fluid_temperature is None or math.isnan(tube_fluid_temperature)
            else tube_fluid_temperature
        )
        tube_density, tube_cp, tube_mass_flow = _transient_phase(
            tube_names, tube_n, tube_z, tube_p, tube_t
        )
        fluid_mass = tube_density * tube_holdup
        if fluid_mass > 0.0 and tube_cp > 0.0:
            mixing = tube_mass_flow * tube_cp * (tube_t - tube_state)
            wall_to_fluid = tube_side_htc * area * (wall - tube_state)
            tube_state += (mixing + wall_to_fluid) / (fluid_mass * tube_cp) * time_step
        tube_for_wall = tube_state
        tube_out_state = tube_state
    else:
        tube_for_wall = tube_t
        tube_out_state = None

    q_shell_to_wall = shell_side_htc * area * (shell_for_wall - wall)
    q_wall_to_tube = tube_side_htc * area * (wall - tube_for_wall)
    wall += (q_shell_to_wall - q_wall_to_tube) / (wall_mass * wall_cp) * time_step

    shell_out_h = shell_h - divide(q_shell_to_wall, shell_n)
    tube_out_h = tube_h + divide(q_wall_to_tube, tube_n)
    if tube_is_hot:
        hot_out_h, cold_out_h = tube_out_h, shell_out_h
    else:
        hot_out_h, cold_out_h = shell_out_h, tube_out_h

    return ExchangerStates(
        hot_in_h=hot_h,
        cold_in_h=cold_h,
        hot_out_h=hot_out_h,
        cold_out_h=cold_out_h,
        hot_out_t=_temperature(hot_mixture, hot_gas, hot_in_p, hot_out_h, hot_in_z),
        cold_out_t=_temperature(cold_mixture, cold_gas, cold_in_p, cold_out_h, cold_in_z),
        hot_out_p=hot_in_p,
        cold_out_p=cold_in_p,
        duty=q_shell_to_wall,
        rating=None,
        ua_value=None,
        delta_t_side=None,
        transient=TransientState(
            wall_temperature=wall,
            shell_fluid_temperature=shell_out_state,
            tube_fluid_temperature=tube_out_state,
        ),
    )


@dataclass(frozen=True, slots=True)
class RatingDiagnostics:
    """The rating calculator's quantities, in SI, as the equipment result publishes them."""

    rating_u: float
    rating_ua: float
    rating_area: float
    tube_side_pressure_drop: float
    shell_side_pressure_drop: float
    tube_side_heat_transfer_coefficient: float
    shell_side_heat_transfer_coefficient: float
    tube_side_reynolds: float
    shell_side_reynolds: float
    tube_side_velocity: float
    shell_side_velocity: float
    jc: float
    jl: float
    jb: float
    js: float
    jr: float


def _phase_properties(
    components: list[str],
    n: float,
    z: list[float],
    pressure: float,
    temperature: float,
) -> dict[str, float | str]:
    """The five properties and the label a single-phase state answers.

    The class's `requireSinglePhaseRatingState` refuses a two-phase inlet, so this does too;
    the phase's composition is the feed's for a single phase, which is what the Rust
    `phase_view` reads.
    """
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flash = pt_flash_solve(mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), z)
    if flash.phase == Phase.TWO_PHASE:
        raise InvalidInputError(
            "stream", "the thermal-hydraulic rating requires single-phase streams"
        )
    root = flash.z_vapour if flash.phase == Phase.ALL_VAPOUR else flash.z_liquid
    composition = list(z)
    reduced = reduced_parameters(mixture, temperature, pressure)
    label = phase_label(mixture, reduced, composition, root)
    molar_mass = sum(
        (component.molar_mass.to("kg/mol").magnitude if component.molar_mass else 0.0) * fraction
        for component, fraction in zip(mixture.components, composition, strict=True)
    )
    volume = (
        pr_molar_volume(root, quantity(temperature, "K"), quantity(pressure, "Pa"))
        .v.to("m**3/mol")
        .magnitude
    )
    density = (
        pr_mass_density(
            quantity(molar_mass, "kg/mol"),
            quantity(volume - mixture.volume_shift(composition), "m**3/mol"),
        )
        .rho.to("kg/m**3")
        .magnitude
    )
    cp_mass = (
        molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            composition,
            root,
        )
        .cp.to("J/(mol*K)")
        .magnitude
        / molar_mass
    )
    transport = phase_transport(
        components, label, quantity(temperature, "K"), quantity(pressure, "Pa"), composition
    )
    return {
        "density": density,
        "cp_mass": cp_mass,
        "mu": transport.mu.to("Pa*s").magnitude,
        "k": transport.k.to("W/(m*K)").magnitude,
        "molar_mass": molar_mass,
        "mass_flow": n * molar_mass,
        "kind": label,
    }


def _rating_states(
    hot_components: list[str],
    hot_in_n: float,
    hot_in_z: list[float],
    hot_in_p: float,
    hot_in_t: float,
    cold_components: list[str],
    cold_in_n: float,
    cold_in_z: list[float],
    cold_in_p: float,
    cold_in_t: float,
    tube_is_hot: bool,
    area: float,
    geometry: dict[str, float | bool | str | None],
    flow_arrangement: str,
) -> tuple[ExchangerStates, RatingDiagnostics]:
    """`runThermalHydraulicRating`: the correlated UA, both drops and the closed outlets."""
    if area <= 0.0 or not math.isfinite(area):
        raise InvalidInputError(
            "rating_area", "the thermal-hydraulic rating needs a positive outside area"
        )
    if tube_is_hot:
        tube_components, tube_n, tube_z, tube_p, tube_t = (
            hot_components,
            hot_in_n,
            hot_in_z,
            hot_in_p,
            hot_in_t,
        )
        shell_components, shell_n, shell_z, shell_p, shell_t = (
            cold_components,
            cold_in_n,
            cold_in_z,
            cold_in_p,
            cold_in_t,
        )
        tube_mixture, tube_gas = _components.mixture_of(hot_components, eos="pr")
        shell_mixture, shell_gas = _components.mixture_of(cold_components, eos="pr")
    else:
        tube_components, tube_n, tube_z, tube_p, tube_t = (
            cold_components,
            cold_in_n,
            cold_in_z,
            cold_in_p,
            cold_in_t,
        )
        shell_components, shell_n, shell_z, shell_p, shell_t = (
            hot_components,
            hot_in_n,
            hot_in_z,
            hot_in_p,
            hot_in_t,
        )
        tube_mixture, tube_gas = _components.mixture_of(cold_components, eos="pr")
        shell_mixture, shell_gas = _components.mixture_of(hot_components, eos="pr")

    hot_mixture, hot_gas = (tube_mixture, tube_gas) if tube_is_hot else (shell_mixture, shell_gas)
    cold_mixture, cold_gas = (shell_mixture, shell_gas) if tube_is_hot else (tube_mixture, tube_gas)
    hot_in_h, _ = enthalpy_at(hot_mixture, hot_gas, hot_in_t, hot_in_p, hot_in_z)
    cold_in_h, _ = enthalpy_at(cold_mixture, cold_gas, cold_in_t, cold_in_p, cold_in_z)
    tube_h, _ = enthalpy_at(tube_mixture, tube_gas, tube_t, tube_p, tube_z)
    shell_h, _ = enthalpy_at(shell_mixture, shell_gas, shell_t, shell_p, shell_z)

    tube = _phase_properties(tube_components, tube_n, tube_z, tube_p, tube_t)
    shell = _phase_properties(shell_components, shell_n, shell_z, shell_p, shell_t)

    method = "kern" if geometry["shell_side_method"] is None else str(geometry["shell_side_method"])
    calculator = _rating_calculator(
        _length(geometry, "tube_outside_diameter", 0.01905),
        _length(geometry, "tube_inside_diameter", 0.01483),
        _length(geometry, "tube_length", 6.096),
        100.0 if geometry["tube_count"] is None else float(geometry["tube_count"]),
        2.0 if geometry["tube_passes"] is None else float(geometry["tube_passes"]),
        _length(geometry, "tube_pitch", 0.02381),
        True if geometry["triangular_pitch"] is None else bool(geometry["triangular_pitch"]),
        _length(geometry, "shell_inside_diameter", 0.5),
        _length(geometry, "baffle_spacing", 0.2),
        10.0 if geometry["baffle_count"] is None else float(geometry["baffle_count"]),
        0.25 if geometry["baffle_cut"] is None else float(geometry["baffle_cut"]),
        from_si(
            52.0
            if geometry["tube_wall_conductivity"] is None
            else float(geometry["tube_wall_conductivity"]),
            "W/(m*K)",
        ),
        from_si(float(tube["density"]), "kg/m**3"),
        from_si(float(tube["mu"]), "Pa*s"),
        from_si(float(tube["cp_mass"]), "J/(kg*K)"),
        from_si(float(tube["k"]), "W/(m*K)"),
        from_si(float(tube["mass_flow"]), "kg/s"),
        from_si(float(shell["density"]), "kg/m**3"),
        from_si(float(shell["mu"]), "Pa*s"),
        from_si(float(shell["cp_mass"]), "J/(kg*K)"),
        from_si(float(shell["k"]), "W/(m*K)"),
        from_si(float(shell["mass_flow"]), "kg/s"),
        0.000176 if geometry["fouling_tube"] is None else float(geometry["fouling_tube"]),
        0.000176 if geometry["fouling_shell"] is None else float(geometry["fouling_shell"]),
        method,
        None
        if geometry["tube_fouling_layer_thickness"] is None
        else from_si(float(geometry["tube_fouling_layer_thickness"]), "m"),
        None
        if geometry["tube_fouling_layer_conductivity"] is None
        else from_si(float(geometry["tube_fouling_layer_conductivity"]), "W/(m*K)"),
        None
        if geometry["shell_viscosity_wall"] is None
        else from_si(float(geometry["shell_viscosity_wall"]), "Pa*s"),
        None if geometry["has_sealing"] is None else bool(geometry["has_sealing"]),
        None if geometry["sealing_pairs"] is None else float(geometry["sealing_pairs"]),
        None
        if geometry["tube_to_baffle_clearance"] is None
        else from_si(float(geometry["tube_to_baffle_clearance"]), "m"),
        None
        if geometry["shell_to_baffle_clearance"] is None
        else from_si(float(geometry["shell_to_baffle_clearance"]), "m"),
        None
        if geometry["bypass_area"] is None
        else from_si(float(geometry["bypass_area"]), "m**2"),
    )
    rating_u = calculator.overall_heat_transfer_coefficient.to("W/(m**2*K)").magnitude
    ua = rating_u * area
    tube_dp = calculator.tube_side_pressure_drop.to("Pa").magnitude
    shell_dp = calculator.shell_side_pressure_drop.to("Pa").magnitude
    tube_p_out = tube_p - tube_dp
    shell_p_out = shell_p - shell_dp
    for field, pressure in (("tube", tube_p_out), ("shell", shell_p_out)):
        if not math.isfinite(pressure) or pressure <= 0.0:
            raise InvalidInputError(
                field, "the correlated pressure drop exhausts the inlet pressure"
            )

    span = tube_t - shell_t
    tube_capacity = shell_capacity = 0.0
    if abs(span) <= 1e-10:
        ntu = effectiveness_value = 0.0
        c_min = c_max = capacity_ratio = heat = 0.0
    else:
        tube_limit = _phase_properties(tube_components, tube_n, tube_z, tube_p, shell_t)
        shell_limit = _phase_properties(shell_components, shell_n, shell_z, shell_p, tube_t)
        for side, limit, original in (
            ("tube", tube_limit, tube),
            ("shell", shell_limit, shell),
        ):
            if limit["kind"] != original["kind"]:
                raise InvalidInputError(
                    side, "the rating temperature interval crosses a phase transition"
                )
        tube_limit_h, _ = enthalpy_at(tube_mixture, tube_gas, shell_t, tube_p, tube_z)
        shell_limit_h, _ = enthalpy_at(shell_mixture, shell_gas, tube_t, shell_p, shell_z)
        tube_capacity = divide(abs(tube_n * (tube_limit_h - tube_h)), abs(span))
        shell_capacity = divide(abs(shell_n * (shell_limit_h - shell_h)), abs(span))
        for field, capacity in (("tube", tube_capacity), ("shell", shell_capacity)):
            if not math.isfinite(capacity) or capacity <= 0.0:
                raise InvalidInputError(
                    field,
                    "the thermal-hydraulic rating requires positive finite heat capacity rates",
                )
        c_min = min(tube_capacity, shell_capacity)
        c_max = max(tube_capacity, shell_capacity)
        capacity_ratio = divide(c_min, c_max)
        ntu = divide(ua, c_min)
        effectiveness_value = (
            ntu / (1.0 + ntu)
            if abs(capacity_ratio - 1.0) < 1e-10 and flow_arrangement != "parallelflow"
            else _effectiveness(ntu, capacity_ratio, flow_arrangement)
        )
        heat = effectiveness_value * c_min * span
        if not math.isfinite(heat) or not math.isfinite(effectiveness_value):
            raise InvalidInputError(
                "rating", "the thermal-hydraulic rating produced an invalid effectiveness or duty"
            )

    if abs(span) <= 1e-10:
        tube_target, shell_target = tube_h, shell_h
    else:
        tube_target = tube_h - divide(heat, tube_n)
        shell_target = shell_h + divide(heat, shell_n)
    tube_solved = ph_flash_solve(
        tube_mixture, tube_gas, from_si(tube_p_out, "Pa"), from_si(tube_target, "J/mol"), tube_z
    )
    shell_solved = ph_flash_solve(
        shell_mixture,
        shell_gas,
        from_si(shell_p_out, "Pa"),
        from_si(shell_target, "J/mol"),
        shell_z,
    )
    for side, solved, mixture, gas, pressure, z, target in (
        ("tube", tube_solved, tube_mixture, tube_gas, tube_p_out, tube_z, tube_target),
        ("shell", shell_solved, shell_mixture, shell_gas, shell_p_out, shell_z, shell_target),
    ):
        reached, _ = enthalpy_at(mixture, gas, solved.T.to("K").magnitude, pressure, z)
        if abs(reached - target) > 1e-5 * max(1.0, abs(target)):
            raise InvalidInputError(side, "the rating outlet PH flash missed its energy closure")

    if tube_is_hot:
        hot_out_h, cold_out_h = tube_target, shell_target
        hot_out_t, cold_out_t = tube_solved.T, shell_solved.T
    else:
        hot_out_h, cold_out_h = shell_target, tube_target
        hot_out_t, cold_out_t = shell_solved.T, tube_solved.T

    diagnostics = RatingDiagnostics(
        rating_u=rating_u,
        rating_ua=ua,
        rating_area=area,
        tube_side_pressure_drop=tube_dp,
        shell_side_pressure_drop=shell_dp,
        tube_side_heat_transfer_coefficient=calculator.tube_side_heat_transfer_coefficient.to(
            "W/(m**2*K)"
        ).magnitude,
        shell_side_heat_transfer_coefficient=calculator.shell_side_heat_transfer_coefficient.to(
            "W/(m**2*K)"
        ).magnitude,
        tube_side_reynolds=calculator.tube_side_reynolds,
        shell_side_reynolds=calculator.shell_side_reynolds,
        tube_side_velocity=calculator.tube_side_velocity.to("m/s").magnitude,
        shell_side_velocity=calculator.shell_side_velocity.to("m/s").magnitude,
        jc=calculator.jc,
        jl=calculator.jl,
        jb=calculator.jb,
        js=calculator.js,
        jr=calculator.jr,
    )
    return (
        ExchangerStates(
            hot_in_h=hot_in_h,
            cold_in_h=cold_in_h,
            hot_out_h=hot_out_h,
            cold_out_h=cold_out_h,
            hot_out_t=hot_out_t,
            cold_out_t=cold_out_t,
            hot_out_p=tube_p_out if tube_is_hot else shell_p_out,
            cold_out_p=shell_p_out if tube_is_hot else tube_p_out,
            duty=heat,
            rating=Rating(
                hot_capacity=tube_capacity if tube_is_hot else shell_capacity,
                cold_capacity=shell_capacity if tube_is_hot else tube_capacity,
                c_min=c_min,
                c_max=c_max,
                capacity_ratio=capacity_ratio,
                ntu=ntu,
                effectiveness=effectiveness_value,
            ),
            ua_value=None,
            delta_t_side=None,
        ),
        diagnostics,
    )


def _length(geometry: dict[str, float | bool | str | None], name: str, default: float) -> Q:
    """One geometry length, as a quantity, with the calculator's own default."""
    value = geometry[name]
    return from_si(default if value is None else float(value), "m")


def _optional_si(spec: dict[str, object], name: str, value: Q | None) -> float | None:
    """An optional declared input in its SI magnitude, or ``None`` when it was not stated."""
    return None if value is None else input_to_si(spec, name, value)


@dataclass(frozen=True, slots=True)
class _Mechanical:
    """The two-stream thermal estimate, in SI."""

    lmtd: float
    approach: float
    used_u: float
    ua: float
    types: tuple[HeatExchangerMechanicalType, ...]
    selected: HeatExchangerMechanicalType


def _as_exchanger_row(row: HeaterMechanicalType) -> HeatExchangerMechanicalType:
    """The shared geometry model's row, under this model's generated record type."""
    return HeatExchangerMechanicalType(
        mechanical_type=row.mechanical_type,
        required_area=row.required_area,
        required_ua=row.required_ua,
        overall_heat_transfer_coefficient=row.overall_heat_transfer_coefficient,
        approach_temperature=row.approach_temperature,
        estimated_length=row.estimated_length,
        inner_diameter=row.inner_diameter,
        outer_diameter=row.outer_diameter,
        wall_thickness=row.wall_thickness,
        tube_count=row.tube_count,
        tube_passes=row.tube_passes,
        fin_surface_area=row.fin_surface_area,
        estimated_pressure_drop=row.estimated_pressure_drop,
        estimated_weight=row.estimated_weight,
        module_length=row.module_length,
        module_width=row.module_width,
        module_height=row.module_height,
    )


def _two_stream_mechanical(
    hot_in: float,
    hot_out: float,
    cold_in: float,
    cold_out: float,
    duty_w: float,
    ua_w_per_k: float,
    effectiveness: float,
    manual: str | None,
    criterion: str,
) -> _Mechanical:
    """``HeatExchangerMechanicalDesign.calcDesign``'s two-stream branch.

    **The sizing coefficient is the type's own typical U**, not the larger of it and a used
    coefficient as the heater branch takes; the used coefficient is the class's 500 and this
    branch never reads the exchanger's ``UAvalue`` for it.
    """
    duty = abs(duty_w)
    used_u = 500.0
    hot_mean = (hot_in + hot_out) / 2.0
    cold_mean = (cold_in + cold_out) / 2.0
    lmtd = None
    approach = None
    if not math.isnan(hot_mean) and not math.isnan(cold_mean):
        ordered = (
            (hot_in, hot_out, cold_in, cold_out)
            if hot_mean >= cold_mean
            else (cold_in, cold_out, hot_in, hot_out)
        )
        delta_t1 = abs(ordered[0] - ordered[3])
        delta_t2 = abs(ordered[1] - ordered[2])
        lmtd = _lmtd(delta_t1, delta_t2)
        approach = min(delta_t1, delta_t2)
    ua = ua_w_per_k if ua_w_per_k > 0.0 else None
    # **Carried although it cannot win through the declared inputs**: the rating that wrote
    # the effectiveness already produced the UA this compares against; a pinned or delta-T
    # run leaves the field at zero.
    if duty > 0.0 and effectiveness > 0.0 and effectiveness < 0.999:
        delta_max = abs(hot_in - cold_in)
        if delta_max > 0.0:
            c_min = duty / (effectiveness * delta_max)
            derived = c_min * -math.log(1.0 - effectiveness)
            if math.isfinite(derived) and derived > 0.0 and (ua is None or derived > ua):
                ua = derived
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
        area = max(ua_for_type / max(typical, 1.0e-6), 0.1)
        rows.append(_as_exchanger_row(_geometry(name, area, ua_for_type, type_approach)))
    if manual is not None:
        selected = next(row for row in rows if row.mechanical_type == manual)
    else:

        def metric(row: HeatExchangerMechanicalType) -> float:
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


def _temperature(
    mixture: object, ideal_gas: object, pressure: float, h: float, z: list[float]
) -> Q:
    """The temperature a side reaches at a shifted enthalpy."""
    solved = ph_flash_solve(
        mixture,  # type: ignore[arg-type]
        ideal_gas,  # type: ignore[arg-type]
        from_si(pressure, "Pa"),
        from_si(h, "J/mol"),
        z,
    )
    return solved.T


def _effectiveness(ntu: float, capacity_ratio: float, arrangement: str) -> float:
    """The effectiveness ``calcThermalEffectivenes`` carries, for the three arrangements."""
    if capacity_ratio == 0.0:
        return 1.0 - math.exp(-ntu)
    if arrangement == "parallelflow":
        return (1.0 - math.exp(-ntu * (1.0 + capacity_ratio))) / (1.0 + capacity_ratio)
    if arrangement == "counterflow" and capacity_ratio == 1.0:
        return ntu / (1.0 + ntu)
    exp = math.exp(-ntu * (1.0 - capacity_ratio))
    return (1.0 - exp) / (1.0 - capacity_ratio * exp)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.heat_exchanger")
