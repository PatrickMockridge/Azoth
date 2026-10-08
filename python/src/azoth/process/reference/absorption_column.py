"""``process.absorption_column`` - the tray absorber, as a registered id.

Spec: ``specs/models/process/absorption_column.toml``.

The Python twin of ``crates/azoth-process/src/models/absorption_column.rs``, whose arithmetic is
the column's own: `AbsorptionColumn extends DistillationColumn` and overrides **no `run`**, so
what this module adds is a shape - no condenser, no reboiler, the gas at stage 0 and the solvent
at the top stage - and the refusals.

**The class's own tests solve an isothermal column, and that is a one-sweep one.** They pin every
stage with `SimpleTray.setOutletTemperature`, which makes the base's own gate - the mean
tray-temperature change - exactly zero, so the solve stops after its first sweep. The capture
holds those rows as evidence and the unpinned columns as the oracles.
"""

from __future__ import annotations

from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AbsorptionColumnResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.process.kernels import Stream
from azoth.process.reference import _unported
from azoth.process.reference._column_stage import StreamRecord
from azoth.process.reference.capacity import (
    DEFAULT_INTERNAL_DIAMETER_M,
    DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER,
    fs_limits,
    gas_load_limits,
)
from azoth.process.reference.designer import (
    DEFAULT_TRAY_SPACING_M,
    DEFAULT_WEIR_HEIGHT_M,
    UNSIZED_COLUMN_DIAMETER_M,
)
from azoth.process.reference.distillation_column import (
    DEFAULT_DOWNCOMMER_AREA_FRACTION,
    DEFAULT_HOLE_AREA_FRACTION,
    DEFAULT_HOLE_DIAMETER_MM,
    UNPORTED_SOLVERS,
    Draws,
    _AbsorberMurphree,
    _feed,
    _Murphree,
    _pumparound_return_specification,
    _pumparound_returns_tear,
    _si,
    _side_draw_flow_specification,
    _side_draw_flow_tear,
    _States,
    _states,
    reactive_section,
)
from azoth.process.reference.mechanical import (
    DEFAULT_CONTACTOR_INTERNALS_TYPE,
    DEFAULT_MATERIAL_GRADE,
    DEFAULT_MAX_FLOODING_FACTOR,
    DEFAULT_MAX_OPERATION_PRESSURE_BARA,
    DEFAULT_TRAY_EFFICIENCY,
    DEFAULT_TRAY_TYPE,
    MechanicalGeometry,
    mechanical_design,
)

#: The class's own `DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR`.
DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR = 0.15


def absorption_column(
    gas_components: list[str],
    gas_n: Q,
    gas_z: list[float],
    gas_p: Q,
    gas_t: Q,
    solvent_components: list[str],
    solvent_n: Q,
    solvent_z: list[float],
    solvent_p: Q,
    solvent_t: Q,
    number_of_stages: int,
    top_pressure: Q,
    bottom_pressure: Q,
    tray_temperatures: list[float] | None = None,
    temperature_tolerance: float = 1.0e-2,
    max_iterations: int = 80,
    murphree_efficiency: float | None = None,
    component_murphree_efficiency: list[float] | None = None,
    max_allowable_gas_load_factor: float | None = None,
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    solver_type: str | None = None,
    tray_murphree_efficiency: list[float] | None = None,
    gas_side_draw_fractions: list[float] | None = None,
    liquid_side_draw_fractions: list[float] | None = None,
    pumparound_fractions: list[float] | None = None,
    side_draw_flow_tray: int | None = None,
    side_draw_flow_phase: str | None = None,
    side_draw_flow_target: Q | None = None,
    side_draw_flow_tolerance: float | None = None,
    side_draw_flow_max_iterations: int | None = None,
    pumparound_return_tray: int | None = None,
    pumparound_draw_tray: int | None = None,
    pumparound_draw_fraction: float | None = None,
    pumparound_temperature_drop: Q | None = None,
    pumparound_tolerance: float | None = None,
    pumparound_max_iterations: int | None = None,
    column_diameter: Q | None = None,
    max_allowable_fs_factor: float | None = None,
    tray_efficiency: float | None = None,
    max_flooding_factor: float | None = None,
    tray_type: str | None = None,
    contactor_internals_type: str | None = None,
    material_grade: str | None = None,
    max_operation_pressure: float | None = None,
) -> AbsorptionColumnResult:
    """Solve a tray absorber.

    Args:
        gas_components: the gas's substances, by name.
        gas_n: the gas feed's molar flow.
        gas_z: the gas feed's composition.
        gas_p: the gas feed's pressure.
        gas_t: the gas feed's temperature.
        solvent_components: the solvent's substances, in the gas's own order.
        solvent_n: the solvent feed's molar flow.
        solvent_z: the solvent feed's composition.
        solvent_p: the solvent feed's pressure.
        solvent_t: the solvent feed's temperature.
        number_of_stages: the trays, numbered 0 at the bottom where the gas enters.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at the bottom stage.
        tray_temperatures: one outlet-temperature pin per tray, `NaN` where a tray has none.
            **A pinned column stops after its first sweep**, because the gate is the
            tray-temperature change it made zero.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        murphree_efficiency: **not ported**; `SimpleTray.setMurphreeEfficiency` would close it.
        component_murphree_efficiency: **not ported**; the per-component override
            `applyMurphreeCorrection` applies.
        max_allowable_gas_load_factor: the Souders-Brown limit the K family is checked against.
            **It does not enter the solve.**
        solver_type: the base's strategy, by `process.distillation_column`'s own names.
        column_diameter: the column's internal diameter, which both capacity families divide the
            gas outlet's volumetric flow by. **The solve is indifferent to it.**
        max_allowable_fs_factor: the `Fs` limit, which this class defaults to `3.0` rather than
            the base's `2.5`.

    Returns:
        The tray profile, the treated gas, the loaded solvent and both capacity families.

    Raises:
        InvalidInputError: for a declared parameter whose arithmetic is not ported, or an inlet
            pair that does not carry the same substances in the same order.
        SolverNotConvergedError: when the solve misses its gate.

    See :func:`azoth.process.reference.distillation_column._states`.
    """
    # **The two efficiency fields, resolved and clamped.** `AbsorptionColumn` inherits the
    # base's two setters and adds `setComponentMurphreeEfficiency`'s map, so the correction reads
    # component, then per-tray, then column-wide; the per-tray-per-component map the class's
    # other overload writes has no palette spelling and is not carried.
    murphree = (
        None
        if murphree_efficiency is None and component_murphree_efficiency is None
        else _AbsorberMurphree(
            base=_Murphree(
                column_wide=_Murphree.clamp(
                    1.0 if murphree_efficiency is None else murphree_efficiency
                ),
                # **The base column's per-stage vector**, folded in where the class's own
                # `getMurphreeEfficiency(tray)` resolves it.
                per_stage=(
                    None
                    if tray_murphree_efficiency is None
                    else tuple(_Murphree.clamp(v) for v in tray_murphree_efficiency)
                ),
            ),
            per_component=(
                None
                if component_murphree_efficiency is None
                else tuple(_Murphree.clamp(v) for v in component_murphree_efficiency)
            ),
        ).checked(len(gas_components))
    )
    # **The same eight `process.distillation_column` refuses, from this id's own declaration**,
    # one literal arm each so `tools/check_unported.py` reads the same shape in both halves.
    if solver_type is not None and solver_type not in ("direct_substitution", "naphtali_sandholm"):
        if solver_type == "damped_substitution":
            raise _unported.refuse("solver_type=damped_substitution")
        if solver_type == "inside_out":
            raise _unported.refuse("solver_type=inside_out")
        if solver_type == "matrix_inside_out":
            raise _unported.refuse("solver_type=matrix_inside_out")
        if solver_type == "wegstein":
            raise _unported.refuse("solver_type=wegstein")
        if solver_type == "sum_rates":
            raise _unported.refuse("solver_type=sum_rates")
        if solver_type == "newton":
            raise _unported.refuse("solver_type=newton")
        if solver_type == "mesh_residual":
            raise _unported.refuse("solver_type=mesh_residual")
        if solver_type == "auto":
            raise _unported.refuse("solver_type=auto")
        raise InvalidInputError(
            "solver_type",
            f"`{solver_type}` is not one of `ColumnSolverFactory`'s ten strategies, which are "
            f"`direct_substitution`, `damped_substitution`, `inside_out`, `matrix_inside_out`, "
            f"`wegstein`, `sum_rates`, `newton`, `naphtali_sandholm`, `mesh_residual` and `auto`",
        )
    if list(gas_components) != list(solvent_components):
        raise InvalidInputError(
            "solvent_components",
            "the two inlets of one column carry the same substances in the same order: the "
            "class's own feed handling indexes both by the same component list, and a solvent "
            "that names different substances in a different order is a different fluid, not an "
            "inlet",
        )

    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    gas_pressure = input_to_si(spec, "gas_p", gas_p)
    gas_temperature = input_to_si(spec, "gas_t", gas_t)
    solvent_pressure = input_to_si(spec, "solvent_p", solvent_p)
    solvent_temperature = input_to_si(spec, "solvent_t", solvent_t)
    top = input_to_si(spec, "top_pressure", top_pressure)
    bottom = input_to_si(spec, "bottom_pressure", bottom_pressure)
    stages = int(_si(spec, "number_of_stages", number_of_stages))
    tolerance = _si(spec, "temperature_tolerance", temperature_tolerance)
    iterations_cap = int(_si(spec, "max_iterations", max_iterations))
    gas_flow = input_to_si(spec, "gas_n", gas_n)
    solvent_flow = input_to_si(spec, "solvent_n", solvent_n)

    # **Only the bounds this solve resolves.** `tray_efficiency` is the vessel's, not the
    # solve's, and it is applied where its value is known - so carrying it here would report a
    # bound as unchecked that the run does check.
    resolvable = {
        "number_of_stages": float(stages),
        "top_pressure": top,
        "bottom_pressure": bottom,
        "temperature_tolerance": tolerance,
        "gas_t": gas_temperature,
        "solvent_t": solvent_temperature,
    }
    apply_checks(
        [check for check in checks.on_input if check.quantity in resolvable],
        resolvable.get,
        warnings,
    )

    # **The gas is the column's own `feed` and the solvent its `top_feed`**, which is the
    # position `addGasInStream` and `addSolventInStream` give them: stage 0 and the top stage.
    #
    # **The draws are the base column's, and the two outer loops are its own too**: a specified
    # draw flow is a search over whole solves and a pumparound return is a fixed point, so this
    # wraps the same tears `process.distillation_column` does rather than re-deriving them.
    def solve_once(
        active_draws: Draws, active_returns: tuple[StreamRecord | None, ...] | None = None
    ) -> _States:
        return _states(
            list(gas_components),
            gas_flow,
            list(gas_z),
            gas_temperature,
            gas_pressure,
            stages,
            0,
            False,
            False,
            top,
            bottom,
            None,
            None,
            tolerance,
            iterations_cap,
            None,
            None,
            solver_type,
            top_feed=_feed(
                list(solvent_components),
                solvent_flow,
                list(solvent_z),
                solvent_temperature,
                solvent_pressure,
            ),
            tray_temperatures=None if tray_temperatures is None else tuple(tray_temperatures),
            reactive=reactive_section(reactive, reactive_start_tray, reactive_end_tray),
            draws=active_draws,
            absorber_murphree=murphree,
            pumparound_inlets=active_returns,
        )

    draws: Draws = (
        (
            None
            if gas_side_draw_fractions is None
            else tuple(float(v) for v in gas_side_draw_fractions),
            None
            if liquid_side_draw_fractions is None
            else tuple(float(v) for v in liquid_side_draw_fractions),
            None if pumparound_fractions is None else tuple(float(v) for v in pumparound_fractions),
        )
        if gas_side_draw_fractions is not None
        or liquid_side_draw_fractions is not None
        or pumparound_fractions is not None
        else None
    )
    # **An absorber has no ends**, so the two end flags are false and a draw on stage 0 or the
    # top stage is legitimate - which is what `_side_draw_flow_tear` is told here.
    flows = _side_draw_flow_specification(
        side_draw_flow_tray,
        side_draw_flow_phase,
        side_draw_flow_target
        if side_draw_flow_target is None
        else input_to_si(spec, "side_draw_flow_target", side_draw_flow_target),
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
    )
    returns = _pumparound_return_specification(
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        None
        if pumparound_temperature_drop is None
        else input_to_si(spec, "pumparound_temperature_drop", pumparound_temperature_drop),
    )
    if returns is not None:
        if flows is not None:
            raise InvalidInputError(
                "pumparound_return_tray",
                "a pumparound with a return is stated beside a side-draw flow specification: "
                "`solveWithColumnTearVariables` solves the two as coordinated tear variables, and "
                "this port carries the pumparound's own fixed point alone",
            )
        states = _pumparound_returns_tear(
            solve_once,
            list(gas_components),
            draws,
            returns,
            stages,
            1.0e-4 if pumparound_tolerance is None else float(pumparound_tolerance),
            12 if pumparound_max_iterations is None else int(pumparound_max_iterations),
        )
    elif flows is not None:
        states = _side_draw_flow_tear(
            solve_once, list(gas_components), draws, flows, stages, False, False
        )
    else:
        states = solve_once(draws)
    warnings.extend(states.warnings)

    # **The vessel around the absorber**, on the same trays. **The six geometry values are the
    # designer's own defaults and not declared inputs**: `AbsorptionColumn` carries no tray
    # spacing, weir height or hole diameter of its own, so this entry declares none either - the
    # capture's row is at exactly those defaults.
    resolved_efficiency = (
        DEFAULT_TRAY_EFFICIENCY if tray_efficiency is None else float(tray_efficiency)
    )
    geometry = MechanicalGeometry(
        tray_type=DEFAULT_TRAY_TYPE if tray_type is None else str(tray_type),
        contactor_internals_type=(
            DEFAULT_CONTACTOR_INTERNALS_TYPE
            if contactor_internals_type is None
            else str(contactor_internals_type)
        ),
        tray_efficiency=resolved_efficiency,
        max_flooding_factor=(
            DEFAULT_MAX_FLOODING_FACTOR
            if max_flooding_factor is None
            else float(max_flooding_factor)
        ),
        material_grade=DEFAULT_MATERIAL_GRADE if material_grade is None else str(material_grade),
        max_operation_pressure_bara=(
            DEFAULT_MAX_OPERATION_PRESSURE_BARA
            if max_operation_pressure is None
            else float(max_operation_pressure)
        ),
        tray_spacing_m=DEFAULT_TRAY_SPACING_M,
        weir_height_m=DEFAULT_WEIR_HEIGHT_M,
        hole_diameter_m=DEFAULT_HOLE_DIAMETER_MM / 1000.0,
        hole_area_fraction=DEFAULT_HOLE_AREA_FRACTION,
        downcommer_area_fraction=DEFAULT_DOWNCOMMER_AREA_FRACTION,
        column_diameter_override_m=UNSIZED_COLUMN_DIAMETER_M,
    )
    # The declared value the vessel cannot take, refused from the id's own `[[unported]]` row.
    if str(geometry.contactor_internals_type).lower() == "packed":
        raise _unported.refuse("contactor_internals_type=packed")
    # **The vessel sizing's own bound, applied where the value is known.** The solve resolves the
    # six *it* takes; this one is the record's.
    apply_checks(
        [check for check in checks_for(_spec()).on_input if check.quantity == "tray_efficiency"],
        {"tray_efficiency": resolved_efficiency}.get,
        warnings,
    )
    mechanical = mechanical_design(states, list(gas_components), geometry)

    return AbsorptionColumnResult(
        tray_temperature=tuple(from_si(value, "K") for value in states.tray_temperature),
        tray_pressure=tuple(from_si(value, "Pa") for value in states.tray_pressure),
        tray_gas_n=tuple(from_si(value, "mol/s") for value in states.tray_gas_n),
        tray_liquid_n=tuple(from_si(value, "mol/s") for value in states.tray_liquid_n),
        gas_out_n=from_si(states.distillate_n, "mol/s"),
        gas_out_z=states.distillate_z,
        gas_out_p=from_si(states.distillate_p, "Pa"),
        gas_out_t=from_si(states.distillate_t, "K"),
        gas_out_h=from_si(states.distillate_h, "J/mol"),
        liquid_out_n=from_si(states.bottoms_n, "mol/s"),
        liquid_out_z=states.bottoms_z,
        liquid_out_p=from_si(states.bottoms_p, "Pa"),
        liquid_out_t=from_si(states.bottoms_t, "K"),
        liquid_out_h=from_si(states.bottoms_h, "J/mol"),
        iterations=states.iterations,
        temperature_residual=from_si(states.temperature_residual, "K"),
        mass_residual=states.mass_residual,
        energy_residual=states.energy_residual,
        **capacity_fields(
            states,
            list(gas_components),
            column_diameter,
            max_allowable_fs_factor,
            max_allowable_gas_load_factor,
        ),
        vessel_diameter=from_si(mechanical.vessel_diameter_m, "m"),
        vessel_height=from_si(mechanical.vessel_height_m, "m"),
        vessel_wall_thickness=from_si(mechanical.vessel_wall_thickness_mm / 1000.0, "mm"),
        actual_trays=mechanical.actual_trays,
        flooding_factor=mechanical.flooding_factor,
        weir_loading=mechanical.weir_loading,
        tray_pressure_drop_mbar=mechanical.tray_pressure_drop_mbar,
        total_pressure_drop_bar=from_si(mechanical.total_pressure_drop_bar * 1.0e5, "bar"),
        reboiler_duty_kw=from_si(mechanical.reboiler_duty_kw * 1000.0, "kW"),
        condenser_duty_kw=from_si(mechanical.condenser_duty_kw * 1000.0, "kW"),
        material_grade=mechanical.material_grade,
        warnings=tuple(warnings),
    )


def capacity_fields(
    states: _States,
    components: list[str],
    column_diameter: Q | None,
    max_allowable_fs_factor: float | None,
    max_allowable_gas_load_factor: float | None,
) -> dict[str, Any]:
    """Both capacity families at the two products, as this class's eight getters answer them.

    **One resolution for the two ids**, because `StrippingColumn` extends this class and inherits
    every one of the eight - so a stripper's limits and an absorber's are the same computation
    under two names, and stating them once is what keeps them one.

    **The gas-load family's second density is the bottom tray's vapour and not the liquid
    product.** A NeqSim `getLiquidOutStream` carries the *tray's* two-phase system, so its phase 0
    is the vapour it leaves with; `process.absorption_column`'s Rust side says why that is a factor
    of `1.26` on the stripper rather than a formality.
    """
    gas = Stream.from_pt(
        components,
        list(states.distillate_z),
        states.distillate_n,
        from_si(states.distillate_p, "Pa"),
        from_si(states.distillate_t, "K"),
    )
    liquid_out_vapour = Stream.from_pt(
        components,
        list(states.tray_gas_z[0]),
        states.tray_gas_n[0],
        from_si(states.tray_pressure[0], "Pa"),
        from_si(states.tray_temperature[0], "K"),
    )
    diameter = (
        quantity(DEFAULT_INTERNAL_DIAMETER_M, "m") if column_diameter is None else column_diameter
    )
    fs = fs_limits(
        gas,
        diameter,
        DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER
        if max_allowable_fs_factor is None
        else float(max_allowable_fs_factor),
    )
    gas_load = gas_load_limits(
        gas,
        liquid_out_vapour,
        diameter,
        DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR
        if max_allowable_gas_load_factor is None
        else float(max_allowable_gas_load_factor),
    )
    return {
        "fs_factor": fs["fs_factor"],
        "fs_factor_utilization": fs["fs_factor_utilization"],
        "fs_factor_within_design_limit": fs["fs_factor_within_design_limit"],
        "minimum_diameter_for_fs_limit": from_si(fs["minimum_diameter_for_fs_limit"], "m"),
        "gas_load_factor": from_si(gas_load["gas_load_factor"], "m/s"),
        "gas_load_factor_utilization": gas_load["gas_load_factor_utilization"],
        "gas_load_factor_within_design_limit": gas_load["gas_load_factor_within_design_limit"],
        "minimum_diameter_for_gas_load_limit": from_si(
            gas_load["minimum_diameter_for_gas_load_limit"], "m"
        ),
    }


def _spec() -> dict[str, object]:
    """This model's generated spec."""
    import azoth._models_gen as models

    return models.model("process.absorption_column")


__all__ = ["DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR", "UNPORTED_SOLVERS", "absorption_column"]
