"""``process.stripping_column`` - the tray stripper, which is the absorber renamed.

Spec: ``specs/models/process/stripping_column.toml``.

The Python twin of ``crates/azoth-process/src/models/stripping_column.rs``. **The arithmetic is
``process.absorption_column``'s**, because the class is: `StrippingColumn extends AbsorptionColumn`
and its ninety lines rename the two inlets and the two products. The class's own javadoc gives
the reason - absorption and stripping are one set of counter-current equilibrium-stage equations,
and the thermodynamic driving force sets the direction of transfer - so the stripping gas enters
stage 0 and the rich liquid the top stage, exactly where an absorber's gas and solvent enter.
"""

from __future__ import annotations

from azoth.core.result import StrippingColumnResult
from azoth.core.units import Q
from azoth.process.reference.absorption_column import absorption_column as absorber


def stripping_column(
    stripping_gas_components: list[str],
    stripping_gas_n: Q,
    stripping_gas_z: list[float],
    stripping_gas_p: Q,
    stripping_gas_t: Q,
    rich_liquid_components: list[str],
    rich_liquid_n: Q,
    rich_liquid_z: list[float],
    rich_liquid_p: Q,
    rich_liquid_t: Q,
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
) -> StrippingColumnResult:
    """Solve a tray stripper.

    Args:
        stripping_gas_components: the stripping gas's substances, by name.
        stripping_gas_n: its molar flow, which enters stage 0.
        stripping_gas_z: its composition.
        stripping_gas_p: its pressure.
        stripping_gas_t: its temperature.
        rich_liquid_components: the rich liquid's substances, in the gas's own order.
        rich_liquid_n: its molar flow, which enters the top stage.
        rich_liquid_z: its composition.
        rich_liquid_p: its pressure.
        rich_liquid_t: its temperature.
        number_of_stages: the trays, numbered 0 at the bottom.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at the bottom stage.
        tray_temperatures: one outlet-temperature pin per tray.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        murphree_efficiency: **not ported**.
        component_murphree_efficiency: **not ported**.
        max_allowable_gas_load_factor: the Souders-Brown limit the K family is checked against.
            **It does not enter the solve.**
        solver_type: the base's strategy, by ``process.distillation_column``'s own names.
        column_diameter: the column's internal diameter, which both capacity families divide the
            gas outlet's volumetric flow by.
        max_allowable_fs_factor: the ``Fs`` limit, which the absorber this class extends
            defaults to ``3.0``.

    Returns:
        The tray profile, the stripped gas, the lean liquid and both capacity families.

    See :func:`azoth.process.reference.absorption_column.absorption_column`.
    """
    out = absorber(
        stripping_gas_components,
        stripping_gas_n,
        stripping_gas_z,
        stripping_gas_p,
        stripping_gas_t,
        rich_liquid_components,
        rich_liquid_n,
        rich_liquid_z,
        rich_liquid_p,
        rich_liquid_t,
        number_of_stages,
        top_pressure,
        bottom_pressure,
        tray_temperatures,
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        component_murphree_efficiency,
        max_allowable_gas_load_factor,
        # **The three follow `absorption_column`'s own signature and not this function's**,
        # because the delegation is positional: the spec declares them before `solver_type`
        # and the absorber's signature carries them there.
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        solver_type,
        tray_murphree_efficiency,
        gas_side_draw_fractions,
        liquid_side_draw_fractions,
        pumparound_fractions,
        side_draw_flow_tray,
        side_draw_flow_phase,
        side_draw_flow_target,
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        pumparound_temperature_drop,
        pumparound_tolerance,
        pumparound_max_iterations,
        column_diameter,
        max_allowable_fs_factor,
    )
    return StrippingColumnResult(
        tray_temperature=out.tray_temperature,
        tray_pressure=out.tray_pressure,
        tray_gas_n=out.tray_gas_n,
        tray_liquid_n=out.tray_liquid_n,
        overhead_gas_n=out.gas_out_n,
        overhead_gas_z=out.gas_out_z,
        overhead_gas_p=out.gas_out_p,
        overhead_gas_t=out.gas_out_t,
        overhead_gas_h=out.gas_out_h,
        lean_liquid_n=out.liquid_out_n,
        lean_liquid_z=out.liquid_out_z,
        lean_liquid_p=out.liquid_out_p,
        lean_liquid_t=out.liquid_out_t,
        lean_liquid_h=out.liquid_out_h,
        iterations=out.iterations,
        temperature_residual=out.temperature_residual,
        mass_residual=out.mass_residual,
        energy_residual=out.energy_residual,
        # **The base's own capacity answers, inherited and overridden in no part**: the
        # absorber's delegation already computed them at this id's two inputs.
        fs_factor=out.fs_factor,
        fs_factor_utilization=out.fs_factor_utilization,
        fs_factor_within_design_limit=out.fs_factor_within_design_limit,
        minimum_diameter_for_fs_limit=out.minimum_diameter_for_fs_limit,
        gas_load_factor=out.gas_load_factor,
        gas_load_factor_utilization=out.gas_load_factor_utilization,
        gas_load_factor_within_design_limit=out.gas_load_factor_within_design_limit,
        minimum_diameter_for_gas_load_limit=out.minimum_diameter_for_gas_load_limit,
        warnings=out.warnings,
    )
