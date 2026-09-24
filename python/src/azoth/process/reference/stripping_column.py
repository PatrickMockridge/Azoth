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
    solver_type: str | None = None,
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
        max_allowable_gas_load_factor: a design limit the solve does not read.
        solver_type: the base's strategy, by ``process.distillation_column``'s own names.

    Returns:
        The tray profile, the stripped gas and the lean liquid.

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
        solver_type,
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
        warnings=out.warnings,
    )
