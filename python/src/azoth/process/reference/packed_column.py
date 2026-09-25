"""``process.packed_column`` - the packed column, which is the tray column at a derived height.

Spec: ``specs/models/process/packed_column.toml``.

The Python twin of ``crates/azoth-process/src/models/packed_column.rs``. **The arithmetic is
``process.distillation_column``'s**, because the class's is: ``PackedColumn extends
DistillationColumn`` in 463 lines and its ``run`` is ``super.run(id)`` followed by
``calcPackingHydraulics()``. What this module adds is the stage count a constructor derives from
a packed height, and the refusal of the one packing parameter the class itself refuses.

**The packing group is the excluded report's**, so four of its five parameters are carried
inertly: ``packing_type``, ``structured_packing``, ``design_flood_fraction`` and
``column_diameter`` are read by ``ColumnInternalsDesigner`` after the column has converged, and
no part of the solve reads them.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.result import PackedColumnResult
from azoth.core.units import Q
from azoth.process.reference.distillation_column import (
    distillation_column as _distillation_column,
)

#: The HETP guess ``PackedColumn``'s constructors divide a packed height by.
HETP_GUESS_M = 0.5


def stage_count(packed_height_m: float) -> int:
    """The stage count ``PackedColumn.estimateStages`` derives from a packed height.

    ``ceil(packed_height / 0.5)`` floored at two. **The floor is reached rather than refused,
    and that is measured**: the capture's stage-count rows report two middle trays for ``0.2`` m,
    for ``0.0`` m and for ``-1.0`` m alike, because neither ``estimateStages`` nor
    ``setPackedHeight`` has a domain check.

    Args:
        packed_height_m: the packed bed height, in metres.

    Returns:
        The middle trays the constructor would make.
    """
    if not math.isfinite(packed_height_m):
        # Java's `(int)` cast maps a `NaN` to zero and saturates an infinity at the type's
        # maximum; the floor decides the first, and the second is an allocation either way.
        return 2
    return max(math.ceil(packed_height_m / HETP_GUESS_M), 2)


def packed_column(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    packed_height: Q,
    feed_stage: int,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float = 1.0e-6,
    max_iterations: int = 200,
    reboiler_temperature: Q | None = None,
    condenser_temperature: Q | None = None,
    packing_type: str | None = None,
    structured_packing: bool | None = None,
    design_flood_fraction: float | None = None,
    packing_hydraulic_capacity_factor: float | None = None,
    column_diameter: Q | None = None,
    murphree_efficiency: float | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
) -> PackedColumnResult:
    """Solve a packed column.

    Args:
        components: the feed's substances, by name.
        feed_n: the feed's molar flow.
        feed_z: the feed composition.
        feed_p: the feed's pressure, which is its own.
        feed_t: the feed's temperature.
        packed_height: the packed bed height, which the constructor divides by a ``0.5`` m HETP
            guess to fix the stage count: ``ceil(packed_height / 0.5)`` floored at two.
        feed_stage: the stage the feed enters, 0-based over the trays including the ends.
        has_reboiler: whether stage 0 is a reboiler.
        has_condenser: whether the top stage is a condenser.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at stage 0.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        reboiler_temperature: the reboiler's temperature, which pins the bottom tray.
        condenser_temperature: the condenser's temperature, which pins the top tray.
        packing_type: the packing's name. **It does not enter the solve.**
        structured_packing: whether the packing is structured. **It does not enter the solve.**
        design_flood_fraction: the fraction of flood the packing is sized to. **It does not enter
            the solve.**
        packing_hydraulic_capacity_factor: the packing's relative hydraulic capacity factor.
            **Refused where it is not positive and finite**, which is
            ``setPackingHydraulicCapacityFactor``'s own check.
        column_diameter: the column's internal diameter. **It does not enter the solve.**
        murphree_efficiency: **not ported**; omitted is the ideal stage.
        solver_type: **only ``direct_substitution`` and ``naphtali_sandholm`` are ported**.
        top_specification_type: the top product's degree of freedom.
        top_specification_target: its target value.
        top_specification_component: the component a purity or a recovery constrains.
        bottom_specification_type: the bottom product's degree of freedom.
        bottom_specification_target: its target value.
        bottom_specification_component: the component a purity or a recovery constrains.

    Returns:
        The base column's own record, which is the id's claim.

    Raises:
        InvalidInputError: for a ``packing_hydraulic_capacity_factor`` that is not positive and
            finite, and for every input ``process.distillation_column`` refuses.

    See :func:`azoth.process.reference.distillation_column.distillation_column`.
    """
    if packing_hydraulic_capacity_factor is not None and (
        not math.isfinite(packing_hydraulic_capacity_factor)
        or packing_hydraulic_capacity_factor <= 0.0
    ):
        raise InvalidInputError(
            "packing_hydraulic_capacity_factor",
            f"a packing hydraulic capacity factor of {packing_hydraulic_capacity_factor} is not "
            "positive and finite, which `PackedColumn.setPackingHydraulicCapacityFactor` refuses "
            "with an `IllegalArgumentException`",
        )
    # The other four are declarations the solve is indifferent to: the class accepts them and
    # `ColumnInternalsDesigner` reads each one after the column has converged.
    del packing_type, structured_packing, design_flood_fraction, column_diameter

    out = _distillation_column(
        components,
        feed_n,
        feed_z,
        feed_p,
        feed_t,
        stage_count(_metres(packed_height)),
        feed_stage,
        has_reboiler,
        has_condenser,
        top_pressure,
        bottom_pressure,
        temperature_tolerance,
        max_iterations,
        reboiler_temperature,
        condenser_temperature,
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
    )
    return PackedColumnResult(
        tray_temperature=out.tray_temperature,
        tray_pressure=out.tray_pressure,
        tray_gas_n=out.tray_gas_n,
        tray_liquid_n=out.tray_liquid_n,
        distillate_n=out.distillate_n,
        distillate_z=out.distillate_z,
        distillate_p=out.distillate_p,
        distillate_t=out.distillate_t,
        distillate_h=out.distillate_h,
        bottoms_n=out.bottoms_n,
        bottoms_z=out.bottoms_z,
        bottoms_p=out.bottoms_p,
        bottoms_t=out.bottoms_t,
        bottoms_h=out.bottoms_h,
        gas_side_draw_n=out.gas_side_draw_n,
        liquid_side_draw_n=out.liquid_side_draw_n,
        pumparound_n=out.pumparound_n,
        condenser_duty=out.condenser_duty,
        reboiler_duty=out.reboiler_duty,
        iterations=out.iterations,
        temperature_residual=out.temperature_residual,
        mass_residual=out.mass_residual,
        energy_residual=out.energy_residual,
        warnings=out.warnings,
    )


def _metres(packed_height: Q) -> float:
    """The packed height in metres, which is the unit the constructor's rule is stated in."""
    return float(packed_height.to("m").magnitude)
