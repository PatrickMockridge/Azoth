"""``process.packed_column`` - the packed column, which is the tray column at a derived height.

Spec: ``specs/models/process/packed_column.toml``.

The Python twin of ``crates/azoth-process/src/models/packed_column.rs``. **The arithmetic is
``process.distillation_column``'s**, because the class's is: ``PackedColumn extends
DistillationColumn`` in 463 lines and its ``run`` is ``super.run(id)`` followed by
``calcPackingHydraulics()``. What this module adds is the stage count a constructor derives from
a packed height, and the refusal of the one packing parameter the class itself refuses.

**The packing group reaches the report and not the solve.** ``packing_type``,
``structured_packing``, ``design_flood_fraction`` and ``column_diameter`` are read by
``ColumnInternalsDesigner`` after the column has converged, so ``PackedColumn.run`` is
``super.run(id)`` then a report - and the last seven fields of the record are that report.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.result import PackedColumnResult
from azoth.core.units import Q, from_si, quantity
from azoth.hydraulics.reference.packing_hydraulics import _packing_hydraulics
from azoth.hydraulics.reference.packing_sizing import packing_sizing as _packing_sizing
from azoth.process.kernels import Stream
from azoth.process.reference.capacity import DEFAULT_MAX_ALLOWABLE_FS_FACTOR, fs_limits
from azoth.process.reference.distillation_column import _distillation_column_states
from azoth.process.reference.rate_based_packed_column import _phase_view

#: The HETP guess ``PackedColumn``'s constructors divide a packed height by.
HETP_GUESS_M = 0.5

#: The packing name ``PackedColumn`` carries when none is stated.
DEFAULT_PACKING_NAME = "Pall-Ring-50"
#: The fraction of flood the class sizes a bed to when none is stated.
DEFAULT_DESIGN_FLOOD_FRACTION = 0.70
#: The relative hydraulic capacity factor the class carries when none is stated.
DEFAULT_HYDRAULIC_CAPACITY_FACTOR = 1.0
#: The internal diameter ``PackedColumn`` carries when none is stated.
#:
#: **At or below zero rather than absent**, which is what makes it the sizing branch.
UNSIZED_COLUMN_DIAMETER_M = -1.0

#: The surface tension ``ColumnInternalsDesigner.getTrayProperties`` falls back to, in N/m.
#:
#: **A third constant, and it is not the rate-based path's.** The designer asks
#: ``fluid.getInterphaseProperties().getSurfaceTension(0, 1)`` and answers ``0.02`` when that is
#: not finite and positive, while ``RateBasedPackedColumn.estimateSurfaceTension`` answers
#: ``0.025`` on the same class of pair. The capture prints both interface readings as ``0.0`` on
#: every tray, so this constant is what answers.
DESIGNER_SURFACE_TENSION_N_PER_M = 0.02

#: The vapour viscosity ``PackingHydraulicsCalculator`` holds when nobody sets one, in Pa*s.
#:
#: **The designer does not set it, and that is read from the jar rather than guessed.**
#: Disassembling the pinned jar shows ``ColumnInternalsDesigner`` calling exactly
#: ``setPackingPreset``, ``setDesignFloodFraction``, ``setVaporDensity``, ``setLiquidDensity``,
#: ``setLiquidViscosity``, ``setSurfaceTension`` and ``setColumnDiameter`` on the calculator, so
#: this field keeps its constructor initialiser. No published quantity reads it.
CALCULATOR_VAPOR_VISCOSITY_PA_S = 1.0e-5
#: The vapour diffusivity the calculator holds when nobody sets one, in m**2/s.
CALCULATOR_VAPOR_DIFFUSIVITY_M2_S = 1.0e-5
#: The liquid diffusivity the calculator holds when nobody sets one, in m**2/s.
CALCULATOR_LIQUID_DIFFUSIVITY_M2_S = 1.0e-9


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
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
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
    max_allowable_fs_factor: float | None = None,
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
        murphree_efficiency: the column-wide Murphree tray efficiency; omitted is the ideal stage.
        tray_murphree_efficiency: one override per *theoretical* stage, ``NaN`` where a stage
            falls through to the column-wide value.
        gas_side_draw_fractions: the vapour each stage withdraws, one entry per stage.
        liquid_side_draw_fractions: the liquid each stage withdraws as a side draw.
        pumparound_fractions: the liquid each stage withdraws as a pumparound.
        side_draw_flow_tray: the stage whose draw flow is specified.
        side_draw_flow_phase: which of that stage's phases the specification controls.
        side_draw_flow_target: the mass flow the draw must deliver.
        side_draw_flow_tolerance: the relative residual the specification's search stops at.
        side_draw_flow_max_iterations: the candidate cap for that search.
        pumparound_return_tray: the stage a pumparound's liquid comes back to.
        pumparound_draw_tray: the stage it leaves.
        pumparound_draw_fraction: the fraction of that stage's liquid withdrawn.
        pumparound_temperature_drop: the cooler's drop on the way back.
        pumparound_tolerance: the relative residual the return's loop stops at.
        pumparound_max_iterations: the cap on that loop.
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
    # **The stages, not the record alone**, because the report is read at a tray's own fluid and
    # the record's profile carries no compositions. `_record` is the same construction
    # `distillation_column` uses, so this is one solve and one mapping.
    states, warnings = _distillation_column_states(
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
        # **A notional stage is still a stage.** The height derives the tray count, so these are
        # one entry per theoretical stage with the same `NaN` fall-through the column's own
        # vector uses.
        tray_murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        # **The section is the base's own**, whose signature carries it after the two
        # specifications - which is the order the spec declares its inputs in, and the order a
        # positional caller has to follow.
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        # **The draws are the base column's and `PackedColumn` inherits every one of them**, so a
        # draw on a notional stage is the split the column makes. The two unit-carrying targets
        # cross as quantities, which is what the base reference's own signature takes.
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
    )
    # ---- The report, at the middle tray: `trays.size() / 2`, the class's own index.
    middle = len(states.tray_temperature) // 2
    vapor = _phase_view(
        list(components),
        "gas",
        states.tray_temperature[middle],
        states.tray_pressure[middle],
        list(states.tray_gas_z[middle]),
        states.tray_gas_n[middle],
    )
    liquid = _phase_view(
        list(components),
        "liquid",
        states.tray_temperature[middle],
        states.tray_pressure[middle],
        list(states.tray_liquid_z[middle]),
        states.tray_liquid_n[middle],
    )
    vapor_mass_flow = states.tray_gas_n[middle] * vapor["molar_mass"]
    liquid_mass_flow = states.tray_liquid_n[middle] * liquid["molar_mass"]
    vapor_density = quantity(vapor["density"], "kg/m**3")
    liquid_density = quantity(liquid["density"], "kg/m**3")

    flood_fraction = (
        DEFAULT_DESIGN_FLOOD_FRACTION
        if design_flood_fraction is None
        else float(design_flood_fraction)
    )
    capacity_factor = (
        DEFAULT_HYDRAULIC_CAPACITY_FACTOR
        if packing_hydraulic_capacity_factor is None
        else float(packing_hydraulic_capacity_factor)
    )
    stated_diameter = (
        UNSIZED_COLUMN_DIAMETER_M
        if column_diameter is None
        else float(column_diameter.to("m").magnitude)
    )
    name = DEFAULT_PACKING_NAME if packing_type is None else str(packing_type)

    # **The sizing branch first, then the hydraulics at the diameter it resolved**, which is the
    # composition `calcPackingHydraulics` performs.
    internal_diameter = (
        quantity(stated_diameter, "m")
        if stated_diameter > 0.0
        else _packing_sizing(
            name,
            flood_fraction,
            quantity(vapor_mass_flow, "kg/s"),
            quantity(liquid_mass_flow, "kg/s"),
            vapor_density,
            liquid_density,
            liquid["transport"].mu,
            capacity_factor,
        ).column_diameter
    )

    # **The Fs family is read at the diameter the sizing just resolved**, which is the class's
    # own order: `calcPackingHydraulics` ends by writing the inherited `internalDiameter` from
    # the stated `columnDiameter` where it is positive and from the designer's sizing otherwise.
    #
    # **And the base's record is not built here.** `PackedColumn.calcPackingHydraulics` builds its
    # designer with `internalsType = "packed"`, so `calculateTrayed` never runs for this class:
    # there is no trayed internals report, and this result carries the base's own quantities
    # without one - the same split the Rust half has.
    limits = fs_limits(
        Stream.from_pt(
            list(components),
            list(states.distillate_z),
            states.distillate_n,
            from_si(states.distillate_p, "Pa"),
            from_si(states.distillate_t, "K"),
        ),
        internal_diameter,
        DEFAULT_MAX_ALLOWABLE_FS_FACTOR
        if max_allowable_fs_factor is None
        else float(max_allowable_fs_factor),
    )

    hydraulics = _packing_hydraulics(
        name,
        internal_diameter,
        packed_height,
        quantity(vapor_mass_flow, "kg/s"),
        quantity(liquid_mass_flow, "kg/s"),
        vapor_density,
        liquid_density,
        # **The tray's liquid viscosity and the calculator's own vapour viscosity**, which is the
        # asymmetry the jar shows: the designer sets the first and not the second.
        quantity(CALCULATOR_VAPOR_VISCOSITY_PA_S, "Pa*s"),
        liquid["transport"].mu,
        quantity(DESIGNER_SURFACE_TENSION_N_PER_M, "N/m"),
        quantity(CALCULATOR_VAPOR_DIFFUSIVITY_M2_S, "m**2/s"),
        quantity(CALCULATOR_LIQUID_DIFFUSIVITY_M2_S, "m**2/s"),
        capacity_factor,
        structured=structured_packing,
    )

    return PackedColumnResult(
        tray_temperature=tuple(from_si(value, "K") for value in states.tray_temperature),
        tray_pressure=tuple(from_si(value, "Pa") for value in states.tray_pressure),
        tray_gas_n=tuple(from_si(value, "mol/s") for value in states.tray_gas_n),
        tray_liquid_n=tuple(from_si(value, "mol/s") for value in states.tray_liquid_n),
        distillate_n=from_si(states.distillate_n, "mol/s"),
        distillate_z=states.distillate_z,
        distillate_p=from_si(states.distillate_p, "Pa"),
        distillate_t=from_si(states.distillate_t, "K"),
        distillate_h=from_si(states.distillate_h, "J/mol"),
        bottoms_n=from_si(states.bottoms_n, "mol/s"),
        bottoms_z=states.bottoms_z,
        bottoms_p=from_si(states.bottoms_p, "Pa"),
        bottoms_t=from_si(states.bottoms_t, "K"),
        bottoms_h=from_si(states.bottoms_h, "J/mol"),
        gas_side_draw_n=tuple(from_si(value, "mol/s") for value in states.gas_side_draw_n),
        liquid_side_draw_n=tuple(from_si(value, "mol/s") for value in states.liquid_side_draw_n),
        pumparound_n=tuple(from_si(value, "mol/s") for value in states.pumparound_n),
        condenser_duty=from_si(states.condenser_duty, "W"),
        reboiler_duty=from_si(states.reboiler_duty, "W"),
        iterations=states.iterations,
        temperature_residual=states.temperature_residual,
        mass_residual=states.mass_residual,
        energy_residual=states.energy_residual,
        fs_factor=limits["fs_factor"],
        fs_factor_utilization=limits["fs_factor_utilization"],
        fs_factor_within_design_limit=limits["fs_factor_within_design_limit"],
        minimum_diameter_for_fs_limit=from_si(limits["minimum_diameter_for_fs_limit"], "m"),
        hetp=hydraulics.hetp,
        theoretical_stages=hydraulics.theoretical_stages,
        percent_flood=hydraulics.percent_flood,
        flooding_velocity=hydraulics.flooding_velocity,
        packing_pressure_drop=hydraulics.total_pressure_drop,
        hydraulics_ok=hydraulics.design_ok,
        internal_diameter=internal_diameter,
        warnings=tuple(warnings),
    )


def _metres(packed_height: Q) -> float:
    """The packed height in metres, which is the unit the constructor's rule is stated in."""
    return float(packed_height.to("m").magnitude)
