"""The process layer: unit operations and flowsheets.

The ``process.*`` ids, one per palette entry under ``specs/unit_ops/``, and the check
that holds a flowsheet to the calculus's rules. Each id is a registered model: a spec
under ``specs/models/process/``, a kernel in Rust and a reference here, compared by
``test_cross_impl`` against a committed NeqSim capture.

The Stream-level arithmetic those models wrap is :mod:`azoth.process.kernels`, which is what
a flowsheet's executor calls and what ``azoth check`` validates the wiring of. To *run* one,
:func:`run_flowsheet` in :mod:`azoth.process.flowsheet`; to hold one and edit it,
:class:`Session` in the same module — the middleware's live document, which is also what
:func:`forms` and :func:`tools` describe.
"""

from __future__ import annotations

import pathlib

from azoth import _core
from azoth._dispatch import resolve
from azoth.core.result import (
    AbsorptionColumnResult,
    ComponentSplitterResult,
    CompressorResult,
    CoolerResult,
    DistillationColumnResult,
    EjectorResult,
    ExpanderResult,
    FilterResult,
    FlareResult,
    GasScrubberResult,
    GibbsReactorResult,
    HeaterResult,
    HeatExchangerResult,
    HemReleaseResult,
    ManifoldResult,
    MapTurboExpanderResult,
    MixerResult,
    PackedColumnResult,
    PipeResult,
    PlugFlowReactorResult,
    PumpResult,
    RateBasedPackedColumnResult,
    SeparatorResult,
    ShortcutDistillationColumnResult,
    SplitterResult,
    StirredTankReactorResult,
    StrippingColumnResult,
    TankResult,
    ThreePhaseSeparatorResult,
    ThrottlingValveResult,
    TurboExpanderResult,
)
from azoth.core.units import Q
from azoth.process.flowsheet import (
    FlowsheetResult,
    FlowsheetStream,
    Residuals,
    Session,
    TearResult,
    catalogue,
    forms,
    run_flowsheet,
    tools,
)
from azoth.process.kernels import Stream

__all__ = [
    "FlowsheetResult",
    "FlowsheetStream",
    "Residuals",
    "Session",
    "Stream",
    "TearResult",
    "absorption_column",
    "catalogue",
    "component_splitter",
    "compressor",
    "cooler",
    "distillation_column",
    "ejector",
    "expander",
    "filter",
    "flare",
    "forms",
    "gas_scrubber",
    "gibbs_reactor",
    "heat_exchanger",
    "heater",
    "hem_release",
    "load_flowsheet",
    "manifold",
    "map_turbo_expander_compressor",
    "mixer",
    "packed_column",
    "pipe",
    "plug_flow_reactor",
    "pump",
    "rate_based_packed_column",
    "run_flowsheet",
    "separator",
    "shortcut_distillation_column",
    "splitter",
    "stirred_tank_reactor",
    "stripping_column",
    "tank",
    "three_phase_separator",
    "throttling_valve",
    "tools",
    "validate",
]

_MANIFOLD = "process.manifold"
_MIXER = "process.mixer"
_PACKED_COLUMN = "process.packed_column"
_RATE_BASED_PACKED_COLUMN = "process.rate_based_packed_column"
_PIPE = "process.pipe"
_GAS_SCRUBBER = "process.gas_scrubber"
_HEAT_EXCHANGER = "process.heat_exchanger"
_COMPONENT_SPLITTER = "process.component_splitter"
_COMPRESSOR = "process.compressor"
_ABSORPTION_COLUMN = "process.absorption_column"
_STRIPPING_COLUMN = "process.stripping_column"
_DISTILLATION_COLUMN = "process.distillation_column"
_COOLER = "process.cooler"
_EJECTOR = "process.ejector"
_EXPANDER = "process.expander"
_TURBO_EXPANDER_COMPRESSOR = "process.turbo_expander_compressor"
_MAP_TURBO_EXPANDER_COMPRESSOR = "process.map_turbo_expander_compressor"
_FILTER = "process.filter"
_FLARE = "process.flare"
_HEATER = "process.heater"
_HEM_RELEASE = "process.hem_release"
_SEPARATOR = "process.separator"
_SHORTCUT_DISTILLATION_COLUMN = "process.shortcut_distillation_column"
_THROTTLING_VALVE = "process.throttling_valve"
_GIBBS_REACTOR = "process.gibbs_reactor"
_PLUG_FLOW_REACTOR = "process.plug_flow_reactor"
_PUMP = "process.pump"
_SPLITTER = "process.splitter"
_STIRRED_TANK_REACTOR = "process.stirred_tank_reactor"
_TANK = "process.tank"
_THREE_PHASE_SEPARATOR = "process.three_phase_separator"


def validate(flowsheet: str, palette_dir: str = "specs/unit_ops") -> list[str]:
    """Validate a flowsheet's TOML text against a palette directory.

    Returns the diagnostic lines; empty when the flowsheet is clean.
    """
    return _core.validate_flowsheet(flowsheet, palette_dir)


def load_flowsheet(path: str, palette_dir: str = "specs/unit_ops") -> list[str]:
    """Read a flowsheet file and validate it; returns the diagnostic lines."""
    return validate(pathlib.Path(path).read_text(), palette_dir)


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
) -> PumpResult:
    """Raise a stream's pressure, adding the pump's work as enthalpy.

    ``components``, ``inlet_n``, ``inlet_z``, ``inlet_p`` and ``inlet_t`` are the inlet's
    record - the fluid, its flow, composition, pressure and temperature. The inlet's molar
    enthalpy is *not* an input: it is a state function of ``(T, P, z)``, so accepting one
    would let a caller hand over a state that does not exist.

    Without a curve the work is the isentropic head - ``(H(P_out, s_in) - H(P_in)) / eta``,
    from a flash - and ``outlet_pressure`` is the pressure it reaches. A non-empty ``curve``
    selects ``PumpChart``'s reduced quadratic instead: the head sets the pressure rise, the
    chart's efficiency sets the shaft power, and ``outlet_pressure`` is ignored. ``speed``,
    ``reference_density``, ``viscosity``, ``use_viscosity_correction``, ``check_npsh`` and
    ``npsh_margin`` are the chart's own surface.

    Raises:
        InvalidInputError: where the shapes disagree, the efficiency is outside ``(0, 1]``, a
            pump with no curve has no outlet pressure, or the curve is not a curve.

    See :func:`azoth.process.reference.pump`.
    """
    return resolve(_PUMP)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_pressure=outlet_pressure,
        isentropic_efficiency=isentropic_efficiency,
        curve=curve,
        speed=speed,
        reference_density=reference_density,
        viscosity=viscosity,
        use_viscosity_correction=use_viscosity_correction,
        check_npsh=check_npsh,
        npsh_margin=npsh_margin,
    )


def distillation_column(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    number_of_stages: int,
    feed_stage: int,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float = 1.0e-6,
    max_iterations: int = 200,
    reboiler_temperature: Q | None = None,
    condenser_temperature: Q | None = None,
    murphree_efficiency: float | None = None,
    tray_murphree_efficiency: list[float] | None = None,
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
    internals_type: str | None = None,
    tray_spacing: Q | None = None,
    weir_height: Q | None = None,
    hole_diameter: Q | None = None,
    hole_area_fraction: float | None = None,
    downcommer_area_fraction: float | None = None,
    design_flood_fraction: float | None = None,
    column_diameter_override: Q | None = None,
    hydraulic_pressure_drop_coupling: bool | None = None,
    hydraulic_pressure_drop_internals_type: str | None = None,
    tray_efficiency: float | None = None,
    max_flooding_factor: float | None = None,
    tray_type: str | None = None,
    contactor_internals_type: str | None = None,
    material_grade: str | None = None,
    max_operation_pressure: float | None = None,
) -> DistillationColumnResult:
    """Solve a distillation column by sequential substitution.

    ``number_of_stages`` counts the trays **between the ends**, so a reboiler or a condenser
    adds one stage each; ``feed_stage`` is 0-based over the trays *including* them, which is
    why stage 0 is the reboiler when there is one. The ends are pinned by temperature -
    ``condenser_temperature`` and ``reboiler_temperature`` reach the trays' own outlet
    specifications, exactly as ``setCondenserTemperature`` does - and a middle tray flashes at
    its own enthalpy.

    **The answer is a profile**, not a scalar: the tray temperatures, pressures and both
    traffic rates, then the two products and the two duties. **Four of the twenty-two declared
    parameters are refused by name** - the Murphree efficiency and the two product
    specifications, whose arithmetic is not ported - as are the nine other solving strategies
    the class carries. A refusal names the class that would close it, so what is owed is
    readable rather than absent.

    Raises:
        InvalidInputError: for a stage or a feed stage outside the column, and for any
            declared parameter whose arithmetic is not ported.
        SolverNotConvergedError: when the solve misses its gate; the message carries the
            residuals.

    See :func:`azoth.process.reference.distillation_column`.
    """
    return resolve(_DISTILLATION_COLUMN)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        number_of_stages=number_of_stages,
        feed_stage=feed_stage,
        has_reboiler=has_reboiler,
        has_condenser=has_condenser,
        top_pressure=top_pressure,
        bottom_pressure=bottom_pressure,
        reboiler_temperature=reboiler_temperature,
        condenser_temperature=condenser_temperature,
        temperature_tolerance=temperature_tolerance,
        max_iterations=max_iterations,
        murphree_efficiency=murphree_efficiency,
        tray_murphree_efficiency=tray_murphree_efficiency,
        reactive=reactive,
        reactive_start_tray=reactive_start_tray,
        reactive_end_tray=reactive_end_tray,
        solver_type=solver_type,
        top_specification_type=top_specification_type,
        top_specification_target=top_specification_target,
        top_specification_component=top_specification_component,
        bottom_specification_type=bottom_specification_type,
        bottom_specification_target=bottom_specification_target,
        bottom_specification_component=bottom_specification_component,
        gas_side_draw_fractions=gas_side_draw_fractions,
        liquid_side_draw_fractions=liquid_side_draw_fractions,
        pumparound_fractions=pumparound_fractions,
        side_draw_flow_tray=side_draw_flow_tray,
        side_draw_flow_phase=side_draw_flow_phase,
        side_draw_flow_target=side_draw_flow_target,
        side_draw_flow_tolerance=side_draw_flow_tolerance,
        side_draw_flow_max_iterations=side_draw_flow_max_iterations,
        pumparound_return_tray=pumparound_return_tray,
        pumparound_draw_tray=pumparound_draw_tray,
        pumparound_draw_fraction=pumparound_draw_fraction,
        pumparound_temperature_drop=pumparound_temperature_drop,
        pumparound_tolerance=pumparound_tolerance,
        pumparound_max_iterations=pumparound_max_iterations,
        column_diameter=column_diameter,
        max_allowable_fs_factor=max_allowable_fs_factor,
        internals_type=internals_type,
        tray_spacing=tray_spacing,
        weir_height=weir_height,
        hole_diameter=hole_diameter,
        hole_area_fraction=hole_area_fraction,
        downcommer_area_fraction=downcommer_area_fraction,
        design_flood_fraction=design_flood_fraction,
        column_diameter_override=column_diameter_override,
        hydraulic_pressure_drop_coupling=hydraulic_pressure_drop_coupling,
        hydraulic_pressure_drop_internals_type=hydraulic_pressure_drop_internals_type,
        tray_efficiency=tray_efficiency,
        max_flooding_factor=max_flooding_factor,
        tray_type=tray_type,
        contactor_internals_type=contactor_internals_type,
        material_grade=material_grade,
        max_operation_pressure=max_operation_pressure,
    )


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
    temperature_tolerance: float,
    max_iterations: int,
    tray_temperatures: list[float] | None = None,
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
    """Solve a tray absorber, or a stripper.

    ``AbsorptionColumn extends DistillationColumn`` and overrides no ``run``, so this is the
    column's counter-current stage solve with **no condenser, no reboiler and two inlets**: the
    gas enters stage 0 and the solvent the top stage, which is what ``addGasInStream`` and
    ``addSolventInStream`` call. A stripper is the same machine with its inlets named for what
    they carry.

    The answer is the profile - each tray's temperature, pressure and both traffic rates - and
    the class's own two products: ``getGasOutStream`` is the treated gas and
    ``getLiquidOutStream`` the loaded solvent.

    **The class's own tests pin every tray, and a pinned column stops after one sweep**:
    ``SimpleTray.setOutletTemperature`` makes the base's own gate - the mean tray-temperature
    change - exactly zero. ``tray_temperatures`` is that mechanism, and the capture's pinned
    rows are kept as evidence rather than as the oracle.

    **Declared and refused by name**: the Murphree efficiency of either kind, and every solving
    strategy this port does not carry. ``max_allowable_gas_load_factor`` is accepted and does
    not enter the solve - ``isGasLoadFactorWithinDesignLimit`` reads it and no part of ``run``
    does.

    See :func:`azoth.process.reference.absorption_column`.
    """
    return resolve(_ABSORPTION_COLUMN)(  # type: ignore[no-any-return]
        gas_components=gas_components,
        gas_n=gas_n,
        gas_z=gas_z,
        gas_p=gas_p,
        gas_t=gas_t,
        solvent_components=solvent_components,
        solvent_n=solvent_n,
        solvent_z=solvent_z,
        solvent_p=solvent_p,
        solvent_t=solvent_t,
        number_of_stages=number_of_stages,
        top_pressure=top_pressure,
        bottom_pressure=bottom_pressure,
        temperature_tolerance=temperature_tolerance,
        max_iterations=max_iterations,
        tray_temperatures=tray_temperatures,
        murphree_efficiency=murphree_efficiency,
        component_murphree_efficiency=component_murphree_efficiency,
        max_allowable_gas_load_factor=max_allowable_gas_load_factor,
        reactive=reactive,
        reactive_start_tray=reactive_start_tray,
        reactive_end_tray=reactive_end_tray,
        solver_type=solver_type,
        tray_murphree_efficiency=tray_murphree_efficiency,
        gas_side_draw_fractions=gas_side_draw_fractions,
        liquid_side_draw_fractions=liquid_side_draw_fractions,
        pumparound_fractions=pumparound_fractions,
        side_draw_flow_tray=side_draw_flow_tray,
        side_draw_flow_phase=side_draw_flow_phase,
        side_draw_flow_target=side_draw_flow_target,
        side_draw_flow_tolerance=side_draw_flow_tolerance,
        side_draw_flow_max_iterations=side_draw_flow_max_iterations,
        pumparound_return_tray=pumparound_return_tray,
        pumparound_draw_tray=pumparound_draw_tray,
        pumparound_draw_fraction=pumparound_draw_fraction,
        pumparound_temperature_drop=pumparound_temperature_drop,
        pumparound_tolerance=pumparound_tolerance,
        pumparound_max_iterations=pumparound_max_iterations,
        column_diameter=column_diameter,
        max_allowable_fs_factor=max_allowable_fs_factor,
        tray_efficiency=tray_efficiency,
        max_flooding_factor=max_flooding_factor,
        tray_type=tray_type,
        contactor_internals_type=contactor_internals_type,
        material_grade=material_grade,
        max_operation_pressure=max_operation_pressure,
    )


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
    tray_murphree_efficiency: list[float] | None = None,
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

    ``PackedColumn extends DistillationColumn`` and its ``run`` is ``super.run(id)`` followed by
    ``calcPackingHydraulics()``, so **the separation is the base column's and the packing does
    not change it**. What the packing does change is the stage count, at construction:
    ``estimateStages(packed_height, 0.5)`` is ``ceil(packed_height / 0.5)`` floored at two, so
    ``packed_height`` replaces the base's ``number_of_stages`` here.

    The five packing parameters are read only by the hydraulics report on the far side of the
    solve, and that report is mechanical design and is not ported, so they are declarations this
    solve is indifferent to - except ``packing_hydraulic_capacity_factor``, which the class's own
    setter refuses where it is not positive and finite.

    See :func:`azoth.process.reference.packed_column`.
    """
    return resolve(_PACKED_COLUMN)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        packed_height=packed_height,
        feed_stage=feed_stage,
        has_reboiler=has_reboiler,
        has_condenser=has_condenser,
        top_pressure=top_pressure,
        bottom_pressure=bottom_pressure,
        temperature_tolerance=temperature_tolerance,
        max_iterations=max_iterations,
        reboiler_temperature=reboiler_temperature,
        condenser_temperature=condenser_temperature,
        packing_type=packing_type,
        structured_packing=structured_packing,
        design_flood_fraction=design_flood_fraction,
        packing_hydraulic_capacity_factor=packing_hydraulic_capacity_factor,
        column_diameter=column_diameter,
        murphree_efficiency=murphree_efficiency,
        tray_murphree_efficiency=tray_murphree_efficiency,
        reactive=reactive,
        reactive_start_tray=reactive_start_tray,
        reactive_end_tray=reactive_end_tray,
        solver_type=solver_type,
        top_specification_type=top_specification_type,
        top_specification_target=top_specification_target,
        top_specification_component=top_specification_component,
        bottom_specification_type=bottom_specification_type,
        bottom_specification_target=bottom_specification_target,
        bottom_specification_component=bottom_specification_component,
        gas_side_draw_fractions=gas_side_draw_fractions,
        liquid_side_draw_fractions=liquid_side_draw_fractions,
        pumparound_fractions=pumparound_fractions,
        side_draw_flow_tray=side_draw_flow_tray,
        side_draw_flow_phase=side_draw_flow_phase,
        side_draw_flow_target=side_draw_flow_target,
        side_draw_flow_tolerance=side_draw_flow_tolerance,
        side_draw_flow_max_iterations=side_draw_flow_max_iterations,
        pumparound_return_tray=pumparound_return_tray,
        pumparound_draw_tray=pumparound_draw_tray,
        pumparound_draw_fraction=pumparound_draw_fraction,
        pumparound_temperature_drop=pumparound_temperature_drop,
        pumparound_tolerance=pumparound_tolerance,
        pumparound_max_iterations=pumparound_max_iterations,
        max_allowable_fs_factor=max_allowable_fs_factor,
    )


def rate_based_packed_column(
    gas_components: list[str],
    gas_n: Q,
    gas_z: list[float],
    gas_p: Q,
    gas_t: Q,
    liquid_components: list[str],
    liquid_n: Q,
    liquid_z: list[float],
    liquid_p: Q,
    liquid_t: Q,
    transfer_components: list[str] | None = None,
    column_diameter: Q | None = None,
    packed_height: Q | None = None,
    number_of_segments: int | None = None,
    packing_type: str | None = None,
    max_iterations: int | None = None,
    convergence_tolerance: Q | None = None,
    mass_transfer_correction: float | None = None,
    mass_transfer_correlation: str | None = None,
    film_model: str | None = None,
    heat_transfer_model: str | None = None,
    segment_solver: str | None = None,
    column_solver: str | None = None,
) -> RateBasedPackedColumnResult:
    """Solve a rate-based packed column.

    **A second physics rather than a packed column**: ``RateBasedPackedColumn`` extends
    ``ProcessEquipmentBaseClass`` and its 4,081 lines are a segment model - axial slices through
    the packing, each with gas and liquid film coefficients, an interface equilibrium and an
    interphase heat balance, coupled counter-currently by a profile solver. There is no stage, no
    condenser, no reboiler and no MESH residual; the four streams are the whole interface.

    **The packing is inside the equations here**, unlike :func:`packed_column`, whose packing is
    read by a report after the solve: the wetted area and the two film coefficients come from
    ``hydraulics.packing_hydraulics`` every segment.

    A profile that reaches its iteration cap is *published* with its last iterate, and
    ``converged`` says so - which is what the class does. A bed of no height is a converged
    solve on its first pass and transfers nothing.

    See :func:`azoth.process.reference.rate_based_packed_column`.
    """
    return resolve(_RATE_BASED_PACKED_COLUMN)(  # type: ignore[no-any-return]
        gas_components=gas_components,
        gas_n=gas_n,
        gas_z=gas_z,
        gas_p=gas_p,
        gas_t=gas_t,
        liquid_components=liquid_components,
        liquid_n=liquid_n,
        liquid_z=liquid_z,
        liquid_p=liquid_p,
        liquid_t=liquid_t,
        transfer_components=transfer_components,
        column_diameter=column_diameter,
        packed_height=packed_height,
        number_of_segments=number_of_segments,
        packing_type=packing_type,
        max_iterations=max_iterations,
        convergence_tolerance=convergence_tolerance,
        mass_transfer_correction=mass_transfer_correction,
        mass_transfer_correlation=mass_transfer_correlation,
        film_model=film_model,
        heat_transfer_model=heat_transfer_model,
        segment_solver=segment_solver,
        column_solver=column_solver,
    )


def stripping_column(
    stripping_gas_components: list[str],
    rich_liquid_components: list[str],
    stripping_gas_n: Q,
    stripping_gas_z: list[float],
    stripping_gas_p: Q,
    stripping_gas_t: Q,
    rich_liquid_n: Q,
    rich_liquid_z: list[float],
    rich_liquid_p: Q,
    rich_liquid_t: Q,
    number_of_stages: int,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float,
    max_iterations: int,
    tray_temperatures: list[float] | None = None,
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
) -> StrippingColumnResult:
    """Strip a rich liquid with a counter-current gas.

    ``StrippingColumn extends AbsorptionColumn`` and its ninety lines rename the two inlets and
    the two products: ``addStrippingGasStream`` is ``addGasInStream`` and
    ``addRichLiquidStream`` is ``addSolventInStream``, so the stripping gas enters stage 0 and
    the rich liquid the top stage. **The class adds no equations** - absorption and stripping
    are one set of counter-current equilibrium-stage equations, and the driving force sets the
    direction of transfer - so this is :func:`absorption_column` under the class's own names, and
    the two products are ``getOverheadGasStream`` and ``getLeanLiquidStream``.

    See :func:`azoth.process.reference.stripping_column`.
    """
    return resolve(_STRIPPING_COLUMN)(  # type: ignore[no-any-return]
        stripping_gas_components=stripping_gas_components,
        rich_liquid_components=rich_liquid_components,
        stripping_gas_n=stripping_gas_n,
        stripping_gas_z=stripping_gas_z,
        stripping_gas_p=stripping_gas_p,
        stripping_gas_t=stripping_gas_t,
        rich_liquid_n=rich_liquid_n,
        rich_liquid_z=rich_liquid_z,
        rich_liquid_p=rich_liquid_p,
        rich_liquid_t=rich_liquid_t,
        number_of_stages=number_of_stages,
        top_pressure=top_pressure,
        bottom_pressure=bottom_pressure,
        temperature_tolerance=temperature_tolerance,
        max_iterations=max_iterations,
        tray_temperatures=tray_temperatures,
        murphree_efficiency=murphree_efficiency,
        component_murphree_efficiency=component_murphree_efficiency,
        max_allowable_gas_load_factor=max_allowable_gas_load_factor,
        reactive=reactive,
        reactive_start_tray=reactive_start_tray,
        reactive_end_tray=reactive_end_tray,
        solver_type=solver_type,
        tray_murphree_efficiency=tray_murphree_efficiency,
        gas_side_draw_fractions=gas_side_draw_fractions,
        liquid_side_draw_fractions=liquid_side_draw_fractions,
        pumparound_fractions=pumparound_fractions,
        side_draw_flow_tray=side_draw_flow_tray,
        side_draw_flow_phase=side_draw_flow_phase,
        side_draw_flow_target=side_draw_flow_target,
        side_draw_flow_tolerance=side_draw_flow_tolerance,
        side_draw_flow_max_iterations=side_draw_flow_max_iterations,
        pumparound_return_tray=pumparound_return_tray,
        pumparound_draw_tray=pumparound_draw_tray,
        pumparound_draw_fraction=pumparound_draw_fraction,
        pumparound_temperature_drop=pumparound_temperature_drop,
        pumparound_tolerance=pumparound_tolerance,
        pumparound_max_iterations=pumparound_max_iterations,
        column_diameter=column_diameter,
        max_allowable_fs_factor=max_allowable_fs_factor,
        tray_efficiency=tray_efficiency,
        max_flooding_factor=max_flooding_factor,
        tray_type=tray_type,
        contactor_internals_type=contactor_internals_type,
        material_grade=material_grade,
        max_operation_pressure=max_operation_pressure,
    )


def stirred_tank_reactor(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    reaction: str,
    limiting_reactant: str,
    conversion: float,
    isothermal: bool,
    reactor_temperature: Q | None = None,
    reactor_pressure: Q | None = None,
    pressure_drop: Q | None = None,
) -> StirredTankReactorResult:
    """React a feed to a conversion of one reactant, then flash the result.

    ``StirredTankReactor.run`` applies its reaction, sets the pressure and holds the
    temperature when it is isothermal, then flashes: a ``TPflash`` at the held temperature
    or a ``PHflash`` at the **feed's own enthalpy**.

    **The reaction is a movement of moles**, not an energy balance: ``react`` takes the
    limiting reactant's moles times ``conversion`` and scales every stoichiometric
    coefficient by the limiting one's, and the heat the reaction releases shows up as the
    temperature the adiabatic flash lands on rather than as a duty. The isothermal branch
    reports the duty it had to supply, which is the outlet's total enthalpy less the feed's.

    **A stoichiometry row naming a substance the feed does not carry is skipped**, which is
    the class's own behaviour, so a reaction whose products are not in the feed moves less
    material than its stoichiometry says.

    See :func:`azoth.process.reference.stirred_tank_reactor`.
    """
    return resolve(_STIRRED_TANK_REACTOR)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        reaction=reaction,
        limiting_reactant=limiting_reactant,
        conversion=conversion,
        isothermal=isothermal,
        reactor_temperature=reactor_temperature,
        reactor_pressure=reactor_pressure,
        pressure_drop=pressure_drop,
    )


def plug_flow_reactor(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    length: float,
    diameter: float,
    number_of_tubes: float,
    energy_mode: str,
    coolant_temperature: Q,
    overall_heat_transfer_coefficient: Q,
    number_of_steps: float,
    integration_method: str,
    property_update_frequency: float,
    thermodynamic_coupling: str,
    reaction: str,
    reaction_orders: list[float],
    rate_type: str,
    pre_exponential_factor: float,
    activation_energy: float,
    temperature_exponent: float,
    heat_of_reaction: float,
    catalyst_bulk_density: Q | None = None,
    catalyst_activity_factor: float | None = None,
    catalyst_particle_diameter: float | None = None,
    catalyst_void_fraction: float | None = None,
    catalyst_molecular_diffusivity: Q | None = None,
    catalyst_effectiveness_enabled: bool | None = None,
    key_component: str | None = None,
) -> PlugFlowReactorResult:
    """March a feed along a tube, reacting as it goes.

    ``PlugFlowReactor.run`` integrates ``[F1..Fn, T, P]`` in a fixed number of steps with a
    scheme it writes out itself - RK4 by default, Euler on request, ``dz = length / steps`` and
    no step-size control - and the answer is the whole profile rather than an outlet row.

    **The properties are frozen by default.** ``thermodynamic_coupling`` is
    ``"frozen_properties"``, so the rate law reads the last re-flash and never one of the RK4
    sub-states; at the default frequency of ten that is why the class's two schemes agree to
    thirteen digits on its own default, and why the step count moves the answer through *where*
    the properties are re-read rather than through the integration order.

    **The pressure row is the empty-tube branch unless a bed is declared.** Supplying
    ``catalyst_bulk_density`` is what sets one, and it is what puts the Ergun equation in the
    pressure row - measured, four orders of magnitude apart on the same gas.

    **The stoichiometry is the reaction data's and the orders are yours.** ``KineticReaction`` is
    a configured object in NeqSim and a declaration carries no coefficient map, so ``reaction``
    resolves against the reaction data exactly as `azoth.process.stirred_tank_reactor`'s does,
    and ``reaction_orders`` supplies the one argument that data does not carry - one entry per
    reactant, in the table's order.

    See :func:`azoth.process.reference.plug_flow_reactor`.
    """
    return resolve(_PLUG_FLOW_REACTOR)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        length=length,
        diameter=diameter,
        number_of_tubes=number_of_tubes,
        energy_mode=energy_mode,
        coolant_temperature=coolant_temperature,
        overall_heat_transfer_coefficient=overall_heat_transfer_coefficient,
        number_of_steps=number_of_steps,
        integration_method=integration_method,
        property_update_frequency=property_update_frequency,
        thermodynamic_coupling=thermodynamic_coupling,
        reaction=reaction,
        reaction_orders=reaction_orders,
        rate_type=rate_type,
        pre_exponential_factor=pre_exponential_factor,
        activation_energy=activation_energy,
        temperature_exponent=temperature_exponent,
        heat_of_reaction=heat_of_reaction,
        catalyst_bulk_density=catalyst_bulk_density,
        catalyst_activity_factor=catalyst_activity_factor,
        catalyst_particle_diameter=catalyst_particle_diameter,
        catalyst_void_fraction=catalyst_void_fraction,
        catalyst_molecular_diffusivity=catalyst_molecular_diffusivity,
        catalyst_effectiveness_enabled=catalyst_effectiveness_enabled,
        key_component=key_component,
    )


def filter(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    pressure_drop: Q,
    pressure_drop_model: str | None = None,
    reference_flow: Q | None = None,
    flow_exponent: float | None = None,
    pressure_drop_curve: list[dict[str, Q]] | None = None,
    media_area: Q | None = None,
    media_bed_depth: Q | None = None,
    media_particle_diameter: Q | None = None,
    media_void_fraction: float | None = None,
    viscosity: Q | None = None,
    performance_curve: list[dict[str, Q]] | None = None,
    particle_size: Q | None = None,
    nominal_removal_efficiency: float | None = None,
    inlet_particle_concentration: float | None = None,
) -> FilterResult:
    """Drop a stream's pressure by a fixed amount at a constant temperature.

    **It holds the temperature**, which is what separates the entry from
    :func:`throttling_valve`: ``Filter.run`` sets the reduced pressure and flashes, so the
    outlet is at the feed's temperature and its enthalpy moves with the pressure - about
    ``24.8`` J/mol per bar on the fluid these cases share.

    A drop larger than the inlet pressure is **clamped rather than refused**, which is what
    ``Filter.run`` does: the outlet lands a millionth of a bar above vacuum and
    ``applied_drop`` reports what was applied, with a warning when it differs from the
    request.

    ``pressure_drop_model`` selects the class's clean-drop model - ``fixed`` by default, then
    ``flow_scaled``, ``tabulated`` and ``ergun`` - and the inputs after it are that model's
    own settings. The beta-ratio ``performance_curve`` (or a scalar
    ``nominal_removal_efficiency``) with an ``inlet_particle_concentration`` enables the
    particle capture.

    See :func:`azoth.process.reference.filter`.
    """
    return resolve(_FILTER)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        pressure_drop=pressure_drop,
        pressure_drop_model=pressure_drop_model,
        reference_flow=reference_flow,
        flow_exponent=flow_exponent,
        pressure_drop_curve=pressure_drop_curve,
        media_area=media_area,
        media_bed_depth=media_bed_depth,
        media_particle_diameter=media_particle_diameter,
        media_void_fraction=media_void_fraction,
        viscosity=viscosity,
        performance_curve=performance_curve,
        particle_size=particle_size,
        nominal_removal_efficiency=nominal_removal_efficiency,
        inlet_particle_concentration=inlet_particle_concentration,
    )


def component_splitter(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: list[float],
) -> ComponentSplitterResult:
    """Divide a stream between two outlets component by component, flashing each.

    **The factor is per component and the outlet count is two.** ``ComponentSplitter.run``
    loops ``for i in 0..2`` and reads ``splitFactor[k]`` as the fraction of component *k* to
    the overhead, so the outlets carry different compositions - which is what separates this
    from :func:`splitter`.

    **The outlet enthalpies are the state's and not the class's**: a fresh NeqSim fluid at
    the first captured row's bottoms gives ``-20309.24832914077`` where the class reports
    ``-14431.31737097011``, the same defect :func:`splitter` records.

    Raises:
        InvalidInputError: where the factors are not one per component, or any is outside
            ``[0, 1]``.

    See :func:`azoth.process.reference.component_splitter`.
    """
    return resolve(_COMPONENT_SPLITTER)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        split_factors=split_factors,
    )


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
    solve_speed: bool | None = None,
    min_speed: float | None = None,
    max_speed: float | None = None,
    limit_speed: bool | None = None,
    chart_type: str | None = None,
    chart_conditions: list[float] | None = None,
    impeller_diameter: Q | None = None,
    surge_curve: list[dict[str, Q]] | None = None,
    surge_control_margin: float | None = None,
    anti_surge_active: bool | None = None,
    surge_control_factor: float | None = None,
) -> CompressorResult:
    """Raise a stream's pressure along an isentrope or the machine's chart.

    **The entropy is derived, not carried.** The record has five fields and ``s`` is not one
    of them: it is a function of ``(T, P, z)``, and the flash at the inlet decides it. At an
    efficiency of one the outlet's entropy equals the inlet's, which is what the captured
    reversible row is for.

    ``h_out = h_in + (h_isentropic - h_in) / eta``. :func:`expander` is the same route with a
    multiplication, because an expansion's isentropic difference is negative.

    ``curve`` selects ``CompressorChart``'s fixed-speed branch instead: the pressure is solved
    from the head, ``speed`` reads the map, and ``chart_type`` names the map class (only the
    default ``fan_law`` is carried). With ``solve_speed`` the speed is the unknown: a Newton
    drives the chart's pressure ratio to ``outlet_pressure``, and ``min_speed``/``max_speed``
    bound it, with ``limit_speed`` pegging an out-of-range answer at the bound.

    ``surge_curve`` is the anti-surge boundary in the same (flow, head) unit as ``curve``.
    With it the result carries the machine's surge read-outs; with ``anti_surge_active`` set,
    the run-time recycle block re-runs the machine above the surge flow and reports the
    fraction and the boosted power.

    Raises:
        InvalidInputError: where the shapes disagree, the efficiency is outside ``(0, 1]``, or
            the chart type is refused.

    See :func:`azoth.process.reference.compressor`.
    """
    return resolve(_COMPRESSOR)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        isentropic_efficiency=isentropic_efficiency,
        outlet_pressure=outlet_pressure,
        curve=curve,
        speed=speed,
        solve_speed=solve_speed,
        min_speed=min_speed,
        max_speed=max_speed,
        limit_speed=limit_speed,
        chart_type=chart_type,
        chart_conditions=chart_conditions,
        impeller_diameter=impeller_diameter,
        surge_curve=surge_curve,
        surge_control_margin=surge_control_margin,
        anti_surge_active=anti_surge_active,
        surge_control_factor=surge_control_factor,
    )


def expander(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
    use_polytropic_calc: bool | None = None,
    polytropic_efficiency: float | None = None,
    number_of_steps: int | None = None,
) -> ExpanderResult:
    """Drop a stream's pressure along an isentrope, multiplying the step by the efficiency.

    ``h_out = h_in + (h_isentropic - h_in) * eta`` - the compressor's rule with the efficiency
    on the other side of the division, because an expansion's isentropic difference is
    negative and dividing would make the machine beat the reversible one at any efficiency
    below one. The captured pair of rows is what holds that: `0.75` of the step is
    ``-618.933`` J/mol and the step divided would be ``-1478.8``.

    Raises:
        InvalidInputError: where the shapes disagree or the efficiency is outside ``(0, 1]``.

    See :func:`azoth.process.reference.expander`.
    """
    return resolve(_EXPANDER)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_pressure=outlet_pressure,
        isentropic_efficiency=isentropic_efficiency,
        use_polytropic_calc=use_polytropic_calc,
        polytropic_efficiency=polytropic_efficiency,
        number_of_steps=number_of_steps,
    )


def map_turbo_expander_compressor(
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
    expander_isentropic_efficiency: float,
    shaft_mode: str | None = None,
    mechanical_efficiency: float | None = None,
    bearing_loss: Q | None = None,
    compressor_outlet_pressure: Q | None = None,
    compressor_isentropic_efficiency: float | None = None,
    compressor_speed: float | None = None,
    min_shaft_speed: float | None = None,
    max_shaft_speed: float | None = None,
    curve: list[dict[str, Q]] | None = None,
    chart_type: str | None = None,
    chart_conditions: list[float] | None = None,
    impeller_diameter: Q | None = None,
) -> MapTurboExpanderResult:
    """Balance a map-driven brake compressor against an expander on one shaft.

    The expander's power, less the mechanical and bearing losses, is what the compressor may
    consume. ``shaft_mode="balanced_speed"`` bisects the shaft speed over the compressor
    map's own measured range until the compressor's power meets the available power - pegged
    at a bound with ``under_power_surge`` or ``over_power_max_speed`` when it cannot.
    ``specified_pressures`` runs the compressor at its own setting and reports the residual.

    See :func:`azoth.process.reference.map_turbo_expander_compressor`.
    """
    return resolve(_MAP_TURBO_EXPANDER_COMPRESSOR)(  # type: ignore[no-any-return]
        expander_components=expander_components,
        compressor_components=compressor_components,
        expander_n=expander_n,
        expander_z=expander_z,
        expander_p=expander_p,
        expander_t=expander_t,
        compressor_n=compressor_n,
        compressor_z=compressor_z,
        compressor_p=compressor_p,
        compressor_t=compressor_t,
        expander_out_pressure=expander_out_pressure,
        expander_isentropic_efficiency=expander_isentropic_efficiency,
        shaft_mode=shaft_mode,
        mechanical_efficiency=mechanical_efficiency,
        bearing_loss=bearing_loss,
        compressor_outlet_pressure=compressor_outlet_pressure,
        compressor_isentropic_efficiency=compressor_isentropic_efficiency,
        compressor_speed=compressor_speed,
        min_shaft_speed=min_shaft_speed,
        max_shaft_speed=max_shaft_speed,
        curve=curve,
        chart_type=chart_type,
        chart_conditions=chart_conditions,
        impeller_diameter=impeller_diameter,
    )


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
    """Balance a turbo-expander and a boost compressor on one shaft.

    **The speed is the only unknown.** The expander's isentropic head comes from one
    ``PSflash`` at ``expander_out_pressure`` and its efficiency from the velocity-ratio
    parabola or the 2-D map; the compressor's absorbed power is its design polytropic head
    scaled by ``(N/N_design)²`` and a Q/N head factor over a Q/N efficiency factor times the
    design efficiency. A damped Newton on ``N`` drives the expander's power to the
    compressor's plus the bearing loss ``10 (N/N_design)²``.

    **``chart`` selects the map branch**, whose rows carry ``(igv, uc, efficiency)``; its
    stage-head grid is not read by the machine. Stating ``expander_out_temperature`` runs
    the inverse: the required efficiency is solved from the target state and the design
    efficiency is back-calculated from the converged ratio.

    Raises:
        InvalidInputError: for a non-positive design or geometry value, an efficiency outside
            ``(0, 1]``, a non-positive isentropic drop or a map that is not a map.

    See :func:`azoth.process.reference.turbo_expander_compressor`.
    """
    return resolve(_TURBO_EXPANDER_COMPRESSOR)(  # type: ignore[no-any-return]
        expander_components=expander_components,
        compressor_components=compressor_components,
        expander_n=expander_n,
        expander_z=expander_z,
        expander_p=expander_p,
        expander_t=expander_t,
        compressor_n=compressor_n,
        compressor_z=compressor_z,
        compressor_p=compressor_p,
        compressor_t=compressor_t,
        expander_out_pressure=expander_out_pressure,
        impeller_diameter=impeller_diameter,
        design_speed=design_speed,
        design_uc=design_uc,
        design_qn=design_qn,
        igv_opening=igv_opening,
        igv_control_mode=igv_control_mode,
        maximum_igv_area=maximum_igv_area,
        igv_area_increase_factor=igv_area_increase_factor,
        bearing_loss=bearing_loss,
        expander_design_isentropic_efficiency=expander_design_isentropic_efficiency,
        compressor_design_polytropic_efficiency=compressor_design_polytropic_efficiency,
        compressor_design_polytropic_head=compressor_design_polytropic_head,
        expander_out_temperature=expander_out_temperature,
        uc_curve=uc_curve,
        qn_efficiency_curve=qn_efficiency_curve,
        qn_head_curve=qn_head_curve,
        chart=chart,
        igv_penalty_curve=igv_penalty_curve,
    )


def cooler(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
) -> CoolerResult:
    """Cool a stream to a stated temperature, or by a stated duty.

    ``unit_ops.cooler``'s three parameters are ``Heater``'s three setters, reached through
    ``Cooler`` - and ``Cooler`` adds no steady-state arithmetic to ``Heater``, so this is the
    same call as :func:`heater` under the entry a palette shows for a machine that removes
    heat. ``duty`` is signed and unchecked: ``-5000`` W is the usual case here, and a
    positive value heats.

    Raises:
        InvalidInputError: for both specifications at once, a pressure drop that leaves a
            non-positive pressure, or a duty on a stream carrying no flow.

    See :func:`azoth.process.reference.cooler`.
    """
    return resolve(_COOLER)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_temperature=outlet_temperature,
        duty=duty,
        pressure_drop=pressure_drop,
    )


def heater(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
) -> HeaterResult:
    """Heat or cool a stream to a stated temperature, or by a stated duty.

    ``components``, ``inlet_n``, ``inlet_z``, ``inlet_p`` and ``inlet_t`` are the inlet's
    record, and ``outlet_temperature``, ``duty`` and ``pressure_drop`` are
    ``unit_ops.heater``'s own three parameters - all optional, and the pressure drop applied
    to the inlet before the flash.

    **A temperature and a duty together are refused**, because ``Heater.setOutletTemperature``
    and ``setDuty`` clear each other's flags: which one the class runs is the order they were
    set in, which these arguments cannot express. With neither, the drop is isothermal,
    because ``run``'s else branch is ``T_in + dT`` with ``dT`` zero.

    Raises:
        InvalidInputError: for both specifications at once, a pressure drop that leaves a
            non-positive pressure, or a duty on a stream carrying no flow.

    See :func:`azoth.process.reference.heater`.
    """
    return resolve(_HEATER)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_temperature=outlet_temperature,
        duty=duty,
        pressure_drop=pressure_drop,
    )


def pipe(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    length: Q,
    diameter: Q,
    roughness: Q,
    fittings: list[dict[str, Q]] | None = None,
    outlet_pressure: Q | None = None,
    correlation: str | None = None,
    number_of_increments: int | None = None,
    elevation: Q | None = None,
) -> PipeResult:
    """Drop a stream's pressure along a line, solving the outlet pressure it implies.

    ``length``, ``diameter`` and ``roughness`` are ``unit_ops.pipe``'s three parameters, and
    the outlet pressure is *not* one of them: a line's drop is a function of its geometry and
    its fluid, which is what makes a pipe a different shape of unit operation from the
    two-port machines around it.

    **Two equations, one per phase**, chosen by the phase label: a gas takes the compressible
    ``P1^2 - P2^2`` form and anything else Darcy-Weisbach. ``fittings`` adds their ``sum(L/D) D``
    to the physical length; ``outlet_pressure`` replaces the solved drop with a solved *flow*,
    ``AdiabaticPipe.calcFlow``'s bisection; ``correlation = "beggs_and_brills"`` runs the ported
    Beggs-and-Brills march over the same geometry instead.

    Raises:
        InvalidInputError: where the geometry is not positive.

    See :func:`azoth.process.reference.pipe`.
    """
    return resolve(_PIPE)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        length=length,
        diameter=diameter,
        roughness=roughness,
        fittings=fittings,
        outlet_pressure=outlet_pressure,
        correlation=correlation,
        number_of_increments=number_of_increments,
        elevation=elevation,
    )


def splitter(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: list[float],
) -> SplitterResult:
    """Split a stream into several with the same state, in proportion to ``split_factors``.

    ``feed_n``, ``feed_z``, ``feed_p`` and ``feed_t`` are the inlet's record and
    ``split_factors`` is ``unit_ops.splitter``'s own parameter, whose length is how many
    outlets there are. A splitter changes no state: every outlet carries the inlet's
    composition, pressure, temperature and molar enthalpy, and the factors divide the
    flow between them. They are normalised, so only their ratios mean anything.

    NeqSim's ``Splitter.run`` returns outlet enthalpies that do not conserve energy, so
    this diverges deliberately rather than reproducing them.

    Raises:
        InvalidInputError: where there is no outlet, or the factors do not sum to a
            positive value.

    See :func:`azoth.process.reference.splitter`.
    """
    return resolve(_SPLITTER)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        split_factors=split_factors,
    )


def mixer(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
    outlet_pressure: Q | None = None,
) -> MixerResult:
    """Join several streams into one, conserving molar flow and enthalpy.

    ``feed_n``, ``feed_z``, ``feed_p`` and ``feed_t`` are the feeds' records - one entry
    per feed, the fluid's ``components`` named once - and ``outlet_pressure`` is
    ``unit_ops.mixer``'s own parameter, optional because a mixer joins at the feeds'
    lowest pressure unless a header's pressure is stated.

    Mixing is isenthalpic, so the outlet's enthalpy is the feeds' weighted by their flows
    and its temperature is the one at which the joined mixture carries it at the outlet
    pressure. That is why the outlet is at neither feed's temperature in general.

    Raises:
        InvalidInputError: where the feeds' shapes disagree, no feed carries anything, or
            the outlet pressure is not positive.

    See :func:`azoth.process.reference.mixer`.
    """
    return resolve(_MIXER)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        outlet_pressure=outlet_pressure,
    )


def gas_scrubber(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    internal_diameter: Q | None = None,
    design_gas_load_factor: Q | None = None,
    heat_input: Q | None = None,
) -> GasScrubberResult:
    """Flash a feed into vapour and liquid, as a scrubber does.

    ``GasScrubber extends Separator`` and does not override ``run``, so this is
    :func:`separator`'s arithmetic under the other entry, and its three parameters are the
    separator's. **The Souders-Brown capacity metric is the class's own addition**, computed
    when ``internal_diameter`` and ``design_gas_load_factor`` are both given: stated one at a
    time it is refused, because a vessel with a size and no design basis has nothing to check.

    Raises:
        InvalidInputError: where the shapes disagree, the pressure drop leaves a
            non-positive pressure, the entrainment fraction is outside ``[0, 1]``, only one
            of the two mechanical inputs is given, or the state makes the metric undefined.

    See :func:`azoth.process.reference.gas_scrubber`.
    """
    return resolve(_GAS_SCRUBBER)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        pressure_drop=pressure_drop,
        gas_in_liquid=gas_in_liquid,
        internal_diameter=internal_diameter,
        design_gas_load_factor=design_gas_load_factor,
        heat_input=heat_input,
    )


def flare(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
) -> FlareResult:
    """A flare's steady state: the record through, and the two numbers beside it.

    ``Flare.run`` clones the inlet into the outlet and changes nothing in it, so the record is
    a pass-through. What the machine computes is its own report: ``heatDuty``, the gas's
    inferior calorific value per normal cubic metre at 0 C times its flow in standard cubic
    metres per second at 15 C, and ``co2Emission``, the carbon the gas carries times
    ``44.01e-3`` kg/mol.

    **The duty multiplies two reference states**: the calorific value is per *normal* cubic
    metre at 0 C and the flow is at 15 C, which the class's own arithmetic does and this
    reproduces - the case records the ratio.

    Raises:
        InvalidInputError: where a component has no row in ISO 6976's table or none in the
            element table.

    See :func:`azoth.process.reference.flare`.
    """
    return resolve(_FLARE)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
    )


def tank(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
) -> TankResult:
    """Join a tank's inlets and split the result into a gas and a liquid outlet.

    **It declares no parameters, and the one it looked like it needed is why.**
    ``Tank.run`` flashes at the *fluid's* own volume and internal energy, so the steady
    state is the flash the feed already carries - measured, the captured two-phase row is
    :func:`separator`'s first row to the last digit. ``setVolume`` is read by the mechanical
    design, the capacity report and the JSON dump and by nothing in ``run``, so a declared
    ``volume`` would be a parameter no steady-state path reads.

    **A tank is not a three-phase machine**: ``run`` never enables ``multiPhaseCheck``, so a
    feed a three-phase separator splits three ways comes back in two, with the water in the
    liquid.

    Raises:
        InvalidInputError: where the feeds' shapes disagree or carry no flow in total.

    See :func:`azoth.process.reference.tank`.
    """
    return resolve(_TANK)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
    )


def ejector(
    motive_components: list[str],
    motive_n: Q,
    motive_z: list[float],
    motive_p: Q,
    motive_t: Q,
    suction_components: list[str],
    suction_n: Q,
    suction_z: list[float],
    suction_p: Q,
    suction_t: Q,
    discharge_pressure: Q,
    motive_nozzle_efficiency: float,
    suction_nozzle_efficiency: float,
    mixing_efficiency: float,
    diffuser_efficiency: float,
) -> EjectorResult:
    """Expand a motive stream, entrain a suction stream with it, and diffuse the mixture.

    **A gas ejector on NeqSim's quasi one-dimensional route**: each stream expands
    isentropically to the mixing pressure, drops to the enthalpy its own nozzle efficiency
    leaves, is given the velocity ``sqrt(2 dh)`` that implies, and the two momenta mix; a
    diffuser then recovers the mixing velocity's kinetic energy up to ``discharge_pressure``.

    **The mixing pressure is the class's own estimate**, not an input - ``setMixingPressure``
    exists and the palette does not declare it - and so are the two design velocities the
    outlet's final ``v^2/2`` is taken at.

    Raises:
        InvalidInputError: where an efficiency is outside ``(0, 1]`` or a stream's shapes
            disagree.

    See :func:`azoth.process.reference.ejector`.
    """
    return resolve(_EJECTOR)(  # type: ignore[no-any-return]
        motive_components=motive_components,
        motive_n=motive_n,
        motive_z=motive_z,
        motive_p=motive_p,
        motive_t=motive_t,
        suction_components=suction_components,
        suction_n=suction_n,
        suction_z=suction_z,
        suction_p=suction_p,
        suction_t=suction_t,
        discharge_pressure=discharge_pressure,
        motive_nozzle_efficiency=motive_nozzle_efficiency,
        suction_nozzle_efficiency=suction_nozzle_efficiency,
        mixing_efficiency=mixing_efficiency,
        diffuser_efficiency=diffuser_efficiency,
    )


def three_phase_separator(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_aqueous: float = 0.0,
    gas_in_oil: float = 0.0,
    oil_in_aqueous: float = 0.0,
    oil_in_gas: float = 0.0,
    aqueous_in_gas: float = 0.0,
    aqueous_in_oil: float = 0.0,
    heat_input: Q | None = None,
) -> ThreePhaseSeparatorResult:
    """Flash a stream into vapour, oil and aqueous outlets.

    **The three-phase vessel turns `multiPhaseCheck` on for its flash**, so a feed that
    separates into a gas, an oil and an aqueous phase is split three ways - a two-phase feed
    into whichever two it has. The vessel holds the feed's temperature, so a pressure drop is
    not a throttling; ``heat_input`` is the one thing that moves the flash.

    **The six entrainment fractions are declared, and their order is part of the answer**:
    each moves a share of the *current* from-phase's moles, component by component, in
    ``run``'s own order, so a later pair reads what an earlier one left.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero or a fraction
            is outside ``[0, 1]``.

    See :func:`azoth.process.reference.three_phase_separator`.
    """
    return resolve(_THREE_PHASE_SEPARATOR)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        pressure_drop=pressure_drop,
        gas_in_aqueous=gas_in_aqueous,
        gas_in_oil=gas_in_oil,
        oil_in_aqueous=oil_in_aqueous,
        oil_in_gas=oil_in_gas,
        aqueous_in_gas=aqueous_in_gas,
        aqueous_in_oil=aqueous_in_oil,
        heat_input=heat_input,
    )


def hem_release(
    components: list[str],
    upstream_z: list[float],
    upstream_t: Q,
    upstream_p: Q,
    diameter: Q,
    discharge_coefficient: float,
    back_pressure: Q,
) -> HemReleaseResult:
    """The release rate of a homogeneous-equilibrium short-opening orifice.

    The upstream stagnation state is expanded isentropically to a 49-pressure grid, and the
    maximum of ``rho sqrt(2 (h0 - h))`` over that grid is the throat. The answer is **four
    stations** - the stagnation state, the throat, the aliased exit and the receiving-pressure
    state - with the rate, the choked flag and, when the flow is choked, the throat's
    equilibrium sound speed.

    A failure is the class's own result and not an exception: ``status`` is ``invalid`` or
    ``unsupported``, the rate and stations are absent, and one diagnostic says why. **No NaN
    crosses the boundary.**

    Raises:
        InvalidInputError: for a non-positive state or opening, a discharge coefficient
            outside ``(0, 1]``, or a composition that is not one entry per component and sums
            to one.

    See :func:`azoth.process.reference.hem_release`.
    """
    return resolve(_HEM_RELEASE)(  # type: ignore[no-any-return]
        components=components,
        upstream_z=upstream_z,
        upstream_t=upstream_t,
        upstream_p=upstream_p,
        diameter=diameter,
        discharge_coefficient=discharge_coefficient,
        back_pressure=back_pressure,
    )


def manifold(
    components: list[str],
    feed_n: list[Q],
    feed_z: list[list[float]],
    feed_p: list[Q],
    feed_t: list[Q],
    split_factors: list[float],
) -> ManifoldResult:
    """Join a manifold's feeds and divide the mixture between its outlets.

    ``Manifold.run`` is a mixer and a splitter composed - ``localmixer.run()``, then
    ``localsplitter.run()`` over the mixture - so every rule :func:`mixer` and
    :func:`splitter` state holds here. The one rule that is the manifold's own is the
    low-flow filter: a feed whose mass flow is at or below ``1e-20`` kg/hr is not mixed, and
    the palette declares no way to change that.

    The outlet count is ``split_factors``' length, which is why ``unit_ops.manifold``'s
    ``outlet`` is a ``many`` one: an entry with a single outlet could not describe the
    machine.

    Raises:
        InvalidInputError: where the feeds' shapes disagree, a factor is negative, or no
            feed is above the low-flow threshold.

    See :func:`azoth.process.reference.manifold`.
    """
    return resolve(_MANIFOLD)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        split_factors=split_factors,
    )


def separator(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    heat_input: Q | None = None,
) -> SeparatorResult:
    """Flash a stream into vapour and liquid outlets.

    ``feed_n``, ``feed_z``, ``feed_p`` and ``feed_t`` are the inlet's record and
    ``pressure_drop``, ``heat_input`` and ``gas_in_liquid`` are ``unit_ops.separator``'s
    own parameters. The flash is at the **feed's** temperature - the vessel holds it - so a
    pressure drop is not a throttling and the outlets do not carry the inlet's enthalpy; a
    ``heat_input`` is the one thing that moves it.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero or the
            entrainment fraction is outside ``[0, 1]``.

    See :func:`azoth.process.reference.separator`.
    """
    return resolve(_SEPARATOR)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        pressure_drop=pressure_drop,
        gas_in_liquid=gas_in_liquid,
        heat_input=heat_input,
    )


def shortcut_distillation_column(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    light_key: str,
    heavy_key: str,
    light_key_recovery_distillate: float,
    heavy_key_recovery_bottoms: float,
    reflux_ratio_multiplier: float,
    condenser_pressure: Q | None = None,
    reboiler_pressure: Q | None = None,
) -> ShortcutDistillationColumnResult:
    """Solve a Fenske-Underwood-Gilliland shortcut column.

    ``light_key`` and ``heavy_key`` name the two components the separation is written
    between, and their recoveries plus ``reflux_ratio_multiplier`` are the four numbers the
    class takes. The answer is **not a profile but a stage count, a feed tray and two
    duties**: Fenske's minimum stages, Underwood's minimum reflux, Molokanov's fit to
    Gilliland for the actual stages, Kirkbride for the feed tray.

    Three things the class does that the numbers alone do not show, and all three are
    reproduced rather than corrected. The **duties are estimates**: a hard-coded 30000 J/mol
    average latent heat for the condenser and one per cent of the feed's enthalpy for the
    reboiler, so the condenser's figure does not depend on what the fluid is. The
    **Kirkbride argument is the class's own**, with the square on the bottoms-to-distillate
    ratio, where the correlation as usually quoted squares the composition ratio. And a
    **reflux multiplier at or below one is refused**, because the class returns an infinite
    stage count at exactly one and a silent fallback below it.

    Raises:
        InvalidInputError: for a key that is not a component, a relative volatility at or
            below one, a reflux multiplier at or below one, a split that empties a product,
            or a feed whose flash is a trivial solution.

    See :func:`azoth.process.reference.shortcut_distillation_column`.
    """
    return resolve(_SHORTCUT_DISTILLATION_COLUMN)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        light_key=light_key,
        heavy_key=heavy_key,
        light_key_recovery_distillate=light_key_recovery_distillate,
        heavy_key_recovery_bottoms=heavy_key_recovery_bottoms,
        reflux_ratio_multiplier=reflux_ratio_multiplier,
        condenser_pressure=condenser_pressure,
        reboiler_pressure=reboiler_pressure,
    )


def throttling_valve(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    valve_opening: float | None = None,
    kv: float | None = None,
    fouling_fraction: float | None = None,
    valve_characterization: str | None = None,
    xt: float | None = None,
    allow_choked: bool | None = None,
    is_calc_pressure: bool | None = None,
    is_isothermal: bool | None = None,
    sizing_standard: str | None = None,
) -> ThrottlingValveResult:
    """Drop a stream to a lower pressure without heat or work.

    ``inlet_n``, ``inlet_z``, ``inlet_p`` and ``inlet_t`` are the inlet's record and
    ``outlet_pressure`` is ``unit_ops.throttling_valve``'s own parameter. The drop is
    isenthalpic, so the outlet's temperature is the one the mixture reaches at the outlet
    pressure and its enthalpy is the inlet's.

    ``outlet_pressure`` is taken at face value even when it is above the inlet: NeqSim's
    ``acceptNegativeDP`` defaults to ``true`` and the capture's third case records it.

    See :func:`azoth.process.reference.throttling_valve`.
    """
    return resolve(_THROTTLING_VALVE)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_pressure=outlet_pressure,
        valve_opening=valve_opening,
        kv=kv,
        fouling_fraction=fouling_fraction,
        valve_characterization=valve_characterization,
        xt=xt,
        allow_choked=allow_choked,
        is_calc_pressure=is_calc_pressure,
        is_isothermal=is_isothermal,
        sizing_standard=sizing_standard,
    )


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
) -> HeatExchangerResult:
    """Exchange heat between two streams.

    **The two sides carry different fluids**, so this is the one unit operation with two
    component lists rather than one. ``hot_in_*`` and ``cold_in_*`` are the two inlets'
    records and ``ua``, ``flow_arrangement``, ``hot_outlet_temperature`` and
    ``cold_outlet_temperature`` are ``unit_ops.heat_exchanger``'s own parameters, of which
    the last three describe two modes.

    The default mode is the **effectiveness-NTU rating** its NeqSim class runs: ``ua`` and
    ``flow_arrangement`` size the exchanger and the duty is what comes out. Pinning one
    outlet temperature is the other mode, and it energy-balances the opposite side.

    Raises:
        InvalidInputError: where neither or both outlet temperatures are given, where the
            rating is asked for without a ``ua``, or where the arrangement is not one of
            ``counterflow``, ``parallelflow`` and ``shell_and_tube``.

    See :func:`azoth.process.reference.heat_exchanger`.
    """
    return resolve(_HEAT_EXCHANGER)(  # type: ignore[no-any-return]
        hot_components=hot_components,
        hot_in_n=hot_in_n,
        hot_in_z=hot_in_z,
        hot_in_p=hot_in_p,
        hot_in_t=hot_in_t,
        cold_components=cold_components,
        cold_in_n=cold_in_n,
        cold_in_z=cold_in_z,
        cold_in_p=cold_in_p,
        cold_in_t=cold_in_t,
        ua=ua,
        flow_arrangement=flow_arrangement,
        hot_outlet_temperature=hot_outlet_temperature,
        cold_outlet_temperature=cold_outlet_temperature,
    )


def gibbs_reactor(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    energy_mode: str,
    damping_composition: float,
    max_iterations: float,
    convergence_tolerance: float,
    min_iterations: float,
    reactor_type: str | None = None,
    catalyst_bulk_density: float | None = None,
    catalyst_particle_diameter: float | None = None,
    bed_void_fraction: float | None = None,
    ghsv: float | None = None,
    lhsv: float | None = None,
    number_of_beds: int | None = None,
    max_operation_pressure: float | None = None,
    max_operation_temperature: float | None = None,
    joint_efficiency: float | None = None,
    corrosion_allowance: float | None = None,
    armijo_line_search: bool | None = None,
    armijo_c1: float | None = None,
    armijo_rho: float | None = None,
    armijo_max_backtracks: int | None = None,
    adaptive_step_size: bool | None = None,
    regularization: bool | None = None,
    regularization_threshold: float | None = None,
    regularization_tau: float | None = None,
) -> GibbsReactorResult:
    """Bring a feed to its Gibbs equilibrium at its own temperature and pressure.

    ``GibbsReactor.run`` carries its **own** Lagrange-multiplier Newton iteration over the
    element balances and its **own** species database, and neither is ``ChemicalEquilibrium``'s -
    so this is not a composition of P10's minimiser and the two would answer differently.

    **The equilibrium temperature is the feed's, and there is no setter.** The class reads
    ``system.getTemperature()`` from the fluid it is handed.

    **The convergence test is on the undamped step and the floor bites.** ``min_iterations`` is
    the pass below which the tolerance is not consulted, and three of the class's own capture
    rows stop there rather than on the tolerance. A run that reaches ``max_iterations`` still
    succeeds, with ``converged`` false.

    **The vessel is ``ReactorMechanicalDesign``, and it is read from the feed.** The eleven
    mechanical inputs default to the class's own field initialisers, so omitting them gives the
    class's design rather than no design; ``ghsv`` and ``lhsv`` are reciprocal hours and cross
    unconverted, the unit vocabulary having no reciprocal-time unit.

    See :func:`azoth.process.reference.gibbs_reactor`.
    """
    return resolve(_GIBBS_REACTOR)(  # type: ignore[no-any-return]
        components=components,
        feed_n=feed_n,
        feed_z=feed_z,
        feed_p=feed_p,
        feed_t=feed_t,
        energy_mode=energy_mode,
        damping_composition=damping_composition,
        max_iterations=max_iterations,
        convergence_tolerance=convergence_tolerance,
        min_iterations=min_iterations,
        reactor_type=reactor_type,
        catalyst_bulk_density=catalyst_bulk_density,
        catalyst_particle_diameter=catalyst_particle_diameter,
        bed_void_fraction=bed_void_fraction,
        ghsv=ghsv,
        lhsv=lhsv,
        number_of_beds=number_of_beds,
        max_operation_pressure=max_operation_pressure,
        max_operation_temperature=max_operation_temperature,
        joint_efficiency=joint_efficiency,
        corrosion_allowance=corrosion_allowance,
        armijo_line_search=armijo_line_search,
        armijo_c1=armijo_c1,
        armijo_rho=armijo_rho,
        armijo_max_backtracks=armijo_max_backtracks,
        adaptive_step_size=adaptive_step_size,
        regularization=regularization,
        regularization_threshold=regularization_threshold,
        regularization_tau=regularization_tau,
    )
