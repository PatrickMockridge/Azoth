"""The process layer: unit operations and flowsheets.

The ``process.*`` ids, one per palette entry under ``specs/unit_ops/``, and the check
that holds a flowsheet to the calculus's rules. Each id is a registered model: a spec
under ``specs/models/process/``, a kernel in Rust and a reference here, compared by
``test_cross_impl`` against a committed NeqSim capture.

The Stream-level arithmetic those models wrap is :mod:`azoth.process.kernels`, which is what
a flowsheet's executor calls and what ``azoth check`` validates the wiring of. To *run* one,
:func:`run_flowsheet` in :mod:`azoth.process.flowsheet`.
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
    ManifoldResult,
    MixerResult,
    PackedColumnResult,
    PipeResult,
    PlugFlowReactorResult,
    PumpResult,
    SeparatorResult,
    ShortcutDistillationColumnResult,
    SplitterResult,
    StirredTankReactorResult,
    StrippingColumnResult,
    TankResult,
    ThreePhaseSeparatorResult,
    ThrottlingValveResult,
)
from azoth.core.units import Q
from azoth.process.flowsheet import (
    FlowsheetResult,
    FlowsheetStream,
    Residuals,
    TearResult,
    run_flowsheet,
)
from azoth.process.kernels import Stream

__all__ = [
    "FlowsheetResult",
    "FlowsheetStream",
    "Residuals",
    "Stream",
    "TearResult",
    "absorption_column",
    "component_splitter",
    "compressor",
    "cooler",
    "distillation_column",
    "ejector",
    "expander",
    "filter",
    "flare",
    "gas_scrubber",
    "gibbs_reactor",
    "heat_exchanger",
    "heater",
    "load_flowsheet",
    "manifold",
    "mixer",
    "packed_column",
    "pipe",
    "plug_flow_reactor",
    "pump",
    "run_flowsheet",
    "separator",
    "shortcut_distillation_column",
    "splitter",
    "stirred_tank_reactor",
    "stripping_column",
    "tank",
    "three_phase_separator",
    "throttling_valve",
    "validate",
]

_MANIFOLD = "process.manifold"
_MIXER = "process.mixer"
_PACKED_COLUMN = "process.packed_column"
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
_FILTER = "process.filter"
_FLARE = "process.flare"
_HEATER = "process.heater"
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
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> PumpResult:
    """Raise a stream's pressure, adding the pump's work as enthalpy.

    ``components``, ``inlet_n``, ``inlet_z``, ``inlet_p`` and ``inlet_t`` are the inlet's
    record - the fluid, its flow, composition, pressure and temperature - and
    ``outlet_pressure`` and ``isentropic_efficiency`` are ``unit_ops.pump``'s own two
    parameters. The inlet's molar enthalpy is *not* an input: it is a state function of
    ``(T, P, z)``, so accepting one would let a caller hand over a state that does not exist.

    The work is the isentropic head - ``(H(P_out, s_in) - H(P_in)) / eta``, from a flash -
    and not the incompressible ``v (P_out - P_in)``, which is its linearisation. NeqSim's
    ``Pump.run`` does the former, because ``calculateAsCompressor`` defaults to ``true``.

    Raises:
        InvalidInputError: where the shapes disagree or the efficiency is outside ``(0, 1]``.

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
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    gas_side_draw_fractions: list[float] | None = None,
    liquid_side_draw_fractions: list[float] | None = None,
    pumparound_fractions: list[float] | None = None,
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
        solver_type=solver_type,
        top_specification_type=top_specification_type,
        top_specification_target=top_specification_target,
        top_specification_component=top_specification_component,
        bottom_specification_type=bottom_specification_type,
        bottom_specification_target=bottom_specification_target,
        bottom_specification_component=bottom_specification_component,
        reactive=reactive,
        reactive_start_tray=reactive_start_tray,
        reactive_end_tray=reactive_end_tray,
        gas_side_draw_fractions=gas_side_draw_fractions,
        liquid_side_draw_fractions=liquid_side_draw_fractions,
        pumparound_fractions=pumparound_fractions,
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
    solver_type: str | None = None,
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
        solver_type=solver_type,
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
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
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
        solver_type=solver_type,
        top_specification_type=top_specification_type,
        top_specification_target=top_specification_target,
        top_specification_component=top_specification_component,
        bottom_specification_type=bottom_specification_type,
        bottom_specification_target=bottom_specification_target,
        bottom_specification_component=bottom_specification_component,
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
    solver_type: str | None = None,
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
        solver_type=solver_type,
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

    See :func:`azoth.process.reference.filter`.
    """
    return resolve(_FILTER)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        pressure_drop=pressure_drop,
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
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> CompressorResult:
    """Raise a stream's pressure along an isentrope, dividing the step by the efficiency.

    **The entropy is derived, not carried.** The record has five fields and ``s`` is not one
    of them: it is a function of ``(T, P, z)``, and the flash at the inlet decides it. At an
    efficiency of one the outlet's entropy equals the inlet's, which is what the captured
    reversible row is for.

    ``h_out = h_in + (h_isentropic - h_in) / eta``. :func:`expander` is the same route with a
    multiplication, because an expansion's isentropic difference is negative.

    Raises:
        InvalidInputError: where the shapes disagree or the efficiency is outside ``(0, 1]``.

    See :func:`azoth.process.reference.compressor`.
    """
    return resolve(_COMPRESSOR)(  # type: ignore[no-any-return]
        components=components,
        inlet_n=inlet_n,
        inlet_z=inlet_z,
        inlet_p=inlet_p,
        inlet_t=inlet_t,
        outlet_pressure=outlet_pressure,
        isentropic_efficiency=isentropic_efficiency,
    )


def expander(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
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
) -> PipeResult:
    """Drop a stream's pressure along a line, solving the outlet pressure it implies.

    ``length``, ``diameter`` and ``roughness`` are ``unit_ops.pipe``'s three parameters, and
    the outlet pressure is *not* one of them: a line's drop is a function of its geometry and
    its fluid, which is what makes a pipe a different shape of unit operation from the
    two-port machines around it.

    **Two equations, one per phase**, chosen by the phase label: a gas takes the compressible
    ``P1^2 - P2^2`` form and anything else Darcy-Weisbach. **The aqueous viscosity branch is
    not ported**, so a water line's Reynolds number is ``1.61`` times NeqSim's; the case
    records it.

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
    heat_input: Q | None = None,
) -> GasScrubberResult:
    """Flash a feed into vapour and liquid, as a scrubber does.

    ``GasScrubber extends Separator`` and does not override ``run``, so this is
    :func:`separator`'s arithmetic under the other entry, and its three parameters are the
    separator's. **The Souders-Brown capacity metric is not ported**: it needs an internal
    diameter and a design gas load factor, neither of which the palette declares, and it asks
    whether the vessel is big enough rather than stating anything about the stream.

    Raises:
        InvalidInputError: where the shapes disagree, the pressure drop leaves a
            non-positive pressure, or the entrainment fraction is outside ``[0, 1]``.

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
    )
