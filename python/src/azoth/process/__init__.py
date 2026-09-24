"""The process layer: unit operations and flowsheets.

The ``process.*`` ids, one per palette entry under ``specs/unit_ops/``, and the check
that holds a flowsheet to the calculus's rules. Each id is a registered model: a spec
under ``specs/models/process/``, a kernel in Rust and a reference here, compared by
``test_cross_impl`` against a committed NeqSim capture.

The Stream-level arithmetic those models wrap is :mod:`azoth.process.kernels`, which is
what a flowsheet's executor will call and what ``azoth check`` validates the wiring of.
"""

from __future__ import annotations

import pathlib

from azoth import _core
from azoth._dispatch import resolve
from azoth.core.result import (
    ComponentSplitterResult,
    CompressorResult,
    CoolerResult,
    DistillationColumnResult,
    ExpanderResult,
    FilterResult,
    GasScrubberResult,
    HeaterResult,
    HeatExchangerResult,
    ManifoldResult,
    MixerResult,
    PipeResult,
    PumpResult,
    SeparatorResult,
    ShortcutDistillationColumnResult,
    SplitterResult,
    TankResult,
    ThrottlingValveResult,
)
from azoth.core.units import Q
from azoth.process.kernels import Stream

__all__ = [
    "Stream",
    "component_splitter",
    "compressor",
    "cooler",
    "distillation_column",
    "expander",
    "filter",
    "gas_scrubber",
    "heat_exchanger",
    "heater",
    "load_flowsheet",
    "manifold",
    "mixer",
    "pipe",
    "pump",
    "separator",
    "shortcut_distillation_column",
    "splitter",
    "tank",
    "throttling_valve",
    "validate",
]

_MANIFOLD = "process.manifold"
_MIXER = "process.mixer"
_PIPE = "process.pipe"
_GAS_SCRUBBER = "process.gas_scrubber"
_HEAT_EXCHANGER = "process.heat_exchanger"
_COMPONENT_SPLITTER = "process.component_splitter"
_COMPRESSOR = "process.compressor"
_DISTILLATION_COLUMN = "process.distillation_column"
_COOLER = "process.cooler"
_EXPANDER = "process.expander"
_FILTER = "process.filter"
_HEATER = "process.heater"
_SEPARATOR = "process.separator"
_SHORTCUT_DISTILLATION_COLUMN = "process.shortcut_distillation_column"
_THROTTLING_VALVE = "process.throttling_valve"
_PUMP = "process.pump"
_SPLITTER = "process.splitter"
_TANK = "process.tank"


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
    reboiler_temperature: Q,
    condenser_temperature: Q,
    temperature_tolerance: float,
    max_iterations: int,
    murphree_efficiency: float | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
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
