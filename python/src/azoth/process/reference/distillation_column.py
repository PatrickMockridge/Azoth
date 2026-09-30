"""``process.distillation_column`` - the tray-by-tray substitution solve.

Spec: ``specs/models/process/distillation_column.toml``

The Python twin of ``crates/azoth-process/src/models/distillation_column.rs`` over
``crates/azoth-process/src/kernels/distillation_column.rs``, written to mirror them.

# Three mechanisms that are not tidiness

**The linear temperature seed enters through the inlets.** ``SimpleTray.init()`` re-states
every inlet stream's temperature to the tray's before its first run, which is how
``DistillationColumn.init``'s linear profile reaches the trays at all; the class then restores
the caller-owned feed, which is why the sweeps see the feed's own state.

**The reboiler's ``init`` inlet is the feed stage's liquid**, not the first tray's - ``init``
links it that way and the sweeps then replace it with the first tray's.

**The ends are pinned by temperature, or by a specification.** ``setCondenserTemperature``
reaches the tray's own ``outTemperature``, so that tray's flash is a ``TPflash`` at the pin
rather than at the mixed enthalpy; a middle tray has no pin and flashes at its enthalpy. A
specification is the end's other route, and ``ColumnSpecification``'s five types split in two:
a ``reflux_ratio`` and a ``duty`` are written onto the end itself, and a purity, a recovery
and a flow rate are driven by an **outer secant on the end's temperature**.

**There is no specification homotopy here, and that is a measurement.**
``getEffectiveSpecificationHomotopySteps`` returns three staged targets only when the solver
is ``AUTO``, and the field it otherwise reads is one - so the continuation the class carries
belongs to the solver this port refuses.

# What this does not carry

The relaxation controller damps on a combined residual of the scaled temperature, mass and
energy terms. Two of those are MESH norms - ``ColumnMeshResidualEvaluator``'s - and this port
carries the temperature term alone, which is why the row NeqSim's own deethanizer is written
on converges here and not there.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DistillationColumnResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.process.kernels import Stream
from azoth.process.reference import _column_stage as _stage
from azoth.process.reference import _unported
from azoth.process.reference._column_stage import (
    StreamRecord,
    reactive_stage,
    split_draws,
    stage,
    validate_draws,
)
from azoth.process.reference.capacity import (
    DEFAULT_INTERNAL_DIAMETER_M,
    DEFAULT_MAX_ALLOWABLE_FS_FACTOR,
    fs_limits,
)
from azoth.process.reference.designer import (
    DEFAULT_TRAY_SPACING_M,
    DEFAULT_WEIR_HEIGHT_M,
    DesignerGeometry,
    DesignerReport,
    designer_report,
)

_DrawVector = tuple[float, ...] | None
_Draws = tuple[_DrawVector, _DrawVector, _DrawVector] | None

#: One draw kind's fractions, or `None` where the case states none of that kind.
DrawVector = tuple[float, ...] | None
#: The three draw kinds, in `DRAW_FIELDS`' order - what `_states` takes.
Draws = tuple[DrawVector, DrawVector, DrawVector] | None
#: The pumparound returns the last outer pass carried, one entry per tray.
PumparoundInlets = tuple[StreamRecord | None, ...]
#: One side-draw flow specification, as the helper returns it: tray, phase, target, tolerance,
#: cap.
SideDrawSpecification = tuple[int, str, float, float, int]

#: `ColumnInternalsDesigner`'s own constructor defaults, one per geometry input.
DEFAULT_INTERNALS_TYPE = "sieve"
DEFAULT_HOLE_DIAMETER_MM = 12.7
DEFAULT_HOLE_AREA_FRACTION = 0.1
DEFAULT_DOWNCOMMER_AREA_FRACTION = 0.1
DEFAULT_DESIGNER_FLOOD_FRACTION = 0.8
UNSIZED_COLUMN_DIAMETER_M = -1.0

#: The class's own adaptive-relaxation constants, from `DistillationColumn`'s initialisers.
MIN_SEQUENTIAL_RELAXATION = 0.5
MAX_ADAPTIVE_RELAXATION = 1.2
RELAXATION_INCREASE_FACTOR = 1.2
RELAXATION_DECREASE_FACTOR = 0.5
#: The floor on the *temperature* update's step, a different clamp from the streams'.
MIN_TEMPERATURE_RELAXATION = 0.2
#: `DistillationColumn`'s own tear settings, which the coupling's outer loop runs at.
COLUMN_TEAR_TOLERANCE = 1.0e-4
MAX_COLUMN_TEAR_ITERATIONS = 12
#: `maxPumparoundIterations`' initialiser, which the tear's limit is taken against.
MAX_PUMPAROUND_ITERATIONS = 8
#: The residual above which the tear's variables count as changed.
CHANGED_RESIDUAL = 1.0e-12

#: `DEFAULT_MASS_BALANCE_TOLERANCE` and `DEFAULT_ENTHALPY_BALANCE_TOLERANCE`.
MASS_BALANCE_TOLERANCE = 1.6e-2
ENTHALPY_BALANCE_TOLERANCE = 1.6e-2

#: `ColumnSpecification`'s five types, by the names the model's enum declares.
SPECIFICATION_KINDS = (
    "product_purity",
    "component_recovery",
    "product_flow_rate",
    "reflux_ratio",
    "duty",
)
#: The two the outer loop does not drive: `needsAdjustment` is false for both, because
#: `applyDirectSpecification` writes them onto the end itself.
DIRECT_KINDS = ("reflux_ratio", "duty")
#: `ColumnSpecification`'s own convergence defaults.
SPECIFICATION_TOLERANCE = 1.0e-4
SPECIFICATION_MAX_ITERATIONS = 20
#: `secantStep`'s step cap, its temperature window, and its first iteration's offset.
SECANT_MAX_STEP = 50.0
SECANT_MIN_TEMPERATURE = 100.0
SECANT_MAX_TEMPERATURE = 1000.0
SECANT_FIRST_OFFSET = 5.0

#: The strategies the class carries and this port does not, read from the declaration.
UNPORTED_SOLVERS = _unported.values_for("solver_type")


class Specification(NamedTuple):
    """One product specification: the type, the target and the component.

    `kind` is one of :data:`SPECIFICATION_KINDS`; the target is in the unit its type implies -
    dimensionless for a purity, a recovery or a ratio, **mol/hr** for a flow rate, W for a duty.

    **Its location is the end it is handed to, and not a field.** NeqSim's
    `ColumnSpecification` carries a `ProductLocation`, but `validateColumnSpecification` refuses
    a `TOP` one given to `setBottomSpecification`, so the field can only agree with its slot -
    which is why the model's inputs name the end and there are two of them.
    """

    #: One of `SPECIFICATION_KINDS`.
    kind: str
    #: The target value.
    target: float
    #: The component a purity or a recovery constrains; unread by the other three.
    component: str | None = None

    def value(self, product: StreamRecord, feed: StreamRecord) -> float:
        """The value the product currently has, which is what the error is measured against."""
        if self.kind == "product_purity":
            index = _component_index(product, self.component)
            return float(product["z"][index])
        if self.kind == "component_recovery":
            index = _component_index(product, self.component)
            supplied = feed["n"] * feed["z"][index]
            if abs(supplied) <= 1.0e-12:
                return 0.0
            return float(product["n"] * product["z"][index] / supplied)
        if self.kind == "product_flow_rate":
            # `getFlowRate("mol/hr")`: the class's own target unit for this type.
            return float(product["n"] * 3600.0)
        return 0.0


def _component_index(product: StreamRecord, name: str | None) -> int:
    """The position of a named component, or a refusal naming what was asked for."""
    if name is None:
        raise InvalidInputError(
            "component", "a purity or a recovery constrains a component, and none was stated"
        )
    if name not in product["components"]:
        raise InvalidInputError("component", f"{name} is not one of this fluid's components")
    return list(product["components"]).index(name)


class _States(NamedTuple):
    """The solved column, in SI magnitudes - the profile and both products."""

    #: **What the stages had to fall back on**, one entry per kind, which the model merges with
    #: its own checks' caveats.
    warnings: tuple[Warning, ...]
    #: The vapour each tray withdrew, one entry per tray and zero where it drew none.
    gas_side_draw_n: tuple[float, ...]
    #: The liquid each tray withdrew as a liquid side draw.
    liquid_side_draw_n: tuple[float, ...]
    #: The liquid each tray withdrew as a pumparound.
    pumparound_n: tuple[float, ...]
    #: Each tray's temperature, K, from the reboiler at stage 0 to the condenser.
    tray_temperature: tuple[float, ...]
    #: Each tray's pressure, Pa.
    tray_pressure: tuple[float, ...]
    #: Each tray's vapour traffic, mol/s.
    tray_gas_n: tuple[float, ...]
    #: Each tray's liquid traffic, mol/s.
    tray_liquid_n: tuple[float, ...]
    #: Each tray's vapour composition.
    tray_gas_z: tuple[tuple[float, ...], ...]
    #: Each tray's liquid composition.
    tray_liquid_z: tuple[tuple[float, ...], ...]
    #: The distillate's molar flow.
    distillate_n: float
    #: The distillate's composition.
    distillate_z: tuple[float, ...]
    #: The distillate's molar enthalpy, J/mol.
    distillate_h: float
    #: The distillate's own temperature, K, which is the end's and not the caller's pin.
    distillate_t: float
    #: The distillate's own pressure, Pa.
    distillate_p: float
    #: The bottoms' molar flow.
    bottoms_n: float
    #: The bottoms' composition.
    bottoms_z: tuple[float, ...]
    #: The bottoms' molar enthalpy, J/mol.
    bottoms_h: float
    #: The bottoms' own temperature, K.
    bottoms_t: float
    #: The bottoms' own pressure, Pa.
    bottoms_p: float
    #: The condenser's duty, W.
    condenser_duty: float
    #: The reboiler's duty, W.
    reboiler_duty: float
    #: Iterations taken.
    iterations: int
    #: The mean tray-temperature change at the last iteration, K.
    temperature_residual: float
    #: The products' worst component imbalance against the feed, relative.
    mass_residual: float
    #: The enthalpy closure.
    energy_residual: float


def _distillation_column_states(
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
) -> tuple[_States, list[Warning]]:
    """Solve a distillation column by sequential substitution, and hand back its stages.

    Args:
        components: the feed's substances, by name.
        feed_n: the feed's molar flow.
        feed_z: the feed composition.
        feed_p: the feed's pressure, which is its own - the tray pressures come from
            ``bottom_pressure`` and ``top_pressure``.
        feed_t: the feed's temperature.
        number_of_stages: the trays between the ends.
        feed_stage: the stage the feed enters, 0-based over the trays including the ends.
        has_reboiler: whether stage 0 is a reboiler.
        has_condenser: whether the top stage is a condenser.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at stage 0.
        reboiler_temperature: the reboiler's temperature, which pins the bottom tray. ``None``
            leaves the end to flash at its own enthalpy, which is what makes a duty
            specification reachable at all.
        condenser_temperature: the condenser's temperature, which pins the top tray. ``None`` as
            the reboiler's is. **A pin wins over a directly-applied duty**, because the tray
            tests its stated outlet temperature before it reads a heat input.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        murphree_efficiency: the column-wide Murphree tray efficiency; omitted is the ideal
            stage. Refused with ``solver_type="naphtali_sandholm"``, whose solver carries a
            different correction.
        tray_murphree_efficiency: one override per stage, ``NaN`` where a stage falls through
            to ``murphree_efficiency``; a length that is not ``number_of_stages`` is refused.
        solver_type: **only ``direct_substitution`` is ported**, which is the class's default.
        top_specification_type: the top product's degree of freedom, one of
            :data:`SPECIFICATION_KINDS`; omitted, the top is pinned by temperature instead.
        top_specification_target: its target, in the unit its type implies: dimensionless for a
            purity or a recovery, **mol/hr** for a flow rate, W for a duty.
        top_specification_component: the component a purity or a recovery constrains, which the
            class refuses to see omitted; the other three types do not read it.
        bottom_specification_type: the bottom product's degree of freedom, as the top's.
        bottom_specification_target: its target, in the unit its type implies.
        bottom_specification_component: the component a purity or a recovery constrains.

    Returns:
        The solved stages, in SI magnitudes, and the warnings their own fallbacks raised.
        **The stages and not only the record, because `process.packed_column`'s reference needs
        the middle tray's compositions** - the profile the record carries is temperatures and
        flows, and the packing's report is read at a tray's own fluid. `distillation_column` is
        this and `_record` together, and a caller wanting the record reads that.

    Raises:
        InvalidInputError: for a stage or a feed stage outside the column, for any declared
            parameter whose arithmetic is not ported, and for a ratio that is not one.
        SolverNotConvergedError: when the solve misses its gate.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = distillation_column(
        ...     ["methane", "n-butane"],
        ...     q(7.490704036290964, "mol/s"),
        ...     [0.5, 0.5],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     4,
        ...     2,
        ...     True,
        ...     True,
        ...     q(19.0, "bar"),
        ...     q(20.0, "bar"),
        ...     1e-06,
        ...     200,
        ...     q(373.15, "K"),
        ...     q(253.15, "K"),
        ... )
        >>> round(r.distillate_n.to("mol/s").magnitude, 4)
        3.7915
    """
    efficiency = _murphree(murphree_efficiency, tray_murphree_efficiency, number_of_stages)
    _refuse_unported(solver_type, efficiency)

    section = reactive_section(reactive, reactive_start_tray, reactive_end_tray)
    if section is not None and solver_type == "naphtali_sandholm":
        raise _unported.refuse("reactive@solver_type=naphtali_sandholm")

    top_specification = _build_specification(
        top_specification_type, top_specification_target, top_specification_component, "top"
    )
    bottom_specification = _build_specification(
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        "bottom",
    )

    draws_stated = (
        gas_side_draw_fractions is not None
        or liquid_side_draw_fractions is not None
        or pumparound_fractions is not None
    )
    draws = None
    if draws_stated:
        # **A kind that states nothing stays `None` rather than becoming an empty vector.** The
        # two mean different things: an absent vector is "no tray draws this", which the kernel
        # skips, and an empty one is a vector of the wrong length, which the kernel refuses. The
        # three are independent - a column may draw vapour on one tray and nothing else anywhere -
        # so encoding them all as tuples would refuse a state the kernel accepts.
        draws = (
            None
            if gas_side_draw_fractions is None
            else tuple(float(v) for v in gas_side_draw_fractions),
            None
            if liquid_side_draw_fractions is None
            else tuple(float(v) for v in liquid_side_draw_fractions),
            None if pumparound_fractions is None else tuple(float(v) for v in pumparound_fractions),
        )

    if draws_stated and solver_type == "naphtali_sandholm":
        raise _unported.refuse("gas_side_draw_fractions@solver_type=naphtali_sandholm")

    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    p = input_to_si(spec, "feed_p", feed_p)
    t = input_to_si(spec, "feed_t", feed_t)
    top = input_to_si(spec, "top_pressure", top_pressure)
    bottom = input_to_si(spec, "bottom_pressure", bottom_pressure)
    # **Absent means no pin**, which is what `setCondenserTemperature` not having been called
    # does - and what makes a duty specification reachable at all.
    t_reb = (
        None
        if reboiler_temperature is None
        else input_to_si(spec, "reboiler_temperature", reboiler_temperature)
    )
    t_cond = (
        None
        if condenser_temperature is None
        else input_to_si(spec, "condenser_temperature", condenser_temperature)
    )

    # **A declared input arrives in one of two shapes.** A dimensionless one comes as the bare
    # number it is, and a unit-carrying one as a quantity; `input_to_si` refuses the first, so
    # these four go through the same laxity the bridge's `_si` has.
    stages = int(_si(spec, "number_of_stages", number_of_stages))
    stage = int(_si(spec, "feed_stage", feed_stage))
    tolerance = _si(spec, "temperature_tolerance", temperature_tolerance)
    iterations_cap = int(_si(spec, "max_iterations", max_iterations))

    apply_checks(
        checks.on_input,
        {
            "number_of_stages": float(stages),
            "top_pressure": top,
            "bottom_pressure": bottom,
            "temperature_tolerance": tolerance,
            "feed_t": t,
        }.get,
        warnings,
    )

    # **A specified draw flow is an outer search over whole column solves**, so the tear owns
    # the call and the block below is its inner solve - see `_side_draw_flows`.
    def solve_once(
        active_draws: Draws, active_returns: tuple[StreamRecord | None, ...] | None = None
    ) -> _States:
        """One inner solve at one set of draws and returns, which is what both outer loops
        iterate: the tear moves the draws and the pumparound return moves the returns."""
        return _states(
            components,
            n,
            list(feed_z),
            t,
            p,
            stages,
            stage,
            has_reboiler,
            has_condenser,
            top,
            bottom,
            t_reb,
            t_cond,
            tolerance,
            iterations_cap,
            top_specification,
            bottom_specification,
            solver_type,
            reactive=section,
            draws=active_draws,
            murphree_efficiency=efficiency,
            pumparound_inlets=active_returns,
        )

    flows = _side_draw_flow_specification(
        side_draw_flow_tray,
        side_draw_flow_phase,
        # **The target arrives as a quantity and the helper takes SI**, which is where every
        # other unit-carrying input of this model crosses too.
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
            list(components),
            draws,
            returns,
            _tray_count(stages, has_reboiler, has_condenser),
            1.0e-4 if pumparound_tolerance is None else float(pumparound_tolerance),
            12 if pumparound_max_iterations is None else int(pumparound_max_iterations),
        )
    else:
        states = (
            solve_once(draws)
            if flows is None
            else _side_draw_flow_tear(
                solve_once, list(components), draws, flows, stages, has_reboiler, has_condenser
            )
        )
    warnings.extend(states.warnings)

    return states, warnings


#: The floor `applyHydraulicPressureDrop` clamps a rewritten top pressure to, Pa (1e-6 bara).
MIN_TOP_PRESSURE_PA = 0.1


def _coupled(
    solve: Callable[[Q, Q], tuple[_States, list[Warning]]],
    components: list[str],
    top_pressure: Q,
    bottom_pressure: Q,
    internals_type: str,
) -> tuple[_States, list[Warning]]:
    """``hydraulicPressureDropCouplingEnabled``'s outer loop, mirroring ``column/coupling.rs``.

    **It is a solver change and not a report**: `hasActiveColumnTearVariables` counts the flag, so
    setting it is what puts the column on `solveWithColumnTearVariables`. Each pass re-solves the
    column and rewrites **one end** so the difference between the two equals the designer's summed
    drop - `bottom = top + drop` where the top is positive, `top = max(1e-6 bara, bottom - drop)`
    otherwise - and converges on the relative change it made.
    """
    from azoth.process.reference.designer import DesignerGeometry, designer_report

    limit = max(MAX_COLUMN_TEAR_ITERATIONS, MAX_PUMPAROUND_ITERATIONS, 1)
    top = float(top_pressure.to("Pa").magnitude)
    bottom = float(bottom_pressure.to("Pa").magnitude)
    residual = float("inf")
    states: _States | None = None
    warnings: list[Warning] = []
    for _ in range(limit):
        states, warnings = solve(quantity(top, "Pa"), quantity(bottom, "Pa"))
        report = designer_report(
            states, components, DesignerGeometry(internals_type=internals_type)
        )
        drop = max(0.0, float(report.total_pressure_drop_pa))
        if math.isfinite(top) and top > 0.0:
            previous = bottom
            bottom = top + drop
            if not (math.isfinite(previous) and previous > 0.0):
                residual = 1.0
            else:
                residual = abs(bottom - previous) / max(CHANGED_RESIDUAL, abs(previous))
        elif math.isfinite(bottom) and bottom > 0.0:
            previous_top = top
            top = max(MIN_TOP_PRESSURE_PA, bottom - drop)
            if not (math.isfinite(previous_top) and previous_top > 0.0):
                residual = 1.0
            else:
                residual = abs(top - previous_top) / max(CHANGED_RESIDUAL, abs(previous_top))
        else:
            residual = 0.0

        if residual <= COLUMN_TEAR_TOLERANCE:
            # **A last solve where the update moved something**, which is the class's own
            # re-solve inside its convergence branch.
            if residual > CHANGED_RESIDUAL:
                return solve(quantity(top, "Pa"), quantity(bottom, "Pa"))
            return states, warnings
        if residual <= CHANGED_RESIDUAL:
            return states, warnings

    return solve(quantity(top, "Pa"), quantity(bottom, "Pa"))


def _record(
    states: _States,
    components: list[str],
    internal_diameter: Q,
    max_allowable_fs_factor: float,
    internals: DesignerReport,
    warnings: list[Warning],
) -> DistillationColumnResult:
    """The base column's record, from the stages that solved it and its two reports.

    **One construction with one caller.** `distillation_column` is this and
    `_distillation_column_states` together.

    **The two capacity inputs and the internals report are the caller's**, because none is read
    by the solve: the diameter is the column's internal one and the geometry reaches nothing on
    the run path. **`process.packed_column` does not come through here** - its class builds the
    designer with `internalsType = "packed"`, so `calculateTrayed` never runs for it and there is
    no trayed report to carry.
    """
    gas_out = Stream.from_pt(
        components,
        list(states.distillate_z),
        states.distillate_n,
        from_si(states.distillate_p, "Pa"),
        from_si(states.distillate_t, "K"),
    )
    limits = fs_limits(gas_out, internal_diameter, max_allowable_fs_factor)
    return DistillationColumnResult(
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
        required_diameter=from_si(internals.required_diameter_m, "m"),
        controlling_tray_index=internals.controlling_tray_index,
        internals_design_ok=internals.design_ok,
        max_percent_flood=internals.max_percent_flood,
        min_percent_flood=internals.min_percent_flood,
        average_tray_efficiency=internals.average_tray_efficiency,
        total_pressure_drop=from_si(internals.total_pressure_drop_pa, "Pa"),
        total_pressure_drop_mbar=internals.total_pressure_drop_pa / 100.0,
        tray_percent_flood=tuple(tray.percent_flood for tray in internals.trays),
        tray_pressure_drop=tuple(
            from_si(tray.total_pressure_drop_pa, "Pa") for tray in internals.trays
        ),
        tray_efficiency=tuple(tray.tray_efficiency for tray in internals.trays),
        warnings=tuple(warnings),
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
) -> DistillationColumnResult:
    """Solve a distillation column by sequential substitution.

    The reference implementation of ``process.distillation_column``, and the record its spec
    declares. **The arithmetic and the argument list are `_distillation_column_states`'s**; this
    adds the mapping from the solved stages onto the record, which is `_record`'s.

    Args:
        components: the feed's substances, by name.
        feed_n: the feed's molar flow.
        feed_z: the feed composition.
        feed_p: the feed's pressure, which is its own.
        feed_t: the feed's temperature.
        number_of_stages: the stage count, the reboiler at stage 0 where there is one.
        feed_stage: the stage the feed enters, 0-based over the trays including the ends.
        has_reboiler: whether stage 0 is a reboiler.
        has_condenser: whether the top stage is a condenser.
        top_pressure: the pressure at the top stage.
        bottom_pressure: the pressure at stage 0.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        reboiler_temperature: the reboiler's temperature, which pins the bottom tray.
        condenser_temperature: the condenser's temperature, which pins the top tray.
        murphree_efficiency: the column-wide Murphree tray efficiency.
        tray_murphree_efficiency: one override per stage, ``NaN`` where a stage falls through.
        solver_type: **only ``direct_substitution`` and ``naphtali_sandholm`` are ported**.
        top_specification_type: the top product's degree of freedom.
        top_specification_target: its target value.
        top_specification_component: the component a purity or a recovery constrains.
        bottom_specification_type: the bottom product's degree of freedom, as the top's.
        bottom_specification_target: its target, in the unit its type implies.
        bottom_specification_component: the component a purity or a recovery constrains.
        reactive: whether the middle trays flash reactively.
        reactive_start_tray: the first reactive middle tray.
        reactive_end_tray: the last reactive middle tray, inclusive.
        gas_side_draw_fractions: the vapour each tray withdraws.
        liquid_side_draw_fractions: the liquid each tray withdraws as a side draw.
        pumparound_fractions: the liquid each tray withdraws as a pumparound.
        side_draw_flow_tray: the tray whose draw flow is specified.
        side_draw_flow_phase: which of that tray's phases the specification controls.
        side_draw_flow_target: the mass flow the draw must deliver.
        side_draw_flow_tolerance: the relative residual that search stops at.
        side_draw_flow_max_iterations: the candidate cap for that search.
        pumparound_return_tray: the tray a pumparound's liquid comes back to.
        pumparound_draw_tray: the tray it leaves.
        pumparound_draw_fraction: the fraction of that tray's liquid withdrawn.
        pumparound_temperature_drop: the cooler's drop on the way back.
        pumparound_tolerance: the relative residual the return's loop stops at.
        pumparound_max_iterations: the cap on that loop.
        column_diameter: the column's internal diameter, which the capacity limits divide the
            gas outlet's volumetric flow by. **The solve is indifferent to it.**
        max_allowable_fs_factor: the ``Fs`` limit the family is checked against.
        internals_type: which tray the internals tree is designed as.
        tray_spacing: the tray spacing the internals are designed at.
        weir_height: the weir height.
        hole_diameter: the sieve hole diameter, **in millimetres**.
        hole_area_fraction: the hole area over the active area.
        downcommer_area_fraction: the downcomer area over the total.
        design_flood_fraction: the fraction of flood the internals are sized to.
        column_diameter_override: a stated diameter for the internals tree, which replaces its
            sizing branch entirely.

    Returns:
        The tray profile, both products, both duties, the three residuals, the Fs family and the
        internals tree.

    Raises:
        InvalidInputError: for a stage or a feed stage outside the column, for any declared
            parameter whose arithmetic is not ported, and for a ratio that is not one.
        SolverNotConvergedError: when the solve misses its gate.
    """

    def solve(top: Q, bottom: Q) -> tuple[_States, list[Warning]]:
        """One pass: the column's own solve at the two pressures the loop currently holds."""
        return _distillation_column_states(
            components,
            feed_n,
            feed_z,
            feed_p,
            feed_t,
            number_of_stages,
            feed_stage,
            has_reboiler,
            has_condenser,
            top,
            bottom,
            temperature_tolerance,
            max_iterations,
            reboiler_temperature,
            condenser_temperature,
            murphree_efficiency,
            tray_murphree_efficiency,
            solver_type,
            top_specification_type,
            top_specification_target,
            top_specification_component,
            bottom_specification_type,
            bottom_specification_target,
            bottom_specification_component,
            reactive,
            reactive_start_tray,
            reactive_end_tray,
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

    # **The coupling is the outer loop and the solve above is its inner one**, which is the
    # class's own split: the flag decides which entry runs rather than adding a step inside
    # either.
    if hydraulic_pressure_drop_coupling:
        states, warnings = _coupled(
            solve,
            components,
            top_pressure,
            bottom_pressure,
            hydraulic_pressure_drop_internals_type
            if hydraulic_pressure_drop_internals_type is not None
            else DEFAULT_INTERNALS_TYPE,
        )
    else:
        states, warnings = solve(top_pressure, bottom_pressure)
    # **The internals tree, read after the solve and from the trays it left.** Its geometry
    # reaches nothing on the run path, which is the class's own split.
    internals = designer_report(
        states,
        components,
        DesignerGeometry(
            internals_type=(
                DEFAULT_INTERNALS_TYPE if internals_type is None else str(internals_type)
            ),
            tray_spacing_m=(
                DEFAULT_TRAY_SPACING_M
                if tray_spacing is None
                else float(tray_spacing.to("m").magnitude)
            ),
            weir_height_m=(
                DEFAULT_WEIR_HEIGHT_M
                if weir_height is None
                else float(weir_height.to("m").magnitude)
            ),
            hole_diameter_m=(
                DEFAULT_HOLE_DIAMETER_MM / 1000.0
                if hole_diameter is None
                else float(hole_diameter.to("mm").magnitude) / 1000.0
            ),
            hole_area_fraction=(
                DEFAULT_HOLE_AREA_FRACTION
                if hole_area_fraction is None
                else float(hole_area_fraction)
            ),
            downcommer_area_fraction=(
                DEFAULT_DOWNCOMMER_AREA_FRACTION
                if downcommer_area_fraction is None
                else float(downcommer_area_fraction)
            ),
            design_flood_fraction=(
                DEFAULT_DESIGNER_FLOOD_FRACTION
                if design_flood_fraction is None
                else float(design_flood_fraction)
            ),
            column_diameter_override_m=(
                UNSIZED_COLUMN_DIAMETER_M
                if column_diameter_override is None
                else float(column_diameter_override.to("m").magnitude)
            ),
        ),
    )

    return _record(
        states,
        components,
        column_diameter
        if column_diameter is not None
        else quantity(DEFAULT_INTERNAL_DIAMETER_M, "m"),
        max_allowable_fs_factor
        if max_allowable_fs_factor is not None
        else DEFAULT_MAX_ALLOWABLE_FS_FACTOR,
        internals,
        warnings,
    )


def reactive_section(
    reactive: bool | None, start: int | None, end: int | None
) -> tuple[int, int] | None:
    """`DistillationColumn.setReactive`'s two forms, resolved from the three declared inputs.

    The class clears both bounds for `setReactive(true)` and sets both for
    `setReactive(true, start, end)`; a single bound is a declaration the class cannot make, so it
    is refused rather than guessed - and so is a section stated without `reactive = true`. **One
    resolution for the four column ids**, because the absorber, the stripper and the packed
    column all inherit this setter and none of them overrides it.
    """
    if (start is None) != (end is None):
        raise InvalidInputError(
            "reactive_start_tray",
            "a reactive section is stated by both bounds or by neither: `setReactive(true)` "
            "covers every middle tray and `setReactive(true, start, end)` a run of them, and the "
            "class has no form that states one end alone",
        )
    if reactive:
        return (-1, -1) if start is None or end is None else (int(start), int(end))
    if start is not None:
        raise InvalidInputError(
            "reactive",
            "a reactive section was stated without `reactive = true`, which is a declaration "
            "that says nothing",
        )
    return None


def _states(
    components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_t: float,
    feed_p: float,
    number_of_stages: int,
    feed_stage: int,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: float,
    bottom_pressure: float,
    reboiler_temperature: float | None,
    condenser_temperature: float | None,
    temperature_tolerance: float,
    max_iterations: int,
    top_specification: Specification | None = None,
    bottom_specification: Specification | None = None,
    solver_type: str | None = None,
    top_feed: StreamRecord | None = None,
    tray_temperatures: tuple[float, ...] | None = None,
    reactive: tuple[int, int] | None = None,
    draws: Draws = None,
    murphree_efficiency: _Murphree | None = None,
    absorber_murphree: _AbsorberMurphree | None = None,
    pumparound_inlets: PumparoundInlets | None = None,
) -> _States:
    """The whole solve, in SI magnitudes: the one arithmetic the kernel and a dump share.

    ``top_feed`` is the class's **second inlet**, which enters the top stage:
    `AbsorptionColumn.addSolventInStream` and `StrippingColumn.addRichLiquidStream` are both
    `addFeedStream(stream, getNumberOfTrays() - 1)`. ``tray_temperatures`` is one outlet pin
    per tray - `SimpleTray.setOutletTemperature` - and a column that pins every tray stops
    after its first sweep, because the base's gate is the tray-temperature change it made zero.
    """
    tray_count = number_of_stages + int(has_reboiler) + int(has_condenser)
    # **What the stages had to fall back on**, one entry per kind: a reactive tray's own
    # fallback, deduplicated the way the Rust kernel deduplicates it.
    collected_warnings: list[Warning] = []
    # What each tray withdrew, one vector per kind, in `DRAW_FIELDS`' order.
    drawn: list[list[StreamRecord | None]] = [
        [None] * tray_count,
        [None] * tray_count,
        [None] * tray_count,
    ]
    if number_of_stages == 0:
        raise InvalidInputError(
            "number_of_stages", "a column with no stages between its ends is not a column"
        )
    if feed_stage >= tray_count:
        raise InvalidInputError(
            "feed_stage",
            f"stage {feed_stage} is outside a column of {tray_count} trays, whose stages run "
            f"0 to {tray_count - 1}",
        )
    if not temperature_tolerance > 0.0:
        raise InvalidInputError(
            "temperature_tolerance",
            f"a convergence tolerance is a positive number, and {temperature_tolerance} is not one",
        )

    # A one-tray column is its own case: the interpolation divides by `tray_count - 1`, so a
    # single tray takes the bottom's own pressure rather than a `0/0` - the mirror of the kernel's
    # own guard.
    span = tray_count - 1
    pressures = (
        [bottom_pressure]
        if span == 0
        else [
            bottom_pressure + (top_pressure - bottom_pressure) * i / span for i in range(tray_count)
        ]
    )

    def pin(i: int) -> float | None:
        if tray_temperatures is not None:
            stated = tray_temperatures[i]
            if stated == stated:  # not NaN
                return stated
        if i == 0 and has_reboiler:
            return reboiler_temperature
        if i == tray_count - 1 and has_condenser:
            return condenser_temperature
        return None

    def end_mode(i: int, specification: Specification | None) -> tuple[str, float]:
        """How a tray is run: its temperature, its own reflux flash, or a duty.

        **A pin wins over a directly-applied duty**, which is NeqSim's own order of tests: the
        tray `run` checks its stated outlet temperature first and takes a `TPflash` there, so a
        heat input a `DUTY` specification wrote is never read. A reflux ratio does win, because
        the end's own `run` tests `refluxIsSet` before anything else.
        """
        stated = pin(i)
        if specification is not None:
            if specification.kind == "reflux_ratio":
                return ("reflux", specification.target)
            if specification.kind == "duty" and stated is None:
                return ("duty", specification.target)
        return ("temperature", stated if stated is not None else float("nan"))

    def specification_at(i: int) -> Specification | None:
        if i == tray_count - 1:
            return top_specification
        if i == 0:
            return bottom_specification
        return None

    modes = [end_mode(i, specification_at(i)) for i in range(tray_count)]

    if solver_type == "naphtali_sandholm":
        # **The mesh solve, warmed from this one.** The class's own `initializeTrayStateFromColumn`
        # maps a column's converged trays onto the MESH variables, and a port that owns both
        # solves can always take that path - so the substitution core runs first and its answer
        # is the seed.
        if top_specification is not None or bottom_specification is not None:
            raise InvalidInputError(
                "solver_type",
                "a product specification re-runs the whole column inside "
                "`solveWithSpecificationTargets`, which is a second integration this port does "
                "not carry: state the ends' temperatures instead",
            )
        from azoth.process.reference._column_mesh import mesh_states

        return mesh_states(
            components,
            feed_n,
            list(feed_z),
            feed_t,
            feed_p,
            feed_stage,
            list(pressures),
            [pinned if (pinned := pin(i)) is not None else float("nan") for i in range(tray_count)],
            seed=_states(
                components,
                feed_n,
                feed_z,
                feed_t,
                feed_p,
                number_of_stages,
                feed_stage,
                has_reboiler,
                has_condenser,
                top_pressure,
                bottom_pressure,
                reboiler_temperature,
                condenser_temperature,
                temperature_tolerance,
                max_iterations,
            ),
        )

    def end_temperature(i: int) -> float | None:
        """The temperature the tray's own mode states, where it states a finite one.

        The seed reads the *tray's* mode rather than the caller's optional pin: an end run on
        a duty or on a reflux ratio states no temperature, and its step then falls back to the
        class's own, which is one kelvin below the feed's.
        """
        kind, value = modes[i]
        if kind != "temperature" or value != value or value in (float("inf"), float("-inf")):
            return None
        return value

    gas: list[StreamRecord | None] = [None] * tray_count
    liquid: list[StreamRecord | None] = [None] * tray_count
    # **A returning pumparound is an inlet like any other**, which is what makes the recycle a
    # recycle: the tray the liquid comes back to mixes it with the rest of its inlets.
    returns: list[StreamRecord | None] = (
        list(pumparound_inlets) if pumparound_inlets is not None else [None] * tray_count
    )
    # **What each tray's own flash found**, before either the draws or the Murphree correction
    # touched it: the vapour is the correction's `y_eq` and `y_in`, and the pair is the phase
    # count its two guards read. The class's correction compares against
    # `trays.get(i - 1).getThermoSystem().getPhase(0)` - the *system*, which
    # `setCachedGasOutStream` never writes - so the stage above receives the corrected vapour
    # while the correction sees the uncorrected one. Both are needed and they differ.
    equilibrium: list[tuple[StreamRecord | None, StreamRecord | None]] = [(None, None)] * tray_count

    def present(stream: StreamRecord | None, what: str) -> StreamRecord:
        """A stream the seed needs, refused rather than carried as an absent phase.

        Both reads mirror a kernel `expect`: by the time `init` seeds a tray, the tray below it
        has run, so an absent phase there is a broken link and not a state.
        """
        if stream is None:
            raise InvalidInputError("column", what)
        return stream

    def is_reactive(i: int) -> bool:
        """**The section is stated over middle trays and the ends are never in it**:
        `replaceMiddleTrays` walks from the reboiler's edge to the condenser's."""
        if reactive is None:
            return False
        if (i == 0 and has_reboiler) or (i + 1 == tray_count and has_condenser):
            return False
        middle = i - int(has_reboiler)
        start, end = reactive
        if start < 0 and end < 0:
            return True
        return start <= middle <= end

    def draws_at(i: int) -> tuple[float, float, float]:
        """The three fractions a tray states, a zero where a vector is absent."""
        if draws is None:
            return (0.0, 0.0, 0.0)
        return tuple(  # type: ignore[return-value]
            0.0 if vector is None or i >= len(vector) else vector[i] for vector in draws
        )

    def run_with(i: int, inlets: list[StreamRecord]) -> None:
        kind, value = modes[i]
        if kind == "temperature":
            # A middle tray carries a NaN here and no pin, so its flash is at its own enthalpy.
            pin = value if value == value else None
            out = (
                reactive_stage(inlets, pressures[i], pin, 0.0)
                if is_reactive(i)
                else stage(inlets, pressures[i], pin, 0.0)
            )
        elif kind == "duty":
            out = (
                reactive_stage(inlets, pressures[i], None, value)
                if is_reactive(i)
                else stage(inlets, pressures[i], None, value)
            )
        else:
            out = _stage.reflux_end(
                inlets, pressures[i], value, "vapour" if i == tray_count - 1 else "liquid"
            )
        stated = draws_at(i)
        if stated != (0.0, 0.0, 0.0):
            # **The draws come off the phase the tray has just formed**, exactly as the Rust
            # stage's own split does.
            (
                out["gas"],
                out["liquid"],
                out["gas_side_draw"],
                out["liquid_side_draw"],
                out["pumparound"],
            ) = split_draws(out["gas"], out["liquid"], stated)
        equilibrium[i] = (out["gas"], out["liquid"])
        if absorber_murphree is not None:
            # **`AbsorptionColumn`'s override**, whose ``below`` is the tray's *gas inlet*:
            # ``gasInStream`` at stage 0 and the stage below's flash above it.
            inlet = (
                _feed(components, feed_n, feed_z, feed_t, feed_p)
                if i == 0
                else equilibrium[i - 1][0]
            )
            corrected = _correct_absorber(out["gas"], out["liquid"], inlet, absorber_murphree, i)
            if corrected is not None:
                out["gas"], out["liquid"] = corrected
        effective = None if murphree_efficiency is None else murphree_efficiency.resolve(i)
        if (
            absorber_murphree is None
            and effective is not None
            and _corrects(i, tray_count, has_condenser, effective)
        ):
            below_gas, below_liquid = equilibrium[i - 1]
            if out["liquid"] is not None and below_liquid is not None:
                out["gas"] = _correct_vapour(
                    equilibrium=out["gas"],
                    inlet=below_gas,
                    efficiency=effective,
                )
        gas[i], liquid[i] = out["gas"], out["liquid"]
        drawn[0][i] = out.get("gas_side_draw")
        drawn[1][i] = out.get("liquid_side_draw")
        drawn[2][i] = out.get("pumparound")
        for warning in out["warnings"]:
            if all(held.code != warning.code for held in collected_warnings):
                collected_warnings.append(warning)

    def at_temperature(stream: StreamRecord, temperature: float) -> StreamRecord:
        if temperature != temperature:  # NaN
            return stream
        return _restate(components, stream, temperature)

    # The three validators, and the ends' refusal: a fraction on a reboiler or a condenser names
    # a stage this port does not have.
    if draws is not None:
        for vector, name in zip(
            draws,
            ("gas_side_draw_fractions", "liquid_side_draw_fractions", "pumparound_fractions"),
            strict=True,
        ):
            if vector is None:
                continue
            if len(vector) != tray_count:
                raise InvalidInputError(
                    name,
                    f"{len(vector)} fraction(s) against {tray_count} tray(s): the vector is one "
                    f"entry per tray, a zero where a tray draws nothing",
                )
            for index, exists, where in (
                (0, has_reboiler, "the reboiler"),
                (tray_count - 1, has_condenser, "the condenser"),
            ):
                if exists and vector[index] != 0.0:
                    raise InvalidInputError(
                        name,
                        f"a draw is stated on {where}: the ends here are `column::reboiler` and "
                        f"`column::condenser` rather than stages",
                    )
        for i in range(tray_count):
            validate_draws(draws_at(i))

    def inlets_of(i: int, seed: list[float] | None) -> list[StreamRecord]:
        inlets: list[StreamRecord] = []
        at = (lambda s: at_temperature(s, seed[i])) if seed is not None else (lambda s: s)
        if i > 0 and gas[i - 1] is not None:
            inlets.append(at(gas[i - 1]))
        if i + 1 < tray_count and liquid[i + 1] is not None:
            inlets.append(at(liquid[i + 1]))
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        if top_feed is not None and i == tray_count - 1:
            inlets.append(top_feed)
        returned = returns[i]
        if returned is not None:
            inlets.append(returned)
        return inlets

    def temperature(i: int) -> float:
        stream = gas[i] if gas[i] is not None else liquid[i]
        return float("nan") if stream is None else stream["t"]

    # ---- `init` ----
    run_with(feed_stage, [_feed(components, feed_n, feed_z, feed_t, feed_p)])
    if has_reboiler:
        run_with(
            0,
            [present(liquid[feed_stage], "the feed stage has a liquid after its flash")],
        )

    feed_temperature = temperature(feed_stage)
    # **The top tray's temperature is the *feed's* and not the condenser's pin**, which is the
    # class's own branch and it is subtle. `init` computes this **before** the upward link, so
    # the top tray has no internal inlet yet and `getNumberOfInputStreams() > 0` is false for
    # every column whose feed is not at the top - which is what makes `feedTrayTemperature -
    # 1.0` the branch a distillation column takes. The other branch is for a tray that carries
    # an *external* feed of its own - an absorber's solvent at the top stage - and there the
    # temperature is that feed's, because the tray has not run.
    # Either external feed the top tray carries, which is the class's
    # `getNumberOfInputStreams() > 0` in this port's terms.
    at_top: StreamRecord | None = None
    if feed_stage == tray_count - 1:
        at_top = _feed(components, feed_n, feed_z, feed_t, feed_p)
    if top_feed is not None:
        at_top = top_feed
    cond_temperature = at_top["t"] if at_top is not None else feed_temperature - 1.0
    reb_temperature = liquid[0]["t"] if liquid[0] is not None else feed_temperature
    temperatures = [float("nan")] * tray_count
    temperatures[feed_stage] = feed_temperature
    # **The mirror of the step below**: a column whose feed tray is its top has no trays above
    # it, so the step is never taken - and on a one-tray column both steps are.
    above = tray_count - feed_stage - 1
    delta_up = (feed_temperature - cond_temperature) / above if above else 0.0
    # **A feed at stage 0 has no trays below it**, so the downward step is never taken and the
    # division would be 0/0. Rust's `0.0/0.0` is a NaN that the empty loop never reads; Python
    # raises, so the step is computed only where there is a tray to move.
    delta_down = (reb_temperature - feed_temperature) / feed_stage if feed_stage else 0.0
    delta = 0.0
    for i in range(feed_stage + 1, tray_count):
        delta += delta_up
        temperatures[i] = feed_temperature - delta
    delta = 0.0
    for i in range(feed_stage - 1, -1, -1):
        delta += delta_down
        temperatures[i] = feed_temperature + delta
    for i in range(tray_count):
        pinned = pin(i)
        if pinned is not None:
            temperatures[i] = pinned

    for i in range(1, tray_count):
        inlets = [at_temperature(present(gas[i - 1], "the tray below has run"), temperatures[i])]
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        if top_feed is not None and i == tray_count - 1:
            inlets.append(top_feed)
        run_with(i, inlets)
    for i in range(tray_count - 2, 0, -1):
        inlets = [
            at_temperature(present(gas[i - 1], "the tray below has run"), temperatures[i]),
            at_temperature(present(liquid[i + 1], "the tray above has run"), temperatures[i]),
        ]
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        if top_feed is not None and i == tray_count - 1:
            inlets.append(top_feed)
        run_with(i, inlets)
    if has_reboiler:
        run_with(
            0,
            [at_temperature(present(liquid[1], "the first tray has run"), temperatures[0])],
        )

    # ---- the sweeps, once, or under the outer specification loop ----
    def sweep(temperatures: list[float]) -> tuple[int, float]:
        """`solveSequential`: the sweeps, once, to the temperature gate."""
        relaxation = 1.0
        previous_combined = float("inf")
        temperature_residual = float("inf")
        took = 0
        for iteration in range(1, max_iterations + 1):
            took = iteration
            old = [temperature(i) for i in range(tray_count)]

            for i in range(feed_stage, 1, -1):
                run_with(i - 1, inlets_of(i - 1, None))
            if has_reboiler:
                run_with(0, inlets_of(0, None))
            for i in range(1, tray_count):
                run_with(i, inlets_of(i, None))
            for i in range(tray_count - 2, feed_stage - 1, -1):
                run_with(i, inlets_of(i, None))

            effective = min(max(relaxation, MIN_TEMPERATURE_RELAXATION), 1.0)
            total = 0.0
            for i in range(tray_count):
                updated = temperature(i)
                if updated != updated or updated in (float("inf"), float("-inf")):
                    updated = old[i]
                total += abs(updated - old[i])
                temperatures[i] = old[i] + effective * (updated - old[i])
            temperature_residual = total / tray_count

            combined = temperature_residual / temperature_tolerance
            if combined > previous_combined * 1.05:
                relaxation = max(MIN_SEQUENTIAL_RELAXATION, relaxation * RELAXATION_DECREASE_FACTOR)
            elif combined < previous_combined * 0.98:
                relaxation = min(MAX_ADAPTIVE_RELAXATION, relaxation * RELAXATION_INCREASE_FACTOR)
            previous_combined = combined

            if temperature_residual <= temperature_tolerance:
                break

        return took, temperature_residual

    def product_of(level: int) -> StreamRecord:
        """The product an end publishes: the top's vapour, the bottom's liquid.

        **The top's vapour falls back to its liquid**, which is how the kernel's own
        `end_product` reads a total condenser's distillate.
        """
        vapour = gas[tray_count - 1]
        if level == 1:
            stream = vapour if vapour is not None else liquid[tray_count - 1]
        else:
            stream = liquid[0]
        if stream is None:
            raise InvalidInputError(
                "column", "an end has no product, which is a column whose ends did not solve"
            )
        return stream

    def adjustable() -> bool:
        """`hasAdjustableSpecifications`: whether any end's specification needs the outer loop."""
        return any(
            spec is not None and spec.kind not in DIRECT_KINDS
            for spec in (top_specification, bottom_specification)
        )

    if adjustable():
        # `solveWithSpecificationTargets`: the outer secant on the ends' temperatures.
        #
        # **There is no homotopy here, and that is a measurement.**
        # `getEffectiveSpecificationHomotopySteps` returns three stages only when
        # `solverType == AUTO`, and the field it otherwise reads is one - so the staged
        # continuation the class carries belongs to the solver this port refuses.
        adjust_top = top_specification is not None and top_specification.kind not in DIRECT_KINDS
        adjust_bottom = (
            bottom_specification is not None and bottom_specification.kind not in DIRECT_KINDS
        )
        top_start = condenser_temperature if condenser_temperature is not None else feed_t - 20.0
        bottom_start = reboiler_temperature if reboiler_temperature is not None else feed_t + 20.0
        top_pair = (top_start, top_start - SECANT_FIRST_OFFSET)
        bottom_pair = (bottom_start, bottom_start + SECANT_FIRST_OFFSET)
        top_errors = (float("nan"), float("nan"))
        bottom_errors = (float("nan"), float("nan"))
        feed_record = _feed(components, feed_n, feed_z, feed_t, feed_p)

        for outer in range(SPECIFICATION_MAX_ITERATIONS):
            if adjust_top:
                modes[-1] = ("temperature", top_pair[0] if outer == 0 else top_pair[1])
            if adjust_bottom:
                modes[0] = ("temperature", bottom_pair[0] if outer == 0 else bottom_pair[1])

            iterations, temperature_residual = sweep(temperatures)

            top_error = 0.0
            bottom_error = 0.0
            if adjust_top and top_specification is not None:
                top_error = (
                    top_specification.value(product_of(1), feed_record) - top_specification.target
                )
            if adjust_bottom and bottom_specification is not None:
                bottom_error = (
                    bottom_specification.value(product_of(0), feed_record)
                    - bottom_specification.target
                )
            if (
                abs(top_error) < SPECIFICATION_TOLERANCE
                and abs(bottom_error) < SPECIFICATION_TOLERANCE
            ):
                break
            if adjust_top:
                if outer == 0:
                    top_errors = (top_error, float("nan"))
                else:
                    top_pair = (top_pair[1], secant_step(top_pair, (top_errors[0], top_error)))
                    top_errors = (top_error, float("nan"))
            if adjust_bottom:
                if outer == 0:
                    bottom_errors = (bottom_error, float("nan"))
                else:
                    bottom_pair = (
                        bottom_pair[1],
                        secant_step(bottom_pair, (bottom_errors[0], bottom_error)),
                    )
                    bottom_errors = (bottom_error, float("nan"))
        else:
            from azoth.core.errors import SolverNotConvergedError

            raise SolverNotConvergedError(iterations, temperature_residual, SPECIFICATION_TOLERANCE)
    else:
        iterations, temperature_residual = sweep(temperatures)

    def outlet_enthalpy(i: int) -> float:
        total = 0.0
        for stream in (gas[i], liquid[i]):
            if stream is not None:
                total += stream["n"] * stream["h"]
        return total

    def inlet_enthalpy(i: int) -> float:
        return sum(s["n"] * s["h"] for s in inlets_of(i, None))

    distillate = product_of(1)
    bottoms = product_of(0)
    reboiler_duty = (outlet_enthalpy(0) - inlet_enthalpy(0)) if has_reboiler else 0.0
    condenser_duty = (
        outlet_enthalpy(tray_count - 1) - inlet_enthalpy(tray_count - 1) if has_condenser else 0.0
    )

    feeds = [_feed(components, feed_n, feed_z, feed_t, feed_p)]
    if top_feed is not None:
        feeds.append(top_feed)
    feed_enthalpy = sum(float(feed["n"]) * float(feed["h"]) for feed in feeds)
    # **A pumparound's return is an inlet and its draw is an outlet**, so the two appear on
    # opposite sides of both closures and cancel - the column's inlets are its feeds *and*
    # whatever comes back to it. The cooler's duty is the enthalpy between them and is not an
    # imbalance.
    returned_enthalpy = sum(
        float(returned["n"]) * float(returned["h"]) for returned in returns if returned is not None
    )
    # **The draws are a third route out of the column**, exactly as they are in the kernel's own
    # closure: a feed leaves as the two products *and* whatever the trays withdrew, so a residual
    # that ignored them would report an imbalance that is not there. Measured on the binary column
    # drawing a quarter of tray 3's vapour, the ignoring form reports a mass residual of 0.248.
    held = [draw for kind in drawn for draw in kind if draw is not None]
    drawn_enthalpy = sum(draw["n"] * draw["h"] for draw in held)
    products_enthalpy = (
        distillate["n"] * distillate["h"] + bottoms["n"] * bottoms["h"] + drawn_enthalpy
    )
    energy_residual = (
        abs(feed_enthalpy + returned_enthalpy + reboiler_duty + condenser_duty - products_enthalpy)
        / abs(feed_enthalpy)
        if abs(feed_enthalpy) > 0.0
        else 0.0
    )

    mass_residual = 0.0
    for c in range(len(feed_z)):
        supplied = sum(float(feed["n"]) * float(feed["z"][c]) for feed in feeds)
        supplied += sum(
            float(returned["n"]) * float(returned["z"][c])
            for returned in returns
            if returned is not None
        )
        withdrawn = sum(draw["n"] * draw["z"][c] for draw in held)
        delivered = (
            distillate["n"] * distillate["z"][c] + bottoms["n"] * bottoms["z"][c] + withdrawn
        )
        if abs(supplied) > 1.0e-12:
            mass_residual = max(mass_residual, abs(supplied - delivered) / abs(supplied))

    # **The two closures gate a conserving column and are only reported on a corrected one.**
    # `getMassBalanceError` - the quantity the class's own gate scales against
    # `baseMassTolerance` - compares a stage's inlet streams with the stage's *system*, which
    # `mixStream` has already made equal, so it is identically zero and never binds. A Murphree
    # efficiency is what makes a column's products stop carrying the feed, and that is the
    # class's own arithmetic rather than a failure of it: what a stage hands up carries the
    # flash's *moles* at the blend's composition, so the interior imbalances do not cancel.
    # `_reconcile_products` below is what the class publishes in that state.
    conserving = murphree_efficiency is None
    if (
        temperature_residual > temperature_tolerance
        or (conserving and mass_residual > MASS_BALANCE_TOLERANCE)
        or (conserving and energy_residual > ENTHALPY_BALANCE_TOLERANCE)
    ):
        from azoth.core.errors import SolverNotConvergedError

        raise SolverNotConvergedError(
            iterations,
            max(temperature_residual, mass_residual, energy_residual),
            min(temperature_tolerance, MASS_BALANCE_TOLERANCE, ENTHALPY_BALANCE_TOLERANCE),
        )

    # **`updateProductsFromExternalComponentBalance`**: the two published products rescaled,
    # per component, so that they carry the feed exactly. It is on the class's run path and a
    # column at a Murphree efficiency cannot be published without it - the correction moves
    # components between a stage's vapour and its liquid without moving its moles, so the
    # column's overall closure cannot close on its own.
    supplied_c = [
        # **The returns join the supply here**, because this function's withdrawal is the draw:
        # the recycled liquid leaves as a draw and comes back, so a reconciliation that
        # subtracted the draw without adding the return would shrink the bottoms by exactly the
        # draw - which is what it did, by `0.516` against the class's unchanged `3.699`.
        sum(float(feed["n"]) * float(feed["z"][c]) for feed in feeds)
        + sum(
            float(returned["n"]) * float(returned["z"][c])
            for returned in returns
            if returned is not None
        )
        for c in range(len(feed_z))
    ]
    withdrawn_c = [sum(draw["n"] * draw["z"][c] for draw in held) for c in range(len(feed_z))]
    top_moles = [0.0] * len(feed_z)
    bottom_moles = [0.0] * len(feed_z)
    for c in range(len(feed_z)):
        here = distillate["n"] * distillate["z"][c]
        below = bottoms["n"] * bottoms["z"][c]
        carried = here + below
        if carried <= 1.0e-20:
            continue
        scale = max(0.0, supplied_c[c] - withdrawn_c[c]) / carried
        top_moles[c] = here * scale
        bottom_moles[c] = below * scale
    distillate = _rebuild(distillate, top_moles)
    bottoms = _rebuild(bottoms, bottom_moles)

    # **`finalizeSolve` writes them back onto the end trays**, and the record has to show it:
    # the capture's `tray0_liquid_n` and `tray5_gas_n` are the *products* and not the ends'
    # flashes, because `finalizeSolve` runs after the solve and before a probe can read a tray.
    gas[tray_count - 1] = distillate
    liquid[0] = bottoms
    if murphree_efficiency is not None:
        reboiler_duty = (outlet_enthalpy(0) - inlet_enthalpy(0)) if has_reboiler else 0.0
        condenser_duty = (
            outlet_enthalpy(tray_count - 1) - inlet_enthalpy(tray_count - 1)
            if has_condenser
            else 0.0
        )

    return _States(
        warnings=tuple(collected_warnings),
        gas_side_draw_n=tuple(draw["n"] if draw is not None else 0.0 for draw in drawn[0]),
        liquid_side_draw_n=tuple(draw["n"] if draw is not None else 0.0 for draw in drawn[1]),
        pumparound_n=tuple(draw["n"] if draw is not None else 0.0 for draw in drawn[2]),
        tray_temperature=tuple(temperature(i) for i in range(tray_count)),
        tray_pressure=tuple(pressures),
        tray_gas_n=tuple(
            stream["n"] if (stream := gas[i]) is not None else 0.0 for i in range(tray_count)
        ),
        tray_liquid_n=tuple(
            stream["n"] if (stream := liquid[i]) is not None else 0.0 for i in range(tray_count)
        ),
        tray_gas_z=tuple(
            tuple(stream["z"]) if (stream := gas[i]) is not None else () for i in range(tray_count)
        ),
        tray_liquid_z=tuple(
            tuple(stream["z"]) if (stream := liquid[i]) is not None else ()
            for i in range(tray_count)
        ),
        distillate_n=distillate["n"],
        distillate_z=tuple(distillate["z"]),
        distillate_h=distillate["h"],
        distillate_t=float(distillate["t"]),
        distillate_p=float(distillate["p"]),
        bottoms_n=bottoms["n"],
        bottoms_z=tuple(bottoms["z"]),
        bottoms_h=bottoms["h"],
        bottoms_t=float(bottoms["t"]),
        bottoms_p=float(bottoms["p"]),
        condenser_duty=condenser_duty,
        reboiler_duty=reboiler_duty,
        iterations=iterations,
        temperature_residual=temperature_residual,
        mass_residual=mass_residual,
        energy_residual=energy_residual,
    )


class _AbsorberMurphree(NamedTuple):
    """`AbsorptionColumn`'s two efficiency maps over the base's two levels.

    `getComponentMurphreeEfficiency(tray, component)` reads
    ``componentMurphreeEfficiencies[component]``, then the base's per-tray and column-wide pair -
    the four levels, of which the per-tray-per-component map the class's other overload writes
    has no palette spelling and is not carried.
    """

    base: _Murphree
    per_component: tuple[float, ...] | None = None

    def resolve(self, tray: int, component: int) -> float:
        """`getComponentMurphreeEfficiency`: the component's own value, else the base's two."""
        if self.per_component is not None and component < len(self.per_component):
            value = self.per_component[component]
            if not math.isnan(value):
                return value
        return self.base.resolve(tray)

    def checked(self, components: int) -> _AbsorberMurphree:
        """The component vector's length held to the components the column carries."""
        if self.per_component is not None and len(self.per_component) != components:
            raise InvalidInputError(
                "component_murphree_efficiency",
                f"{len(self.per_component)} component efficiency(ies) against {components} "
                f"component(s): the vector is one entry per component in the column's own "
                f"order, and a component that states none is written `NaN` rather than left out",
            )
        return self


MOLE_TOLERANCE = 1.0e-15


def _allocate_vapour_moles(
    mole_fraction: Sequence[float], available_moles: Sequence[float], vapour_moles: float
) -> list[float]:
    """`AbsorptionColumn.allocateVaporMoles`: the flash's vapour moles across the components.

    The loop is the class's, including its two fallbacks: while a component's trial share is
    more than the moles available to it, that component is fixed at its availability and the
    remainder is redistributed. The ``remainingWeight`` branch divides by the *available* moles
    instead, reached only when every unfixed component's corrected fraction has collapsed.
    """
    count = len(mole_fraction)
    allocated = [0.0] * count
    fixed = [False] * count
    remaining = vapour_moles
    for _ in range(count):
        if remaining <= MOLE_TOLERANCE:
            break
        weight = sum(
            fraction for fraction, held in zip(mole_fraction, fixed, strict=True) if not held
        )
        limited = False
        for component in range(count):
            if fixed[component]:
                continue
            trial = (
                remaining * mole_fraction[component] / weight
                if weight > MOLE_TOLERANCE
                else remaining
                * available_moles[component]
                / max(
                    sum(
                        moles
                        for moles, held in zip(available_moles, fixed, strict=True)
                        if not held
                    ),
                    MOLE_TOLERANCE,
                )
            )
            if trial > available_moles[component] + MOLE_TOLERANCE:
                allocated[component] = available_moles[component]
                remaining -= allocated[component]
                fixed[component] = True
                limited = True
        if not limited:
            for component in range(count):
                if not fixed[component]:
                    allocated[component] = (
                        remaining * mole_fraction[component] / weight
                        if weight > MOLE_TOLERANCE
                        else remaining
                        * available_moles[component]
                        / max(
                            sum(
                                moles
                                for moles, held in zip(available_moles, fixed, strict=True)
                                if not held
                            ),
                            MOLE_TOLERANCE,
                        )
                    )
            remaining = 0.0
    return allocated


def _correct_absorber(
    equilibrium_gas: StreamRecord | None,
    equilibrium_liquid: StreamRecord | None,
    inlet: StreamRecord | None,
    efficiency: _AbsorberMurphree,
    tray: int,
) -> tuple[StreamRecord, StreamRecord] | None:
    """`AbsorptionColumn.applyMurphreeCorrection`, which replaces the base column's.

    The base blends one phase; this blends **both** and re-allocates the flash's vapour moles
    across the components, so what a stage hands up carries a different composition at the same
    total moles.
    """
    if equilibrium_gas is None or equilibrium_liquid is None or inlet is None:
        return None
    count = len(equilibrium_gas["z"])
    corrected = [0.0] * count
    total = [0.0] * count
    required = False
    total_fraction = 0.0
    for component in range(count):
        density = efficiency.resolve(tray, component)
        required |= density < 1.0 - 1.0e-10
        equilibrium_fraction = equilibrium_gas["z"][component]
        inlet_fraction = inlet["z"][component]
        blended = max(0.0, inlet_fraction + density * (equilibrium_fraction - inlet_fraction))
        corrected[component] = blended
        total_fraction += blended
        total[component] = max(
            0.0,
            equilibrium_gas["z"][component] * equilibrium_gas["n"]
            + equilibrium_liquid["z"][component] * equilibrium_liquid["n"],
        )
    if not required or total_fraction <= MOLE_TOLERANCE:
        return None
    corrected = [value / total_fraction for value in corrected]

    vapour_moles = equilibrium_gas["n"]
    liquid_moles = equilibrium_liquid["n"]
    allocated = _allocate_vapour_moles(corrected, total, vapour_moles)
    gas_z = tuple(
        moles / vapour_moles if vapour_moles > MOLE_TOLERANCE else 0.0 for moles in allocated
    )
    liquid_z = tuple(
        (max(0.0, available - gas) / liquid_moles) if liquid_moles > MOLE_TOLERANCE else 0.0
        for gas, available in zip(allocated, total, strict=True)
    )
    return (
        _stage.stream_at(
            list(equilibrium_gas["components"]),
            vapour_moles,
            list(gas_z),
            float(equilibrium_gas["t"]),
            float(equilibrium_gas["p"]),
        ),
        _stage.stream_at(
            list(equilibrium_liquid["components"]),
            liquid_moles,
            list(liquid_z),
            float(equilibrium_liquid["t"]),
            float(equilibrium_liquid["p"]),
        ),
    )


class _Murphree(NamedTuple):
    """`DistillationColumn`'s two efficiency fields and the resolution between them.

    The class carries ``murphreeEfficiency``, a column-wide scalar defaulting to ``1.0``, and
    ``perStageMurphreeEfficiency``, a nullable array; ``getEffectiveMurphreeEfficiency(stage)``
    is the whole rule. **Two levels, not four** - the per-component ones are
    ``AbsorptionColumn``'s. ``NaN`` is the fall-through sentinel rather than an absence.
    """

    column_wide: float
    per_stage: tuple[float, ...] | None = None

    def resolve(self, index: int) -> float:
        """`getEffectiveMurphreeEfficiency`: the override, else the column-wide value."""
        if self.per_stage is not None and index < len(self.per_stage):
            value = self.per_stage[index]
            if not math.isnan(value):
                return value
        return self.column_wide

    @staticmethod
    def clamp(efficiency: float) -> float:
        """`clampMurphreeEfficiency`: the class clamps a request into ``[0, 1]``."""
        return max(0.0, min(1.0, efficiency))


def _murphree(
    column_wide: float | None,
    per_stage: Sequence[float] | None,
    tray_count: int,
) -> _Murphree | None:
    """The two fields resolved and clamped, or ``None`` where neither is stated.

    The per-stage **length** is the one thing refused rather than clamped, because a length is
    a statement about the column and not a value - `setMurphreeEfficiencies` refuses it in its
    own words.
    """
    if column_wide is None and per_stage is None:
        return None
    if per_stage is not None and len(per_stage) != tray_count:
        raise InvalidInputError(
            "tray_murphree_efficiency",
            f"{len(per_stage)} override(s) for a column of {tray_count} stage(s): "
            "`DistillationColumn.setMurphreeEfficiencies` refuses an array whose length is not "
            "the stage count, and a stage that states no override is written `NaN` rather than "
            "left out",
        )
    return _Murphree(
        column_wide=_Murphree.clamp(1.0 if column_wide is None else column_wide),
        per_stage=None
        if per_stage is None
        else tuple(_Murphree.clamp(value) for value in per_stage),
    )


def _corrects(index: int, tray_count: int, has_condenser: bool, efficiency: float) -> bool:
    """Whether `DistillationColumn.applyMurphreeCorrection` would correct the stage at ``index``.

    The class's three tests, in its order: an ideal stage, the reboiler at zero, the condenser
    when there is one. ``tray_count`` counts the ends.
    """
    return (
        efficiency < 1.0 - 1.0e-10 and index > 0 and not (has_condenser and index >= tray_count - 1)
    )


def _correct_vapour(
    equilibrium: StreamRecord | None, inlet: StreamRecord | None, efficiency: float
) -> StreamRecord | None:
    """`DistillationColumn.applyMurphreeCorrection`: the vapour leaving one stage, blended.

    ``None`` where the class would return without touching anything. The two guards that are
    *phases* - `getNumberOfPhases() < 2` on the stage and on the one below - are the caller's,
    because an absent phase is the only way this library says a phase is not there.
    """
    if equilibrium is None or inlet is None:
        return equilibrium
    if list(equilibrium["components"]) != list(inlet["components"]):
        raise InvalidInputError(
            "murphree_efficiency",
            "the stage and the stage below it do not carry the same substances in the same "
            "order, so the correction has no basis for pairing one's components with the "
            "other's: `DistillationColumn.applyMurphreeCorrection` indexes both by position",
        )
    actual = [
        max(0.0, y_in + efficiency * (y_eq - y_in))
        for y_eq, y_in in zip(equilibrium["z"], inlet["z"], strict=True)
    ]
    total = sum(actual)
    if total > 1.0e-15:
        actual = [value / total for value in actual]
    # **The total moles are the flash's own and the composition is the blended one.** The class
    # builds the corrected stream from `fluid.phaseToSystem(0)`, replaces every component's
    # fraction and mole count, and re-initialises - so what the stage above *receives* is a
    # stream whose overall composition is the blend, at the stage's own temperature and pressure.
    return _stage.stream_at(
        list(equilibrium["components"]),
        equilibrium["n"],
        actual,
        float(equilibrium["t"]),
        float(equilibrium["p"]),
    )


def _rebuild(source: StreamRecord, moles: list[float]) -> StreamRecord:
    """A stream restated over stated component moles, at its own temperature and pressure."""
    total = sum(moles)
    return _stage.stream_at(
        list(source["components"]),
        total,
        [value / total for value in moles],
        float(source["t"]),
        float(source["p"]),
    )


def secant_step(guess: tuple[float, float], error: tuple[float, float]) -> float:
    """`secantStep`, with its two guards: a step cap and a physically reasonable window."""
    denominator = error[1] - error[0]
    if abs(denominator) < 1.0e-15:
        next_guess = guess[1] + 2.0
    else:
        next_guess = guess[1] - error[1] * (guess[1] - guess[0]) / denominator
    if next_guess - guess[1] > SECANT_MAX_STEP:
        next_guess = guess[1] + SECANT_MAX_STEP
    elif next_guess - guess[1] < -SECANT_MAX_STEP:
        next_guess = guess[1] - SECANT_MAX_STEP
    return min(max(next_guess, SECANT_MIN_TEMPERATURE), SECANT_MAX_TEMPERATURE)


def _feed(components: list[str], n: float, z: list[float], t: float, p: float) -> StreamRecord:
    """The external feed, rebuilt so a seed can never reach it."""
    return _stage.stream_at(components, n, list(z), t, p)


def _restate(components: list[str], stream: StreamRecord, temperature: float) -> StreamRecord:
    """A stream at another temperature, at its own pressure and composition."""
    return _stage.stream_at(components, stream["n"], list(stream["z"]), temperature, stream["p"])


def _refuse_unported(solver_type: str | None, murphree_efficiency: _Murphree | None) -> None:
    """Refuse every value this id's spec declares as not carried.

    **Declared and refused, rather than withdrawn.** The palette declares them because they
    are the machine's own form fields, and each refusal reads the row the spec holds. The key
    is a literal, one arm per value, in this half and in the Rust half - which is what makes
    the pair countable by `tools/check_unported.py` instead of comparable by reading.
    """
    # **The mesh solve carries a different efficiency arithmetic**, refused rather than
    # silently ignored: `NaphtaliSandholmSolver.applyMurphreeEfficiencyToK` corrects a tray's
    # K-values by an Edmister `K^eta` proxy, where the sequential core's
    # `applyMurphreeCorrection` - the one this port carries - blends the vapour leaving a stage
    # against the vapour entering it.
    if solver_type == "naphtali_sandholm" and murphree_efficiency is not None:
        raise _unported.refuse("murphree_efficiency@solver_type=naphtali_sandholm")
    if solver_type is None or solver_type in ("direct_substitution", "naphtali_sandholm"):
        return
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


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.distillation_column")


#: `SIDE_DRAW_CANDIDATE_SCAN_STEP`, the grid the bounded scan walks.
_SIDE_DRAW_SCAN_STEP = 5.0e-3
#: `addSideDrawFlowSpecification`'s seed for an uncontrolled draw.
_SIDE_DRAW_SEED_FRACTION = 0.05
#: `wasSideDrawFractionAttempted`'s identity tolerance.
_SIDE_DRAW_IDENTITY = 1.0e-12


def _pumparound_return_specification(
    return_tray: int | None,
    draw_tray: int | None,
    draw_fraction: float | None,
    temperature_drop: float | None,
) -> tuple[int, int, float, float] | None:
    """**One pumparound with a return**, from its four declared scalars.

    One and not a list, for the reason the side-draw flow specification gives:
    `addLiquidPumparound` owns one draw tray each and refuses a second pumparound on the same
    tray. A partial statement is refused, which is the judgement the end specifications make.
    """
    if (
        return_tray is None
        or draw_tray is None
        or draw_fraction is None
        or temperature_drop is None
    ):
        if (
            return_tray is not None
            or draw_tray is not None
            or draw_fraction is not None
            or temperature_drop is not None
        ):
            raise InvalidInputError(
                "pumparound_return_tray",
                "a pumparound with a return is stated by its draw tray, its return tray, its "
                "fraction and its temperature drop together: `addLiquidPumparound(name, "
                "drawTray, returnTray, fraction, drop)` takes all four, and a declaration that "
                "states some of "
                "them says nothing",
            )
        return None
    return (
        int(draw_tray),
        int(return_tray),
        float(draw_fraction),
        float(temperature_drop),
    )


def _pumparound_returns_tear(
    solve: Callable[[Draws, tuple[StreamRecord | None, ...] | None], _States],
    components: list[str],
    draws: Draws,
    specification: tuple[int, int, float, float],
    tray_count: int,
    tolerance: float,
    max_iterations: int,
) -> _States:
    """`ColumnPumparound`'s return, iterated to a fixed point on its own flow.

    **The class's own loop and not a search**: each pass runs the whole column, rebuilds every
    return from the draw it has just produced - `updateReturnStream`'s own clone, restate at
    `T - drop`, re-flash - and stops when the return's *flow* stops moving. There are no
    candidates here and no acceptance rule, because a return is a stream the column already
    knows how to take.
    """
    draw_tray, return_tray, fraction, temperature_drop = specification
    for name, tray in (
        ("pumparound_draw_tray", draw_tray),
        ("pumparound_return_tray", return_tray),
    ):
        if tray >= tray_count:
            raise InvalidInputError(
                name,
                f"tray {tray} of a column with {tray_count} tray(s), which is not a tray it has",
            )
    if not 0.0 <= fraction <= 1.0:
        raise InvalidInputError(
            "pumparound_draw_fraction",
            f"a draw fraction of {fraction} is outside [0, 1], which `addLiquidPumparound` refuses",
        )

    # **The draw fraction is stated on the tray**, which is the half this model already carries.
    active_draws = _pumparound_with_fraction(draws, draw_tray, fraction, tray_count)
    returns: list[StreamRecord | None] = [None] * tray_count
    outcome = solve(active_draws, None)
    for _ in range(max_iterations):
        drawn = outcome.pumparound_n[draw_tray]
        if drawn <= 0.0:
            raise InvalidInputError(
                "pumparound_draw_tray",
                f"the pumparound drawing from tray {draw_tray} withdrew nothing, so its return "
                f"has no stream to carry to tray {return_tray}",
            )
        returned = _cooled_pumparound(components, outcome, draw_tray, temperature_drop)
        change = _return_flow_change(returns[return_tray], returned)
        returns[return_tray] = returned
        outcome = solve(active_draws, tuple(returns))
        if change <= tolerance:
            break
    return outcome


def _pumparound_with_fraction(
    draws: Draws, draw_tray: int, fraction: float, tray_count: int
) -> Draws:
    """The same draws with one tray's pumparound fraction stated."""
    gas, liquid, pumparound = draws if draws is not None else (None, None, None)
    vector: list[float] = [0.0] * tray_count if pumparound is None else list(pumparound)
    if len(vector) < tray_count:
        vector = vector + [0.0] * (tray_count - len(vector))
    vector[draw_tray] = fraction
    stated: tuple[float, ...] = tuple(vector)
    return (gas, liquid, stated)


def _cooled_pumparound(
    components: list[str], states: _States, draw_tray: int, temperature_drop: float
) -> StreamRecord:
    """**`updateReturnStream`**: the draw, cooled by the stated drop and re-flashed.

    The draw is the tray's own liquid phase scaled by the fraction, so its composition is the
    tray's liquid composition and its flow is `pumparound_n` - which is the whole of what
    `getLiquidPumparoundDrawStream` carries.
    """
    moles = float(states.pumparound_n[draw_tray])
    composition = list(states.tray_liquid_z[draw_tray])
    temperature = float(states.tray_temperature[draw_tray]) - temperature_drop
    if not temperature > 0.0:
        raise InvalidInputError(
            "pumparound_temperature_drop",
            f"cooling the draw by {temperature_drop} K puts the return at {temperature} K, "
            f"which is not a temperature: `updateReturnStream` raises on the same state",
        )
    return _stage.stream_at(
        list(components),
        moles,
        composition,
        temperature,
        float(states.tray_pressure[draw_tray]),
    )


def _return_flow_change(previous: StreamRecord | None, current: StreamRecord) -> float:
    """`updateReturnStream`'s own relative change: the two flows against the larger of them."""
    prior = 0.0 if previous is None else float(previous["n"])
    current_n = float(current["n"])
    return abs(prior - current_n) / max(prior, current_n, 1.0e-12)


def _side_draw_flow_specification(
    tray: int | None,
    phase: str | None,
    target: float | None,
    tolerance: float | None,
    max_iterations: int | None,
) -> SideDrawSpecification | None:
    """The one side-draw flow specification this model declares, or ``None``.

    **One and not a list, because one is what the port implements**: the class solves a list of
    them as *coordinated* tear variables, which is a different search. A target without a tray or
    a phase is refused, as the two end specifications are.
    """
    if tray is None and phase is None and target is None:
        return None
    if tray is None or phase is None or target is None:
        raise InvalidInputError(
            "side_draw_flow_tray",
            "a side-draw flow specification is stated by its tray, its phase and its target "
            "together: `addSideDrawFlowSpecification(tray, phase, flow, unit)` takes all three, "
            "and a declaration that states one or two of them says nothing",
        )
    if not (target == target and abs(target) != float("inf")) or target < 0.0:
        raise InvalidInputError(
            "side_draw_flow_target",
            f"a side-draw target flow of {target} is not finite and non-negative, which "
            f"`ColumnSideDrawSpecification`'s own constructor refuses",
        )
    if phase not in ("gas", "liquid"):
        raise InvalidInputError(
            "side_draw_flow_phase",
            f"`{phase}` is not a side-draw phase: `SideDrawPhase` carries `gas` and `liquid`",
        )
    return (
        int(tray),
        phase,
        float(target),
        1.0e-4 if tolerance is None else float(tolerance),
        12 if max_iterations is None else int(max_iterations),
    )


def _side_draw_flows_refused(draws: Draws) -> None:
    """A specification beside a pumparound is the class's *coordinated* tear."""
    if draws is not None and draws[2] is not None:
        raise InvalidInputError(
            "side_draw_flow_target",
            "a side-draw flow specification is stated beside a pumparound: the class solves the "
            "two together as coordinated tear variables, and this port carries the independent "
            "single-variable search",
        )


def _tray_count(number_of_stages: int, has_reboiler: bool, has_condenser: bool) -> int:
    return int(number_of_stages) + int(has_reboiler) + int(has_condenser)


def _side_draw_fraction(draws: Draws, tray: int, phase: str) -> float:
    vector = None if draws is None else draws[0 if phase == "gas" else 1]
    if vector is None or tray >= len(vector):
        return 0.0
    return float(vector[tray])


def _side_draw_maximum(draws: Draws, tray: int, phase: str) -> float:
    if phase == "gas":
        return 1.0
    pumparound = 0.0 if draws is None or draws[2] is None else float(draws[2][tray])
    return max(0.0, 1.0 - pumparound)


def _side_draw_with_fraction(
    draws: Draws, tray: int, phase: str, fraction: float, tray_count: int
) -> Draws:
    """The same draws at one tray's stated fraction, the vector grown where the caller stated
    none - which is what `setSideDrawFraction` does to the tray."""
    gas, liquid, pumparound = draws if draws is not None else (None, None, None)
    current = gas if phase == "gas" else liquid
    vector: list[float] = [0.0] * tray_count if current is None else list(current)
    if len(vector) < tray_count:
        vector = vector + [0.0] * (tray_count - len(vector))
    vector[tray] = fraction
    stated: tuple[float, ...] = tuple(vector)
    return (stated, liquid, pumparound) if phase == "gas" else (gas, stated, pumparound)


def _side_draw_mass_flow(states: object, components: list[str], tray: int, phase: str) -> float:
    """The draw's mass flow, kg/s - `getSideDrawStream(...).getFlowRate(unit)`.

    **The draw is the tray's own phase scaled**, so its composition is the tray's and its mass
    flow is the withdrawn moles times that composition's molar mass - which is exactly what
    `Stream::mass_flow` computes on the Rust side.
    """
    from azoth.eos import components as databank

    moles = (states.gas_side_draw_n if phase == "gas" else states.liquid_side_draw_n)[tray]  # type: ignore[attr-defined]
    composition = (states.tray_gas_z if phase == "gas" else states.tray_liquid_z)[tray]  # type: ignore[attr-defined]
    if moles <= 0.0 or not composition:
        return 0.0
    fluid = databank.mixture_of(components, eos="pr")[0]
    mass = sum(
        zi * component.molar_mass.to("kg/mol").magnitude  # type: ignore[union-attr]
        for zi, component in zip(composition, fluid.components, strict=True)
    )
    return float(moles) * float(mass)


def _side_draw_attempted(attempted: list[float], fraction: float) -> bool:
    if fraction != fraction or abs(fraction) == float("inf"):
        return True
    return any(abs(prior - fraction) <= _SIDE_DRAW_IDENTITY for prior in attempted)


def _side_draw_next(
    target: float,
    accepted_fractions: list[float],
    accepted_flows: list[float],
    attempted: list[float],
    maximum: float,
) -> float:
    """`selectNextSingleSideDrawCandidate`: accepted probes only can propose."""
    if not accepted_fractions:
        origin = attempted[-1] if attempted else 0.0
        return _side_draw_grid(origin, attempted, maximum)

    best_index = 0
    best_residual = float("inf")
    lower: int | None = None
    lower_residual = float("inf")
    upper: int | None = None
    upper_residual = float("inf")
    for index, flow in enumerate(accepted_flows):
        residual = abs(flow - target)
        if residual < best_residual:
            best_residual = residual
            best_index = index
        if flow <= target and residual < lower_residual:
            lower_residual = residual
            lower = index
        if flow >= target and residual < upper_residual:
            upper_residual = residual
            upper = index

    if lower is not None and upper is not None and lower != upper:
        denominator = accepted_flows[upper] - accepted_flows[lower]
        interpolated = (
            accepted_fractions[lower]
            + (target - accepted_flows[lower])
            * (accepted_fractions[upper] - accepted_fractions[lower])
            / denominator
            if abs(denominator) > 1.0e-12
            else 0.5 * (accepted_fractions[lower] + accepted_fractions[upper])
        )
        interpolated = min(max(0.0, interpolated), maximum)
        if not _side_draw_attempted(attempted, interpolated):
            return interpolated

    best_fraction = accepted_fractions[best_index]
    actual = accepted_flows[best_index]
    multiplicative = (
        0.0
        if target <= 1.0e-12
        else min(
            max(
                0.0,
                best_fraction + _SIDE_DRAW_SCAN_STEP
                if abs(actual) <= 1.0e-12
                else best_fraction * target / actual,
            ),
            maximum,
        )
    )
    if not _side_draw_attempted(attempted, multiplicative):
        return multiplicative
    nearest = _side_draw_nearest_rejected(best_fraction, attempted, accepted_fractions)
    if nearest == nearest and abs(nearest - best_fraction) < _SIDE_DRAW_SCAN_STEP:
        opposite = best_fraction + (5.0 * _SIDE_DRAW_SCAN_STEP) * (
            1.0 if best_fraction - nearest >= 0.0 else -1.0
        )
        opposite = min(max(0.0, opposite), maximum)
        if not _side_draw_attempted(attempted, opposite):
            return opposite
    return _side_draw_grid(best_fraction, attempted, maximum)


def _side_draw_nearest_rejected(
    origin: float, attempted: list[float], accepted: list[float]
) -> float:
    nearest = float("nan")
    nearest_distance = float("inf")
    for fraction in attempted:
        if _side_draw_attempted(accepted, fraction):
            continue
        distance = abs(fraction - origin)
        if distance < nearest_distance:
            nearest_distance = distance
            nearest = fraction
    return nearest


def _side_draw_grid(origin: float, attempted: list[float], maximum: float) -> float:
    """`nextUntriedSideDrawGridFractionAround`: the deterministic scan, above then below."""
    first_upper = -(-(origin + _SIDE_DRAW_IDENTITY) // _SIDE_DRAW_SCAN_STEP) * _SIDE_DRAW_SCAN_STEP
    first_lower = (origin - _SIDE_DRAW_IDENTITY) // _SIDE_DRAW_SCAN_STEP * _SIDE_DRAW_SCAN_STEP
    maximum_steps = int(-(-maximum // _SIDE_DRAW_SCAN_STEP)) + 1
    for step in range(maximum_steps + 1):
        offset = step * _SIDE_DRAW_SCAN_STEP
        upper = first_upper + offset
        if upper <= maximum + _SIDE_DRAW_IDENTITY:
            upper = min(max(0.0, upper), maximum)
            if not _side_draw_attempted(attempted, upper):
                return upper
        lower = first_lower - offset
        if lower >= -_SIDE_DRAW_IDENTITY:
            lower = min(max(0.0, lower), maximum)
            if not _side_draw_attempted(attempted, lower):
                return lower
    return float("nan")


def _side_draw_flow_tear(
    solve: Callable[[Draws], _States],
    components: list[str],
    draws: Draws,
    specification: SideDrawSpecification,
    number_of_stages: int,
    has_reboiler: bool,
    has_condenser: bool,
) -> _States:
    """`solveSingleSideDrawFlowSpecification`: the class's own candidate search.

    **The copy is the part a pure function does not need** - each candidate is another solve -
    and acceptance is ``Ok``, because this implementation refuses a state it cannot reach
    rather than publishing a fallback for it. **The one-shot continuation retry is not ported**:
    it re-solves a rejected fraction *warm*, and `_states` always seeds from its own `init`.
    """
    _side_draw_flows_refused(draws)
    tray, phase, target, tolerance, max_iterations = specification
    tray_count = _tray_count(number_of_stages, has_reboiler, has_condenser)
    if tray >= tray_count:
        raise InvalidInputError(
            "side_draw_flow_tray",
            f"a side-draw flow specification names tray {tray} of a column with {tray_count} "
            f"tray(s), which is not a tray it has",
        )
    if (tray == 0 and has_reboiler) or (tray + 1 == tray_count and has_condenser):
        raise InvalidInputError(
            "side_draw_flow_tray",
            f"a side-draw flow specification names {tray}, which is this port's reboiler or "
            f"condenser: the ends are `column::reboiler` and `column::condenser` rather than "
            f"stages, so there is no draw for the specification to move",
        )

    maximum = _side_draw_maximum(draws, tray, phase)
    candidate_fraction = _side_draw_fraction(draws, tray, phase)
    if candidate_fraction <= 0.0:
        candidate_fraction = _SIDE_DRAW_SEED_FRACTION
    attempted: list[float] = []
    accepted_fractions: list[float] = []
    accepted_flows: list[float] = []
    best_residual = float("inf")
    best_states: _States | None = None

    for _ in range(max_iterations):
        if candidate_fraction != candidate_fraction or abs(candidate_fraction) == float("inf"):
            break
        candidate_fraction = min(max(0.0, candidate_fraction), maximum)
        if _side_draw_attempted(attempted, candidate_fraction):
            candidate_fraction = _side_draw_next(
                target, accepted_fractions, accepted_flows, attempted, maximum
            )
            if candidate_fraction != candidate_fraction:
                break
        attempted.append(candidate_fraction)
        candidate_draws = _side_draw_with_fraction(
            draws, tray, phase, candidate_fraction, tray_count
        )
        try:
            states = solve(candidate_draws)
        except Exception:  # a rejected candidate is the class's own state
            candidate_fraction = _side_draw_next(
                target, accepted_fractions, accepted_flows, attempted, maximum
            )
            continue
        flow = _side_draw_mass_flow(states, components, tray, phase)
        accepted_fractions.append(candidate_fraction)
        accepted_flows.append(flow)
        residual = abs(flow - target) / max(1.0e-12, abs(target))
        if residual < best_residual:
            best_residual = residual
            best_states = states
        if residual <= tolerance:
            return states
        candidate_fraction = _side_draw_next(
            target, accepted_fractions, accepted_flows, attempted, maximum
        )

    if best_states is None:
        raise InvalidInputError(
            "side_draw_flow_target",
            f"every candidate for a {phase.upper()} draw of {target} kg/s on tray {tray} was "
            f"rejected by the column solve, so the class's own "
            f"`finalizeColumnTearConvergenceStatus` has nothing to publish",
        )
    return best_states


def _build_specification(
    kind: str | None, target: float | None, component: str | None, which: str
) -> Specification | None:
    """One end's specification from its three declared parameters.

    `None` where no type was stated, which is the temperature-pinned route the model takes by
    default. A target without a type, or a type without a target, is a declaration that cannot
    be read rather than a default worth guessing.
    """
    if kind is None:
        if target is not None or component is not None:
            raise InvalidInputError(
                "specification",
                f"the {which} has a target or a component and no type, so nothing says what it "
                f"constrains",
            )
        return None
    if kind not in SPECIFICATION_KINDS:
        raise InvalidInputError(
            "specification",
            f"{kind} is not one of `ColumnSpecification`'s five types: "
            f"{', '.join(SPECIFICATION_KINDS)}",
        )
    if target is None:
        raise InvalidInputError(
            "specification", f"the {which} specification states a type and no target"
        )
    if kind in ("product_purity", "component_recovery") and component is None:
        raise InvalidInputError(
            "specification",
            f"a {which} purity or recovery constrains a component, and none was stated",
        )
    return Specification(kind=kind, target=float(target), component=component)
