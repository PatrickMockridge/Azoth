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

from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DistillationColumnResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.process.reference import _column_stage as _stage
from azoth.process.reference._column_stage import (
    StreamRecord,
    reactive_stage,
    split_draws,
    stage,
    validate_draws,
)

#: The class's own adaptive-relaxation constants, from `DistillationColumn`'s initialisers.
MIN_SEQUENTIAL_RELAXATION = 0.5
MAX_ADAPTIVE_RELAXATION = 1.2
RELAXATION_INCREASE_FACTOR = 1.2
RELAXATION_DECREASE_FACTOR = 0.5
#: The floor on the *temperature* update's step, a different clamp from the streams'.
MIN_TEMPERATURE_RELAXATION = 0.2
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

#: The nine strategies the class carries and this port does not.
UNPORTED_SOLVERS = (
    "damped_substitution",
    "inside_out",
    "matrix_inside_out",
    "wegstein",
    "sum_rates",
    "newton",
    "naphtali_sandholm",
    "mesh_residual",
    "auto",
)


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
        murphree_efficiency: **not ported**; omitted is the ideal stage.
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
        The tray profile, both products, both duties and the three residuals.

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
    _refuse_unported(murphree_efficiency, solver_type)

    # **`setReactive`'s two forms, and neither states half a section** - the mirror of the Rust
    # model's own resolution.
    if (reactive_start_tray is None) != (reactive_end_tray is None):
        raise InvalidInputError(
            "reactive_start_tray",
            "a reactive section is stated by both bounds or by neither: `setReactive(true)` "
            "covers every middle tray and `setReactive(true, start, end)` a run of them, and the "
            "class has no form that states one end alone",
        )
    section: tuple[int, int] | None = None
    if reactive:
        start, end = reactive_start_tray, reactive_end_tray
        section = (-1, -1) if start is None or end is None else (int(start), int(end))
    elif reactive_start_tray is not None:
        raise InvalidInputError(
            "reactive",
            "a reactive section was stated without `reactive = true`, which is a declaration "
            "that says nothing",
        )
    else:
        section = None
    if section is not None and solver_type == "naphtali_sandholm":
        raise InvalidInputError(
            "reactive",
            "a reactive section under `naphtali_sandholm` is not ported: `NaphtaliSandholmSolver` "
            "reads its fugacities from the MESH equations, and this port's mesh does not route a "
            "tray's flash at all, so the flag would be ignored rather than honoured",
        )

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
        raise InvalidInputError(
            "gas_side_draw_fractions",
            "side draws under `naphtali_sandholm` are not ported: this port's mesh solves the "
            "MESH equations together and never forms a tray's own outlet streams",
        )

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

    states = _states(
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
        draws=draws,
    )
    warnings.extend(states.warnings)

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
        warnings=tuple(warnings),
    )


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
    draws: tuple[tuple[float, ...] | None, tuple[float, ...] | None, tuple[float, ...] | None]
    | None = None,
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
                        f"a draw is stated on {where}, which is not ported: the ends here are "
                        f"`column::reboiler` and `column::condenser` rather than stages",
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
    cond_temperature = end_temperature(tray_count - 1)
    if cond_temperature is None:
        cond_temperature = feed_temperature - 1.0
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
        abs(feed_enthalpy + reboiler_duty + condenser_duty - products_enthalpy) / abs(feed_enthalpy)
        if abs(feed_enthalpy) > 0.0
        else 0.0
    )

    mass_residual = 0.0
    for c in range(len(feed_z)):
        supplied = sum(float(feed["n"]) * float(feed["z"][c]) for feed in feeds)
        withdrawn = sum(draw["n"] * draw["z"][c] for draw in held)
        delivered = (
            distillate["n"] * distillate["z"][c] + bottoms["n"] * bottoms["z"][c] + withdrawn
        )
        if abs(supplied) > 1.0e-12:
            mass_residual = max(mass_residual, abs(supplied - delivered) / abs(supplied))

    if (
        temperature_residual > temperature_tolerance
        or mass_residual > MASS_BALANCE_TOLERANCE
        or energy_residual > ENTHALPY_BALANCE_TOLERANCE
    ):
        from azoth.core.errors import SolverNotConvergedError

        raise SolverNotConvergedError(
            iterations,
            max(temperature_residual, mass_residual, energy_residual),
            min(temperature_tolerance, MASS_BALANCE_TOLERANCE, ENTHALPY_BALANCE_TOLERANCE),
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


def _refuse_unported(murphree_efficiency: float | None, solver_type: str | None) -> None:
    """Refuse every parameter the palette declares and this tranche does not implement.

    **Declared and refused, rather than withdrawn.** The palette declares them because they
    are the machine's own form fields; each refusal names the class that would close it, so a
    reader learns what is owed rather than meeting an absence.
    """
    if murphree_efficiency is not None:
        raise InvalidInputError(
            "murphree_efficiency",
            f"a Murphree efficiency of {murphree_efficiency} is not ported: "
            f"`SimpleTray.setMurphreeEfficiency` and the per-tray correction the column solver "
            f"applies after each run are the classes that would close it. Omitted means the "
            f"ideal stage, which is the class's own default of one",
        )
    if solver_type is None or solver_type in ("direct_substitution", "naphtali_sandholm"):
        return
    raise InvalidInputError(
        "solver_type",
        f"solver_type = {solver_type} is not ported: `ColumnSolverFactory."
        f"{unported_solver_class(solver_type)}` is the class that would close it. **The "
        f"capture measures why it is refused rather than ported**: "
        f"`validation/neqsim/captures/process_column_solvers.tsv` puts every one of the ten "
        f"strategies within `2.5e-6` K of every other on the binary column's tray 1 and within "
        f"`1.1e-7` relative on its distillate, so they are path variants rather than different "
        f"physics",
    )


def unported_solver_class(strategy: str) -> str:
    """The `ColumnSolverFactory` class behind a strategy this port does not carry.

    **Named per strategy, because that is what a refusal owes a caller.** `columnSolver` hands
    back one of these for each `SolverType`, and `AutoSolver` is the ladder rather than a
    method: `candidateSolvers` returns `NAPHTALI_SANDHOLM` first, then `MATRIX_INSIDE_OUT`,
    `INSIDE_OUT` and `DAMPED_SUBSTITUTION`, and this port has the first and the last.
    """
    return {
        "damped_substitution": "DampedSubstitutionSolver",
        "inside_out": "InsideOutSolver",
        "matrix_inside_out": "MatrixInsideOutSolver",
        "wegstein": "WegsteinSolver",
        "sum_rates": "SumRatesSolver",
        "newton": "TemperatureNewtonSolver",
        "mesh_residual": "MeshResidualSolver",
    }.get(strategy, "AutoSolver")


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.distillation_column")


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
