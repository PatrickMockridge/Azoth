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

**The ends are pinned by temperature.** ``setCondenserTemperature`` reaches the tray's own
``outTemperature``, so that tray's flash is a ``TPflash`` at the pin rather than at the mixed
enthalpy; a middle tray has no pin and flashes at its enthalpy.

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
from azoth.process.reference._column_stage import stage

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


class _States(NamedTuple):
    """The solved column, in SI magnitudes - the profile and both products."""

    #: Each tray's temperature, K, from the reboiler at stage 0 to the condenser.
    tray_temperature: tuple[float, ...]
    #: Each tray's pressure, Pa.
    tray_pressure: tuple[float, ...]
    #: Each tray's vapour traffic, mol/s.
    tray_gas_n: tuple[float, ...]
    #: Each tray's liquid traffic, mol/s.
    tray_liquid_n: tuple[float, ...]
    #: The distillate's molar flow.
    distillate_n: float
    #: The distillate's composition.
    distillate_z: tuple[float, ...]
    #: The distillate's molar enthalpy, J/mol.
    distillate_h: float
    #: The bottoms' molar flow.
    bottoms_n: float
    #: The bottoms' composition.
    bottoms_z: tuple[float, ...]
    #: The bottoms' molar enthalpy, J/mol.
    bottoms_h: float
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
        reboiler_temperature: the reboiler's temperature, which pins the bottom tray.
        condenser_temperature: the condenser's temperature, which pins the top tray.
        temperature_tolerance: the gate on the mean tray-temperature change.
        max_iterations: the iteration cap.
        murphree_efficiency: **not ported**; omitted is the ideal stage.
        solver_type: **only ``direct_substitution`` is ported**, which is the class's default.
        top_specification_type: **not ported**; the top is pinned by temperature instead.
        top_specification_target: **not ported**.
        top_specification_component: **not ported**.
        bottom_specification_type: **not ported**.
        bottom_specification_target: **not ported**.
        bottom_specification_component: **not ported**.

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
        ...     q(373.15, "K"),
        ...     q(253.15, "K"),
        ...     1e-06,
        ...     200,
        ... )
        >>> round(r.distillate_n.to("mol/s").magnitude, 4)
        3.7915
    """
    _refuse_unported(
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
    )

    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    p = input_to_si(spec, "feed_p", feed_p)
    t = input_to_si(spec, "feed_t", feed_t)
    top = input_to_si(spec, "top_pressure", top_pressure)
    bottom = input_to_si(spec, "bottom_pressure", bottom_pressure)
    t_reb = input_to_si(spec, "reboiler_temperature", reboiler_temperature)
    t_cond = input_to_si(spec, "condenser_temperature", condenser_temperature)

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
    )

    return DistillationColumnResult(
        tray_temperature=[from_si(value, "K") for value in states.tray_temperature],
        tray_pressure=[from_si(value, "Pa") for value in states.tray_pressure],
        tray_gas_n=[from_si(value, "mol/s") for value in states.tray_gas_n],
        tray_liquid_n=[from_si(value, "mol/s") for value in states.tray_liquid_n],
        distillate_n=from_si(states.distillate_n, "mol/s"),
        distillate_z=states.distillate_z,
        distillate_p=from_si(top if has_condenser else p, "Pa"),
        distillate_t=from_si(t_cond if has_condenser else t, "K"),
        distillate_h=from_si(states.distillate_h, "J/mol"),
        bottoms_n=from_si(states.bottoms_n, "mol/s"),
        bottoms_z=states.bottoms_z,
        bottoms_p=from_si(bottom if has_reboiler else p, "Pa"),
        bottoms_t=from_si(t_reb if has_reboiler else t, "K"),
        bottoms_h=from_si(states.bottoms_h, "J/mol"),
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
    reboiler_temperature: float,
    condenser_temperature: float,
    temperature_tolerance: float,
    max_iterations: int,
) -> _States:
    """The whole solve, in SI magnitudes: the one arithmetic the kernel and a dump share."""
    tray_count = number_of_stages + int(has_reboiler) + int(has_condenser)
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

    pressures = [
        bottom_pressure
        + (top_pressure - bottom_pressure) * i / (tray_count - 1)
        for i in range(tray_count)
    ]

    def pin(i: int) -> float | None:
        if i == 0 and has_reboiler:
            return reboiler_temperature
        if i == tray_count - 1 and has_condenser:
            return condenser_temperature
        return None

    gas: list[dict | None] = [None] * tray_count
    liquid: list[dict | None] = [None] * tray_count

    def run_with(i: int, inlets: list[dict]) -> None:
        out = stage(inlets, pressures[i], pin(i), 0.0)
        gas[i], liquid[i] = out["gas"], out["liquid"]

    def at_temperature(stream: dict, temperature: float) -> dict:
        if temperature != temperature:  # NaN
            return stream
        return _restate(components, stream, temperature)

    def inlets_of(i: int, seed: list[float] | None) -> list[dict]:
        inlets: list[dict] = []
        at = (lambda s: at_temperature(s, seed[i])) if seed is not None else (lambda s: s)
        if i > 0 and gas[i - 1] is not None:
            inlets.append(at(gas[i - 1]))
        if i + 1 < tray_count and liquid[i + 1] is not None:
            inlets.append(at(liquid[i + 1]))
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        return inlets

    def temperature(i: int) -> float:
        stream = gas[i] if gas[i] is not None else liquid[i]
        return float("nan") if stream is None else stream["t"]

    # ---- `init` ----
    run_with(feed_stage, [_feed(components, feed_n, feed_z, feed_t, feed_p)])
    if has_reboiler:
        run_with(0, [liquid[feed_stage]])

    feed_temperature = temperature(feed_stage)
    cond_temperature = (
        condenser_temperature if has_condenser else feed_temperature - 1.0
    )
    reb_temperature = liquid[0]["t"] if liquid[0] is not None else feed_temperature
    temperatures = [float("nan")] * tray_count
    temperatures[feed_stage] = feed_temperature
    delta_up = (feed_temperature - cond_temperature) / (tray_count - feed_stage - 1.0)
    delta_down = (reb_temperature - feed_temperature) / feed_stage
    delta = 0.0
    for i in range(feed_stage + 1, tray_count):
        delta += delta_up
        temperatures[i] = feed_temperature - delta
    delta = 0.0
    for i in range(feed_stage - 1, -1, -1):
        delta += delta_down
        temperatures[i] = feed_temperature + delta
    for i in range(tray_count):
        if pin(i) is not None:
            temperatures[i] = pin(i)

    for i in range(1, tray_count):
        inlets = [at_temperature(gas[i - 1], temperatures[i])]
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        run_with(i, inlets)
    for i in range(tray_count - 2, 0, -1):
        inlets = [
            at_temperature(gas[i - 1], temperatures[i]),
            at_temperature(liquid[i + 1], temperatures[i]),
        ]
        if i == feed_stage:
            inlets.append(_feed(components, feed_n, feed_z, feed_t, feed_p))
        run_with(i, inlets)
    if has_reboiler:
        run_with(0, [at_temperature(liquid[1], temperatures[0])])

    # ---- the sweeps ----
    relaxation = 1.0
    previous_combined = float("inf")
    temperature_residual = float("inf")
    iterations = 0
    for iteration in range(1, max_iterations + 1):
        iterations = iteration
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
            relaxation = max(
                MIN_SEQUENTIAL_RELAXATION, relaxation * RELAXATION_DECREASE_FACTOR
            )
        elif combined < previous_combined * 0.98:
            relaxation = min(
                MAX_ADAPTIVE_RELAXATION, relaxation * RELAXATION_INCREASE_FACTOR
            )
        previous_combined = combined

        if temperature_residual <= temperature_tolerance:
            break

    def outlet_enthalpy(i: int) -> float:
        total = 0.0
        for stream in (gas[i], liquid[i]):
            if stream is not None:
                total += stream["n"] * stream["h"]
        return total

    def inlet_enthalpy(i: int) -> float:
        return sum(s["n"] * s["h"] for s in inlets_of(i, None))

    distillate = gas[tray_count - 1] if has_condenser else liquid[tray_count - 1]
    bottoms = liquid[0] if has_reboiler else gas[0]
    reboiler_duty = (outlet_enthalpy(0) - inlet_enthalpy(0)) if has_reboiler else 0.0
    condenser_duty = (
        outlet_enthalpy(tray_count - 1) - inlet_enthalpy(tray_count - 1)
        if has_condenser
        else 0.0
    )

    feed_enthalpy = feed_n * _feed(components, feed_n, feed_z, feed_t, feed_p)["h"]
    products_enthalpy = distillate["n"] * distillate["h"] + bottoms["n"] * bottoms["h"]
    energy_residual = (
        abs(feed_enthalpy + reboiler_duty + condenser_duty - products_enthalpy)
        / abs(feed_enthalpy)
        if abs(feed_enthalpy) > 0.0
        else 0.0
    )

    mass_residual = 0.0
    for c, zi in enumerate(feed_z):
        supplied = feed_n * zi
        delivered = distillate["n"] * distillate["z"][c] + bottoms["n"] * bottoms["z"][c]
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
        tray_temperature=tuple(temperature(i) for i in range(tray_count)),
        tray_pressure=tuple(pressures),
        tray_gas_n=tuple(gas[i]["n"] if gas[i] is not None else 0.0 for i in range(tray_count)),
        tray_liquid_n=tuple(
            liquid[i]["n"] if liquid[i] is not None else 0.0 for i in range(tray_count)
        ),
        distillate_n=distillate["n"],
        distillate_z=tuple(distillate["z"]),
        distillate_h=distillate["h"],
        bottoms_n=bottoms["n"],
        bottoms_z=tuple(bottoms["z"]),
        bottoms_h=bottoms["h"],
        condenser_duty=condenser_duty,
        reboiler_duty=reboiler_duty,
        iterations=iterations,
        temperature_residual=temperature_residual,
        mass_residual=mass_residual,
        energy_residual=energy_residual,
    )


def _feed(
    components: list[str], n: float, z: list[float], t: float, p: float
) -> dict[str, object]:
    """The external feed, rebuilt so a seed can never reach it."""
    return _stage.stream_at(components, n, list(z), t, p)


def _restate(components: list[str], stream: dict, temperature: float) -> dict[str, object]:
    """A stream at another temperature, at its own pressure and composition."""
    return _stage.stream_at(components, stream["n"], list(stream["z"]), temperature, stream["p"])


def _refuse_unported(
    murphree_efficiency: float | None,
    solver_type: str | None,
    top_specification_type: str | None,
    top_specification_target: float | None,
    top_specification_component: str | None,
    bottom_specification_type: str | None,
    bottom_specification_target: float | None,
    bottom_specification_component: str | None,
) -> None:
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
    if solver_type is not None and solver_type != "direct_substitution":
        raise InvalidInputError(
            "solver_type",
            f"solver_type = {solver_type} is not ported. Only `direct_substitution` is, which "
            f"is the class's own default; `naphtali_sandholm` is `NaphtaliSandholmSolver` and "
            f"the rest are `ColumnSolverFactory`'s inside-out family, where `auto` is a ladder "
            f"rather than one method",
        )
    for which, end_temperature, spec_type, target, component in (
        ("top", "condenser_temperature", top_specification_type, top_specification_target,
         top_specification_component),
        ("bottom", "reboiler_temperature", bottom_specification_type, bottom_specification_target,
         bottom_specification_component),
    ):
        if spec_type is not None or target is not None or component is not None:
            raise InvalidInputError(
                "specification",
                f"a {which} product specification is not ported: `ColumnSpecification`'s five "
                f"types and two locations, and the outer loop that drives one to its target, "
                f"are what would close it. The {which} is pinned by `{end_temperature}` "
                f"instead, which is the class's other route",
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
