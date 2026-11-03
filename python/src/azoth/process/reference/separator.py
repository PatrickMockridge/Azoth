"""``process.separator`` - the separator's kernel.

Spec: ``specs/models/process/separator.toml``

The Python twin of ``crates/azoth-process/src/models/separator.rs`` over
``crates/azoth-process/src/kernels/separator.rs``, written to mirror them.

# The vessel holds the temperature, and that is not a throttling

``Separator.run`` sets the pressure to ``P - pressureDrop`` and runs a ``TPflash``, so the
outlets come out at the **feed's** temperature - not at a temperature a caller states, and
not at the one an isenthalpic drop would give. Measured on the captured fluid, a 2 bar drop
takes -2705.35 J/mol in to a weighted -2364.95 out. That is the class's own model of the
vessel, and reproducing it means the balance is open whenever a pressure drop is asked for.
A ``heat_input`` is the one thing that moves the flash, and it moves it by ``Q / n`` J/mol -
``Heater.run``'s one-second convention, where the duty is W and the flow is mol/s.

# Entrainment moves material and not energy

``gas_in_liquid`` carries that fraction of the vapour's moles into the liquid outlet,
component by component, which is NeqSim's ``addPhaseFractionToPhase`` on its
``feed``/``mole`` basis. Both outlets are then at the same temperature, so the transfer
does not close the balance either - and it acts only where both phases exist, because a
from-phase or a to-phase that is not there has nothing to move.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.numerics import divide
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, SeparatorResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.pt_flash import pt_flash
from azoth.process.reference._phase_outlet import phase_enthalpy

#: The class's `DEFAULT_LIQUID_DENSITY_FOR_SIZING`, which the sizing velocity substitutes
#: where the system carries no liquid.
DEFAULT_LIQUID_DENSITY_FOR_SIZING = 1000.0
#: The class's `DEFAULT_LIQUID_DENSITY`, which the gas-load factor substitutes twice: for a
#: dry system and for a second phase no denser than the gas.
DEFAULT_LIQUID_DENSITY = 1000.0
#: The class's `MIN_LIQUID_GAS_DENSITY_DIFFERENCE`. **The load factor's rule and not the
#: sizing family's**: the two are deliberately separate.
MIN_LIQUID_GAS_DENSITY_DIFFERENCE = 10.0
#: The class's constructor defaults for the two inputs that shape the gas area.
DEFAULT_ORIENTATION = "horizontal"
DEFAULT_DESIGN_LIQUID_LEVEL_FRACTION = 0.8


def separator(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    heat_input: Q | None = None,
    internal_diameter: Q | None = None,
    design_gas_load_factor: Q | None = None,
    orientation: str | None = None,
    design_liquid_level_fraction: float | None = None,
    separator_length: Q | None = None,
    liquid_level_fraction: float | None = None,
    safety_factor: float | None = None,
    mechanical_design: bool | None = None,
    auto_size: bool | None = None,
    fg: float | None = None,
    volume_safety_factor: float | None = None,
    retention_time: Q | None = None,
    max_operation_pressure: Q | None = None,
    allowable_stress: Q | None = None,
    joint_efficiency: float | None = None,
    corrosion_allowance: Q | None = None,
) -> SeparatorResult:
    """Flash a stream into vapour and liquid outlets.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature, which the vessel holds.
        pressure_drop: the pressure drop across the vessel.
        heat_input: heat added before the flash, or ``None`` for a flash at the feed's
            temperature.
        gas_in_liquid: the fraction of the vapour's moles carried into the liquid outlet.
        internal_diameter: the vessel's internal diameter, or ``None`` for no capacity
            check. Stated with ``design_gas_load_factor`` or not at all.
        design_gas_load_factor: the Souders-Brown ``K`` the vessel is designed to.
        orientation: ``"horizontal"`` or ``"vertical"``; absent is the class's
            ``"horizontal"``.
        design_liquid_level_fraction: the horizontal gas area's derating; absent is the
            class's ``0.8``.
        separator_length: the vessel's tangent-to-tangent length; absent is the class's
            ``5.0`` m. Requires the diameter and design load factor.
        liquid_level_fraction: the operating level as a fraction of the maximum liquid
            height, which is what ``setLiquidLevel`` takes; absent is the class's ``0.5``.
        safety_factor: the factor ``sizeFromFlow`` multiplies the gas flow by; absent
            leaves ``required_diameter`` absent.

    Returns:
        Both outlets' records, and the five capacity getters where a vessel was stated.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero, the
            entrainment fraction is outside ``[0, 1]``, exactly one of the diameter and
            the design load factor is stated, an orientation or liquid level is stated
            without them, the orientation is not one of the two, or a density the
            Souders-Brown velocity divides by is not positive.
        OutOfRangeError: for a pressure drop or a fraction the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = separator(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.7, 0.3],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(0.0, "Pa"),
        ... )
        >>> round(r.vapour_n.to("mol/s").magnitude, 6)
        0.818221
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    # Converted for the unit check alone: the route below takes the quantity itself.
    input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)

    diameter = (
        None
        if internal_diameter is None
        else input_to_si(spec, "internal_diameter", internal_diameter)
    )
    load_factor = (
        None
        if design_gas_load_factor is None
        else input_to_si(spec, "design_gas_load_factor", design_gas_load_factor)
    )
    # A dimensionless input crosses as its number, not as a quantity: there is no unit to
    # convert, and `input_to_si` refuses a bare float where the spec's dimension is
    # dimensionless - the same rule every other fraction here follows.
    level_fraction = (
        None if design_liquid_level_fraction is None else float(design_liquid_level_fraction)
    )
    length = (
        None
        if separator_length is None
        else input_to_si(spec, "separator_length", separator_length)
    )
    operating_level = None if liquid_level_fraction is None else float(liquid_level_fraction)
    safety = None if safety_factor is None else float(safety_factor)
    design = bool(mechanical_design) if mechanical_design is not None else False
    auto = bool(auto_size) if auto_size is not None else False
    design_fg = None if fg is None else float(fg)
    volume_factor = None if volume_safety_factor is None else float(volume_safety_factor)
    retention = (
        None if retention_time is None else input_to_si(spec, "retention_time", retention_time)
    )
    max_pressure = (
        None
        if max_operation_pressure is None
        else input_to_si(spec, "max_operation_pressure", max_operation_pressure)
    )
    stress = (
        None
        if allowable_stress is None
        else input_to_si(spec, "allowable_stress", allowable_stress)
    )
    joint = None if joint_efficiency is None else float(joint_efficiency)
    corrosion = (
        None
        if corrosion_allowance is None
        else input_to_si(spec, "corrosion_allowance", corrosion_allowance) * 1000.0
    )
    apply_checks(
        checks.on_input,
        {
            "pressure_drop": drop,
            "gas_in_liquid": gas_in_liquid,
            "internal_diameter": diameter,
            "design_gas_load_factor": load_factor,
            "design_liquid_level_fraction": level_fraction,
            "separator_length": length,
            "liquid_level_fraction": operating_level,
            "safety_factor": safety,
            "fg": design_fg,
            "volume_safety_factor": volume_factor,
            "retention_time": retention,
            "max_operation_pressure": max_pressure,
            "allowable_stress": stress,
            "joint_efficiency": joint,
            "corrosion_allowance": corrosion,
        }.get,
        warnings,
    )

    # **A stated diameter needs its design basis; a design basis does not need a diameter.**
    # The mechanical design sizes its own diameter from the flow and the design `K`.
    if diameter is not None and load_factor is None:
        raise InvalidInputError(
            "design_gas_load_factor",
            "`internal_diameter` was stated and `design_gas_load_factor` was not, and the "
            "Souders-Brown capacity metric needs both - a vessel with a size but no design "
            "basis has nothing to check",
        )
    if design and auto:
        raise InvalidInputError(
            "auto_size",
            "`mechanical_design` and `auto_size` are two ways to design the same vessel; "
            "state one",
        )
    if auto:
        if diameter is not None:
            raise InvalidInputError(
                "internal_diameter",
                "`auto_size` sizes the vessel itself, so a stated diameter would be overwritten",
            )
        if volume_safety_factor is not None:
            raise InvalidInputError(
                "volume_safety_factor",
                "`auto_size` sets the volume factor from `safety_factor`; state it there",
            )
    if not design and not auto and load_factor is not None and diameter is None:
        raise InvalidInputError(
            "internal_diameter",
            "`design_gas_load_factor` was stated without `internal_diameter` and without "
            "`mechanical_design`, and the design basis has no vessel to be used on",
        )
    if design and load_factor is None:
        raise InvalidInputError(
            "design_gas_load_factor",
            "`mechanical_design` was stated without `design_gas_load_factor`, and the design's "
            "Souders-Brown basis is that number",
        )
    # The two that shape the gas area have no meaning without the vessel they shape.
    if diameter is None:
        if orientation is not None and not design and not auto:
            raise InvalidInputError(
                "orientation",
                "`orientation` was stated without `internal_diameter`, "
                "`design_gas_load_factor` and `mechanical_design`, and the gas area it shapes "
                "is not computed until the vessel is",
            )
        if design_liquid_level_fraction is not None:
            raise InvalidInputError(
                "design_liquid_level_fraction",
                "`design_liquid_level_fraction` was stated without `internal_diameter` and "
                "`design_gas_load_factor`, and the gas area it shapes is not computed until "
                "the vessel is",
            )
        for stated, value in (
            ("separator_length", separator_length),
            ("liquid_level_fraction", liquid_level_fraction),
            ("safety_factor", None if auto else safety_factor),
        ):
            if value is not None:
                raise InvalidInputError(
                    stated,
                    f"`{stated}` was stated without `internal_diameter` and "
                    f"`design_gas_load_factor`, and the vessel it sizes is not computed until "
                    f"both are stated",
                )
    if not design and not auto:
        for stated, value in (
            ("fg", fg),
            ("volume_safety_factor", volume_safety_factor),
            ("retention_time", retention_time),
            ("max_operation_pressure", max_operation_pressure),
            ("allowable_stress", allowable_stress),
            ("joint_efficiency", joint_efficiency),
            ("corrosion_allowance", corrosion_allowance),
        ):
            if value is not None:
                raise InvalidInputError(
                    stated, f"`{stated}` was stated without `mechanical_design`"
                )
    kind = DEFAULT_ORIENTATION if orientation is None else orientation
    if kind not in ("horizontal", "vertical"):
        raise InvalidInputError(
            "orientation",
            f"a vessel is `horizontal` or `vertical`, and {kind!r} is neither",
        )
    level = (
        DEFAULT_DESIGN_LIQUID_LEVEL_FRACTION if level_fraction is None else float(level_fraction)
    )

    duty = None if heat_input is None else input_to_si(spec, "heat_input", heat_input)
    # The temperature crosses as the quantity, not as its magnitude: the flash takes one
    # and `enthalpy_at` takes the other, and the route is where those two meet.
    states = _route(components, n, feed_z, p, feed_t, drop, gas_in_liquid, duty)

    capacity = (
        None
        if diameter is None or load_factor is None
        else _capacity_states(
            components,
            states.temperature.to("K").magnitude,
            states.p_out,
            list(feed_z),
            states.n_vapour,
            list(states.liquid_z),
            diameter,
            load_factor,
            kind,
            level,
        )
    )
    vessel = (
        None
        if not (design or auto) or (not auto and load_factor is None)
        else _mechanical_design(
            components,
            states.temperature.to("K").magnitude,
            states.p_out,
            list(states.y),
            states.n_vapour,
            list(states.liquid_z),
            states.n_liquid,
            kind == "horizontal",
            load_factor if load_factor is not None else 0.107,
            (0.5 if auto else 0.8) if design_fg is None else design_fg,
            (
                (1.2 if safety_factor is None else float(safety_factor))
                if auto
                else (1.0 if volume_factor is None else volume_factor)
            ),
            (120.0 if auto else 60.0) if retention is None else retention,
            100.0e5 if max_pressure is None else max_pressure,
            (425.0e6 if stress is None else stress) / 1.0e6,
            1.0 if joint is None else joint,
            0.0 if corrosion is None else corrosion,
            perform_sizing=auto,
        )
    )
    sizing = (
        None
        if diameter is None or load_factor is None
        else _sizing_getters(
            components,
            states.temperature.to("K").magnitude,
            states.p_out,
            list(feed_z),
            states.n_vapour,
            list(states.liquid_z),
            diameter,
            load_factor,
            kind,
            level,
            5.0 if length is None else length,
            0.5 if operating_level is None else operating_level,
            safety,
        )
    )

    return SeparatorResult(
        vapour_n=from_si(states.n_vapour, "mol/s"),
        vapour_z=states.y,
        vapour_p=from_si(states.p_out, "Pa"),
        vapour_t=states.temperature,
        vapour_h=from_si(states.vapour_h, "J/mol"),
        liquid_n=from_si(states.n_liquid, "mol/s"),
        liquid_z=states.liquid_z,
        liquid_p=from_si(states.p_out, "Pa"),
        liquid_t=states.temperature,
        liquid_h=from_si(states.liquid_h, "J/mol"),
        capacity_utilization=None if capacity is None else capacity.utilization,
        max_allowable_gas_velocity=(
            None
            if capacity is None or capacity.max_velocity is None
            else from_si(capacity.max_velocity, "m/s")
        ),
        gas_superficial_velocity=(
            None
            if capacity is None or capacity.superficial_velocity is None
            else from_si(capacity.superficial_velocity, "m/s")
        ),
        gas_load_factor=(
            None
            if capacity is None or capacity.load_factor is None
            else from_si(capacity.load_factor, "m/s")
        ),
        is_overloaded=None if capacity is None else capacity.overloaded,
        gas_area_fraction=None if sizing is None else sizing.gas_area_fraction,
        rated_gas_load_factor_by_phase=(
            None
            if sizing is None or not sizing.rated_by_phase
            else tuple(from_si(value, "m/s") for value in sizing.rated_by_phase)
        ),
        required_diameter=(
            None
            if sizing is None or sizing.required_diameter is None
            else from_si(sizing.required_diameter, "m")
        ),
        design_inner_diameter=(None if vessel is None else from_si(vessel.inner_diameter, "m")),
        design_tantan_length=(None if vessel is None else from_si(vessel.tantan_length, "m")),
        design_outer_diameter=(None if vessel is None else from_si(vessel.outer_diameter, "m")),
        design_wall_thickness=(None if vessel is None else from_si(vessel.wall_thickness, "m")),
        max_design_gas_volume_flow=(
            None
            if vessel is None
            else from_si(vessel.max_design_gas_volume_flow_m3_hr / 3600.0, "m**3/h")
        ),
        vessel_shell_weight=(
            None if vessel is None else from_si(vessel.vessel_shell_weight_kg, "kg")
        ),
        piping_weight=None if vessel is None else from_si(vessel.piping_weight_kg, "kg"),
        structural_weight=(None if vessel is None else from_si(vessel.structural_weight_kg, "kg")),
        electrical_weight=(None if vessel is None else from_si(vessel.electrical_weight_kg, "kg")),
        internals_weight=(None if vessel is None else from_si(vessel.internals_weight_kg, "kg")),
        total_weight=None if vessel is None else from_si(vessel.total_weight_kg, "kg"),
        module_width=None if vessel is None else from_si(vessel.module_width, "m"),
        module_length=None if vessel is None else from_si(vessel.module_length, "m"),
        module_height=None if vessel is None else from_si(vessel.module_height, "m"),
        vessel_volume=(None if vessel is None else from_si(vessel.vessel_volume_m3, "m**3")),
        hhll=None if vessel is None else from_si(vessel.hhll, "m"),
        hll=None if vessel is None else from_si(vessel.hll, "m"),
        nll=None if vessel is None else from_si(vessel.nll, "m"),
        lll=None if vessel is None else from_si(vessel.lll, "m"),
        llll=None if vessel is None else from_si(vessel.llll, "m"),
        weir_height=None if vessel is None else from_si(vessel.weir_height, "m"),
        nil=None if vessel is None else from_si(vessel.nil, "m"),
        effective_length_liquid=(
            None if vessel is None else from_si(vessel.effective_length_liquid, "m")
        ),
        effective_length_gas=(
            None if vessel is None else from_si(vessel.effective_length_gas, "m")
        ),
        inlet_nozzle_id=(None if vessel is None else from_si(vessel.inlet_nozzle_id, "m")),
        gas_outlet_nozzle_id=(
            None if vessel is None else from_si(vessel.gas_outlet_nozzle_id, "m")
        ),
        oil_outlet_nozzle_id=(
            None if vessel is None else from_si(vessel.oil_outlet_nozzle_id, "m")
        ),
        surge_volume=(None if vessel is None else from_si(vessel.surge_volume_m3, "m**3")),
        surge_time=None if vessel is None else from_si(vessel.surge_time_s, "s"),
        liquid_retention_time=(
            None if vessel is None else from_si(vessel.liquid_retention_time_s, "s")
        ),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class SepStates:
    """The separator's intermediates, in the order it computes them, in SI.

    ``beta`` is the one number here the record does not carry: the vapour fraction is
    implied by the two flows, and reading it back off them would be a ratio of two
    arithmetic results rather than the split the flash reported.
    """

    p_out: float
    n_vapour: float
    n_liquid: float
    y: tuple[float, ...]
    liquid_z: tuple[float, ...]
    temperature: Q
    vapour_h: float
    liquid_h: float


def _route(
    components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_p: float,
    feed_t: Q,
    pressure_drop: float,
    gas_in_liquid: float,
    heat_input: float | None,
) -> SepStates:
    """The separator's arithmetic, in SI, in the order it computes them.

    **One arithmetic, two consumers.** :func:`separator` builds the result from this and
    :func:`azoth.process.reference.tank` reaches the same steps for the vessel whose steady
    state is this one at zero drop and no entrainment - so the two models cannot drift into
    two answers. The kernel-level refusals are here rather than in the model because
    ``kernels/separator.rs`` is where Rust puts them.
    """
    # **A zero flow is not a state**, and the heat-per-mole divides by it.
    if feed_n <= 0.0:
        raise OutOfRangeError("feed_n", feed_n, "a unit operation needs a positive molar flow")
    if not 0.0 <= gas_in_liquid <= 1.0:
        raise InvalidInputError(
            "gas_in_liquid",
            f"an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]",
        )
    if feed_p <= 0.0:
        # The inlet is flashed before the vessel's arithmetic in Rust, so this is where the
        # two backends meet: the same refusal the flash's own pressure bound makes.
        raise OutOfRangeError("P", feed_p, "an absolute pressure; zero and below are not states")
    p_out = feed_p - pressure_drop
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {pressure_drop} Pa takes a {feed_p} Pa inlet to {p_out} Pa, "
            f"which is not a pressure",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    if heat_input is None:
        # The vessel holds the feed's temperature, so the flash is a TP one and the
        # outlets share the temperature they came in at.
        flash = pt_flash(mixture, feed_t, from_si(p_out, "Pa"), feed_z)
        temperature = feed_t
        x, y = _phase_compositions(flash.phase, feed_z, list(flash.x), list(flash.y))
        beta = _vapour_fraction(flash.phase, flash.vapour_fraction)
    else:
        # W over (mol/s) is J/mol, for the reason `Heater.run` gives.
        h_in, _ = enthalpy_at(mixture, ideal_gas, feed_t.to("K").magnitude, feed_p, feed_z)
        moved = ph_flash_solve(
            mixture,
            ideal_gas,
            from_si(p_out, "Pa"),
            from_si(h_in + divide(heat_input, feed_n), "J/mol"),
            feed_z,
        )
        temperature = moved.T
        x, y = _phase_compositions(moved.phase, feed_z, list(moved.x), list(moved.y))
        beta = _vapour_fraction(moved.phase, moved.vapour_fraction)

    n_vapour = feed_n * beta
    n_liquid = feed_n * (1.0 - beta)
    # Only where both phases exist: a from-phase or a to-phase that is not there has
    # nothing to move, which is NeqSim's own guard.
    carried = n_vapour * gas_in_liquid if n_vapour > 0.0 and n_liquid > 0.0 else 0.0

    liquid_n = n_liquid + carried
    liquid_z = (
        [(n_liquid * x[i] + carried * y[i]) / liquid_n for i in range(len(x))]
        if carried > 0.0
        else x
    )

    # **The outlet is a phase, except where the class re-runs it as a stream.** `Separator.run`
    # calls `liquidOutStream.run(id)` under `gasInLiquid != 0.0`, and a stream `run` is a
    # `TPflash` of that composition at the outlet's own temperature - so the liquid is a phase
    # root while nothing was carried into it and a re-flash once something was. The vapour has
    # no such branch here: the two fields that re-run it, `oilInGas` and `waterInGas`, are the
    # three-phase entry's. See `_phase_outlet`, which is where the difference is measured.
    t_out = temperature.to("K").magnitude
    vapour_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, y, liquid=False)
    if gas_in_liquid != 0.0:
        liquid_h, _ = enthalpy_at(mixture, ideal_gas, t_out, p_out, list(liquid_z))
    else:
        liquid_h = phase_enthalpy(mixture, ideal_gas, t_out, p_out, list(liquid_z), liquid=True)

    return SepStates(
        p_out=p_out,
        n_vapour=n_vapour - carried,
        n_liquid=liquid_n,
        y=tuple(y),
        liquid_z=tuple(liquid_z),
        temperature=temperature,
        vapour_h=vapour_h,
        liquid_h=liquid_h,
    )


@dataclass(frozen=True, slots=True)
class CapacityStates:
    """The five capacity getters, in the order the class computes them.

    The three velocities are absent where the system carries no gas: the class's getters
    would read phase 0, which is then the liquid, and a metric formed from it is not a gas
    metric. `gas_density`, `liquid_density`, `gas_area`, `volumetric_flow` and `max_flow`
    are the intermediates on the way there, carried because `layers` dumps them under the
    capture's own names rather than recomputing them.
    """

    utilization: float
    max_velocity: float | None
    superficial_velocity: float | None
    load_factor: float | None
    overloaded: bool
    gas_density: float | None
    liquid_density: float | None
    gas_area: float | None
    volumetric_flow: float | None
    max_flow: float | None
    #: Whether the system carries a condensed phase. **Not `liquid_density is not None`**:
    #: the sizing fallback is a density too, and a dry system has none to print.
    has_liquid: bool


@dataclass(frozen=True, slots=True)
class SizingStates:
    """The sizing family's evidence, before the gas-load factor is formed.

    **The assessment reads this and not the load factor.** `SeparatorCapacityAssessment` folds
    the sizing velocity and the two flows; the load factor's near-dry floor is a getter it
    never reports, so a state the floor refuses is not a state the assessment cannot evidence.
    """

    gas_density: float
    liquid_density: float
    gas_area: float
    volumetric_flow: float
    max_velocity: float
    superficial_velocity: float
    max_flow: float
    utilization: float
    has_liquid: bool


def _sizing_states(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    gas_n: float,
    liquid_z: list[float],
    diameter_m: float,
    design_load_factor: float,
    orientation: str,
    level_fraction: float,
) -> SizingStates | None:
    """The sizing half: densities, area, flows, velocity and utilisation.

    ``None`` is the class's own first answer - a stream with no gas phase uses none of its
    capacity - and is not an error. ``liquid_z`` is the entrained outlet's composition, whose
    density is taken at the root that composition's own flash puts it on.
    """
    from azoth.eos.reference.hydrate_inhibitor_wt import GAS
    from azoth.process.reference.capacity import _labelled_phases, _mass_density, _molar_mass

    mixture, labelled, _shares = _labelled_phases(components, t, p, feed_z, from_feed=True)
    gas = next((entry for entry in labelled if entry[0] == GAS), None)
    if gas is None:
        return None
    _label, gas_z, gas_root = gas
    gas_density = _mass_density(mixture, t, p, gas_z, gas_root)
    gas_molar_mass = _molar_mass(mixture, gas_z)
    if not math.isfinite(gas_density) or gas_density <= 0.0:
        raise InvalidInputError(
            "capacity_utilization",
            f"the gas density is {gas_density} kg/m3, and the Souders-Brown velocity divides "
            f"by its square root",
        )

    has_liquid = any(entry[0] != GAS for entry in labelled)
    # The liquid's own density, or the sizing fallback where the system carries none.
    if has_liquid:
        entrained = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(liquid_z))
        liquid_density = _mass_density(mixture, t, p, list(liquid_z), entrained.z_liquid)
    else:
        liquid_density = DEFAULT_LIQUID_DENSITY_FOR_SIZING

    area = math.pi * diameter_m * diameter_m / 4.0
    gas_area = area * (1.0 - level_fraction) if orientation == "horizontal" else area
    volumetric_flow = gas_n * gas_molar_mass / gas_density

    # **The sizing rule has no floor**: a liquid denser than the gas is all it needs, and the
    # class's `getCapacityUtilization` would answer NaN where this refuses.
    if not math.isfinite(liquid_density) or liquid_density <= gas_density:
        raise InvalidInputError(
            "capacity_utilization",
            f"the liquid phase is {liquid_density} kg/m3 against the gas's {gas_density}, so "
            f"`rho_l - rho_g` is not positive and the Souders-Brown velocity is the square "
            f"root of a negative number. NeqSim answers NaN there; this refuses the state",
        )
    max_velocity = design_load_factor * math.sqrt((liquid_density - gas_density) / gas_density)
    max_flow = max_velocity * gas_area
    if not math.isfinite(max_flow) or max_flow <= 0.0:
        raise InvalidInputError(
            "capacity_utilization",
            f"the maximum allowable gas flow is {max_flow} m3/s, and utilisation divides by it",
        )
    return SizingStates(
        gas_density=gas_density,
        liquid_density=liquid_density,
        gas_area=gas_area,
        volumetric_flow=volumetric_flow,
        max_velocity=max_velocity,
        superficial_velocity=volumetric_flow / gas_area,
        max_flow=max_flow,
        utilization=volumetric_flow / max_flow,
        has_liquid=has_liquid,
    )


def _capacity_states(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    gas_n: float,
    liquid_z: list[float],
    diameter_m: float,
    design_load_factor: float,
    orientation: str,
    level_fraction: float,
) -> CapacityStates:
    """`Separator`'s steady capacity getters, on the system ``run`` leaves behind.

    **Two density rules, and they are not one.** The sizing family - the maximum velocity
    and the utilisation over it - reads the liquid's own density and substitutes `1000.0`
    only where there is no liquid, with **no floor**. The gas-load factor applies the
    class's `10.0` floor to the same density first. A near-dry state makes the two disagree
    by construction.
    """
    sizing = _sizing_states(
        components,
        t,
        p,
        feed_z,
        gas_n,
        liquid_z,
        diameter_m,
        design_load_factor,
        orientation,
        level_fraction,
    )
    if sizing is None:
        # The class's own first answer: a stream with no gas phase uses none of its capacity.
        return CapacityStates(0.0, None, None, None, False, None, None, None, None, None, False)

    # **The load factor's own rule, applied after the sizing one.** The floor substitutes
    # `1000.0` for a second phase no denser than the gas; the positivity check is repeated
    # against the substituted value, which is a state the class would also answer NaN for.
    floored = _load_factor_liquid_density(sizing.liquid_density, sizing.gas_density)
    if not math.isfinite(floored) or floored <= sizing.gas_density:
        raise InvalidInputError(
            "gas_load_factor",
            f"the liquid phase is {floored} kg/m3 after the near-dry floor against the gas's "
            f"{sizing.gas_density}, and the load factor divides by their difference",
        )
    term = (floored - sizing.gas_density) / sizing.gas_density
    load_factor = sizing.superficial_velocity * math.sqrt(1.0 / term)
    if not math.isfinite(load_factor):
        raise InvalidInputError("gas_load_factor", f"the gas load factor is {load_factor} m/s")
    return CapacityStates(
        utilization=sizing.utilization,
        max_velocity=sizing.max_velocity,
        superficial_velocity=sizing.superficial_velocity,
        load_factor=load_factor,
        overloaded=sizing.utilization > 1.0,
        gas_density=sizing.gas_density,
        liquid_density=sizing.liquid_density,
        gas_area=sizing.gas_area,
        volumetric_flow=sizing.volumetric_flow,
        max_flow=sizing.max_flow,
        has_liquid=sizing.has_liquid,
    )


@dataclass(frozen=True, slots=True)
class SizingGetters:
    """`getGasLoadFactor(int)` and `sizeFromFlow(safety)`, in SI.

    ``rated_by_phase`` is one entry per liquid phase, in the class's phase order; the
    phase-indexed overload divides its radicand by the **operating** gas area fraction while
    the superficial velocity itself is the design-area one.
    """

    gas_area_fraction: float
    liquid_volume: float
    separator_volume: float
    rated_by_phase: tuple[float, ...]
    required_diameter: float | None


def _sizing_getters(
    components: list[str],
    t: float,
    p: float,
    feed_z: list[float],
    gas_n: float,
    liquid_z: list[float],
    diameter_m: float,
    design_load_factor: float,
    orientation: str,
    design_level_fraction: float,
    separator_length_m: float,
    liquid_level_fraction: float,
    safety_factor: float | None,
) -> SizingGetters | None:
    """The phase-indexed load factors and the closed-form sizing.

    ``None`` is the class's own first answer for a system with no gas phase: the sizing
    returns without a diameter and the capacity family publishes nothing.
    """
    from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, GAS, OIL
    from azoth.process.reference.capacity import _labelled_phases, _mass_density

    sizing = _sizing_states(
        components,
        t,
        p,
        feed_z,
        gas_n,
        liquid_z,
        diameter_m,
        design_load_factor,
        orientation,
        design_level_fraction,
    )
    if sizing is None:
        return None

    horizontal = orientation == "horizontal"
    cross_area = math.pi * diameter_m * diameter_m / 4.0
    separator_volume = cross_area * separator_length_m
    level = liquid_level_fraction * (diameter_m if horizontal else separator_length_m)
    liquid_volume = (
        _liquid_area(level, diameter_m) * separator_length_m if horizontal else cross_area * level
    )
    gas_area_fraction = 1.0 - liquid_volume / separator_volume if horizontal else 1.0
    if not math.isfinite(gas_area_fraction) or gas_area_fraction <= 0.0:
        raise InvalidInputError(
            "liquid_level_fraction",
            f"the operating liquid level of {liquid_level_fraction} leaves a gas area "
            f"fraction of {gas_area_fraction}, and the phase-indexed load factor divides by it",
        )

    mixture, labelled, _shares = _labelled_phases(components, t, p, feed_z, from_feed=True)
    rated: list[float] = []
    oil_density: float | None = None
    aqueous_density: float | None = None
    for label, phase_z, root in labelled:
        if label == GAS:
            continue
        density = _mass_density(mixture, t, p, list(phase_z), root)
        if label == OIL:
            oil_density = density
        elif label == AQUEOUS:
            aqueous_density = density
        floored = _load_factor_liquid_density(density, sizing.gas_density)
        if not math.isfinite(floored) or floored <= sizing.gas_density:
            raise InvalidInputError(
                "gas_load_factor",
                f"the liquid phase is {floored} kg/m3 after the near-dry floor against the "
                f"gas's {sizing.gas_density}, and the indexed load factor divides by their "
                f"difference",
            )
        term = (1.0 / gas_area_fraction) * (floored - sizing.gas_density) / sizing.gas_density
        if not math.isfinite(term) or term <= 0.0:
            raise InvalidInputError(
                "gas_load_factor",
                f"the phase-indexed load factor's radicand is {term}",
            )
        rated.append(sizing.superficial_velocity * math.sqrt(1.0 / term))

    required_diameter: float | None = None
    if safety_factor is not None:
        liquid_density = (
            oil_density
            if oil_density is not None
            else aqueous_density
            if aqueous_density is not None
            else DEFAULT_LIQUID_DENSITY_FOR_SIZING
        )
        if not math.isfinite(liquid_density) or liquid_density <= sizing.gas_density:
            raise InvalidInputError(
                "safety_factor",
                f"the liquid phase is {liquid_density} kg/m3 against the gas's "
                f"{sizing.gas_density}, so the Souders-Brown velocity `sizeFromFlow` divides "
                f"by is not positive",
            )
        max_velocity = design_load_factor * math.sqrt(
            (liquid_density - sizing.gas_density) / sizing.gas_density
        )
        gas_area_share = 1.0 - design_level_fraction if horizontal else 1.0
        required_gas_area = sizing.volumetric_flow * safety_factor / max_velocity
        total_area = required_gas_area / gas_area_share
        if not math.isfinite(total_area) or total_area <= 0.0:
            raise InvalidInputError(
                "safety_factor",
                f"`sizeFromFlow` answers a total area of {total_area} m2",
            )
        required_diameter = math.sqrt(4.0 * total_area / math.pi)

    return SizingGetters(
        gas_area_fraction=gas_area_fraction,
        liquid_volume=liquid_volume,
        separator_volume=separator_volume,
        rated_by_phase=tuple(rated),
        required_diameter=required_diameter,
    )


def _liquid_area(level: float, diameter: float) -> float:
    """`Separator.liquidArea` on a horizontal vessel: the circular segment below ``level``."""
    cross_area = math.pi * diameter * diameter / 4.0
    if level <= 0.0:
        return 0.0
    if level >= diameter:
        return cross_area
    radius = diameter / 2.0
    if level < radius:
        d = radius - level
        theta = math.acos(d / radius)
        a = radius * math.sin(theta)
        return theta * radius * radius - a * d
    if level > radius:
        d = level - radius
        theta = math.acos(d / radius)
        a = radius * math.sin(theta)
        return (math.pi - theta) * radius * radius + a * d
    return 0.5 * math.pi * radius * radius


@dataclass(frozen=True, slots=True)
class MechanicalDesign:
    """`SeparatorMechanicalDesign.calcDesign`'s vessel, in SI."""

    inner_diameter: float
    tantan_length: float
    outer_diameter: float
    wall_thickness: float
    max_design_volume_flow_m3_hr: float
    max_design_gas_volume_flow_m3_hr: float
    vessel_shell_weight_kg: float
    piping_weight_kg: float
    structural_weight_kg: float
    electrical_weight_kg: float
    internals_weight_kg: float
    total_weight_kg: float
    module_width: float
    module_length: float
    module_height: float
    vessel_volume_m3: float
    hhll: float
    hll: float
    nll: float
    lll: float
    llll: float
    weir_height: float
    nil: float
    effective_length_liquid: float
    effective_length_gas: float
    inlet_nozzle_id: float
    gas_outlet_nozzle_id: float
    oil_outlet_nozzle_id: float
    surge_volume_m3: float
    surge_time_s: float
    liquid_retention_time_s: float


#: NeqSim's gas constant, which its phase volume uses.
_NEQSIM_R = 8.3144621


def _mechanical_design(
    components: list[str],
    t: float,
    p: float,
    gas_z: list[float],
    gas_n: float,
    liquid_z: list[float],
    liquid_n: float,
    horizontal: bool,
    gas_load_factor: float,
    fg: float,
    volume_safety_factor: float,
    retention_time_s: float,
    max_operation_pressure_pa: float,
    allowable_stress_mpa: float,
    joint_efficiency: float,
    corrosion_allowance_mm: float,
    vessel: str = "separator",
    perform_sizing: bool = False,
) -> MechanicalDesign:
    """`calcDesign` on the two outlets, in the order the class computes it.

    **The volume flows are the cubic's own, not mass over the physical-properties density.**
    `Phase.getVolume()` is untranslated while the density carries the Peneloux shift, which the
    capture's PR row measures as `0.11035858` m3 against `0.110877`.
    """
    from azoth.core.units import quantity
    from azoth.eos.reference.phase_transport import phase_transport
    from azoth.process.reference.capacity import _mass_density

    mixture, _ideal_gas = _components.mixture_of(components, eos="pr")
    gas_flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(gas_z))
    gas_root = gas_flash.z_vapour if gas_flash.phase == "two_phase" else gas_flash.z_liquid
    liquid_flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), list(liquid_z))
    gas_density = _mass_density(mixture, t, p, list(gas_z), gas_root)
    liquid_density = _mass_density(mixture, t, p, list(liquid_z), liquid_flash.z_liquid)
    liquid_viscosity = (
        phase_transport(components, "oil", quantity(t, "K"), quantity(p, "Pa"), list(liquid_z))
        .mu.to("Pa*s")
        .magnitude
    )
    if not (
        math.isfinite(gas_density)
        and gas_density > 0.0
        and math.isfinite(liquid_density)
        and liquid_density > gas_density
        and math.isfinite(liquid_viscosity)
        and liquid_viscosity > 0.0
    ):
        raise InvalidInputError(
            "mechanical_design",
            f"the design needs a gas denser than zero, a liquid denser than the gas and a "
            f"positive liquid viscosity; it has {gas_density}, {liquid_density} and "
            f"{liquid_viscosity}",
        )

    gas_volume_flow = gas_n * gas_root * _NEQSIM_R * t / p
    liquid_volume_flow = liquid_n * liquid_flash.z_liquid * _NEQSIM_R * t / p

    max_design_volume_flow = volume_safety_factor * gas_volume_flow * 3600.0
    max_gas_velocity = gas_load_factor * math.sqrt((liquid_density - gas_density) / gas_density)
    gas_area_fraction = fg if horizontal else 1.0
    liquid_area_fraction = 1.0 - fg if horizontal else 1.0
    inner_diameter = math.sqrt(
        4.0 * (max_design_volume_flow / 3600.0) / (math.pi * max_gas_velocity * gas_area_fraction)
    )
    cross_area = math.pi * (inner_diameter / 2.0) ** 2
    max_design_gas_volume_flow = max_gas_velocity * cross_area * gas_area_fraction * 3600.0
    if vessel == "scrubber":
        # **The scrubber's length is five diameters and nothing else sizes it.** Its class
        # never reads the retention time in `calcDesign`.
        tantan_length = inner_diameter * 5.0
    elif perform_sizing and liquid_volume_flow <= 1.0e-10:
        # `performSizingCalculations`' no-liquid branch.
        tantan_length = inner_diameter * 4.0
    else:
        tantan_length = math.sqrt(
            4.0
            * retention_time_s
            * liquid_volume_flow
            * volume_safety_factor
            / (math.pi * inner_diameter * inner_diameter * liquid_area_fraction)
        )
        total = tantan_length + inner_diameter
        if total / inner_diameter > 6.0 or total / inner_diameter < 3.0:
            # `calcDesign` resets to five diameters; `performSizingCalculations` to four.
            tantan_length = inner_diameter * (4.0 if perform_sizing else 5.0)

        if not perform_sizing:
            bubble_velocity = (
                9.82 * (250.0e-6) ** 2 * (liquid_density - gas_density) / (18.0 * liquid_viscosity)
            )
            ar = liquid_volume_flow / bubble_velocity
            daim = math.sqrt(ar / 4.0)
            if daim > inner_diameter:
                inner_diameter = daim
                tantan_length = 4.0 * daim

    if perform_sizing:
        # `performSizingCalculations` on the design-standard state `autoSize` read, which on
        # the capture has no pressure-vessel code: wall, weights and outer diameter stay at
        # their zero defaults and `volumeTotal` stays at the class's 2.0.
        wall_thickness = 0.0
        outer_diameter = inner_diameter
        vessel_shell_weight = 0.0
        piping_weight = 0.0
        structural_weight = 0.0
        electrical_weight = 0.0
        total_weight = 0.0
        vessel_volume = 2.0
    else:
        pressure_bar = max_operation_pressure_pa / 1.0e5
        wall_thickness_mm = (
            pressure_bar
            / 10.0
            * inner_diameter
            * 1.0e3
            / (2.0 * allowable_stress_mpa * joint_efficiency - 1.2 * pressure_bar / 10.0)
            + corrosion_allowance_mm
        )
        wall_thickness = wall_thickness_mm / 1000.0
        outer_diameter = inner_diameter + 2.0 * wall_thickness
        vessel_shell_weight = (
            0.032 * (wall_thickness * 1.0e3) * (inner_diameter * 1.0e3) * tantan_length
        )
        piping_weight = vessel_shell_weight * 0.4
        structural_weight = vessel_shell_weight * 0.1
        electrical_weight = vessel_shell_weight * 0.08
        total_weight = (
            vessel_shell_weight + piping_weight + structural_weight + electrical_weight
        )
        vessel_volume = math.pi * (inner_diameter / 2.0) ** 2 * tantan_length
    module_width = inner_diameter * 2.0
    module_length = tantan_length * 1.5
    module_height = inner_diameter * 2.0 if vessel == "scrubber" else inner_diameter * 2.0 + 1.0

    effective_length_liquid = tantan_length * 0.82
    effective_length_gas = tantan_length * 0.64

    total_volume_flow_m3_hr = (gas_volume_flow + liquid_volume_flow) * 3600.0
    total_mass_flow = (gas_volume_flow * gas_density + liquid_volume_flow * liquid_density) * 3600.0
    mix_density = total_mass_flow / total_volume_flow_m3_hr
    nozzle_velocity = math.sqrt(8000.0 / mix_density)
    inlet_nozzle_id = max(
        0.05,
        math.ceil(
            math.sqrt(4.0 * ((total_volume_flow_m3_hr / 3600.0) / nozzle_velocity) / math.pi) / 0.05
        )
        * 0.05,
    )
    gas_outlet_nozzle_id = (
        math.ceil(math.sqrt(4.0 * (gas_volume_flow / 20.0) / math.pi) / 0.05) * 0.05
    )
    oil_outlet_nozzle_id = max(
        0.05,
        math.ceil(math.sqrt(4.0 * (liquid_volume_flow / 1.5) / math.pi) / 0.025) * 0.025,
    )

    hhll = 0.80 * inner_diameter
    hll = 0.70 * inner_diameter
    nll = 0.50 * inner_diameter
    lll = 0.30 * inner_diameter
    llll = 0.15 * inner_diameter
    weir_height = 0.25 * inner_diameter
    nil = 0.20 * inner_diameter

    def segment(level: float) -> float:
        radius = inner_diameter / 2.0
        height = min(level, inner_diameter)
        theta = 2.0 * math.acos((radius - height) / radius)
        return (radius * radius / 2.0) * (theta - math.sin(theta)) * effective_length_liquid

    surge_volume = segment(hll) - segment(nll)
    liquid_hourly = liquid_volume_flow / 3600.0
    liquid_retention_time = segment(nll) / liquid_hourly
    surge_time = surge_volume / liquid_hourly

    return MechanicalDesign(
        inner_diameter=inner_diameter,
        tantan_length=tantan_length,
        outer_diameter=outer_diameter,
        wall_thickness=wall_thickness,
        max_design_volume_flow_m3_hr=max_design_volume_flow,
        max_design_gas_volume_flow_m3_hr=max_design_gas_volume_flow,
        vessel_shell_weight_kg=vessel_shell_weight,
        piping_weight_kg=piping_weight,
        structural_weight_kg=structural_weight,
        electrical_weight_kg=electrical_weight,
        internals_weight_kg=0.0,
        total_weight_kg=total_weight,
        module_width=module_width,
        module_length=module_length,
        module_height=module_height,
        vessel_volume_m3=vessel_volume,
        hhll=hhll,
        hll=hll,
        nll=nll,
        lll=lll,
        llll=llll,
        weir_height=weir_height,
        nil=nil,
        effective_length_liquid=effective_length_liquid,
        effective_length_gas=effective_length_gas,
        inlet_nozzle_id=inlet_nozzle_id,
        gas_outlet_nozzle_id=gas_outlet_nozzle_id,
        oil_outlet_nozzle_id=oil_outlet_nozzle_id,
        surge_volume_m3=surge_volume,
        surge_time_s=surge_time,
        liquid_retention_time_s=liquid_retention_time,
    )


def _phase_compositions(
    phase: Phase, z: list[float], x: list[float], y: list[float]
) -> tuple[list[float], list[float]]:
    """The two compositions a flash's phase implies.

    **A single phase has the feed's composition, and the flash's own ``x`` or ``y`` is not
    it.** The flash answers a trial there - and on its negative-flash branch that trial is the
    other side of the cubic, so a 95/5 methane/ethane gas at 50 bar would leave the vessel with
    the trial composition and a density 2.6 times the feed's. The two-phase branch keeps the
    flash's split.
    """
    if phase is Phase.TWO_PHASE:
        return x, y
    return list(z), list(z)


def _load_factor_liquid_density(liquid_density: float, gas_density: float) -> float:
    """The gas-load factor's near-dry floor, kept apart from the sizing rule.

    `Separator.getGasLoadFactor()` substitutes `1000.0` whenever `rho_l - rho_g` is under
    `MIN_LIQUID_GAS_DENSITY_DIFFERENCE`; the sizing family does **not**. A state two
    kilograms apart therefore has two liquid densities by design, and the pair is the
    reason this is a function rather than an inline expression.
    """
    if liquid_density - gas_density < MIN_LIQUID_GAS_DENSITY_DIFFERENCE:
        return DEFAULT_LIQUID_DENSITY
    return liquid_density


def _vapour_fraction(phase: Phase, beta: float | None) -> float:
    """The vapour fraction a flash's phase and ``beta`` imply.

    ``beta`` is absent for a single-phase answer - the flash reports ``None`` rather than a
    number outside ``[0, 1]``, because there is no split to report - so the phase is what
    says which side of the interval the outlet is on.
    """
    if phase in (Phase.ALL_VAPOUR, Phase.TRIVIAL):
        return 1.0
    return 0.0 if beta is None else float(beta)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.separator")
