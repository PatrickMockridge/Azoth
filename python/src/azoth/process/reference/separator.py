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
    apply_checks(
        checks.on_input,
        {
            "pressure_drop": drop,
            "gas_in_liquid": gas_in_liquid,
            "internal_diameter": diameter,
            "design_gas_load_factor": load_factor,
            "design_liquid_level_fraction": level_fraction,
        }.get,
        warnings,
    )

    # **Both or neither.** The capacity is a function of the vessel, and a caller who stated
    # half of it has stated nothing: a vessel with a size but no design basis, or the reverse,
    # has no capacity to be using.
    if (diameter is None) != (load_factor is None):
        named, missing = (
            ("internal_diameter", "design_gas_load_factor")
            if diameter is not None
            else ("design_gas_load_factor", "internal_diameter")
        )
        raise InvalidInputError(
            missing,
            f"`{named}` was stated and `{missing}` was not, and the Souders-Brown capacity "
            f"metric needs both - a vessel with a size but no design basis, or the reverse, "
            f"has nothing to check",
        )
    # The two that shape the gas area have no meaning without the vessel they shape.
    if diameter is None:
        for stated, value in (
            ("orientation", orientation),
            ("design_liquid_level_fraction", design_liquid_level_fraction),
        ):
            if value is not None:
                raise InvalidInputError(
                    stated,
                    f"`{stated}` was stated without `internal_diameter` and "
                    f"`design_gas_load_factor`, and the gas area it shapes is not computed "
                    f"until the vessel is",
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
        x, y = list(flash.x), list(flash.y)
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
        x, y = list(moved.x), list(moved.y)
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

    ``liquid_z`` is the outlet's liquid, which is the entrained composition after
    `gas_in_liquid`; its density is taken at the root that composition's own flash puts it
    on, because `initPhysicalProperties` recomputes it after the transfer.
    """
    from azoth.eos.reference.hydrate_inhibitor_wt import GAS
    from azoth.process.reference.capacity import _labelled_phases, _mass_density, _molar_mass

    mixture, labelled, _shares = _labelled_phases(components, t, p, feed_z)
    gas = next((entry for entry in labelled if entry[0] == GAS), None)
    if gas is None:
        # The class's own first answer: a stream with no gas phase uses none of its capacity.
        return CapacityStates(0.0, None, None, None, False, None, None, None, None, None, False)
    _label, gas_z, gas_root = gas
    gas_density = _mass_density(mixture, t, p, gas_z, gas_root)
    gas_molar_mass = _molar_mass(mixture, gas_z)
    if not math.isfinite(gas_density) or gas_density <= 0.0:
        raise InvalidInputError(
            "capacity_utilization",
            f"the gas density is {gas_density} kg/m3, and the Souders-Brown velocity divides "
            f"by its square root",
        )

    # The liquid's own density, or the sizing fallback where the system carries none.
    if any(entry[0] != GAS for entry in labelled):
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
    utilization = volumetric_flow / max_flow
    superficial = volumetric_flow / gas_area

    # **The load factor's own rule, applied after the sizing one.** The floor substitutes
    # `1000.0` for a second phase no denser than the gas; the positivity check is repeated
    # against the substituted value, which is a state the class would also answer NaN for.
    floored = _load_factor_liquid_density(liquid_density, gas_density)
    if not math.isfinite(floored) or floored <= gas_density:
        raise InvalidInputError(
            "gas_load_factor",
            f"the liquid phase is {floored} kg/m3 after the near-dry floor against the gas's "
            f"{gas_density}, and the load factor divides by their difference",
        )
    term = (floored - gas_density) / gas_density
    load_factor = superficial * math.sqrt(1.0 / term)
    if not math.isfinite(load_factor):
        raise InvalidInputError("gas_load_factor", f"the gas load factor is {load_factor} m/s")
    return CapacityStates(
        utilization=utilization,
        max_velocity=max_velocity,
        superficial_velocity=superficial,
        load_factor=load_factor,
        overloaded=utilization > 1.0,
        gas_density=gas_density,
        liquid_density=liquid_density,
        gas_area=gas_area,
        volumetric_flow=volumetric_flow,
        max_flow=max_flow,
        has_liquid=any(entry[0] != GAS for entry in labelled),
    )


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
