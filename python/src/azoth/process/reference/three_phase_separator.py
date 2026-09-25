"""``process.three_phase_separator`` - the separator's kernel.

Spec: ``specs/models/process/three_phase_separator.toml``

The Python twin of ``crates/azoth-process/src/models/three_phase_separator.rs`` over
``crates/azoth-process/src/kernels/three_phase_separator.rs``, written to mirror them.

# Three phases, and the check `run` turns on for them

``run`` sets ``multiPhaseCheck`` before the flash and clears it after, so the split is the
multiphase one: a feed that separates into gas, oil and aqueous is split three ways, and a
two-phase feed into whichever two it has. The flash is at the feed's temperature, or at the
temperature the enthalpy a ``heat_input`` implies is carried at - the vessel holds the
temperature, as ``Separator.run`` does, so a pressure drop is not a throttling.

The phases are named by ``PhaseEos.init``'s label rule rather than by the order the solve
reports them in: ``V/b > 1.75`` is the gas, and of the two liquids the one whose hydrocarbons
outweigh its aqueous components is the oil.

# What the outlets are, and the two branches

The **vapour** is the flash's gas phase unless ``oil_in_gas`` or ``aqueous_in_gas`` is set,
which are the two fractions that carry material out of it - then `run` re-runs the outlet as
a stream. The **oil and aqueous** outlets are re-run whenever those phases exist, without a
condition. A stream ``run`` is a ``TPflash`` of that phase's composition at its own
temperature and pressure, which is the state the composition settles on.

**And that state is not always one this library's two-phase flash can name.** Measured on
this model's own aqueous outlet - ``0.999999986`` water with ``1e-8`` methane and ``1e-15``
n-butane at 300 K and 20 bar - ``pt_flash`` reports a two-phase split whose vapour
composition sums to ``0.9602200424212375``, which is not a composition, and the enthalpy
path refuses it. ``tp_multiflash`` answers the same question with one liquid phase, which is
what NeqSim reports: its outlet's enthalpy is ``-44702.486`` J/mol against the liquid root's
``-44702.552``. ``_settled`` asks the general flash and takes the state it names.

# The entrainment moves moles and no energy

The six fractions each move a share of the *current* from-phase's moles, component by
component, in ``run``'s order - so a later pair reads what an earlier one left. Both outlets
end at the same temperature, so the balance across the transfer is open by what it moved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ThreePhaseSeparatorResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.hydrate_inhibitor_wt import AQUEOUS, GAS, OIL, phase_label
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.tp_multiflash import tp_multiflash


def three_phase_separator(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_aqueous: float,
    gas_in_oil: float,
    oil_in_aqueous: float,
    oil_in_gas: float,
    aqueous_in_gas: float,
    aqueous_in_oil: float,
    heat_input: Q | None = None,
) -> ThreePhaseSeparatorResult:
    """Flash a stream into vapour, oil and aqueous outlets.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature, which the vessel holds.
        pressure_drop: the pressure drop across the vessel.
        gas_in_aqueous: the fraction of the vapour's moles carried into the aqueous phase.
        gas_in_oil: the fraction of the vapour's moles carried into the oil.
        oil_in_aqueous: the fraction of the oil's moles carried into the aqueous phase.
        oil_in_gas: the fraction of the oil's moles carried into the vapour.
        aqueous_in_gas: the fraction of the aqueous phase's moles carried into the vapour.
        aqueous_in_oil: the fraction of the aqueous phase's moles carried into the oil.
        heat_input: heat added before the flash, or ``None`` for a flash at the feed's
            temperature.

    Returns:
        All three outlets' records.

    Raises:
        InvalidInputError: where the pressure drop takes the outlet below zero or a
            fraction is outside ``[0, 1]``.
        OutOfRangeError: for a state the spec bounds, or an enthalpy no temperature
            in the bracket carries.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = three_phase_separator(
        ...     ["methane", "n-butane", "water"],
        ...     q(1.0, "mol/s"),
        ...     [0.5, 0.3, 0.2],
        ...     q(20.0, "bar"),
        ...     q(300.0, "K"),
        ...     q(0.0, "Pa"),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ... )
        >>> round(r.light_liquid_n.to("mol/s").magnitude, 6)
        0.226705
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "feed_n", feed_n)
    t = input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)

    fractions = {
        "gas_in_aqueous": gas_in_aqueous,
        "gas_in_oil": gas_in_oil,
        "oil_in_aqueous": oil_in_aqueous,
        "oil_in_gas": oil_in_gas,
        "aqueous_in_gas": aqueous_in_gas,
        "aqueous_in_oil": aqueous_in_oil,
    }
    apply_checks(
        checks.on_input,
        {"feed_t": t, "pressure_drop": drop, **fractions}.get,
        warnings,
    )

    duty = None if heat_input is None else input_to_si(spec, "heat_input", heat_input)
    states = _route(components, n, feed_z, p, feed_t, drop, fractions, duty)

    return ThreePhaseSeparatorResult(
        vapour_n=from_si(states.vapour.n, "mol/s"),
        vapour_z=states.vapour.z,
        vapour_p=from_si(states.p_out, "Pa"),
        vapour_t=states.temperature,
        vapour_h=from_si(states.vapour.h, "J/mol"),
        light_liquid_n=from_si(states.oil.n, "mol/s"),
        light_liquid_z=states.oil.z,
        light_liquid_p=from_si(states.p_out, "Pa"),
        light_liquid_t=states.temperature,
        light_liquid_h=from_si(states.oil.h, "J/mol"),
        heavy_liquid_n=from_si(states.aqueous.n, "mol/s"),
        heavy_liquid_z=states.aqueous.z,
        heavy_liquid_p=from_si(states.p_out, "Pa"),
        heavy_liquid_t=states.temperature,
        heavy_liquid_h=from_si(states.aqueous.h, "J/mol"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class Outlet:
    """One outlet, in SI: its molar flow, its composition and its molar enthalpy."""

    n: float
    z: tuple[float, ...]
    h: float


@dataclass(frozen=True, slots=True)
class SepStates:
    """The separator's intermediates, in the order it computes them, in SI."""

    p_out: float
    temperature: Q
    vapour: Outlet
    oil: Outlet
    aqueous: Outlet


def _route(
    components: list[str],
    feed_n: float,
    feed_z: list[float],
    feed_p: float,
    feed_t: Q,
    pressure_drop: float,
    fractions: dict[str, float],
    heat_input: float | None,
) -> SepStates:
    """The separator's arithmetic, in SI, in the order it computes it.

    **One arithmetic, two consumers.** :func:`three_phase_separator` builds the result from
    this, and it is the twin of ``kernels/three_phase_separator.rs`` - so the two cannot
    drift into two answers. The kernel-level refusals are here rather than in the model
    because Rust puts them in the kernel.
    """
    for name, fraction in fractions.items():
        if not 0.0 <= fraction <= 1.0:
            raise InvalidInputError(
                name, f"an entrainment fraction is a fraction, and {fraction} is not in [0, 1]"
            )
    p_out = feed_p - pressure_drop
    if p_out <= 0.0:
        raise InvalidInputError(
            "pressure_drop",
            f"a pressure drop of {pressure_drop} Pa takes a {feed_p} Pa inlet to {p_out} Pa, "
            f"which is not a pressure",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    if heat_input is None:
        # The vessel holds the feed's temperature, so the flash is a TP one.
        temperature = feed_t
        t_si = feed_t.to("K").magnitude
    else:
        # W over (mol/s) is J/mol, for the reason `Heater.run` gives.
        h_in, _ = enthalpy_at(mixture, ideal_gas, feed_t.to("K").magnitude, feed_p, feed_z)
        t_si = _solve_temperature(mixture, ideal_gas, p_out, h_in + heat_input / feed_n, feed_z)
        temperature = from_si(t_si, "K")

    phases = _split(mixture, t_si, p_out, feed_z)
    _carry(phases, fractions)

    reflashed_vapour = fractions["oil_in_gas"] != 0.0 or fractions["aqueous_in_gas"] != 0.0
    # **The outlet is the phase, except where `run` re-runs it - and here that is most of
    # them.** The vapour is the gas phase unless `oil_in_gas` or `aqueous_in_gas` carries
    # material out of it; the oil and the aqueous outlets are re-run whenever those phases
    # exist. A stream `run` is a `TPflash` of that phase's composition at its own
    # temperature and pressure, which is a different state from the phase's own root
    # wherever the composition separates once its parent's other phases are gone.
    return SepStates(
        p_out=p_out,
        temperature=temperature,
        vapour=_outlet(
            feed_n, feed_z, phases.get(GAS), GAS, p_out, t_si, reflashed_vapour, mixture, ideal_gas
        ),
        oil=_outlet(feed_n, feed_z, phases.get(OIL), OIL, p_out, t_si, True, mixture, ideal_gas),
        aqueous=_outlet(
            feed_n, feed_z, phases.get(AQUEOUS), AQUEOUS, p_out, t_si, True, mixture, ideal_gas
        ),
    )


def _split(mixture: Any, t_si: float, p_out: float, feed_z: list[float]) -> dict[str, list[float]]:
    """The phases a multiphase flash at ``(t, p_out)`` reports, by the label each carries.

    Each phase is its component moles **per mole of feed**, which is the basis the
    entrainment transfers move.
    """
    flash = tp_multiflash(mixture, from_si(t_si, "K"), from_si(p_out, "Pa"), feed_z)
    reduced = reduced_parameters(mixture, t_si, p_out)
    phases: dict[str, list[float]] = {}
    for index in range(flash.phase_count):
        label = phase_label(mixture, reduced, list(flash.x[index]), flash.z_factor[index])
        phases[label] = [fraction * flash.beta[index] for fraction in flash.x[index]]
    return phases


def _carry(phases: dict[str, list[float]], fractions: dict[str, float]) -> None:
    """Move each entrainment fraction, component by component, in ``run``'s order.

    ``addPhaseFractionToPhase`` computes ``change_i = moles_i(from) * fraction`` and moves it
    between the two phases, so the transfer follows the from-phase's own composition and
    needs no equilibrium. It returns unchanged unless both phases are present and the
    fraction is above ``1e-30``, which is the guard here.
    """
    for from_phase, to_phase, fraction in (
        (GAS, AQUEOUS, fractions["gas_in_aqueous"]),
        (GAS, OIL, fractions["gas_in_oil"]),
        (OIL, AQUEOUS, fractions["oil_in_aqueous"]),
        (OIL, GAS, fractions["oil_in_gas"]),
        (AQUEOUS, GAS, fractions["aqueous_in_gas"]),
        (AQUEOUS, OIL, fractions["aqueous_in_oil"]),
    ):
        if fraction < 1.0e-30:
            continue
        source = phases.get(from_phase)
        target = phases.get(to_phase)
        if source is None or target is None:
            continue
        moved = [amount * fraction for amount in source]
        for index, change in enumerate(moved):
            source[index] -= change
            target[index] += change


def _outlet(
    feed_n: float,
    feed_z: list[float],
    amounts: list[float] | None,
    label: str,
    p_out: float,
    t_si: float,
    reflashed: bool,
    mixture: Any,
    ideal_gas: Any,
) -> Outlet:
    """One outlet: the phase, or the state its composition settles on where it is re-run.

    An absent phase is a **zero-flow outlet** at the outlet state, carrying the feed's
    composition - see the spec's assumptions for why NeqSim's empty clone is not reproduced.
    """
    if amounts is None:
        h, _ = enthalpy_at(mixture, ideal_gas, t_si, p_out, list(feed_z))
        return Outlet(n=0.0, z=tuple(feed_z), h=float(h))
    beta = sum(amounts)
    composition = [amount / beta for amount in amounts] if beta > 0.0 else list(amounts)
    if reflashed:
        h = _settled(mixture, ideal_gas, t_si, p_out, composition)
    else:
        h = _phase_h(mixture, ideal_gas, t_si, p_out, composition, label != GAS)
    return Outlet(n=feed_n * beta, z=tuple(composition), h=h)


def _phase_h(
    mixture: Any, ideal_gas: Any, t_si: float, p_out: float, composition: list[float], liquid: bool
) -> float:
    """One phase's molar enthalpy at a state, on its own side of the cubic."""
    reduced = reduced_parameters(mixture, t_si, p_out)
    root = phase_state(reduced, mixture.kij, list(composition), liquid=liquid).z
    state = molar_enthalpy_entropy(
        mixture, ideal_gas, from_si(t_si, "K"), from_si(p_out, "Pa"), list(composition), root
    )
    return float(state.h.to_base_units().magnitude)


def _settled(
    mixture: Any, ideal_gas: Any, t_si: float, p_out: float, composition: list[float]
) -> float:
    """The state a composition settles on at ``(p, t)`` - what a stream ``run`` reports.

    Two phases or more is the two-phase flash's answer, which is ``enthalpy_at``. One phase
    is *that* phase, and its root is the answer - the branch the module's docstring explains.
    """
    flash = tp_multiflash(mixture, from_si(t_si, "K"), from_si(p_out, "Pa"), composition)
    if flash.phase_count > 1:
        h, _ = enthalpy_at(mixture, ideal_gas, t_si, p_out, composition)
        return float(h)
    labelled = phase_label(
        mixture, reduced_parameters(mixture, t_si, p_out), composition, flash.z_factor[0]
    )
    return _phase_h(mixture, ideal_gas, t_si, p_out, composition, labelled != GAS)


def _solve_temperature(
    mixture: Any, ideal_gas: Any, p_out: float, target: float, z: list[float]
) -> float:
    """The temperature at which a multiphase flash carries the target molar enthalpy.

    The two-phase ``PH`` flash's own temperature is the seed, and the iteration is a
    bisection on a bracket that widens from there: a mixture's enthalpy is monotone in the
    temperature, so a bracket that straddles the target holds the root.
    """

    def at(t_si: float) -> float:
        flash = tp_multiflash(mixture, from_si(t_si, "K"), from_si(p_out, "Pa"), z)
        total = 0.0
        for index in range(flash.phase_count):
            state = molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                from_si(t_si, "K"),
                from_si(p_out, "Pa"),
                list(flash.x[index]),
                flash.z_factor[index],
            )
            total += flash.beta[index] * float(state.h.to_base_units().magnitude)
        return total

    seed = ph_flash_solve(mixture, ideal_gas, from_si(p_out, "Pa"), from_si(target, "J/mol"), z)
    low, high = seed.T.to("K").magnitude - 25.0, seed.T.to("K").magnitude + 25.0
    h_low, h_high = at(low), at(high)
    for _ in range(40):
        if h_low <= target <= h_high:
            break
        step = 0.5 * (high - low)
        if h_low > target:
            if low - step <= 0.0:
                break
            low, h_low = low - step, at(low - step)
        else:
            high, h_high = high + step, at(high + step)
    if not h_low <= target <= h_high:
        raise InvalidInputError(
            "heat_input", f"no temperature in the bracket carries {target} J/mol"
        )
    for _ in range(200):
        middle = 0.5 * (low + high)
        h_middle = at(middle)
        if abs(h_middle - target) < 1.0e-9:
            return middle
        if h_middle < target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.three_phase_separator")
