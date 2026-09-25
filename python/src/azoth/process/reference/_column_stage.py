"""The column's stage, in Python: ``SimpleTray``.

The twin of ``crates/azoth-process/src/column/tray.rs``, and internal to the column for the
same reason: a stage has no palette entry and no model id, because a stage is not a unit
operation a flowsheet connects.

A **stage** is a mixer, a duty and a flash. Its inlets are mixed - the same arithmetic
``kernels::mixer`` runs - and then the mixture flashes either at a stated outlet temperature or
at the enthalpy the mix and the duty imply. Its two outlets are the flash's phases, and an
absent phase is ``None``: NeqSim publishes a zero-flow stream carrying the tray's composition
whose enthalpy is ``-Infinity``, which is not a state a stream can hold.
"""

from __future__ import annotations

from typing import TypedDict

from azoth.core.errors import InvalidInputError
from azoth.core.result import Phase, PhFlashResult, PtFlashResult, ReactiveTpFlashResult
from azoth.core.units import quantity
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash
from azoth.eos.reference.pt_flash import pt_flash
from azoth.reactions.reference import _tables
from azoth.reactions.reference._reactive_flash import MIN_PHASE_FRACTION as _MIN_PHASE_FRACTION
from azoth.reactions.reference.reactive_ph_flash import reactive_ph_flash
from azoth.reactions.reference.reactive_tp_flash import reactive_tp_flash


class StreamRecord(TypedDict):
    """A stream, as the column's internals carry one: the port record plus the names."""

    components: list[str]
    n: float
    z: list[float]
    p: float
    t: float
    h: float


class StageOutcome(TypedDict):
    """A stage's or an end's run: its state and its two phases, either of which may be absent."""

    #: **What the stage had to fall back on, one entry per kind.** A reactive tray's route falls
    #: back the way the class does and says so; every other stage has nothing to report.
    warnings: list[Warning]
    t: float
    p: float
    gas: StreamRecord | None
    liquid: StreamRecord | None


def stage(
    inlets: list[StreamRecord],
    pressure: float,
    out_temperature: float | None,
    heat_input: float,
) -> StageOutcome:
    """One equilibrium stage: mix the inlets, add the duty, flash.

    ``out_temperature`` is the row's pin - the tray the column seeds its ends from flashes at
    that temperature rather than at its enthalpy, which is what ``setCondenserTemperature``
    reaches through the tray's own ``outTemperature``.
    """
    if not inlets:
        raise InvalidInputError("inlets", "a tray needs at least one inlet")

    components = list(inlets[0]["components"])
    n_total = sum(stream["n"] for stream in inlets)
    if n_total <= 0.0:
        raise InvalidInputError(
            "inlets",
            f"a tray's inlets carry {n_total} mol/s in total, so there is no mixture to flash",
        )
    z = [0.0] * len(components)
    h_total = 0.0
    for stream in inlets:
        for i, zi in enumerate(stream["z"]):
            z[i] += stream["n"] * zi
        h_total += stream["n"] * stream["h"]
    z = [value / n_total for value in z]

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    flash: PhFlashResult | PtFlashResult
    if out_temperature is None:
        target = h_total / n_total + heat_input / n_total
        flash = ph_flash(mixture, ideal_gas, quantity(pressure, "Pa"), quantity(target, "J/mol"), z)
        temperature = flash.T.to("K").magnitude
    else:
        flash = pt_flash(mixture, quantity(out_temperature, "K"), quantity(pressure, "Pa"), z)
        temperature = out_temperature
    x, y, phase, beta = list(flash.x), list(flash.y), flash.phase, flash.beta

    gas_share, liquid_share = _phase_fractions(phase, beta)
    # A single phase has the tray's composition, and the flash's trial `x` or `y` is not it.
    single = phase in (Phase.ALL_VAPOUR, Phase.ALL_LIQUID)
    gas_z = z if single else y
    liquid_z = z if single else x

    return {
        "t": temperature,
        "p": pressure,
        "gas": _phase(components, gas_z, n_total * gas_share, pressure, temperature),
        "liquid": _phase(components, liquid_z, n_total * liquid_share, pressure, temperature),
        "warnings": [],
    }


def _phase(
    components: list[str],
    z: list[float],
    n: float,
    pressure: float,
    temperature: float,
) -> StreamRecord | None:
    """One phase as a stream, or ``None`` where the flash found none."""
    if n <= 0.0:
        return None
    return stream_at(components, n, z, temperature, pressure)


def stream_at(
    components: list[str],
    n: float,
    z: list[float],
    temperature: float,
    pressure: float,
) -> StreamRecord:
    """A stream record at a state, with its enthalpy the state's."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    h, _ = enthalpy_at(mixture, ideal_gas, temperature, pressure, list(z))
    return {
        "components": list(components),
        "n": n,
        "z": list(z),
        "p": pressure,
        "t": temperature,
        "h": h,
    }


def _phase_fractions(phase: Phase, beta: float | None) -> tuple[float, float]:
    """The vapour and liquid shares of a flash.

    **The phase decides, and `beta` does not.** NeqSim's outlet getters look for a phase
    *type*: a subcooled stage has no gas phase and reports a vapour of zero flow while its
    ``getBeta()`` still answers one, so reading the split off ``beta`` would give that stage all
    of its flow as vapour.
    """
    if phase == Phase.TWO_PHASE:
        if beta is None:
            raise InvalidInputError("flash", "a two-phase answer with no vapour fraction")
        return float(beta), 1.0 - float(beta)
    if phase == Phase.ALL_VAPOUR:
        return 1.0, 0.0
    if phase == Phase.ALL_LIQUID:
        return 0.0, 1.0
    raise InvalidInputError(
        "flash",
        "the flash is a trivial solution - every K-value straddled one, so nothing here proves "
        "which single phase is present",
    )


def reflux_end(
    inlets: list[StreamRecord],
    pressure: float,
    ratio: float,
    phase: str,
) -> StageOutcome:
    """An end solved at a phase ratio, as `Condenser.run` and `Reboiler.run` call `PVrefluxflash`.

    ``phase`` is ``"vapour"`` for a condenser - the ratio is the liquid it returns over the
    vapour it takes - and ``"liquid"`` for a reboiler, where the same search makes the ratio the
    boilup ``V/B``. The two ratios are reciprocals, so the choice is not cosmetic.

    **The starting temperature is the mixed stream's**, which is where `mixer` puts it - the
    `ph_flash` of the joined enthalpy - and where `prepareMixedStreamForRefluxFlash` leaves it.
    The Rust end reads the same one off its mixed stream rather than choosing a start, so a
    temperature invented here would be a second beginning for one search.
    """
    from azoth.eos.reference.pv_reflux_flash import pv_reflux_flash

    if not inlets:
        raise InvalidInputError("inlets", "an end needs at least one inlet")
    components = list(inlets[0]["components"])
    n_total = sum(stream["n"] for stream in inlets)
    z = [0.0] * len(components)
    h_total = 0.0
    for stream in inlets:
        for i, zi in enumerate(stream["z"]):
            z[i] += stream["n"] * zi
        h_total += stream["n"] * stream["h"]
    z = [value / n_total for value in z]

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    mixed = ph_flash(
        mixture, ideal_gas, quantity(pressure, "Pa"), quantity(h_total / n_total, "J/mol"), z
    )
    flash = pv_reflux_flash(mixture, quantity(pressure, "Pa"), ratio, phase, mixed.T, z)
    temperature = flash.T.to("K").magnitude

    if flash.beta is not None:
        vapour_share = float(flash.beta)
    else:
        # No split to report, so the ratio itself is the share: `1/(1 + R)` for a vapour-ratio
        # search and `R/(1 + R)` for a liquid one.
        vapour_share = 1.0 / (ratio + 1.0) if phase == "vapour" else ratio / (ratio + 1.0)

    return {
        "t": temperature,
        "p": pressure,
        "gas": _phase(components, list(flash.y), n_total * vapour_share, pressure, temperature),
        "liquid": _phase(
            components, list(flash.x), n_total * (1.0 - vapour_share), pressure, temperature
        ),
        "warnings": [],
    }


def reactive_stage(
    inlets: list[StreamRecord],
    pressure: float,
    out_temperature: float | None,
    heat_input: float,
) -> StageOutcome:
    """**A reactive stage: ``SimpleTray.run`` with ``useReactiveFlash`` set.**

    The class's own route, in its own order. A tray whose outlet temperature is stated takes a
    reactive TP flash there; a tray that flashes at its own enthalpy takes a reactive PH flash,
    and - because that answer is a temperature and not a phase list - a reactive TP flash at it,
    which is what NeqSim's own system holds afterwards. A reactive PH flash that fails falls back
    the way the class falls back, and the stage says which step it took.

    **The fluid is the process layer's cubic**, which is why ``reactions.reactive_tp_flash``
    takes one: this column is PR, so a reactive tray that flashed it with SRK would be a
    different machine from the trays beside it.

    The mirror of ``crates/azoth-process/src/column/tray.rs``'s ``reactive_stage``.
    """
    if not inlets:
        raise InvalidInputError("inlets", "a tray needs at least one inlet")
    components = list(inlets[0]["components"])
    n_total = sum(stream["n"] for stream in inlets)
    if n_total <= 0.0:
        raise InvalidInputError(
            "inlets",
            f"a tray's inlets carry {n_total} mol/s in total, so there is no mixture to flash",
        )
    z = [0.0] * len(components)
    h_total = 0.0
    for stream in inlets:
        for i, zi in enumerate(stream["z"]):
            z[i] += stream["n"] * zi
        h_total += stream["n"] * stream["h"]
    z = [value / n_total for value in z]
    moles = [n_total * zi for zi in z]
    warnings: list[Warning] = []

    if out_temperature is None:
        # ``calcMixStreamEnthalpy``: the inlets' total enthalpy plus the duty, and then the
        # *formation inventory* the reactive specification is stated against.
        target = h_total + heat_input + _formation_inventory(components, moles)
        try:
            ph = reactive_ph_flash(
                components,
                quantity(inlets[0]["t"], "K"),
                quantity(pressure, "Pa"),
                [quantity(value, "mol") for value in moles],
                quantity(target, "J"),
                2.0,
                "pr",
            )
            temperature = ph.temperature.to("K").magnitude
        except Exception as error:  # the class's own fallback, and it says so
            temperature = inlets[0]["t"]
            warnings.append(
                Warning(
                    code=WarningCode.SOLVER_NOT_CONVERGED,
                    message=(
                        f"the reactive PH flash did not answer ({error}), so the tray flashed "
                        f"reactively at {temperature} K instead"
                    ),
                )
            )
        outcome = reactive_tp_flash(
            components,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            [quantity(value, "mol") for value in moles],
            2.0,
            "pr",
        )
    else:
        temperature = out_temperature
        outcome = reactive_tp_flash(
            components,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            [quantity(value, "mol") for value in moles],
            2.0,
            "pr",
        )

    gas, liquid = _reactive_outlets(components, pressure, temperature, outcome, warnings)
    return {"t": temperature, "p": pressure, "gas": gas, "liquid": liquid, "warnings": warnings}


def _formation_inventory(components: list[str], moles: list[float]) -> float:
    """``sum_i n_i dHf_i``: the formation inventory the reactive specification carries."""
    total = 0.0
    for name, amount in zip(components, moles, strict=True):
        properties = _tables.formation_properties(name)
        if properties is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no formation row, so its thermochemical enthalpy is unknown",
            )
        total += amount * properties.enthalpy_of_formation
    return total


def _reactive_outlets(
    components: list[str],
    pressure: float,
    temperature: float,
    outcome: ReactiveTpFlashResult,
    warnings: list[Warning],
) -> tuple[StreamRecord | None, StreamRecord | None]:
    """The two outlets a reactive answer holds, or a refusal where it holds no way to tell them
    apart - the mirror of the Rust stage's own rule.

    **A phase the class calls negligible is no outlet**: on a fluid that turns out to be
    single-phase, the delegation's answer carries a row of `1e-15` mol beside the whole feed's,
    because the class's constructor leaves two phase objects each holding it.

    **One phase, and this library's own flash says which it is** - the delegation's indices are
    its bookkeeping rather than a type. **Two, and the label decides**, which is only asked for
    where both are real.
    """
    rows = [[float(value.to("mol").magnitude) for value in row] for row in outcome.phase_moles]
    totals = [sum(row) for row in rows]
    grand = sum(totals)
    held = [
        (index, total)
        for index, total in enumerate(totals)
        if grand > 0.0 and total / grand >= _MIN_PHASE_FRACTION
    ]
    if not held:
        raise InvalidInputError(
            "flash", "the reactive answer holds no phase above the class's own negligible fraction"
        )
    if len(held) == 1:
        index, amount = held[0]
        z = [value / amount for value in rows[index]]
        mixture, _ = _components.mixture_of(components, eos="pr")
        flash = pt_flash(mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), z)
        stream = _phase(components, z, amount, pressure, temperature)
        if flash.phase == Phase.ALL_VAPOUR:
            return stream, None
        if flash.phase == Phase.ALL_LIQUID:
            return None, stream
        raise InvalidInputError(
            "flash",
            "a reactive answer of one phase that this flash does not hold as one, so nothing "
            "here says which outlet it leaves by",
        )
    if len(held) != 2:
        raise InvalidInputError("flash", "the reactive answer holds more phases than two")
    labels = [outcome.phase_type[index] for index, _ in held]
    if labels == ["vapour", "liquid"]:
        (gas_index, gas_moles), (liquid_index, liquid_moles) = held
    elif labels == ["liquid", "vapour"]:
        (liquid_index, liquid_moles), (gas_index, gas_moles) = held
    else:
        warnings.append(
            Warning(
                code=WarningCode.SOLVER_NOT_CONVERGED,
                message=(
                    "the reactive answer's two phases carry no vapour/liquid label, so this "
                    "stage cannot form two outlets from it"
                ),
            )
        )
        raise InvalidInputError(
            "flash",
            "the reactive flash returned two phases without a vapour/liquid label: where its "
            "two phases converge to the same composition, naming one of them the vapour would "
            "be picking a row rather than measuring one",
        )
    gas_row = rows[gas_index]
    liquid_row = rows[liquid_index]
    gas_z = [value / gas_moles for value in gas_row]
    liquid_z = [value / liquid_moles for value in liquid_row]
    return (
        _phase(components, gas_z, gas_moles, pressure, temperature),
        _phase(components, liquid_z, liquid_moles, pressure, temperature),
    )
