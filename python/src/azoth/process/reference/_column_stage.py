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
from azoth.core.result import Phase, PhFlashResult, PtFlashResult
from azoth.core.units import quantity
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at, ph_flash
from azoth.eos.reference.pt_flash import pt_flash


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
    }
