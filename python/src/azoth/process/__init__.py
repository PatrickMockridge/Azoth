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
    CoolerResult,
    HeaterResult,
    HeatExchangerResult,
    MixerResult,
    PumpResult,
    SeparatorResult,
    SplitterResult,
    ThrottlingValveResult,
)
from azoth.core.units import Q
from azoth.process.kernels import Stream

__all__ = [
    "Stream",
    "cooler",
    "heat_exchanger",
    "heater",
    "load_flowsheet",
    "mixer",
    "pump",
    "separator",
    "splitter",
    "throttling_valve",
    "validate",
]

_MIXER = "process.mixer"
_HEAT_EXCHANGER = "process.heat_exchanger"
_COOLER = "process.cooler"
_HEATER = "process.heater"
_SEPARATOR = "process.separator"
_THROTTLING_VALVE = "process.throttling_valve"
_PUMP = "process.pump"
_SPLITTER = "process.splitter"


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
