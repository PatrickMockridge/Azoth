"""Unit operations: the process layer.

A **unit operation** takes one or more streams and returns one or more streams,
changing their state by a physical rule: a separator splits a feed into a gas and a
liquid, a valve drops the pressure at constant enthalpy, a mixer blends several feeds
into one. This namespace is the fourth level of the composition every page in this
book describes - above a calculation and a model, below a flowsheet.

A stream here is not a type of its own. It is the arguments a unit operation takes: a
:class:`~azoth.eos.mixture.Mixture`, the :class:`~azoth.eos.IdealGasModel` an enthalpy
is measured from, a temperature, a pressure, a molar flow and a composition. That is
deliberate rather than unfinished - a stream *type* would be a second way to state the
same six things, and the six are what the specs, both implementations and the tests
already agree about. Where a unit has several inlets the state becomes one entry per
inlet, and the compositions become one row per inlet.

**Every unit operation in this namespace is ported from `NeqSim
<https://github.com/equinor/neqsim>`_**, which is developed at NTNU and maintained by
Equinor and is Apache-2.0. See ``NOTICE`` at the repository root for the attribution,
``docs/src/roadmap.md`` for what is in scope and what is deliberately not, and each
function's own module in :mod:`azoth.process.reference` for what its port took from
which line of the original.

The Pareto argument for the scope, and it is read from NeqSim's source rather than
asserted: nearly every unit operation is a flash call plus arithmetic. NeqSim's
``Separator.run`` is 109 lines of a 4,511-line file and does one ``TPflash`` and a phase
split; its ``ThrottlingValve.run`` is 91 lines of 1,894 and does one isenthalpic
``PHflash``. The rest of each file is performance charts, entrainment models, geometry
and mechanical design, and `docs/src/roadmap.md` records that boundary.

As in every other namespace, each function dispatches to the Rust extension when it is
built and to :mod:`azoth.process.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import (
    CompressorResult,
    ExpanderResult,
    HeaterResult,
    MixerResult,
    PumpResult,
    SeparatorResult,
    SplitterResult,
    ThrottlingValveResult,
)
from azoth.core.units import Q
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

__all__ = [
    "compressor",
    "expander",
    "heater",
    "mixer",
    "pump",
    "separator",
    "splitter",
    "throttling_valve",
]

_COMPRESSOR = "process.compressor"
_EXPANDER = "process.expander"
_HEATER = "process.heater"
_MIXER = "process.mixer"
_PUMP = "process.pump"
_SEPARATOR = "process.separator"
_SPLITTER = "process.splitter"
_THROTTLING_VALVE = "process.throttling_valve"


def separator(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> SeparatorResult:
    """One feed split into a gas and a liquid at a single temperature and pressure.

    Both outlets leave at the flash's temperature and pressure, because that is what a
    vessel held at one state does. They differ in how much of the feed they carry and
    what it is made of::

        from azoth import q
        r = azoth.process.separator(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], pressure_drop=q(0, "bar"), heat_duty=q(0, "W"),
        )
        r.gas_flow, r.liquid_flow, r.phase

    **Read ``phase``, not ``beta``.** ``beta`` is ``None`` for a single-phase feed, and
    the number the flash would report for one is the *extrapolated* split - values
    outside ``[0, 1]`` are ordinary.

    **A zero-flow outlet carries the feed's composition.** The gas composition of an
    all-liquid feed is not a physical quantity, and a row of zeros is not a composition
    - it does not sum to one, and the first unit to read it would produce a plausible
    wrong answer. ``gas_flow = 0`` is the field that says the outlet is empty.

    See :func:`azoth.process.reference.separator`.
    """
    return resolve(_SEPARATOR)(  # type: ignore[no-any-return]
        mixture=mixture,
        ideal_gas=ideal_gas,
        T=T,
        P=P,
        n=n,
        z=z,
        pressure_drop=pressure_drop,
        heat_duty=heat_duty,
    )


def mixer(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: list[Q],
    P: list[Q],
    n: list[Q],
    z: list[list[float]],
) -> MixerResult:
    """Several feeds blended into one, at the lowest inlet pressure.

    The outlet temperature is solved for rather than averaged: the blend happens at
    constant enthalpy, and the temperature of a mixture is not the mean of its parts'
    temperatures::

        r = azoth.process.mixer(
            fluid, ig, T=[q(300, "K"), q(350, "K")], P=[q(20, "bar"), q(15, "bar")],
            n=[q(6, "mol/s"), q(4, "mol/s")], z=[[0.6, 0.4], [0.4, 0.6]],
        )
        r.T, r.P, r.flow, r.z_out

    The outlet pressure is the **lowest** inlet pressure - a mixer is a vessel, and
    nothing in it can be above the lowest pressure any feed arrives at.

    **Every inlet shares one component set**, which is why there is one ``mixture``
    argument. NeqSim accumulates component by component and handles a pseudo-fraction
    the others do not have; azoth does not need to.

    See :func:`azoth.process.reference.mixer`.
    """
    return resolve(_MIXER)(  # type: ignore[no-any-return]
        mixture=mixture, ideal_gas=ideal_gas, T=T, P=P, n=n, z=z
    )


def splitter(
    mixture: Mixture,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    fractions: list[float],
) -> SplitterResult:
    """One feed divided into branches at the feed's own temperature and pressure.

    The branches differ only in how much of the feed each carries::

        r = azoth.process.splitter(
            fluid, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], fractions=[0.3, 0.7],
        )
        r.flows   # (3.0, 7.0)

    **The fractions must sum to one**, and that is checked rather than corrected. NeqSim
    normalises its split factors; a splitter's entire output *is* those numbers, so
    silently rescaling them would make a caller's arithmetic error invisible while
    changing every number downstream.

    **There is no ideal-gas argument.** A splitter does no energy balance, so it needs no
    datum - and it is the only unit operation here that takes none.

    See :func:`azoth.process.reference.splitter`.
    """
    return resolve(_SPLITTER)(  # type: ignore[no-any-return]
        mixture=mixture, T=T, P=P, n=n, z=z, fractions=fractions
    )


def throttling_valve(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    z: list[float],
    pressure_drop: Q,
) -> ThrottlingValveResult:
    """A stream's pressure dropped at constant enthalpy.

    The classic Joule-Thomson device, and the outlet temperature is **below** the
    inlet's for a real gas::

        r = azoth.process.throttling_valve(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), z=[0.6, 0.4],
            pressure_drop=q(5, "bar"),
        )
        r.T   # 293.91 K

    **There is no flow rate in this signature**, and that is a statement about the
    physics rather than a saving: an isenthalpic flash is a molar property, so the outlet
    state does not depend on the flow at all.

    Whether the valve is *large enough* for the flow is a different question, and
    ``hydraulics.control_valve_cv`` is where it is asked.

    See :func:`azoth.process.reference.throttling_valve`.
    """
    return resolve(_THROTTLING_VALVE)(  # type: ignore[no-any-return]
        mixture=mixture, ideal_gas=ideal_gas, T=T, P=P, z=z, pressure_drop=pressure_drop
    )


def heater(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> HeaterResult:
    """A duty applied to a stream at a fixed pressure - a heater or a cooler.

    The duty is the specification and the temperature is the answer. **A negative duty is
    a cooler**: NeqSim's ``Cooler`` has no ``run()`` of its own and inherits ``Heater``'s
    steady state, so this is one model rather than two::

        r = azoth.process.heater(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], pressure_drop=q(0, "bar"), heat_duty=q(50, "kW"),
        )
        r.T   # 328.56 K

    **The flow matters.** A duty is an extensive quantity and the flash inverts a molar
    enthalpy, so ``Q / n`` is part of the model: the same 50 kW on 1 mol/s is ten times
    the temperature rise.

    The other ways NeqSim lets a heater be specified - an outlet temperature, a
    ``deltaT`` - are not here, because neither is a unit operation: both reduce to
    ``eos.pt_flash``, where the temperature is already an input.

    See :func:`azoth.process.reference.heater`.
    """
    return resolve(_HEATER)(  # type: ignore[no-any-return]
        mixture=mixture,
        ideal_gas=ideal_gas,
        T=T,
        P=P,
        n=n,
        z=z,
        pressure_drop=pressure_drop,
        heat_duty=heat_duty,
    )


def compressor(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> CompressorResult:
    """A pressure rise at a stated isentropic efficiency.

    The fluid is taken isentropically to the outlet pressure to find the ideal outlet
    enthalpy, and the real one is further from the inlet by exactly the efficiency::

        r = azoth.process.compressor(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], outlet_pressure=q(40, "bar"), efficiency=0.75,
        )
        r.T, r.isentropic_temperature, r.power

    ``isentropic_temperature`` is reported rather than left to be recomputed, because it
    is what the efficiency is measured against and a caller sizing an intercooler needs
    it.

    **The efficiency is required rather than defaulted.** NeqSim defaults it to ``1.0``
    and clamps rather than refusing; a compressor assumed ideal is one that understates
    every duty it is asked for.

    ``process.pump`` and ``process.expander`` run the same procedure with the efficiency
    scaling the other way.

    See :func:`azoth.process.reference.compressor`.
    """
    return resolve(_COMPRESSOR)(  # type: ignore[no-any-return]
        mixture=mixture,
        ideal_gas=ideal_gas,
        T=T,
        P=P,
        n=n,
        z=z,
        outlet_pressure=outlet_pressure,
        efficiency=efficiency,
    )


def pump(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> PumpResult:
    """A pressure rise in a liquid, at a stated isentropic efficiency.

    The same procedure as :func:`compressor`, applied to a stream that is liquid rather
    than gas. What differs is the answer rather than the arithmetic: a liquid is nearly
    incompressible, so the temperature barely moves and the power is not small at all::

        r = azoth.process.pump(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], outlet_pressure=q(30, "bar"), efficiency=0.8,
        )
        r.T, r.power   # 315.60 K, 8.29 kW

    **Nothing checks that the inlet is a liquid.** A pump run on a vapour is
    arithmetically fine and returns a pressure rise nobody can achieve with a pump.

    See :func:`azoth.process.reference.pump`.
    """
    return resolve(_PUMP)(  # type: ignore[no-any-return]
        mixture=mixture,
        ideal_gas=ideal_gas,
        T=T,
        P=P,
        n=n,
        z=z,
        outlet_pressure=outlet_pressure,
        efficiency=efficiency,
    )


def expander(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    P: Q,
    n: Q,
    z: list[float],
    outlet_pressure: Q,
    efficiency: float,
) -> ExpanderResult:
    """A pressure drop that produces work.

    The same procedure again, with the efficiency **multiplying** the ideal enthalpy
    change rather than dividing it - which is the whole difference from
    :func:`compressor`, and it is one line in NeqSim's source::

        r = azoth.process.expander(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], outlet_pressure=q(5, "bar"), efficiency=0.9,
        )
        r.T, r.isentropic_temperature, r.power   # 259.84 K, 258.33 K, -18.91 kW

    ``power`` is **negative** - the fluid is doing the work - on the same convention the
    compressor uses, so the three machines can be summed in a flowsheet without a rule
    about which sign means what.

    ``T`` is always **above** ``isentropic_temperature``: an expander cannot deliver more
    work than the isentropic drop contains. A model that divided instead of multiplying
    would break that, which is what its test asserts.

    See :func:`azoth.process.reference.expander`.
    """
    return resolve(_EXPANDER)(  # type: ignore[no-any-return]
        mixture=mixture,
        ideal_gas=ideal_gas,
        T=T,
        P=P,
        n=n,
        z=z,
        outlet_pressure=outlet_pressure,
        efficiency=efficiency,
    )
