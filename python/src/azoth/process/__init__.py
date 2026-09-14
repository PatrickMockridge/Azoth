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
already agree about.

Every unit operation in this namespace is ported from `NeqSim
<https://github.com/equinor/neqsim>`_, which is developed at NTNU and maintained by
Equinor and is Apache-2.0. See ``NOTICE`` at the repository root for the attribution,
and ``docs/src/roadmap.md`` for what is in scope and what is deliberately not.

As in every other namespace, each function dispatches to the Rust extension when it is
built and to :mod:`azoth.process.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import SeparatorResult
from azoth.core.units import Q
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

__all__ = [
    "separator",
]

_SEPARATOR = "process.separator"


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

    The first unit operation ported from NeqSim. Its physics is one flash and the
    arithmetic that follows from it, which is why it is a *model* with a spec rather
    than a new kind of thing::

        from azoth import q
        r = azoth.process.separator(
            fluid, ig, T=q(300, "K"), P=q(20, "bar"), n=q(10, "mol/s"),
            z=[0.6, 0.4], pressure_drop=q(0, "bar"), heat_duty=q(0, "W"),
        )
        r.gas_flow, r.liquid_flow, r.phase

    Both outlets leave at the flash's temperature and pressure, because that is what a
    vessel held at one state does. They differ in how much of the feed they carry and
    what it is made of.

    **Read ``phase``, not ``beta``.** ``beta`` is ``None`` for a single-phase feed, and
    the number the flash would report for one is the *extrapolated* split - values
    outside ``[0, 1]`` are ordinary. The split here is decided on ``phase``, and so
    should anything built on top of it.

    **A zero-flow outlet carries the feed's composition.** The gas composition of an
    all-liquid feed is not a physical quantity, and a row of zeros is not a composition
    - it does not sum to one, and the first unit to read it would produce a plausible
    wrong answer. ``gas_flow = 0`` is the field that says the outlet is empty.

    ``heat_duty`` is what makes this more than a flash: zero gives an isothermal
    separator at its feed temperature, and a non-zero duty runs an isenthalpic flash
    instead, so the temperature becomes an answer rather than an input.

    Raises:
        InvalidInputError: if ``z`` is not a composition, or if the flash reports a
            two-phase split with no vapour fraction.
        OutOfRangeError: if an input is outside the spec's declared range - a zero
            flow, a negative pressure drop, a non-positive absolute state.
        SolverNotConvergedError: if the flash hits its cap, or a duty puts the answer
            outside the bracket ``eos.ph_flash`` searches.

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
