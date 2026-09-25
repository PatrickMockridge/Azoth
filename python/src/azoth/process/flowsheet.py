"""Running a flowsheet, and reading what it reached.

The executor is Rust - :mod:`azoth._core`'s ``run_flowsheet``, which is
``azoth_process::executor`` - and this module is the Python side of it, exactly as
:func:`azoth.process.validate` is the Python side of the checker. There is no second
implementation here and there is no second codec: the result is read from the JSON document
``executor::json`` wrote, and :attr:`FlowsheetResult.document` hands that document back
unchanged rather than writing a Python one.

**A flowsheet is self-contained.** Its ``[[inputs]]`` declare each boundary inlet's fluid
and state — ``components``, ``n``, ``z``, ``P`` and ``T``, with ``h`` calculated from the
other three — so ``run_flowsheet`` takes no feed argument for the ordinary case. The
``feeds`` parameter is an *override*, for a sweep or a form that fills a boundary in
without editing the document.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

from azoth import _core
from azoth.core.units import Q, from_si
from azoth.process.kernels import Stream

__all__ = [
    "FlowsheetResult",
    "FlowsheetStream",
    "Residuals",
    "TearResult",
    "run_flowsheet",
]


def _quantity(field: Mapping[str, Any]) -> Q:
    """Rebuild one scalar from the codec's ``{magnitude_si, unit}``.

    The unit name travels with the number, so this converts through the same vocabulary
    the inputs do rather than assuming a unit the codec happened to choose.
    """
    return from_si(float(field["magnitude_si"]), str(field["unit"]))


@dataclass(frozen=True, slots=True, eq=False)
class FlowsheetStream:
    """A stream the run produced, on the five fields a port declaration names.

    ``p`` and ``t`` are this class's spelling and ``P`` and ``T`` are the declaration's:
    the JSON document carries the declaration's, because a front-end reads it against the
    palette, and Python spells them the way :class:`azoth.process.kernels.Stream` does.

    **The enthalpy is the one the executor carried and not a re-derived one.** A stream
    here is what the upstream unit wrote - a mixer's is a flow-weighted mix of its inlets -
    so this is not the same reading as ``Stream.from_pt`` gives for the same ``(T, P, z)``.
    """

    #: Molar flow, ``mol/s``.
    n: Q
    #: Composition, one entry per component.
    z: tuple[float, ...]
    #: Pressure.
    p: Q
    #: Temperature.
    t: Q
    #: Molar enthalpy.
    h: Q


@dataclass(frozen=True, slots=True, eq=False)
class Residuals:
    """A tear's four residuals, in the units ``Recycle.solved()`` compares them.

    **They are not in the same units**, which is the class's own doing: ``flow`` is an
    absolute ``kg/s`` difference below ``1 kg/s`` and a *percentage* at or above it, so one
    tolerance means two things, while temperature and pressure are always percentages and
    composition is a sum of absolute mole-fraction differences.
    """

    #: ``flowBalanceCheck``, in its own mixed unit.
    flow: float
    #: ``compositionBalanceCheck``.
    composition: float
    #: ``temperatureBalanceCheck``, as a percentage.
    temperature: float
    #: ``pressureBalanceCheck``, as a percentage.
    pressure: float


@dataclass(frozen=True, slots=True, eq=False)
class TearResult:
    """One declared recycle, and what its fixed point cost."""

    #: The recycle's own name.
    stream: str
    #: Passes the tear was evaluated on. **Its own count and not the loop's**: a tear switched
    #: off by the low-flow cutoff stops counting, because the class skips an inactive unit
    #: rather than running it - so a loop the outer loop ran twice can report one pass here.
    iterations: int
    #: Whether the last pass's ``solved()`` held.
    solved: bool
    #: Whether the tear was evaluated at all; ``False`` is ``deactivateOnLowFlow``.
    #:
    #: **Read beside ``solved``**, because a tear that is not active is solved by being
    #: *absent* rather than by closing - its four residuals are declared zero rather than
    #: measured. The shipped ``demo.toml`` is that case: its only product is the vapour, so the
    #: recycled liquid is nil and the cutoff switches the loop off on the first pass.
    active: bool
    #: Its four residuals at the last pass, or ``None`` when no pass measured them.
    residuals: Residuals | None


@dataclass(frozen=True, slots=True, eq=False)
class FlowsheetResult:
    """A run's steady state: every named stream, and every tear's convergence.

    **``document`` is the JSON ``executor::json`` wrote**, held rather than re-encoded, so
    that a front-end reading it and a caller reading the fields below cannot disagree about
    what a quantity looks like on the wire. Writing a second document here would be the
    second codec that module exists to avoid.

    ``iterations`` is the outer loop's pass count, which is **not** any one tear's: a
    flowsheet with no recycle runs one pass, and one whose tear stops at the zero-flow
    floor still takes two, because ``Recycle.solved()`` needs a second observation before it
    will call a state converged.
    """

    #: The flowsheet's id.
    flowsheet: str
    #: Whether every tear solved within the cap.
    converged: bool
    #: Passes the outer loop took.
    iterations: int
    #: Every stream, by the endpoint that produced it.
    streams: Mapping[str, FlowsheetStream]
    #: Every declared tear, in declaration order.
    tears: tuple[TearResult, ...]
    #: The JSON document this result is a reading of.
    document: str


def _read(document: str) -> FlowsheetResult:
    """Build a result from the codec's own document."""
    raw = json.loads(document)
    streams = {
        endpoint: FlowsheetStream(
            n=_quantity(record["n"]),
            z=tuple(float(entry) for entry in record["z"]),
            p=_quantity(record["P"]),
            t=_quantity(record["T"]),
            h=_quantity(record["h"]),
        )
        for endpoint, record in raw["streams"].items()
    }
    tears = tuple(
        TearResult(
            stream=tear["stream"],
            iterations=int(tear["iterations"]),
            solved=bool(tear["solved"]),
            active=bool(tear["active"]),
            residuals=None
            if tear["residuals"] is None
            else Residuals(
                flow=float(tear["residuals"]["flow"]),
                composition=float(tear["residuals"]["composition"]),
                temperature=float(tear["residuals"]["temperature"]),
                pressure=float(tear["residuals"]["pressure"]),
            ),
        )
        for tear in raw["tears"]
    )
    return FlowsheetResult(
        flowsheet=str(raw["flowsheet"]),
        converged=bool(raw["converged"]),
        iterations=int(raw["iterations"]),
        streams=MappingProxyType(streams),
        tears=tears,
        document=document,
    )


def run_flowsheet(
    flowsheet: str,
    feeds: Mapping[str, Stream] | None = None,
    palette_dir: str = "specs/unit_ops",
    execution_order: str = "insertion",
) -> FlowsheetResult:
    """Run a flowsheet's TOML text to its steady state.

    **The document declares its own inputs**, so ``feeds`` is omitted in the ordinary case
    and the run is the document's. When it is given, each entry replaces the value of a
    boundary the document already declares, keyed by the input's ``name`` — which is how a
    sweep or a front-end form varies a feed without editing the file. A name the document
    does not declare is **refused**, because a boundary nothing consumes is a wiring error
    rather than an extra feed.

    The palette is what an instance dispatches through, so it is the same directory
    :func:`azoth.process.validate` is given.

    **``execution_order`` is a flag and not a decision to make lightly.**
    ``ProcessSystem.useGraphBasedExecution`` is ``false``, so ``insertion`` - the order the
    instances are declared in - is the class's default and what this defaults to.
    ``topological`` is a depth-first post-order over the connection graph that *drops* the
    declared recycle edges, which is coherent: a torn stream's value is the previous pass's,
    so its reader is not waiting on this pass's producer. A name that is neither is refused
    with the two that are.

    Raises:
        ValueError: where the palette cannot be read or the flowsheet does not parse -
            the same two failures :func:`azoth.process.validate` raises.
        InvalidInputError: where an input's fluid does not resolve, a supplied feed is not
            one the document declares, a connection names an instance or a port the palette
            does not declare, an inlet of single multiplicity receives more than one stream,
            or a unit refuses its arguments. A cycle that no ``[[recycles]]`` entry declares
            is refused rather than hung on.

    See :func:`azoth.process.validate` for the same document checked without running it,
    including the rule that an input's composition is one mole fraction per substance.
    """
    return _read(
        _core.run_flowsheet(
            flowsheet,
            None if feeds is None else {name: s._inner for name, s in feeds.items()},
            palette_dir,
            execution_order,
        )
    )
