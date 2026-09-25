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
import pathlib
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

from azoth import _core
from azoth.core.units import Q, from_si
from azoth.process.kernels import Stream

__all__ = [
    "EMBEDDED_PALETTE",
    "FlowsheetResult",
    "FlowsheetStream",
    "Residuals",
    "Session",
    "TearResult",
    "catalogue",
    "forms",
    "run_flowsheet",
    "tools",
]


def _document(text: str) -> dict[str, Any]:
    """The wire's JSON as a mapping, which is what every call to the session answers with."""
    document: dict[str, Any] = json.loads(text)
    return document


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


#: What a document is opened against when a caller names no palette. ``None`` is the bundle
#: the extension carries, which is what makes a notebook need no path: the same 29 specs are
#: compiled in, and a test in the Rust crate holds the bundle to the tree byte for byte.
EMBEDDED_PALETTE: None = None


def forms(palette_dir: str | None = EMBEDDED_PALETTE) -> list[dict[str, Any]]:
    """The unit-operation palette as one form per entry.

    This is the declaration a unit-op window *is*: each entry's ports and parameters, the
    unit and dimension of each, and the model's own bounds with the rationale a field shows
    on a violation. ``kind`` is what picks a widget — ``quantity``, ``vector``, ``boolean``,
    ``enum``, ``string``, or ``unknown`` for the two entries no model declares.

    The document is the middleware's, not this module's: :func:`azoth.process.catalogue`
    returns what the CLI prints and the browser renders.
    """
    return list(json.loads(catalogue(palette_dir))["unit_ops"])


def tools(palette_dir: str | None = EMBEDDED_PALETTE) -> list[dict[str, Any]]:
    """The agent's tool schema: one tool per command, with the envelope as every call's answer.

    Derived from the command model rather than written beside it, so what an agent may do and
    what an editor may do are the same set.
    """
    return list(json.loads(catalogue(palette_dir, with_tools=True))["tools"])


def catalogue(palette_dir: str | None = EMBEDDED_PALETTE, *, with_tools: bool = False) -> str:
    """The palette as a front-end reads it, as the wire's JSON.

    The raw document, for a caller that would rather parse it itself; :func:`forms` and
    :func:`tools` are the two halves of it.
    """
    return _core.catalogue(palette_dir, with_tools)


class Session:
    """A live flowsheet: the one stateful object the middleware has.

    :func:`run_flowsheet` is one call and an immutable result — right for a script, and not
    enough for a form. A session holds the document, so an edit re-checks it in the same call
    and a run's values belong to the document that produced them.

    **Two states, and ``dirty`` is the difference.** The diagnostics are always current
    because a check is cheap and local; the values are not, because a run is the physics. An
    edit marks them stale rather than leaving them to be read as the current document's, and a
    run that *failed* discards them entirely — they described the document as it was.

    Every call that changes something returns the whole envelope: the document, the graph, the
    diagnostics, the paths and the run. That is the same object a browser is handed, so a
    notebook and an editor read one surface rather than two.
    """

    __slots__ = ("_inner",)

    def __init__(self, document: str, palette_dir: str | None = EMBEDDED_PALETTE) -> None:
        """Open a document from its TOML.

        Raises:
            ValueError: on a document this schema cannot read. A document that reads but does
                not check opens with its diagnostics rather than failing — a form has to be
                able to show a broken flowsheet, because showing it is how a user fixes it.
        """
        self._inner = _core.Session(document, palette_dir)

    @classmethod
    def open(cls, path: str, palette_dir: str | None = EMBEDDED_PALETTE) -> Session:
        """Open a document from a file."""
        return cls(pathlib.Path(path).read_text(encoding="utf-8"), palette_dir)

    def apply(self, command: Mapping[str, Any] | str) -> dict[str, Any]:
        """Apply one command, and return the envelope it left.

        The command is the command model's, as a mapping or as its JSON:
        ``{"command": "remove_instance", "id": "hx1"}``. A command that names something
        absent is not refused — it leaves the document as it was and the check that follows
        is what reports it, so the envelope's diagnostics are the answer.
        """
        text = command if isinstance(command, str) else json.dumps(command)
        return _document(self._inner.apply(text))

    def set_order(self, order: str) -> dict[str, Any]:
        """Set which of the class's two orders the next run takes.

        ``insertion`` is the class's default — ``ProcessSystem.useGraphBasedExecution`` is
        ``false`` — and ``topological`` is the depth-first walk that flag reaches. Setting it does
        not re-run: an order is a property of the next run, and the envelope that comes back
        carries the name, so a caller reads the order from the same object as everything else.
        """
        return _document(self._inner.set_order(order))

    def run(self) -> dict[str, Any]:
        """Run the document, and return the envelope it left.

        A run the document cannot survive is reported in the envelope's ``run_error`` rather
        than raised: a failed run is a state a caller shows, not an exception.
        """
        return _document(self._inner.run())

    def envelope(self) -> dict[str, Any]:
        """Everything as it stands, without running anything."""
        return _document(self._inner.envelope())

    @property
    def document(self) -> str:
        """The document as TOML, which is what :meth:`save` writes."""
        return self._inner.document()

    @property
    def graph(self) -> dict[str, Any]:
        """The connection graph, as the editor's own node/edge document."""
        return _document(self._inner.graph())

    @property
    def diagnostics(self) -> list[dict[str, Any]]:
        """What the checker last said, each with a code, a severity and a target."""
        envelope: dict[str, Any] = self.envelope()
        return list(envelope["diagnostics"])

    @property
    def paths(self) -> list[str]:
        """Every value's path, empty until the document has run."""
        return list(self._inner.paths)

    @property
    def ok(self) -> bool:
        """Whether the document can run: no diagnostic of error severity."""
        return bool(self._inner.ok)

    @property
    def dirty(self) -> bool:
        """Whether the values are older than the document."""
        return bool(self._inner.dirty)

    @property
    def run_error(self) -> str | None:
        """Why the last run failed, where one did."""
        error = self._inner.run_error
        return None if error is None else str(error)

    def value(self, path: str) -> float:
        """One value by path, e.g. ``p1.outlet.P``, in the unit the path's field declares."""
        return float(self._inner.value(path))

    def save(self, path: str) -> None:
        """Write the document to a file, as the TOML a :meth:`open` reads back."""
        pathlib.Path(path).write_text(self.document, encoding="utf-8")
