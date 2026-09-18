"""The process layer: unit operations and flowsheets.

A thin binding to the Rust process layer (`crates/azoth-process`), not a second
implementation. A :class:`Stream` is built here and passed to the kernels, which
run in Rust; units cross as `pint` quantities on this side and SI magnitudes on
the other.

The schema these run against is `specs/unit_ops/` and `specs/flowsheets/`, and
:func:`validate` holds a flowsheet to it — the same check `azoth check` runs.
"""

from __future__ import annotations

import importlib
import pathlib
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from azoth import _core

from azoth.core.units import Q, from_si, to_si

__all__ = [
    "Stream",
    "heat_exchanger",
    "load_flowsheet",
    "mixer",
    "pump",
    "separator",
    "splitter",
    "throttling_valve",
    "validate",
]


def __getattr__(name: str) -> Any:
    """Load the compiled extension on first use, so `import azoth` needs no build."""
    if name == "_core":
        module = importlib.import_module("azoth._core")
        globals()[name] = module
        return module
    raise AttributeError(name)


class Stream:
    """A material stream: composition, molar flow, pressure, temperature and enthalpy.

    ``components`` is a list of substance names, ``z`` their mole fractions, and
    ``n`` the molar flow in mol/s. ``p``, ``t`` and ``h`` are pint quantities.
    """

    __slots__ = ("_inner",)

    def __init__(self, inner: _core.Stream) -> None:
        self._inner = inner

    @classmethod
    def from_pt(cls, components: list[str], z: list[float], n: float, p: Q, t: Q) -> Stream:
        """A stream at a known pressure and temperature."""
        inner = _core.Stream(
            list(components),
            [float(x) for x in z],
            float(n),
            to_si(p, "Pa", "p"),
            to_si(t, "K", "t"),
        )
        return cls(inner)

    @classmethod
    def from_ph(cls, components: list[str], z: list[float], n: float, p: Q, h: Q) -> Stream:
        """A stream at a known pressure and molar enthalpy."""
        inner = _core.Stream.from_ph(
            list(components),
            [float(x) for x in z],
            float(n),
            to_si(p, "Pa", "p"),
            to_si(h, "J/mol", "h"),
        )
        return cls(inner)

    @property
    def components(self) -> tuple[str, ...]:
        return tuple(self._inner.components)

    @property
    def z(self) -> tuple[float, ...]:
        return tuple(self._inner.z)

    @property
    def n(self) -> float:
        return self._inner.n

    @property
    def p(self) -> Q:
        return from_si(self._inner.p, "Pa")

    @property
    def t(self) -> Q:
        return from_si(self._inner.t, "K")

    @property
    def h(self) -> Q:
        return from_si(self._inner.h, "J/mol")

    def __repr__(self) -> str:
        return (
            f"Stream(components={self.components!r}, z={self.z!r}, n={self.n!r}, "
            f"p={self.p!r}, t={self.t!r}, h={self.h!r})"
        )


def splitter(stream: Stream, fractions: Sequence[float]) -> list[Stream]:
    """Split a stream into several with the same state, scaled by ``fractions``."""
    out = _core.splitter(stream._inner, [float(f) for f in fractions])
    return [Stream(s) for s in out]


def mixer(inlets: Sequence[Stream], outlet_pressure: Q | None = None) -> Stream:
    """Join several inlets into one, conserving molar flow and enthalpy."""
    p = None if outlet_pressure is None else to_si(outlet_pressure, "Pa", "outlet_pressure")
    return Stream(_core.mixer([s._inner for s in inlets], p))


def separator(stream: Stream, temperature: Q) -> tuple[Stream, Stream]:
    """Flash a stream into vapour and liquid outlets at ``temperature``."""
    vapour, liquid = _core.separator(stream._inner, to_si(temperature, "K", "temperature"))
    return Stream(vapour), Stream(liquid)


def throttling_valve(stream: Stream, outlet_pressure: Q) -> Stream:
    """Drop a stream to ``outlet_pressure`` without heat or work."""
    return Stream(
        _core.throttling_valve(stream._inner, to_si(outlet_pressure, "Pa", "outlet_pressure"))
    )


def heat_exchanger(hot: Stream, cold: Stream, duty: Q) -> tuple[Stream, Stream]:
    """Move ``duty`` from the hot stream to the cold stream."""
    hot_out, cold_out = _core.heat_exchanger(hot._inner, cold._inner, to_si(duty, "W", "duty"))
    return Stream(hot_out), Stream(cold_out)


def pump(stream: Stream, outlet_pressure: Q, efficiency: float) -> Stream:
    """Raise a liquid stream to ``outlet_pressure`` at ``efficiency``."""
    return Stream(
        _core.pump(
            stream._inner, to_si(outlet_pressure, "Pa", "outlet_pressure"), float(efficiency)
        )
    )


def validate(flowsheet: str, palette_dir: str = "specs/unit_ops") -> list[str]:
    """Validate a flowsheet's TOML text against a palette directory.

    Returns the diagnostic lines; empty when the flowsheet is clean.
    """
    return _core.validate_flowsheet(flowsheet, palette_dir)


def load_flowsheet(path: str, palette_dir: str = "specs/unit_ops") -> list[str]:
    """Read a flowsheet file and validate it; returns the diagnostic lines."""
    return validate(pathlib.Path(path).read_text(), palette_dir)
