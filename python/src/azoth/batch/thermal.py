"""Thermal calculations, over arrays.

Same arrangement as :mod:`azoth.batch.hydraulics`: plain numbers in the spec's canonical
units in, SI base magnitudes out, and a loop over the scalar reference rather than a
second implementation of it.
"""

from __future__ import annotations

from array import array
from collections.abc import Sequence
from dataclasses import dataclass

from azoth.batch._core import run, sequence
from azoth.batch._result import BatchResult
from azoth.core.warnings import Warning

__all__ = ["ConductionPlaneWallBatch", "conduction_plane_wall"]

_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ConductionPlaneWallBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.conduction_plane_wall`."""

    #: Heat flow rate per element, in watts. Signed, following the sign of `dT`.
    q: array[float]


def _build(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ConductionPlaneWallBatch:
    return ConductionPlaneWallBatch(warnings=warnings, units=units, q=columns["q"])  # type: ignore[arg-type]


def conduction_plane_wall(
    *,
    k: Sequence[float],
    A: Sequence[float],
    dT: Sequence[float],
    L: Sequence[float],
) -> ConductionPlaneWallBatch:
    """Steady conduction through a plane wall, over arrays.

    Inputs in W/(m*K), m**2, K and m. ``dT`` is a temperature **difference**, not an
    absolute temperature: ``dT=[30.0]`` means a 30 kelvin difference, and no 273.15 is
    added. The scalar API takes a `pint` quantity here, where passing ``30 degC`` would
    mean an absolute 303.15 K - see the note in ``docs/src/batch.md``.

    See :func:`azoth.thermal.conduction_plane_wall` for the calculation itself.
    """
    result: ConductionPlaneWallBatch = run(
        _CONDUCTION_PLANE_WALL,
        {
            "k": sequence(k, "k"),
            "A": sequence(A, "A"),
            "dT": sequence(dT, "dT"),
            "L": sequence(L, "L"),
        },
        _build,
    )
    return result
