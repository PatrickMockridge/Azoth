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

    Inputs in W/(m*K), m**2, K and m. ``dT`` is a temperature **difference**: a bare
    number here means kelvin of difference and no 273.15 is added.

    The hazard this avoids one layer up is worth naming. The scalar API takes a `pint`
    quantity for `dT`, and a difference and an absolute temperature share a dimension,
    so ``Q(30, "degC")`` is the tempting mistake - it converts to 303.15 K and returns a
    plausible answer ten times too large. The spec marks that input ``interval: true``
    and the scalar boundary refuses the offset unit. A batch call cannot make the
    mistake at all, because it takes no quantities: there is no unit to be wrong about.
    See the note in ``docs/src/batch.md``.

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
