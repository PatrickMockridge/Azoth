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

__all__ = [
    "ConductionPlaneWallBatch",
    "LmtdCorrectionFactorBatch",
    "conduction_plane_wall",
    "lmtd_correction_factor",
]

_CONDUCTION_PLANE_WALL = "thermal.conduction_plane_wall"
_LMTD_CORRECTION_FACTOR = "thermal.lmtd_correction_factor"


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
    Inputs are plain numbers in the spec's canonical unit.

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


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class LmtdCorrectionFactorBatch(BatchResult):
    """Result of a batch :func:`azoth.thermal.lmtd_correction_factor`."""

    #: The capacity ratio per element.
    r: array[float]
    #: The thermal effectiveness per element.
    p: array[float]
    #: The correction factor for the stated shell passes.
    ft: array[float]
    #: The correction factor for one shell pass.
    ft_1_shell_pass: array[float]
    #: The correction factor for two shell passes.
    ft_2_shell_passes: array[float]
    #: The first pass count reaching 0.75, or -1.
    required_shell_passes: array[float]


def _build_lmtd(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> LmtdCorrectionFactorBatch:
    return LmtdCorrectionFactorBatch(
        warnings=warnings,
        units=units,
        r=columns["r"],  # type: ignore[arg-type]
        p=columns["p"],  # type: ignore[arg-type]
        ft=columns["ft"],  # type: ignore[arg-type]
        ft_1_shell_pass=columns["ft_1_shell_pass"],  # type: ignore[arg-type]
        ft_2_shell_passes=columns["ft_2_shell_passes"],  # type: ignore[arg-type]
        required_shell_passes=columns["required_shell_passes"],  # type: ignore[arg-type]
    )


def lmtd_correction_factor(
    *,
    hot_inlet_temperature: Sequence[float],
    hot_outlet_temperature: Sequence[float],
    cold_inlet_temperature: Sequence[float],
    cold_outlet_temperature: Sequence[float],
    shell_passes: Sequence[float] | None = None,
) -> LmtdCorrectionFactorBatch:
    """The LMTD correction factor, over arrays. Inputs in K.

    ``shell_passes`` is optional **per batch, not per element**: present for every
    element or for none. Omit it and every element is one shell pass.

    See :func:`azoth.thermal.lmtd_correction_factor` for the calculation itself.
    """
    arrays: dict[str, object] = {
        "hot_inlet_temperature": sequence(hot_inlet_temperature, "hot_inlet_temperature"),
        "hot_outlet_temperature": sequence(hot_outlet_temperature, "hot_outlet_temperature"),
        "cold_inlet_temperature": sequence(cold_inlet_temperature, "cold_inlet_temperature"),
        "cold_outlet_temperature": sequence(cold_outlet_temperature, "cold_outlet_temperature"),
    }
    if shell_passes is not None:
        arrays["shell_passes"] = sequence(shell_passes, "shell_passes")
    result: LmtdCorrectionFactorBatch = run(_LMTD_CORRECTION_FACTOR, arrays, _build_lmtd)
    return result
