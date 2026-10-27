"""``ExpanderChartKhader`` - a turbo-expander map over the velocity ratio and IGV opening.

The Python twin of ``crates/azoth-process/src/expander_chart.rs``, written to mirror it.

# The abscissa is ``U/C``, not flow

Each measured IGV opening is one curve of isentropic efficiency and stage head against the
blade-speed/spouting-velocity ratio. The chart interpolates **bilinearly** - linearly along
each curve, then linearly between the two bracketing openings - and the U/C direction is
**linearly extrapolated** beyond the measured ends where the IGV direction is **clamped** to
the first or last curve.

# The head is a coefficient over the reference sound speed

``setCurves`` divides every stage head by ``cs_ref²``; a read multiplies the stored
coefficient by the *local* sound speed squared. Azoth has no real-fluid sound-speed model, so
both speeds are stated inputs.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Final

from azoth.core.errors import InvalidInputError

#: ``getOptimumVelocityRatio``'s default when no map is defined.
DEFAULT_OPTIMUM_UC: Final = 0.7


@dataclass(frozen=True, slots=True)
class CurvePoint:
    """One measured point: an IGV opening, a velocity ratio, and the efficiency and head.

    ``head`` is a specific energy in J/kg.
    """

    igv: float
    uc: float
    efficiency: float
    head: float


@dataclass
class _Curve:
    uc: list[float]
    efficiency: list[float]
    head_coefficient: list[float]


class ExpanderChart:
    """The bilinear map and the reference sound speed."""

    def __init__(self, points: list[CurvePoint], reference_sound_speed: float) -> None:
        """Build the map from a flattened point list and the reference sound speed.

        Raises:
            InvalidInputError: for no points, a non-finite value, a non-positive reference
                sound speed, or IGV openings that do not strictly increase in order.
        """
        if not points:
            raise InvalidInputError("chart", "an expander map needs at least one measured point")
        if not math.isfinite(reference_sound_speed) or reference_sound_speed <= 0.0:
            raise InvalidInputError(
                "reference_sound_speed",
                f"the map's head coefficient divides by the reference sound speed squared, "
                f"and {reference_sound_speed} m/s is not positive",
            )
        for index, point in enumerate(points):
            for name, value in (
                ("igv", point.igv),
                ("uc", point.uc),
                ("efficiency", point.efficiency),
                ("head", point.head),
            ):
                if not math.isfinite(value):
                    raise InvalidInputError(
                        "chart", f"point {index}'s {name} is {value}, not a finite number"
                    )

        igv_positions: list[float] = []
        groups: list[list[CurvePoint]] = []
        for point in points:
            if igv_positions and igv_positions[-1] == point.igv:
                groups[-1].append(point)
            else:
                if igv_positions and point.igv <= igv_positions[-1]:
                    raise InvalidInputError(
                        "chart",
                        f"the IGV openings are not strictly increasing: {point.igv} follows "
                        f"{igv_positions[-1]}",
                    )
                igv_positions.append(point.igv)
                groups.append([point])

        curves: list[_Curve] = []
        for group in groups:
            rows = sorted(group, key=lambda row: row.uc)
            curves.append(
                _Curve(
                    uc=[row.uc for row in rows],
                    efficiency=[row.efficiency for row in rows],
                    head_coefficient=[
                        row.head / (reference_sound_speed * reference_sound_speed) for row in rows
                    ],
                )
            )

        self._reference_sound_speed = reference_sound_speed
        self._igv_positions = igv_positions
        self._curves = curves

    @property
    def reference_sound_speed(self) -> float:
        """The reference sound speed the map was normalised with, m/s."""
        return self._reference_sound_speed

    def efficiency(self, uc: float, igv: float) -> float:
        """``getEfficiency(uc, igv)``, as a fraction."""
        return self._interpolate(uc, igv, lambda curve: curve.efficiency)

    def stage_head_drop(self, uc: float, igv: float, local_sound_speed: float) -> float:
        """``getStageHeadDrop(uc, igv, localCs)``, J/kg."""
        coefficient = self._interpolate(uc, igv, lambda curve: curve.head_coefficient)
        return coefficient * local_sound_speed * local_sound_speed

    def optimum_velocity_ratio(self, igv: float) -> float:
        """``getOptimumVelocityRatio(igv)``: the U/C at peak efficiency on the nearest curve."""
        if not self._curves:
            return DEFAULT_OPTIMUM_UC
        nearest = 0
        difference = abs(self._igv_positions[0] - igv)
        for index, position in enumerate(self._igv_positions[1:], start=1):
            candidate = abs(position - igv)
            if candidate < difference:
                difference = candidate
                nearest = index
        curve = self._curves[nearest]
        best = 0
        for index in range(1, len(curve.efficiency)):
            if curve.efficiency[index] > curve.efficiency[best]:
                best = index
        return curve.uc[best]

    def _interpolate(self, uc: float, igv: float, grid: Callable[[_Curve], list[float]]) -> float:
        if not self._curves:
            return -3.56 * (uc - 1.0) * (uc - 1.0) + 1.0
        if len(self._curves) == 1 or igv <= self._igv_positions[0]:
            return _interpolate_curve(uc, self._curves[0].uc, grid(self._curves[0]))
        last = len(self._curves) - 1
        if igv >= self._igv_positions[last]:
            return _interpolate_curve(uc, self._curves[last].uc, grid(self._curves[last]))
        lower = 0
        for index in range(last):
            if self._igv_positions[index] <= igv <= self._igv_positions[index + 1]:
                lower = index
                break
        v0 = _interpolate_curve(uc, self._curves[lower].uc, grid(self._curves[lower]))
        v1 = _interpolate_curve(uc, self._curves[lower + 1].uc, grid(self._curves[lower + 1]))
        t = (igv - self._igv_positions[lower]) / (
            self._igv_positions[lower + 1] - self._igv_positions[lower]
        )
        return v0 + t * (v1 - v0)


def _interpolate_curve(uc: float, x: list[float], y: list[float]) -> float:
    """``interpolateCurve``: piecewise linear, **extrapolating** beyond the U/C ends."""
    n = len(x)
    if n == 1:
        return y[0]
    if uc <= x[0]:
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return y[0] + slope * (uc - x[0])
    if uc >= x[n - 1]:
        slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
        return y[n - 1] + slope * (uc - x[n - 1])
    for index in range(n - 1):
        if x[index] <= uc <= x[index + 1]:
            t = (uc - x[index]) / (x[index + 1] - x[index])
            return y[index] + t * (y[index + 1] - y[index])
    return y[n - 1]
