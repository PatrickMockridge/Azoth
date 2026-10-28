"""``CompressorChart``, the default fan-law map.

The Python twin of ``crates/azoth-process/src/compressor_chart.rs``, written to mirror it.

# One quadratic in ``Q/N`` per table

Every measured flow divides by its speed and every head by the speed squared; one degree-2
least-squares fit through all the reduced points is the head, read back as
``head = f(Q/N) N²``. The polytropic efficiency is fitted against ``Q/N`` with no speed
reduction. **Nothing is clamped**: the class's getters extrapolate the quadratic wherever they
are asked, and the compressor clamps the efficiency downstream in ``run``.

# The inverses

``speed_for_head`` is the class's own damped Newton with its bisection fallback; ``flow_for_head``
reproduces the transient solver's residual, which divides the head by the efficiency fraction
rather than multiplying it. Both are kept verbatim because they are what the class does, not
because the algebra is the obvious one.

# The surge and stonewall rows

The per-speed lookups take the measured curve closest to the requested speed and read its
minimum or maximum flow. No boundary spline is built here; those belong to the anti-surge
machinery.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

from azoth.core.errors import InvalidInputError
from azoth.process.reference.surge_curve import (
    StoneWallCurve,
    SurgeCurve,
    _Spline,
)

#: The class's degree-2 polynomial fit has three coefficients.
_COEFFICIENTS: Final = 3
_MAX_NEWTON: Final = 50
_SPEED_TOLERANCE: Final = 1.0e-6
_DAMPING: Final = 0.7
_MAX_STEP_FRACTION: Final = 0.3
_MAX_BISECTION: Final = 50
_MIN_SPEED_INTERVAL: Final = 1.0
_MAX_FLOW_ITERATIONS: Final = 100
_FLOW_TOLERANCE: Final = 1.0e-6


@dataclass(frozen=True, slots=True)
class CurvePoint:
    """One measured point: a flow and the head and polytropic efficiency there.

    ``head`` is a specific energy in J/kg and ``efficiency`` a fraction.
    """

    speed: float
    flow: float
    head: float
    efficiency: float


@dataclass
class _SpeedCurve:
    speed: float
    flow: list[float] = field(default_factory=list)
    head: list[float] = field(default_factory=list)
    efficiency: list[float] = field(default_factory=list)


#: Which map the chart interpolates: the reduced quadratics, or the measured lines read
#: directly (``CompressorChartAlternativeMapLookup`` and its extrapolating subclass).
FAN_LAW = "fan_law"
MAP_LOOKUP = "map_lookup"
MAP_LOOKUP_EXTRAPOLATE = "map_lookup_extrapolate"
KHADER_2015 = "khader_2015"


class CompressorChart:
    """The reduced fits and the speed lines, as ``CompressorChart`` holds them."""

    @staticmethod
    def from_khader(
        points: list[CurvePoint],
        reference_sound_speed: float,
        process_sound_speed: float,
        impeller_diameter: float,
    ) -> CompressorChart:
        """`CompressorChartKhader2015`: the extrapolating lookup in similarity coordinates.

        The map is fitted in the *reference* fluid's phase sound speed and read in the
        *process* fluid's system value; the head's corrected flow divides by ``3600`` and the
        efficiency's does not - the class's own asymmetries.
        """
        for name, value in (
            ("reference_sound_speed", reference_sound_speed),
            ("process_sound_speed", process_sound_speed),
            ("impeller_diameter", impeller_diameter),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise InvalidInputError(
                    name, f"the Khader map divides by {name}, and {value} is not positive"
                )
        chart = CompressorChart(points, KHADER_2015)
        area = impeller_diameter * impeller_diameter
        lines: list[_KhaderLine] = []
        for curve in chart._curves:
            lines.append(
                _KhaderLine(
                    mach=curve.speed / 60.0 * impeller_diameter / reference_sound_speed,
                    head_flow=[flow / 3600.0 / reference_sound_speed / area for flow in curve.flow],
                    head=[
                        head / (reference_sound_speed * reference_sound_speed)
                        for head in curve.head
                    ],
                    efficiency=list(curve.efficiency),
                )
            )
        chart._khader = _KhaderMap(
            reference_sound_speed=reference_sound_speed,
            process_sound_speed=process_sound_speed,
            impeller_diameter=impeller_diameter,
            lines=lines,
        )
        return chart

    def __init__(self, points: list[CurvePoint], kind: str = FAN_LAW) -> None:
        """Build the fits and the speed lines from a flattened map.

        Raises:
            InvalidInputError: for fewer than three points, a non-finite value, a
                non-positive speed, or an unknown ``kind``.
        """
        if kind not in (FAN_LAW, MAP_LOOKUP, MAP_LOOKUP_EXTRAPOLATE, KHADER_2015):
            raise InvalidInputError("curve", f"'{kind}' is not a chart interpolation")
        if len(points) < _COEFFICIENTS:
            raise InvalidInputError(
                "curve",
                f"a degree-2 reduced curve needs at least {_COEFFICIENTS} measured points, "
                f"and this chart has {len(points)}",
            )
        for index, point in enumerate(points):
            for name, value in (
                ("speed", point.speed),
                ("flow", point.flow),
                ("head", point.head),
                ("efficiency", point.efficiency),
            ):
                if not math.isfinite(value):
                    raise InvalidInputError(
                        "curve", f"point {index}'s {name} is {value}, not a finite number"
                    )
            if point.speed <= 0.0:
                raise InvalidInputError(
                    "curve",
                    f"point {index}'s speed is {point.speed} rpm, and the reduced variables "
                    f"divide by it",
                )

        reduced = [point.flow / point.speed for point in points]
        self._head = _fit(reduced, [point.head / (point.speed * point.speed) for point in points])
        self._efficiency = _fit(reduced, [point.efficiency for point in points])

        curves: list[_SpeedCurve] = []
        for point in points:
            if curves and curves[-1].speed == point.speed:
                curves[-1].flow.append(point.flow)
                curves[-1].head.append(point.head)
                curves[-1].efficiency.append(point.efficiency)
            else:
                curves.append(
                    _SpeedCurve(point.speed, [point.flow], [point.head], [point.efficiency])
                )
        speeds = [point.speed for point in points]
        self._curves = curves
        self._kind = kind
        self._khader: _KhaderMap | None = None
        self._surge: SurgeCurve | None = None
        self._stonewall: StoneWallCurve | None = None
        self._min_speed = min(speeds)
        self._max_speed = max(speeds)
        self._reference_speed = 0.5 * (self._max_speed + self._min_speed)

    @property
    def reference_speed(self) -> float:
        """The midpoint of the measured speed range, rpm."""
        return self._reference_speed

    def polytropic_head(self, flow: float, speed: float) -> float:
        """``getPolytropicHead(flow, speed)``, J/kg. No clamp.

        Raises:
            InvalidInputError: a non-extrapolating lookup's flow outside a speed line's knots -
                the class's spline throws there and the lookup does not catch it.
        """
        if self._kind == FAN_LAW:
            return _evaluate(self._head, flow / speed) * speed * speed
        if self._kind == KHADER_2015:
            map_ = self._khader_map()
            area = map_.impeller_diameter * map_.impeller_diameter
            corrected = flow / 3600.0 / map_.process_sound_speed / area
            mach = speed / 60.0 * map_.impeller_diameter / map_.process_sound_speed
            value = _lookup_lines(
                map_.map_lines(),
                corrected,
                mach,
                head=True,
                extrapolate=True,
                mean=False,
            )
            return value * map_.process_sound_speed * map_.process_sound_speed
        return self._lookup(flow, speed, head=True)

    def polytropic_efficiency(self, flow: float, speed: float) -> float:
        """``getPolytropicEfficiency(flow, speed)``, as a fraction. No clamp.

        Raises:
            InvalidInputError: the same refusal as :meth:`polytropic_head`.
        """
        if self._kind == FAN_LAW:
            return _evaluate(self._efficiency, flow / speed)
        if self._kind == KHADER_2015:
            map_ = self._khader_map()
            area = map_.impeller_diameter * map_.impeller_diameter
            corrected = flow / 3600.0 / map_.process_sound_speed / area
            mach = speed / 60.0 * map_.impeller_diameter / map_.process_sound_speed
            return _lookup_lines(
                map_.map_lines(),
                corrected,
                mach,
                head=False,
                extrapolate=True,
                mean=False,
            )
        return self._lookup(flow, speed, head=False)

    @property
    def kind(self) -> str:
        """The interpolation scheme this chart uses."""
        return self._kind

    def _khader_map(self) -> _KhaderMap:
        if self._khader is None:
            raise InvalidInputError("chart_type", "the Khader kind has no map built")
        return self._khader

    def _lookup(self, flow: float, speed: float, *, head: bool) -> float:
        """The alternative lookups: the nearest line(s), splined along ``flow``."""
        lines = [
            _MapLine(
                speed=curve.speed,
                head_flow=curve.flow,
                efficiency_flow=curve.flow,
                head=curve.head,
                efficiency=curve.efficiency,
            )
            for curve in self._curves
        ]
        return _lookup_lines(
            lines,
            flow,
            speed,
            head=head,
            extrapolate=self._kind == MAP_LOOKUP_EXTRAPOLATE,
            mean=self._kind == MAP_LOOKUP,
        )

    def _closest_indices(self, speed: float) -> list[int]:
        """``getClosestRefSpeeds``: the line itself, the two it sits between, or the end line."""
        order = sorted(range(len(self._curves)), key=lambda index: self._curves[index].speed)
        index = 0
        while index < len(order) and self._curves[order[index]].speed < speed:
            index += 1
        if index < len(order) and self._curves[order[index]].speed == speed:
            return [order[index]]
        if index == 0:
            return [order[0]]
        if index == len(order):
            return [order[len(order) - 1]]
        return [order[index - 1], order[index]]

    def speed_for_head(self, flow: float, head: float) -> float:
        """``getSpeedValue(flow, head)``: the speed that reaches a head, rpm."""
        reference_head = self.polytropic_head(flow, self._reference_speed)
        if reference_head > 0.0 and head > 0.0:
            initial = self._reference_speed * math.sqrt(head / reference_head)
        else:
            initial = self._reference_speed

        lower = self._min_speed if self._min_speed > 0.0 else 100.0
        upper = self._max_speed
        speed = min(upper, max(lower, initial))

        old_speed = speed * 1.01
        old_error = self.polytropic_head(flow, old_speed) - head
        for _ in range(_MAX_NEWTON):
            new_head = self.polytropic_head(flow, speed)
            error = new_head - head
            if abs(error) < _SPEED_TOLERANCE:
                return speed
            derivative = (error - old_error) / (speed - old_speed)
            if abs(derivative) < 1.0e-10:
                derivative = 2.0 * new_head / speed
                if abs(derivative) < 1.0e-10:
                    break
            update = _DAMPING * error / derivative
            update = max(-_MAX_STEP_FRACTION * speed, min(_MAX_STEP_FRACTION * speed, update))
            old_speed = speed
            old_error = error
            speed = min(upper, max(lower, speed - update))
            if abs(speed - old_speed) < 1.0e-10:
                break

        low, high = lower, upper
        head_at_low = self.polytropic_head(flow, low)
        head_at_high = self.polytropic_head(flow, high)
        if head < head_at_low and head < head_at_high:
            low *= 0.5
            high = low * 2.0 if head_at_low < head_at_high else high
        elif head > head_at_low and head > head_at_high:
            high *= 2.0

        for _ in range(_MAX_BISECTION):
            mid = 0.5 * (low + high)
            mid_head = self.polytropic_head(flow, mid)
            if abs(mid_head - head) < _SPEED_TOLERANCE:
                return mid
            if mid_head < head:
                low = mid
            else:
                high = mid
            if high - low < _MIN_SPEED_INTERVAL:
                return mid
        return 0.5 * (low + high)

    def speed(self, flow: float, head: float) -> int:
        """``getSpeed(flow, head)``: the rounded speed, rpm."""
        return round(self.speed_for_head(flow, head))

    def flow_for_head(self, head: float, speed: float, guess: float) -> float:
        """``getFlow(head, speed, guessFlow)``: the flow where ``head/eta`` hits the target."""
        flow = guess
        old_flow = flow * 1.1
        old_error = self.polytropic_head(old_flow, speed) - head
        error = 1.0
        iteration = 1
        while True:
            iteration += 1
            new_head = self.polytropic_head(flow, speed) / self.polytropic_efficiency(flow, speed)
            error = new_head - head
            derivative = (error - old_error) / (flow - old_flow)
            if abs(derivative) < 1.0e-10:
                break
            old_flow = flow
            old_error = error
            flow -= error / derivative
            if flow < 0.0:
                flow = guess * 0.1
            if not (abs(error) > _FLOW_TOLERANCE and iteration < _MAX_FLOW_ITERATIONS):
                break
        return max(0.0, flow)

    def surge_flow_at_speed(self, speed: float) -> float:
        """``getSurgeFlowAtSpeed``: the minimum flow of the closest speed line."""
        curve = self._closest(speed)
        return min(curve.flow) if curve is not None else math.nan

    def stonewall_flow_at_speed(self, speed: float) -> float:
        """``getStonewallFlowAtSpeed``: the maximum flow of the closest speed line."""
        curve = self._closest(speed)
        return max(curve.flow) if curve is not None else math.nan

    def surge_head_at_speed(self, speed: float) -> float:
        """``getSurgeHeadAtSpeed``: the head at the minimum-flow point."""
        curve = self._closest(speed)
        if curve is None:
            return math.nan
        index = curve.flow.index(min(curve.flow))
        return curve.head[index]

    def stonewall_head_at_speed(self, speed: float) -> float:
        """``getStoneWallHeadAtSpeed``: the head at the maximum-flow point."""
        curve = self._closest(speed)
        if curve is None:
            return math.nan
        index = curve.flow.index(max(curve.flow))
        return curve.head[index]

    def set_surge_curve(self, curve: SurgeCurve) -> None:
        """``setSurgeCurve``: the explicit surge boundary, in the chart's head unit."""
        self._surge = curve

    def set_stonewall_curve(self, curve: StoneWallCurve) -> None:
        """``setStoneWallCurve``: the explicit stonewall boundary."""
        self._stonewall = curve

    @property
    def surge_curve(self) -> SurgeCurve | None:
        """The surge curve, where one is set."""
        return self._surge

    @property
    def stonewall_curve(self) -> StoneWallCurve | None:
        """The stonewall curve, where one is set."""
        return self._stonewall

    def surge_flow(self, head: float) -> float | None:
        """``getSurgeCurve().getSurgeFlow(head)``, or ``None`` without a curve."""
        if self._surge is None:
            return None
        return self._surge.flow_at_head(head)

    def stonewall_flow(self, head: float) -> float | None:
        """``getStoneWallCurve().getStoneWallFlow(head)``, or ``None`` without a curve."""
        if self._stonewall is None:
            return None
        return self._stonewall.flow_at_head(head)

    def surge_head(self, flow: float) -> float | None:
        """``getSurgeCurve().getSurgeHead(flow)``, or ``None`` without a curve."""
        if self._surge is None:
            return None
        return self._surge.head_at_flow(flow)

    def stonewall_head(self, flow: float) -> float | None:
        """``getStoneWallCurve().getStoneWallHead(flow)``, or ``None`` without a curve."""
        if self._stonewall is None:
            return None
        return self._stonewall.head_at_flow(flow)

    def is_surge(self, head: float, flow: float) -> bool:
        """``isSurge(head, flow)``: false where no surge curve is set."""
        if self._surge is None:
            return False
        return self._surge.is_limit(head, flow)

    def is_stonewall(self, head: float, flow: float) -> bool:
        """``isStoneWall(head, flow)``: false where no stonewall curve is set."""
        if self._stonewall is None:
            return False
        return self._stonewall.is_limit(head, flow)

    def generate_surge_curve(self) -> None:
        """``generateSurgeCurve``: one point per speed row - its lowest flow - deduped."""
        flow, head = self._extreme_points(lowest=True)
        self._surge = SurgeCurve.from_points(flow, head)

    def generate_stonewall_curve(self) -> None:
        """``generateStoneWallCurve``: one point per speed row - its highest flow - deduped."""
        flow, head = self._extreme_points(lowest=False)
        self._stonewall = StoneWallCurve.from_points(flow, head)

    def _extreme_points(self, *, lowest: bool) -> tuple[list[float], list[float]]:
        """The extreme point of every speed row, deduped by flow, first wins."""
        if not self._curves:
            raise InvalidInputError(
                "curve", "the chart has no speed rows, so it has no boundary to generate"
            )
        seen: list[float] = []
        flow: list[float] = []
        head: list[float] = []
        for curve in self._curves:
            target = min(curve.flow) if lowest else max(curve.flow)
            if target in seen:
                continue
            seen.append(target)
            flow.append(target)
            head.append(curve.head[curve.flow.index(target)])
        return flow, head

    def flow_range_status(self, flow: float, speed: float) -> str:
        """``getFlowRangeStatus``: ``IN_RANGE``, an extrapolation name, or ``NO_CHART``."""
        if not self._curves:
            return "NO_CHART"
        low = self.surge_flow_at_speed(speed)
        high = self.stonewall_flow_at_speed(speed)
        if math.isnan(low) or math.isnan(high):
            return "NO_CHART"
        tolerance = 1.0e-6 * max(1.0, high)
        if flow < low - tolerance:
            return "EXTRAPOLATED_LOW_FLOW"
        if flow > high + tolerance:
            return "EXTRAPOLATED_HIGH_FLOW"
        return "IN_RANGE"

    def _closest(self, speed: float) -> _SpeedCurve | None:
        best: _SpeedCurve | None = None
        best_difference = math.inf
        for curve in self._curves:
            difference = abs(curve.speed - speed)
            if difference < best_difference:
                best_difference = difference
                best = curve
        return best


def _fit(reduced_flow: list[float], values: list[float]) -> tuple[float, float, float]:
    """One degree-2 least-squares fit, by the normal equations."""
    gram = [
        [math.fsum(x ** (row + column) for x in reduced_flow) for column in range(_COEFFICIENTS)]
        for row in range(_COEFFICIENTS)
    ]
    rhs = [
        math.fsum(value * x**row for x, value in zip(reduced_flow, values, strict=True))
        for row in range(_COEFFICIENTS)
    ]
    solution = _solve3(gram, rhs)
    return solution[0], solution[1], solution[2]


def _solve3(a: list[list[float]], b: list[float]) -> list[float]:
    """A 3x3 linear solve by Gauss-Jordan elimination with partial pivoting."""
    matrix = [[*row, rhs] for row, rhs in zip(a, b, strict=True)]
    for column in range(_COEFFICIENTS):
        pivot = max(range(column, _COEFFICIENTS), key=lambda row: abs(matrix[row][column]))
        if matrix[pivot][column] == 0.0:
            raise InvalidInputError("curve", "the reduced fit's normal equations are singular")
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(_COEFFICIENTS):
            if row == column:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[column], strict=True)
            ]
    return [matrix[row][_COEFFICIENTS] for row in range(_COEFFICIENTS)]


@dataclass(frozen=True, slots=True)
class _MapLine:
    """One speed line, as the lookup reads it."""

    speed: float
    head_flow: list[float]
    efficiency_flow: list[float]
    head: list[float]
    efficiency: list[float]


@dataclass(frozen=True, slots=True)
class _KhaderLine:
    """One speed line in the similarity coordinates."""

    mach: float
    head_flow: list[float]
    head: list[float]
    efficiency: list[float]


@dataclass(frozen=True, slots=True)
class _KhaderMap:
    """`CompressorChartKhader2015`'s map, in similarity coordinates."""

    reference_sound_speed: float
    process_sound_speed: float
    impeller_diameter: float
    lines: list[_KhaderLine]

    def map_lines(self) -> list[_MapLine]:
        """**Both getters read the head's corrected flow**; the class's separate efficiency
        column is stored and never read on this path."""
        return [
            _MapLine(
                speed=line.mach,
                head_flow=line.head_flow,
                efficiency_flow=line.head_flow,
                head=line.head,
                efficiency=line.efficiency,
            )
            for line in self.lines
        ]


def _lookup_lines(
    lines: list[_MapLine],
    flow: float,
    speed: float,
    *,
    head: bool,
    extrapolate: bool,
    mean: bool,
) -> float:
    """The mean-over-lines rule or the extrapolating blend, over any set of lines."""
    if not lines:
        raise InvalidInputError(
            "curve", "the chart has no speed lines, so a map lookup has nothing to read"
        )
    order = sorted(range(len(lines)), key=lambda index: lines[index].speed)
    index = 0
    while index < len(order) and lines[order[index]].speed < speed:
        index += 1
    if index < len(order) and lines[order[index]].speed == speed:
        closest = [order[index]]
    elif index == 0:
        closest = [order[0]]
    elif index == len(order):
        closest = [order[len(order) - 1]]
    else:
        closest = [order[index - 1], order[index]]

    values: list[tuple[float, float]] = []
    for position in closest:
        line = lines[position]
        knots, ordinates = (
            (line.head_flow, line.head) if head else (line.efficiency_flow, line.efficiency)
        )
        values.append(
            (
                line.speed,
                _line_value(knots, ordinates, flow, extrapolate=extrapolate, speed=line.speed),
            )
        )
    if mean:
        return sum(value for _, value in values) / len(values)
    if len(values) == 1:
        reference, value = values[0]
        return value * speed / reference if head else value
    (first, first_value), (second, second_value) = values
    slope = (second_value - first_value) / (second - first)
    blended = first_value + slope * (speed - first)
    # **The pinned class scales the blended head by ``speed / max(s1, s2)``**, and leaves the
    # efficiency unscaled.
    return blended * speed / max(first, second) if head else blended


def _line_value(
    knots: list[float],
    ordinates: list[float],
    flow: float,
    *,
    extrapolate: bool,
    speed: float,
) -> float:
    """One speed line's head or efficiency: the spline inside, a refusal or a line outside."""
    n = len(knots)
    if n < 2 or any(knots[i + 1] <= knots[i] for i in range(n - 1)):
        raise InvalidInputError(
            "curve", f"the speed line at {speed} rpm is not strictly increasing in flow"
        )
    if flow < knots[0] or flow > knots[n - 1]:
        if not extrapolate:
            raise InvalidInputError(
                "curve",
                f"the map lookup's spline at {speed} rpm covers {knots[0]} to {knots[n - 1]} "
                f"m3/hr, and {flow} is outside it; only the extrapolating chart reads past the "
                f"knots",
            )
        low, high = (0, 1) if flow < knots[0] else (n - 2, n - 1)
        slope = (ordinates[high] - ordinates[low]) / (knots[high] - knots[low])
        return ordinates[low] + slope * (flow - knots[low])
    return _Spline(knots, ordinates).at(flow)


def _evaluate(coefficients: tuple[float, float, float], x: float) -> float:
    """Horner's rule over the fitted coefficients, ascending."""
    return (coefficients[2] * x + coefficients[1]) * x + coefficients[0]
