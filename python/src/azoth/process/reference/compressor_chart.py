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


class CompressorChart:
    """The reduced fits and the speed lines, as ``CompressorChart`` holds them."""

    def __init__(self, points: list[CurvePoint]) -> None:
        """Build the fits and the speed lines from a flattened map.

        Raises:
            InvalidInputError: for fewer than three points, a non-finite value, or a
                non-positive speed.
        """
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
            else:
                curves.append(_SpeedCurve(point.speed, [point.flow], [point.head]))
        speeds = [point.speed for point in points]
        self._curves = curves
        self._min_speed = min(speeds)
        self._max_speed = max(speeds)
        self._reference_speed = 0.5 * (self._max_speed + self._min_speed)

    @property
    def reference_speed(self) -> float:
        """The midpoint of the measured speed range, rpm."""
        return self._reference_speed

    def polytropic_head(self, flow: float, speed: float) -> float:
        """``getPolytropicHead(flow, speed)``, J/kg. No clamp."""
        return _evaluate(self._head, flow / speed) * speed * speed

    def polytropic_efficiency(self, flow: float, speed: float) -> float:
        """``getPolytropicEfficiency(flow, speed)``, as a fraction. No clamp."""
        return _evaluate(self._efficiency, flow / speed)

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


def _evaluate(coefficients: tuple[float, float, float], x: float) -> float:
    """Horner's rule over the fitted coefficients, ascending."""
    return (coefficients[2] * x + coefficients[1]) * x + coefficients[0]
