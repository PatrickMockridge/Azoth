"""``SafeSplineSurgeCurve`` / ``SafeSplineStoneWallCurve`` - the compressor map's boundaries.

The Python twin of ``crates/azoth-process/src/surge_curve.rs``, written to mirror it.

# Two splines, one per direction

``setCurve`` sorts the measured ``(flow, head)`` pairs by flow and fits a natural cubic spline
``head = f(flow)``; it then dedupes the pairs by head (ascending, the last duplicate winning, as
the class's ``TreeMap`` does) and fits a second natural spline ``flow = f(head)``. A read inside
the measured range is that spline; outside it is a **linear extrapolation** from the two end
points, clamped at zero, and a single measured point is a constant.

# The two curves differ in one comparison

``SurgeCurve`` calls a point a limit when the curve's flow at its head exceeds the operating
flow (``getFlow(head) > flow``); ``StoneWallCurve`` calls it a limit when the flow is greater
than the curve's.

# A natural spline, not a not-a-knot one

The class fits with Apache Commons Math's ``SplineInterpolator``, whose second derivative is
zero at both ends; ``_Spline`` solves that same tridiagonal system.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError


class _Spline:
    """A natural cubic spline through strictly increasing abscissae."""

    __slots__ = ("second", "x", "y")

    def __init__(self, x: list[float], y: list[float]) -> None:
        n = len(x)
        if n < 2:
            raise InvalidInputError("surge_curve", "a spline needs at least two points")
        for i in range(1, n):
            if not x[i] > x[i - 1]:
                raise InvalidInputError(
                    "surge_curve", f"the spline's abscissae are not strictly increasing at {x[i]}"
                )
        h = [x[i + 1] - x[i] for i in range(n - 1)]
        second = [0.0] * n
        diagonal = [0.0] * n
        upper = [0.0] * n
        rhs = [0.0] * n
        for i in range(1, n - 1):
            diagonal[i] = 2.0 * (h[i - 1] + h[i])
            upper[i] = h[i]
            rhs[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
        for i in range(2, n - 1):
            factor = h[i - 1] / diagonal[i - 1]
            diagonal[i] -= factor * upper[i - 1]
            rhs[i] -= factor * rhs[i - 1]
        for i in range(n - 2, 0, -1):
            below = h[i] * second[i + 1] if i + 1 < n - 1 else 0.0
            second[i] = (rhs[i] - below) / diagonal[i]
        self.x = x
        self.y = y
        self.second = second

    def at(self, value: float) -> float:
        n = len(self.x)
        index = 0
        for i in range(n - 1):
            if self.x[i] <= value <= self.x[i + 1]:
                index = i
                break
        h = self.x[index + 1] - self.x[index]
        a = (self.x[index + 1] - value) / h
        b = (value - self.x[index]) / h
        return (
            a * self.y[index]
            + b * self.y[index + 1]
            + ((a * a * a - a) * self.second[index] + (b * b * b - b) * self.second[index + 1])
            * h
            * h
            / 6.0
        )


class _Boundary:
    """The shared curve: the flow-sorted pairs, the head-deduped pairs, and the two splines."""

    __slots__ = (
        "_flow",
        "_flow_from_head",
        "_head",
        "_head_from_flow",
        "_single",
        "_single_flow",
        "_single_head",
        "_sorted_flow",
        "_sorted_head",
    )

    def __init__(self, flow: list[float], head: list[float]) -> None:
        if len(flow) != len(head):
            raise InvalidInputError(
                "surge_curve", "the curve's flow and head arrays must have the same length"
            )
        if not flow:
            raise InvalidInputError("surge_curve", "an anti-surge curve needs at least one point")
        for index, (flow_value, head_value) in enumerate(zip(flow, head, strict=True)):
            if not math.isfinite(flow_value) or not math.isfinite(head_value):
                raise InvalidInputError(
                    "surge_curve",
                    f"point {index} is not a finite pair, and a curve cannot be fitted",
                )
        if len(flow) == 1:
            self._flow = [flow[0]]
            self._head = [head[0]]
            self._sorted_head = [head[0]]
            self._sorted_flow = [flow[0]]
            self._head_from_flow = None
            self._flow_from_head = None
            self._single = True
            self._single_flow = flow[0]
            self._single_head = head[0]
            return

        pairs = sorted(zip(flow, head, strict=True), key=lambda pair: pair[0])
        if len(pairs) == 2:
            raise InvalidInputError(
                "surge_curve",
                "a two-point boundary cannot be splined - the class's interpolator needs three "
                "points, and one is its constant",
            )
        for first, second in itertools.pairwise(pairs):
            if first[0] == second[0]:
                raise InvalidInputError(
                    "surge_curve",
                    f"two points share the flow {first[0]}, and a curve needs distinct flows",
                )
        self._flow = [pair[0] for pair in pairs]
        self._head = [pair[1] for pair in pairs]
        by_head: dict[float, tuple[float, float]] = {}
        for pair in pairs:
            by_head[pair[1]] = pair
        ordered = [by_head[key] for key in sorted(by_head)]
        self._sorted_head = [pair[1] for pair in ordered]
        self._sorted_flow = [pair[0] for pair in ordered]
        if len(self._sorted_head) < 2:
            raise InvalidInputError(
                "surge_curve", "not enough distinct head values for spline interpolation"
            )
        self._head_from_flow = _Spline(self._flow, self._head)
        self._flow_from_head = _Spline(self._sorted_head, self._sorted_flow)
        self._single = False
        self._single_flow = 0.0
        self._single_head = 0.0

    def flow_at_head(self, head_value: float) -> float:
        if self._single:
            return self._single_flow
        n = len(self._sorted_head)
        if self._sorted_head[0] <= head_value <= self._sorted_head[n - 1]:
            return max(0.0, self._flow_from_head.at(head_value))  # type: ignore[union-attr]
        if head_value < self._sorted_head[0]:
            slope = (self._sorted_flow[1] - self._sorted_flow[0]) / (
                self._sorted_head[1] - self._sorted_head[0]
            )
            base_flow, base_head = self._sorted_flow[0], self._sorted_head[0]
        else:
            slope = (self._sorted_flow[n - 1] - self._sorted_flow[n - 2]) / (
                self._sorted_head[n - 1] - self._sorted_head[n - 2]
            )
            base_flow, base_head = self._sorted_flow[n - 1], self._sorted_head[n - 1]
        return max(0.0, base_flow + slope * (head_value - base_head))

    def head_at_flow(self, flow_value: float) -> float:
        if self._single:
            return self._single_head
        n = len(self._flow)
        if self._flow[0] <= flow_value <= self._flow[n - 1]:
            return max(0.0, self._head_from_flow.at(flow_value))  # type: ignore[union-attr]
        if flow_value < self._flow[0]:
            slope = (self._head[1] - self._head[0]) / (self._flow[1] - self._flow[0])
            base_flow, base_head = self._flow[0], self._head[0]
        else:
            slope = (self._head[n - 1] - self._head[n - 2]) / (
                self._flow[n - 1] - self._flow[n - 2]
            )
            base_flow, base_head = self._flow[n - 1], self._head[n - 1]
        return max(0.0, base_head + slope * (flow_value - base_flow))


@dataclass(frozen=True, slots=True)
class SurgeCurve:
    """``SafeSplineSurgeCurve``: the surge boundary, ``flow(head) > flow`` is a limit."""

    _boundary: _Boundary

    @staticmethod
    def from_points(flow: list[float], head: list[float]) -> SurgeCurve:
        """Fit the curve from measured flow and head arrays, in any order."""
        return SurgeCurve(_Boundary(flow, head))

    def flow_at_head(self, head: float) -> float:
        """``getFlow(head)`` / ``getSurgeFlow(head)``: the curve's flow at a head, m3/h."""
        return self._boundary.flow_at_head(head)

    def head_at_flow(self, flow: float) -> float:
        """``getSurgeHead(flow)``: the curve's head at a flow."""
        return self._boundary.head_at_flow(flow)

    def is_limit(self, head: float, flow: float) -> bool:
        """``isSurge(head, flow)``: the curve's flow at the head exceeds the operating flow."""
        return self.flow_at_head(head) > flow


@dataclass(frozen=True, slots=True)
class StoneWallCurve:
    """``SafeSplineStoneWallCurve``: the stonewall boundary, ``flow(head) < flow`` is a limit."""

    _boundary: _Boundary

    @staticmethod
    def from_points(flow: list[float], head: list[float]) -> StoneWallCurve:
        """Fit the curve from measured flow and head arrays, in any order."""
        return StoneWallCurve(_Boundary(flow, head))

    def flow_at_head(self, head: float) -> float:
        """``getFlow(head)``: the curve's flow at a head, m3/h."""
        return self._boundary.flow_at_head(head)

    def head_at_flow(self, flow: float) -> float:
        """The curve's head at a flow."""
        return self._boundary.head_at_flow(flow)

    def is_limit(self, head: float, flow: float) -> bool:
        """``isStoneWall(head, flow)``: the curve's flow at the head is below the flow."""
        return self.flow_at_head(head) < flow
