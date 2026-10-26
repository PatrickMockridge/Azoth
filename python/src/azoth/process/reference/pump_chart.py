"""``PumpChart`` - a manufacturer's head and efficiency curve, reduced by the affinity laws.

The Python twin of ``crates/azoth-process/src/pump_chart.rs``, written to mirror it step for
step. **The two halves are held to the same capture**: the Rust module's own tests and
``python/tests/test_cross_impl.py`` compare the rows of
``validation/neqsim/captures/process_pump_curve.tsv``.

# One quadratic in ``Q/N``

Every measured flow divides by its speed and every head by the speed squared; one degree-2
least-squares fit through all the reduced points is the head, and the efficiency is fitted and
read at the same reduced flow without a speed reduction. The reduced flow is clamped to the
measured interval in ``getHead`` and ``getEfficiency``; the NPSH curve is evaluated wherever it
is asked and only a negative reduced NPSH is floored.

# The corrections

``getCorrectedHead`` multiplies by ``referenceDensity/actualDensity``.
``getFullyCorrectedHead`` converts the flow through the Hydraulic Institute's flow factor,
evaluates the water curve there and multiplies by the head factor, then applies the density
correction. The factors are :func:`viscosity_factors`, the class's
``calculateViscosityCorrection``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

from azoth.core.errors import InvalidInputError

#: ``ThermodynamicConstantsInterface.gravity``, m/s^2.
GRAVITY: Final = 9.80665
#: The class's ``m3/hr -> gpm`` factor in the Hydraulic Institute correlation.
GPM_PER_CUBIC_METRE_PER_HOUR: Final = 4.40287
#: The class's ``m -> ft`` factor in the same correlation.
FEET_PER_METRE: Final = 3.28084
#: ``getEfficiency``'s lower clamp, as a fraction (the class clamps at 0.1 percent).
MIN_EFFICIENCY: Final = 0.001
#: ``getEfficiency``'s upper clamp, as a fraction (the class clamps at 99 percent).
MAX_EFFICIENCY: Final = 0.99
#: ``getNPSHRequired``'s floor for a negative reduced NPSH, m per rpm squared.
MIN_REDUCED_NPSH: Final = 0.01

#: The class's degree-2 polynomial fit has three coefficients.
_COEFFICIENTS: Final = 3
#: ``getBestEfficiencyFlowRate``'s own scan resolution.
_BEP_INTERVALS: Final = 50


@dataclass(frozen=True, slots=True)
class CurvePoint:
    """One measured point: a flow and the head and efficiency there, at one speed.

    ``efficiency`` is a fraction, not the class's percentage.
    """

    speed: float
    flow: float
    head: float
    efficiency: float
    npsh: float | None = None


@dataclass(frozen=True, slots=True)
class ViscosityFactors:
    """The Hydraulic Institute's three correction factors for a viscous fluid."""

    flow: float
    head: float
    efficiency: float


class PumpChart:
    """The reduced quadratic and the corrections, as ``PumpChart`` holds them."""

    def __init__(
        self,
        points: list[CurvePoint],
        reference_density: float | None = None,
        use_viscosity_correction: bool = False,
    ) -> None:
        """Build the reduced fits from a flattened chart table.

        Raises:
            InvalidInputError: for fewer than three points, a non-finite value, a non-positive
                speed, a partial NPSH curve, or points that share one reduced flow.
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
            if point.npsh is not None and not math.isfinite(point.npsh):
                raise InvalidInputError(
                    "curve", f"point {index}'s npsh is {point.npsh}, not a finite number"
                )
        stated = sum(point.npsh is not None for point in points)
        if stated not in (0, len(points)):
            raise InvalidInputError(
                "curve",
                "the chart states a required NPSH at some points and not at others, and the "
                "class's NPSH curve is one value per flow point",
            )

        reduced_flow = [point.flow / point.speed for point in points]
        self._min_reduced_flow = min(reduced_flow)
        self._max_reduced_flow = max(reduced_flow)
        if self._max_reduced_flow == self._min_reduced_flow:
            raise InvalidInputError(
                "curve",
                f"every point has the same reduced flow {self._min_reduced_flow}, so the "
                f"reduced curve has no interval to fit",
            )

        self._head = _fit(
            reduced_flow, [point.head / (point.speed * point.speed) for point in points]
        )
        self._efficiency = _fit(reduced_flow, [point.efficiency for point in points])
        self._npsh = (
            _fit(
                reduced_flow,
                [
                    (point.npsh if point.npsh is not None else 0.0) / (point.speed * point.speed)
                    for point in points
                ],
            )
            if stated
            else None
        )
        speeds = [point.speed for point in points]
        self._reference_speed = 0.5 * (max(speeds) + min(speeds))
        self._reference_density = (
            reference_density if reference_density is not None and reference_density > 0.0 else None
        )
        self._use_viscosity_correction = use_viscosity_correction

    @property
    def reference_speed(self) -> float:
        """The fitted curve's own midpoint speed, rpm."""
        return self._reference_speed

    @property
    def reduced_flow_range(self) -> tuple[float, float]:
        """The measured reduced-flow interval, m3/hr per rpm."""
        return self._min_reduced_flow, self._max_reduced_flow

    def has_npsh_curve(self) -> bool:
        """Whether the chart carries an NPSH curve."""
        return self._npsh is not None

    def has_density_correction(self) -> bool:
        """Whether a reference density makes the density correction active."""
        return self._reference_density is not None

    def uses_viscosity_correction(self) -> bool:
        """``isUseViscosityCorrection``."""
        return self._use_viscosity_correction

    def head(self, flow: float, speed: float) -> float:
        """``getHead(flow, speed)``, m."""
        reduced = self._clamped_reduced_flow(flow, speed)
        return _evaluate(self._head, reduced) * speed * speed

    def efficiency(self, flow: float, speed: float) -> float:
        """``getEfficiency(flow, speed)``, as a fraction."""
        reduced = self._clamped_reduced_flow(flow, speed)
        return min(MAX_EFFICIENCY, max(MIN_EFFICIENCY, _evaluate(self._efficiency, reduced)))

    def npsh_required(self, flow: float, speed: float) -> float:
        """``getNPSHRequired(flow, speed)``, m; zero where the chart has no NPSH curve."""
        if self._npsh is None:
            return 0.0
        reduced = _evaluate(self._npsh, flow / speed)
        # **Only a negative reduced NPSH is floored**, which is the class's own branch.
        reduced = MIN_REDUCED_NPSH if reduced < 0.0 else reduced
        return reduced * speed * speed

    def best_efficiency_flow_rate(self, speed: float) -> float:
        """``getBestEfficiencyFlowRate(speed)``, m3/hr."""
        if speed <= 0.0:
            return 0.0
        best_flow = 0.0
        best_efficiency = 0.0
        span = self._max_reduced_flow - self._min_reduced_flow
        for index in range(_BEP_INTERVALS + 1):
            reduced = self._min_reduced_flow + span * index / _BEP_INTERVALS
            efficiency = _evaluate(self._efficiency, reduced)
            if efficiency > best_efficiency:
                best_efficiency = efficiency
                best_flow = reduced * speed
        return best_flow

    def corrected_head(self, flow: float, speed: float, actual_density: float) -> float:
        """``getCorrectedHead(flow, speed, actualDensity)``, m."""
        head = self.head(flow, speed)
        if self._reference_density is not None and actual_density > 0.0:
            return head * self._reference_density / actual_density
        return head

    def fully_corrected_head(
        self, flow: float, speed: float, actual_density: float, actual_viscosity: float
    ) -> float:
        """``getFullyCorrectedHead(flow, speed, actualDensity, actualViscosity)``, m."""
        if not (self._use_viscosity_correction and actual_viscosity > 1.0):
            return self.corrected_head(flow, speed, actual_density)
        factors = self._factors_at(actual_viscosity, speed)
        water_head = self.head(flow / factors.flow, speed) * factors.head
        if self._reference_density is not None and actual_density > 0.0:
            return water_head * self._reference_density / actual_density
        return water_head

    def corrected_efficiency(self, flow: float, speed: float, actual_viscosity: float) -> float:
        """``getCorrectedEfficiency(flow, speed, actualViscosity)``, as a fraction."""
        base = self.efficiency(flow, speed)
        if self._use_viscosity_correction and actual_viscosity > 1.0:
            return base * self._factors_at(actual_viscosity, speed).efficiency
        return base

    def _factors_at(self, viscosity: float, speed: float) -> ViscosityFactors:
        """The class's ``calculateViscosityCorrection`` at this chart's best-efficiency point."""
        flow_bep = self.best_efficiency_flow_rate(self._reference_speed)
        head_bep = self.head(flow_bep, self._reference_speed)
        return viscosity_factors(viscosity, flow_bep, head_bep, speed)

    def _clamped_reduced_flow(self, flow: float, speed: float) -> float:
        return min(self._max_reduced_flow, max(self._min_reduced_flow, flow / speed))


def viscosity_factors(
    viscosity: float, flow_bep: float, head_bep: float, speed: float
) -> ViscosityFactors:
    """``calculateViscosityCorrection``: the Hydraulic Institute factors, in its own units."""
    unity = ViscosityFactors(flow=1.0, head=1.0, efficiency=1.0)
    if viscosity <= 1.0:
        return unity
    flow_gpm = flow_bep * GPM_PER_CUBIC_METRE_PER_HOUR
    head_ft = head_bep * FEET_PER_METRE
    if flow_gpm <= 0.0 or head_ft <= 0.0 or speed <= 0.0:
        return unity
    b = 26.6 * math.sqrt(viscosity) * head_ft**0.0625 / (flow_gpm**0.375 * speed**0.25)
    if b <= 1.0:
        return unity
    if b <= 40.0:
        return ViscosityFactors(
            flow=max(0.6, 1.0 - 0.01 * b**0.9),
            head=max(0.7, 1.0 - 0.008 * b),
            efficiency=max(0.4, 1.0 - 0.015 * b**0.85),
        )
    return ViscosityFactors(flow=0.6, head=0.7, efficiency=0.4)


def _fit(reduced_flow: list[float], values: list[float]) -> tuple[float, float, float]:
    """One degree-2 least-squares fit, by the normal equations.

    **The Rust half solves the design matrix by SVD pseudo-inverse and this solves its
    normal equations.** The three columns are ``1, x, x^2`` over a measured reduced-flow
    interval, so the Gram matrix is well conditioned; the two agree to the last few digits,
    which the cross-impl's tolerance measures rather than assumes.
    """
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
