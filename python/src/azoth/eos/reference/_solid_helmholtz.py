"""Second-order forward automatic differentiation, the Python twin of ``azoth_eos::dual``.

The solid reference equations compute their Helmholtz derivatives to second order with a
dual number carrying ``a, a_T, a_V, a_TT, a_TV, a_VV``. This is that arithmetic in Python,
shared by the argon and para-hydrogen solid references rather than re-transcribed per
equation. The Debye oscillator machinery (64-point Gauss-Legendre quadrature) is here for
the same reason.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Dual:
    """A value and its derivatives up to second order in temperature and molar volume."""

    value: float
    dt: float
    dv: float
    dtt: float
    dtv: float
    dvv: float

    @staticmethod
    def constant(value: float) -> Dual:
        return Dual(value, 0.0, 0.0, 0.0, 0.0, 0.0)

    @staticmethod
    def temperature(value: float) -> Dual:
        return Dual(value, 1.0, 0.0, 0.0, 0.0, 0.0)

    @staticmethod
    def volume(value: float) -> Dual:
        return Dual(value, 0.0, 1.0, 0.0, 0.0, 0.0)

    def apply_unary(self, f: float, first: float, second: float) -> Dual:
        return Dual(
            f,
            first * self.dt,
            first * self.dv,
            second * self.dt * self.dt + first * self.dtt,
            second * self.dt * self.dv + first * self.dtv,
            second * self.dv * self.dv + first * self.dvv,
        )

    def exp(self) -> Dual:
        result = math.exp(self.value)
        return self.apply_unary(result, result, result)

    def powf(self, exponent: float) -> Dual:
        return self.apply_unary(
            math.pow(self.value, exponent),
            exponent * math.pow(self.value, exponent - 1.0),
            exponent * (exponent - 1.0) * math.pow(self.value, exponent - 2.0),
        )

    def reciprocal(self) -> Dual:
        return self.apply_unary(
            1.0 / self.value,
            -1.0 / (self.value * self.value),
            2.0 / (self.value * self.value * self.value),
        )

    def log_one_minus_exp_negative(self) -> Dual:
        value = self.value
        if value > 50.0:
            exponential = math.exp(-value)
            return self.apply_unary(math.log1p(-exponential), exponential, -exponential)
        denominator = math.expm1(value)
        return self.apply_unary(
            math.log(-math.expm1(-value)),
            1.0 / denominator,
            -math.exp(value) / (denominator * denominator),
        )

    def debye_free_energy(self) -> Dual:
        z = self.value
        debye = debye_function3(z)
        if z < 1.0e-3:
            debye_first = -1.0 / 8.0 + z / 30.0 - z**3 / 1260.0 + z**5 / 45360.0
            debye_second = 1.0 / 30.0 - z**2 / 420.0 + z**4 / 9072.0
        else:
            if z > 50.0:
                exponential = math.exp(-z)
                q = exponential
                q_first = -exponential
            else:
                denominator = math.expm1(z)
                exponential = math.exp(z)
                q = 1.0 / denominator
                q_first = -exponential / (denominator * denominator)
            debye_first = q - 3.0 * debye / z
            debye_second = q_first - 3.0 * debye_first / z + 3.0 * debye / (z * z)

        if z > 50.0:
            exponential = math.exp(-z)
            logarithm = math.log1p(-exponential)
            logarithm_first = exponential
            logarithm_second = -exponential
        else:
            denominator = math.expm1(z)
            logarithm = math.log(-math.expm1(-z))
            logarithm_first = 1.0 / denominator
            logarithm_second = -math.exp(z) / (denominator * denominator)

        return self.apply_unary(
            logarithm - debye,
            logarithm_first - debye_first,
            logarithm_second - debye_second,
        )

    def __add__(self, other: Dual | float) -> Dual:
        other = _coerce(other)
        return Dual(
            self.value + other.value,
            self.dt + other.dt,
            self.dv + other.dv,
            self.dtt + other.dtt,
            self.dtv + other.dtv,
            self.dvv + other.dvv,
        )

    def __radd__(self, other: float) -> Dual:
        return _coerce(other).__add__(self)

    def __sub__(self, other: Dual | float) -> Dual:
        other = _coerce(other)
        return Dual(
            self.value - other.value,
            self.dt - other.dt,
            self.dv - other.dv,
            self.dtt - other.dtt,
            self.dtv - other.dtv,
            self.dvv - other.dvv,
        )

    def __rsub__(self, other: float) -> Dual:
        return _coerce(other).__sub__(self)

    def __mul__(self, other: Dual | float) -> Dual:
        other = _coerce(other)
        return Dual(
            self.value * other.value,
            self.dt * other.value + self.value * other.dt,
            self.dv * other.value + self.value * other.dv,
            self.dtt * other.value + 2.0 * self.dt * other.dt + self.value * other.dtt,
            self.dtv * other.value
            + self.dt * other.dv
            + self.dv * other.dt
            + self.value * other.dtv,
            self.dvv * other.value + 2.0 * self.dv * other.dv + self.value * other.dvv,
        )

    def __rmul__(self, other: float) -> Dual:
        return _coerce(other).__mul__(self)

    def __truediv__(self, other: Dual | float) -> Dual:
        if isinstance(other, Dual):
            return self * other.reciprocal()
        return Dual(
            self.value / other,
            self.dt / other,
            self.dv / other,
            self.dtt / other,
            self.dtv / other,
            self.dvv / other,
        )

    def __rtruediv__(self, other: float) -> Dual:
        return _coerce(other).__truediv__(self)


def _coerce(value: Dual | float) -> Dual:
    return value if isinstance(value, Dual) else Dual.constant(value)


def debye_function3(z: float) -> float:
    """The normalized third-order Debye function, ``D3(0) = 1/3``."""
    if z < 1.0e-3:
        return 1.0 / 3.0 - z / 8.0 + z * z / 60.0 - z**4 / 5040.0 + z**6 / 272160.0
    if z > 50.0:
        return math.pi**4 / (15.0 * z**3)
    integral = 0.0
    for node, weight in _gauss_legendre():
        x = 0.5 * z * (node + 1.0)
        integral += weight * x**3 / math.expm1(x)
    return 0.5 * z * integral / z**3


_QUADRATURE_ORDER = 64


def _gauss_legendre() -> list[tuple[float, float]]:
    nodes: list[tuple[float, float]] = [(0.0, 0.0)] * _QUADRATURE_ORDER
    midpoint = (_QUADRATURE_ORDER + 1) // 2
    for i in range(midpoint):
        root = math.cos(math.pi * (i + 0.75) / (_QUADRATURE_ORDER + 0.5))
        while True:
            previous = root
            p0 = 1.0
            p1 = root
            for order in range(2, _QUADRATURE_ORDER + 1):
                polynomial = ((2.0 * order - 1.0) * root * p1 - (order - 1.0) * p0) / order
                p0 = p1
                p1 = polynomial
            derivative = _QUADRATURE_ORDER * (root * p1 - p0) / (root * root - 1.0)
            root = previous - p1 / derivative
            if abs(root - previous) <= 1.0e-15:
                break
        weight = 2.0 / ((1.0 - root * root) * derivative * derivative)
        nodes[i] = (-root, weight)
        nodes[_QUADRATURE_ORDER - 1 - i] = (root, weight)
    return nodes
