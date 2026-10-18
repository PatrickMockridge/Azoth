"""Number semantics the two implementations share."""

from __future__ import annotations

import math


def natural_log(value: float) -> float:
    """``ln`` with ``f64::ln``'s total semantics: ``-inf`` at zero, NaN below it.

    ``math.log`` raises on a non-positive argument; Rust's ``f64::ln`` returns ``-inf`` at
    ``+0`` and NaN at a negative. Every kernel here is written twice, so a composition or a
    phase fraction that reaches an exact zero must take the same path in both languages.
    """
    if value > 0.0:
        return math.log(value)
    if value == 0.0:
        return -math.inf
    return math.nan


def square_root(value: float) -> float:
    """``sqrt`` with ``f64::sqrt``'s semantics: NaN below zero rather than an exception."""
    if value < 0.0:
        return math.nan
    return math.sqrt(value)


def cube_root(value: float) -> float:
    """``cbrt`` with ``f64::cbrt``'s semantics: a negative base has a negative root."""
    return math.copysign(abs(value) ** (1.0 / 3.0), value)


def ceil(value: float) -> float:
    """``f64::ceil``: total, so an infinity stays one and a NaN stays NaN."""
    if not math.isfinite(value):
        return value
    return float(math.ceil(value))


def divide(numerator: float, denominator: float) -> float:
    """``a / b`` with IEEE semantics: a zero divisor gives an infinity or a NaN.

    Java's and Rust's floating division is total; Python's raises. A kernel written
    twice must take the same branch on a state that reaches the divisor.
    """
    if denominator == 0.0:
        if numerator == 0.0 or math.isnan(numerator):
            return math.nan
        sign = math.copysign(1.0, numerator) * math.copysign(1.0, denominator)
        return math.copysign(math.inf, sign)
    return numerator / denominator


def powf(base: float, exponent: float) -> float:
    """``f64::powf``'s semantics where Python's ``**`` differs.

    A negative base with a non-integral exponent is NaN rather than a complex number,
    and zero to a negative power is an infinity rather than an exception.
    """
    if base < 0.0 and not float(exponent).is_integer():
        return math.nan
    if base == 0.0 and exponent < 0.0:
        return math.copysign(math.inf, base) if exponent % 2.0 == 1.0 else math.inf
    return float(base**exponent)
