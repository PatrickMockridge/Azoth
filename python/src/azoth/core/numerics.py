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
