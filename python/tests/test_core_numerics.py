"""``natural_log``: the logarithm semantics the two implementations share."""

from __future__ import annotations

import math

import pytest

from azoth.core.numerics import natural_log


def test_it_matches_f64_ln_at_and_below_zero() -> None:
    """``math.log`` raises here; Rust's ``f64::ln`` is total on the same arguments."""
    assert natural_log(0.0) == -math.inf
    assert natural_log(-0.0) == -math.inf
    assert math.isnan(natural_log(-1.0))
    assert math.isnan(natural_log(math.nan))


def test_it_is_the_ordinary_logarithm_above_zero() -> None:
    assert natural_log(math.e) == pytest.approx(1.0)
    assert natural_log(math.inf) == math.inf
