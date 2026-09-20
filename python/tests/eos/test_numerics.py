"""Numerical safety: the rule in `docs/src/calculus/numerics.md`, applied.

The first application, and the one that found the defect: the Wagner vapour-pressure form
is `x = 1 - T/Tc` raised to `1.5`, and above the critical temperature `x` is negative.
Rust's `powf` gives `NaN` there and Python's `**` promotes to a *complex* - so the two
kernels disagreed about what kind of thing the function is, and neither refusal was
reported. `the_conventions_disagree` in `lean/Azoth/Pow.lean` is the general statement;
this is the call site.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import OutOfRangeError
from azoth.core.units import quantity
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure

WAGNER = (-7.0, 1.5, -2.0, 0.5, 0.0)
TC = 190.56
PC = 4.5992e6


def test_the_wagner_form_is_refused_above_the_critical_temperature() -> None:
    """A refusal, not a clamp: there is no saturation pressure above `Tc` to report."""
    below = antoine_vapor_pressure(
        *WAGNER, "wagner", quantity(TC, "K"), quantity(PC, "Pa"), quantity(150.0, "K")
    )
    assert below.p_sat.magnitude > 0.0

    for t in (190.57, 200.0, 300.0):
        with pytest.raises(OutOfRangeError) as excinfo:
            antoine_vapor_pressure(
                *WAGNER, "wagner", quantity(TC, "K"), quantity(PC, "Pa"), quantity(t, "K")
            )
        assert "critical temperature" in str(excinfo.value)


def test_the_wagner_form_is_exact_at_the_critical_point() -> None:
    """`x = 0` is the boundary, and it is *inside* the domain rather than outside it.

    The numerator vanishes and the exponential is one, so the form returns `Pc` - which is
    the right answer at the critical point and the reason the guard is `x < 0` and not
    `x <= 0`. A guard one comparison too strict would refuse a valid state.
    """
    at_tc = antoine_vapor_pressure(
        *WAGNER, "wagner", quantity(TC, "K"), quantity(PC, "Pa"), quantity(TC, "K")
    )
    assert at_tc.p_sat.magnitude == pytest.approx(PC)
