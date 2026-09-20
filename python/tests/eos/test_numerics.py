"""Numerical safety: the rule in `docs/src/calculus/numerics.md`, applied.

The first application, and the one that found the defect: the Wagner vapour-pressure form
is `x = 1 - T/Tc` raised to `1.5`, and above the critical temperature `x` is negative.
Rust's `powf` gives `NaN` there and Python's `**` promotes to a *complex* - so the two
kernels disagreed about what kind of thing the function is, and neither refusal was
reported. `the_conventions_disagree` in `lean/Azoth/Pow.lean` is the general statement;
this is the call site.
"""

from __future__ import annotations

import math

import pytest

from azoth.core.errors import OutOfRangeError
from azoth.core.units import quantity
from azoth.eos import from_names
from azoth.eos.reference.antoine_vapor_pressure import antoine_vapor_pressure
from azoth.eos.reference.wilson_activity_coefficients import wilson_activity_coefficients

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


def test_an_even_root_has_two_real_values() -> None:
    """**The page's table, checked.** `x^(p/q)` is `q`-valued, and for even `q` with a
    positive radicand there are **two** real roots, not one.

    `Azoth.Pow.even_power_has_two_roots` is the fact underneath - an even power is not
    injective, so its inverse is not a function - and `sign(x)^p |x|^(p/q)` is a *selection*
    among them rather than a formula.
    """
    # Two values, and both are roots.
    for root in (2.0, -2.0):
        assert root**2 == 4.0
    assert math.sqrt(4.0) == 2.0, "the single-valued form picks one, which is the choice"

    # **An odd denominator is the easy case**: one real root whatever the sign, so there is
    # nothing for a formula to choose. Python's `**` will not give it - it promotes the
    # negative base, which is the next test - so the real odd root is written out.
    assert math.cbrt(-8.0) == -2.0
    assert pytest.approx(2.0, abs=1e-12) == 32.0 ** (1 / 5)


def test_python_promotes_a_negative_base_and_lean_does_not() -> None:
    """The *fourth* convention, and the reason the two azoth kernels disagreed in kind.

    Python's `**` with a negative base and a fractional exponent returns a **complex** -
    `(1.0 + 1.732j)` for `(-8)^(1/3)`, whose real part is `1.0`. That is Lean's
    `Real.rpow` value to the digit, and neither is the real cube root, which is `-2`. Rust's
    `powf` returns `NaN`.

    The Wagner form's refusal exists because of this: one kernel would have returned `NaN`
    as a result and the other would have raised a `TypeError`, so they disagreed about what
    the function *is*, and no case reached the state to say so.
    """
    value = (-8.0) ** (1 / 3)
    assert isinstance(value, complex)
    assert value == pytest.approx(complex(1.0, 1.7320508075688772), abs=1e-12)
    assert value.real == pytest.approx(1.0), "Lean's Real.rpow value"
    assert math.cbrt(-8.0) == -2.0, "and the real odd root is neither"


def test_a_supercritical_component_propagates_nan_in_both_kernels() -> None:
    """**The second site the rule found, and the same defect class as the Wagner form.**

    `eos.wilson_activity_coefficients`' correlation is a power of `1 - T/Tc`, so a
    component above its critical temperature has a negative base and fractional
    exponents. Its spec says a supercritical component "yields `NaN` exactly as NeqSim
    does, which is not refused here" - and that was true of the Rust kernel and false of
    the Python one, which raised `ValueError: math domain error` from `math.pow`. The two
    disagreed in *kind* about a state neither refuses.

    The policy here is **propagate** rather than refuse, because the spec states NaN as
    the model's own domain ending and NeqSim does the same. So the fix is Python returning
    NaN, not both refusing.
    """
    mixture = from_names(["n-heptane", "n-octane"])
    # n-heptane's Tc is 540.2 K, so 560 K is supercritical for one component.
    below = wilson_activity_coefficients(mixture, quantity(400.0, "K"), [0.5, 0.5])
    assert all(math.isfinite(v) for v in below.ln_gamma)

    above = wilson_activity_coefficients(mixture, quantity(560.0, "K"), [0.5, 0.5])
    assert all(math.isnan(v) for v in above.ln_gamma), (
        "a supercritical component must propagate NaN, which is what NeqSim and the Rust "
        "kernel both do - and what the spec says"
    )
