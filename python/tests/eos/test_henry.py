"""The Henry reference state, as `azoth.eos.reference._henry` carries it.

The Python mirror of the Rust tests in `crates/azoth-eos/src/henry.rs`, asserting the
same constants - so the two kernels are held to the same numbers rather than to each
other, which is the rule a *shared* function follows rather than a registered model.
"""

from __future__ import annotations

import math

import pytest

from azoth.core.errors import PropertyUnavailableError
from azoth.eos import components
from azoth.eos.reference import _henry as henry


def test_the_co2_correlation_is_neqsims_expression() -> None:
    """NeqSim's own arithmetic, evaluated by hand at one state.

    CO2's set - the only real correlation in the table - with the `* 0.01802 * 100` kept
    as written rather than folded to `1.802`, so the test fails if either factor is
    dropped. `H(298.15) = 51.6545 bar` is what NeqSim computes; the literature value for
    CO2 in the mole-fraction convention is about 0.5 bar, and the difference is NeqSim's,
    not this port's - the task is to reproduce it.
    """
    co2 = components.entry("co2")
    value = henry.coefficient(co2.henry, 298.15)
    assert value == pytest.approx(51.6545, abs=1.0e-3)
    assert henry.coefficient(co2.henry, 323.15) > value
    assert henry.coefficient(co2.henry, 273.15) < value


def test_the_sentinel_rows_cap_rather_than_overflow() -> None:
    """**The cap is load-bearing for most of the table, not an edge case.**

    Forty rows carry `900` as their constant. `exp(900)` overflows, and the cap is what
    turns that into a finite insoluble coefficient - so a caller that skipped it would
    propagate an infinity into a fugacity coefficient rather than refusing. Methane is
    the case, and it is a substance every model here carries.

    `math.exp` raises where Rust's `.exp()` returns infinity, so the two languages reach
    the cap by different routes; this is what says they agree.
    """
    methane = components.entry("methane")
    assert methane.henry.is_fitted(), "the row lists a correlation"
    assert not math.isfinite(henry.coefficient(methane.henry, 298.15))
    assert henry.effective_coefficient(methane, 298.15) == henry.INSOLUBLE_HENRY_COEFFICIENT

    # An ion is capped whatever its row says, which is what made the cap a model
    # statement rather than a numerical guard.
    sodium = components.entry("na+")
    assert sodium.component_type == components.ION
    assert henry.is_capped(sodium, 1.0)


def test_an_unfitted_substance_is_refused() -> None:
    """A substance the table fits nothing for is refused rather than evaluated."""
    ethane = components.entry("ethane")
    assert not ethane.henry.is_fitted()
    with pytest.raises(PropertyUnavailableError):
        henry.effective_coefficient(ethane, 298.15)
    # And the zeros really would have given a number, which is why this is a refusal.
    assert henry.coefficient(ethane.henry, 298.15) == pytest.approx(1.802, abs=1.0e-12)


@pytest.mark.parametrize("temperature", [273.15, 298.15, 323.15, 373.15])
def test_the_temperature_derivative_matches_a_finite_difference(temperature: float) -> None:
    """`getHenryCoefdT` is `H` times the logarithmic derivative, so a dropped factor shows."""
    co2 = components.entry("co2").henry
    step = 1.0e-4
    numeric = (
        henry.coefficient(co2, temperature + step) - henry.coefficient(co2, temperature - step)
    ) / (2.0 * step)
    assert henry.coefficient_dt(co2, temperature) == pytest.approx(numeric, rel=1.0e-6)
