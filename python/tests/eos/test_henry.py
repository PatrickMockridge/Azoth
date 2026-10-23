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

    CO2's sourced set, with the `* 0.01802 * 100` kept as written rather than folded to
    `1.802`, so the test fails if either factor is dropped. `H(298.15) = 28.95 bar` is what
    the table's coefficients give after upstream PR #4045 replaced the placeholder
    correlations with sourced data.
    """
    co2 = components.entry("co2")
    value = henry.coefficient(co2.henry, 298.15)
    assert value == pytest.approx(28.95, abs=1.0e-3)
    assert henry.coefficient(co2.henry, 323.15) > value
    assert henry.coefficient(co2.henry, 273.15) < value


def test_the_cap_survives_the_sourced_data() -> None:
    """**The `exp(900)` sentinel is gone and the cap is not.**

    Upstream PR #4045 replaced the four placeholder sets - forty rows of `900`, the
    constant `2.92`, the `0.059565` quartet and the one CO2 set copied onto `co`,
    `cos` - with sourced two-term coefficients, so methane's row now evaluates to a finite
    `723.75 bar` rather than overflowing. The cap remains the model's rule for an ion, and
    it still catches a value that would overflow.

    `math.exp` raises where Rust's `.exp()` returns infinity, so the two languages reach
    the cap by different routes; `is_capped_value` is what says they agree.
    """
    methane = components.entry("methane")
    assert methane.henry.is_fitted(), "the row lists a correlation"
    assert henry.coefficient(methane.henry, 298.15) == pytest.approx(723.75, abs=1.0e-6)

    sodium = components.entry("na+")
    assert sodium.component_type == components.ION
    assert henry.is_capped(sodium, 1.0)
    assert henry.is_capped_value(math.inf, ion=False)


def test_an_unfitted_substance_is_refused() -> None:
    """A substance the table fits nothing for is refused rather than evaluated."""
    mdea = components.entry("mdea")
    assert not mdea.henry.is_fitted()
    with pytest.raises(PropertyUnavailableError):
        henry.effective_coefficient(mdea, 298.15)
    # And the zeros really would have given a number, which is why this is a refusal.
    assert henry.coefficient(mdea.henry, 298.15) == pytest.approx(1.802, abs=1.0e-12)


def test_the_ion_clause_is_asked_before_the_fittedness_check() -> None:
    """`ComponentGE` caps an ion before it asks whether the row was fitted.

    An ion's four zeros would otherwise reach the unfitted refusal, and NeqSim never
    refuses there: the cap decides first, so the ion is `1e12 bar` whatever its row says.
    """
    assert henry.is_capped_value(1.802, ion=True)
    assert not henry.is_capped_value(1.802, ion=False)
    assert henry.is_capped_value(math.inf, ion=False)


@pytest.mark.parametrize("temperature", [273.15, 298.15, 323.15, 373.15])
def test_the_temperature_derivative_matches_a_finite_difference(temperature: float) -> None:
    """`getHenryCoefdT` is `H` times the logarithmic derivative, so a dropped factor shows."""
    co2 = components.entry("co2").henry
    step = 1.0e-4
    numeric = (
        henry.coefficient(co2, temperature + step) - henry.coefficient(co2, temperature - step)
    ) / (2.0 * step)
    assert henry.coefficient_dt(co2, temperature) == pytest.approx(numeric, rel=1.0e-6)
