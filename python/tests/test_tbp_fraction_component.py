"""`characterization.tbp_fraction_component`: the card, and the two acentric factors in it.

The cases themselves are run by `python/tests/models/test_model_cases.py` and compared across the
two kernels by `test_cross_impl.py`; this file pins the two facts a case states but does not
explain - that the card's acentric factor is the fitted coefficient's root rather than the
model's Edmister value, and that the Rackett compressibility is a reference flash rather than a
correlation.
"""

from __future__ import annotations

import pytest

from azoth import characterization, ureg
from azoth.core.errors import OutOfRangeError

q = ureg.Quantity


def test_the_card_states_the_roots_acentric_factor_not_the_models() -> None:
    """`setm` rewrites the acentric factor, and `init` re-derives the coefficient from it.

    The Edmister value the TBP model computes (0.45870775) and the root of the polynomial the
    fitted `m` belongs to (0.46493310) are different numbers, and the card states the root.
    """
    cut = characterization.tbp_fraction_component(q(0.095, "kg/mol"), q(730.0, "kg/m**3"))

    assert cut.acentric_factor == pytest.approx(0.4649331015568601, rel=1e-12)
    assert cut.acentric_factor != pytest.approx(0.4587077521283205, rel=1e-6)
    assert cut.attraction_exponent == pytest.approx(1.1737600510000004, rel=1e-12)


def test_the_racket_compressibility_is_the_reference_flash() -> None:
    """The recorded Rackett value for C7, which a correlation on the two numbers would miss."""
    cut = characterization.tbp_fraction_component(q(0.095, "kg/mol"), q(730.0, "kg/m**3"))

    assert cut.racket_z == pytest.approx(0.27686767147391916, rel=1e-12)
    assert cut.racket_z_cpa == cut.racket_z


def test_the_heat_capacity_coefficient_states_the_getter() -> None:
    """`cp_d` is stored as zero and its getter adds `1e-10`, which is what the card states."""
    cut = characterization.tbp_fraction_component(q(0.095, "kg/mol"), q(730.0, "kg/m**3"))

    assert cut.cp_d.to("J/(mol*K**4)").magnitude == pytest.approx(1.0e-10, rel=1e-12)


def test_a_cut_without_a_positive_mass_or_density_is_refused() -> None:
    with pytest.raises(OutOfRangeError):
        characterization.tbp_fraction_component(q(0.0, "kg/mol"), q(730.0, "kg/m**3"))
    with pytest.raises(OutOfRangeError):
        characterization.tbp_fraction_component(q(0.095, "kg/mol"), q(0.0, "kg/m**3"))
