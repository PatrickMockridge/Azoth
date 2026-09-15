"""A keycard overrides parameter by parameter and keeps the rest."""

import pytest

import azoth


def test_keycard_overrides_one_parameter_and_keeps_the_rest() -> None:
    card = azoth.keycard.use(
        {
            "schema_version": 2,
            "keyholder": {"name": "smoke"},
            "components": {"methane": {"omega": {"value": 0.1111, "unit": "dimensionless"}}},
        }
    )
    methane = azoth.eos.component("methane", card=card)
    assert methane.omega == pytest.approx(0.1111)  # the card's
    assert methane.Tc.magnitude == pytest.approx(190.56, rel=1e-3)  # the databank's, kept
