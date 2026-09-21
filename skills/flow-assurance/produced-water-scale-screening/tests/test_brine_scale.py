"""The brine screens `SKILL.md` documents."""

import pytest

import azoth

Q = azoth.ureg.Quantity

BRINE = ["water", "Na+", "Cl-", "CO3--", "HCO3-"]
Z = [
    0.802568218298555,
    0.0963081861958266,
    0.0963081861958266,
    0.00321027287319422,
    0.00160513643659711,
]


def test_a_supersaturated_brine_takes_halite_and_lands_at_one() -> None:
    taken = azoth.eos.salt_precipitation(
        BRINE, salt="NaCl", T=Q(298.15, "K"), P=Q(10.0, "bar"), z=Z
    )
    assert taken.initial_saturation_ratio == pytest.approx(1.2722, rel=1e-4)
    assert taken.precipitated_moles > 0.0
    # The answer is the extent at which the ratio reaches one, so it is not an exhaustion.
    assert taken.final_saturation_ratio == pytest.approx(1.0, abs=1e-6)
    assert taken.extent_of_maximum < 1.0


def test_a_brine_that_cannot_form_the_mineral_is_refused() -> None:
    with pytest.raises(azoth.InvalidInputError, match="is built from"):
        azoth.eos.salt_precipitation(
            ["water", "Na+", "Cl-"],
            salt="BaSO4",
            T=Q(298.15, "K"),
            P=Q(10.0, "bar"),
            z=[0.9, 0.05, 0.05],
        )
