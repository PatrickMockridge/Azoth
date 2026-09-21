"""The wax amount and its appearance temperature, as `SKILL.md` documents."""

import pytest

import azoth

Q = azoth.ureg.Quantity

COMPONENTS = ["methane", "n-heptane", "nc14", "nc20"]
Z = [0.7, 0.1, 0.1, 0.1]


def test_the_wax_fraction_rises_as_the_temperature_falls() -> None:
    amounts = [
        azoth.eos.tp_multiflash_wax(
            COMPONENTS, T=Q(t, "K"), P=Q(5.0, "bar"), z=Z, eos="srk"
        ).wax_fraction
        for t in (275.0, 265.0, 255.0)
    ]
    assert amounts[0] == 0.0
    assert amounts[1] == pytest.approx(0.0707, rel=1e-2)
    assert amounts[0] < amounts[1] < amounts[2]


def test_nothing_forms_above_the_appearance_temperature() -> None:
    absent = azoth.eos.tp_multiflash_wax(
        COMPONENTS, T=Q(275.0, "K"), P=Q(5.0, "bar"), z=Z, eos="srk"
    )
    assert absent.wax_fraction == 0.0
    assert absent.phase_count == 2
    # The flag that says which of the two flashes answered: a WAT search bisects on it.
    assert absent.converged is False
