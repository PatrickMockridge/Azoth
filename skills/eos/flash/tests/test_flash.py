"""The temperature-pressure flash `SKILL.md` documents."""

import pytest

import azoth


def test_pt_flash_splits_a_binary() -> None:
    fluid = azoth.eos.from_names(["methane", "n-butane"])
    r = azoth.eos.pt_flash(
        fluid, T=azoth.ureg.Quantity(330.0, "K"), P=azoth.ureg.Quantity(2.5e6, "Pa"), z=[0.6, 0.4]
    )
    assert r.beta == pytest.approx(0.842, rel=1e-3)
    assert str(r.phase) == "two_phase"
