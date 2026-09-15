"""The composition lookup and mixture build `SKILL.md` documents."""

import pytest

import azoth


def test_from_names_builds_a_mixture() -> None:
    fluid = azoth.eos.from_names(["methane", "n-butane"])
    assert len(fluid) == 2


def test_component_lookup_returns_critical_constants() -> None:
    methane = azoth.eos.component("methane")
    assert methane.Tc.magnitude == pytest.approx(190.56, rel=1e-3)
