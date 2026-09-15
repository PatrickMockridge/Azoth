"""The unit constructor and the bare-float rule `SKILL.md` states."""

import pytest

import azoth
from azoth import UnitMismatchError


def test_q_constructs_and_converts() -> None:
    assert azoth.ureg.Quantity(1000.0, "Pa").to("kPa").magnitude == pytest.approx(1.0)


def test_a_bare_float_is_refused_not_misread() -> None:
    q = azoth.ureg.Quantity
    with pytest.raises(UnitMismatchError):
        azoth.hydraulics.darcy_weisbach(
            0.02, 100.0, q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
        )
