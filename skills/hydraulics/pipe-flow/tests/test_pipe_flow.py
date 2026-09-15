"""The pipe-flow result contract `SKILL.md` documents."""

import pytest

import azoth


def test_darcy_weisbach_contract() -> None:
    q = azoth.ureg.Quantity
    r = azoth.hydraulics.darcy_weisbach(
        0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
    )
    assert r.dp.magnitude == pytest.approx(22455.0, rel=1e-3)
    assert azoth.WarningCode.RANGE_CHECK_SKIPPED in {w.code for w in r.warnings}
