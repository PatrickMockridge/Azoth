"""The ideal-gas heat capacity `SKILL.md` documents."""

import pytest

import azoth


def test_ideal_gas_cp_constant_polynomial() -> None:
    q = azoth.ureg.Quantity
    r = azoth.eos.ideal_gas_cp(
        q(30.0, "J/(mol*K)"),
        q(0.0, "J/(mol*K**2)"),
        q(0.0, "J/(mol*K**3)"),
        q(0.0, "J/(mol*K**4)"),
        q(0.0, "J/(mol*K**5)"),
        T=q(300.0, "K"),
    )
    assert r.cp.magnitude == pytest.approx(30.0)
