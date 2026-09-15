"""The API surface `SKILL.md` documents: it must still exist and behave as stated.

These run against whatever backend is installed. The point is not the numbers — the
library's own suite covers those — but that the imports, result fields and error
types the skill tells an agent to use have not drifted.
"""

import pytest

import azoth
from azoth import UnitMismatchError

q = azoth.ureg.Quantity


def test_backends_are_available() -> None:
    assert azoth.backends()


def test_darcy_weisbach_result_contract() -> None:
    r = azoth.hydraulics.darcy_weisbach(
        0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
    )
    assert r.dp.magnitude == pytest.approx(22455.0, rel=1e-3)
    assert azoth.WarningCode.RANGE_CHECK_SKIPPED in {w.code for w in r.warnings}
    assert not r.is_clean


def test_a_bare_number_is_rejected_not_misread() -> None:
    with pytest.raises(UnitMismatchError):
        azoth.hydraulics.darcy_weisbach(
            0.02, 100.0, q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
        )


def test_keycard_is_a_value_not_a_setting() -> None:
    card = azoth.keycard.use({"schema_version": 2, "keyholder": {"name": "smoke"}})
    assert card.keyholder == "smoke"
