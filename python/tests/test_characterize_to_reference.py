"""The three flags on `characterize_to_reference`, on whichever backends are built.

The cases are run by `python/tests/models/test_model_cases.py` on the selected backend and by
`test_cross_impl.py` on both. What this file adds is what no case can hold: **a refusal has no
case**, and "the flag rebins the source without moving its total" is a claim about the arithmetic
rather than about one slate's numbers.

Both backends, because the Python kernel is the reference the Rust one is read against - a refusal
either of them forgot is a divergence no case would show.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from azoth import characterization
from azoth._dispatch import available, use_backend
from azoth.core.errors import InvalidInputError
from azoth.core.units import ureg

q = ureg.Quantity

BACKENDS = sorted(available())


@pytest.fixture(params=BACKENDS)
def backend(request: pytest.FixtureRequest) -> Iterator[str]:
    """Run the test on each backend that is built, so neither is left unexercised."""
    with use_backend(request.param):
        yield request.param


def row(
    moles: float,
    molar_mass: float,
    density: float | None = None,
    boiling_point: float | None = None,
) -> dict[str, object]:
    """One source row, with only the fields the test is about stated."""
    out: dict[str, object] = {"moles": q(moles, "mol"), "molar_mass": q(molar_mass, "kg/mol")}
    if density is not None:
        out["density"] = q(density, "kg/m**3")
    if boiling_point is not None:
        out["boiling_point"] = q(boiling_point, "K")
    return out


def cut(molar_mass: float, boiling_point: float, moles: float | None = None) -> dict[str, object]:
    """One reference row. Its mole amount is stated only where a test reads it."""
    out: dict[str, object] = {
        "molar_mass": q(molar_mass, "kg/mol"),
        "boiling_point": q(boiling_point, "K"),
    }
    if moles is not None:
        out["moles"] = q(moles, "mol")
    return out


SOURCE = [
    row(1.0, 0.10, 700.0, 378.77),
    row(1.5, 0.14, 760.0, 444.67),
    row(2.0, 0.18, 800.0, 502.058),
    row(1.0, 0.24, 840.0, 573.86),
    row(0.5, 0.32, 880.0, 646.972),
]
REFERENCE = [cut(0.12, 412.832), cut(0.20, 527.8), cut(0.30, 630.83)]


def test_the_delumping_flag_rebins_the_source_without_moving_its_total(backend: str) -> None:
    """Set, the flag moves the group **moles** - and the source's total stays what it was.

    The split decides which group each sub-fraction falls in, so the distribution moves; the mole
    weights conserve each row's moles, so the total does not. A port that dropped the weights
    would move one and not the other, whichever way it dropped them.
    """
    plain = characterization.characterize_to_reference(SOURCE, REFERENCE)
    split = characterization.characterize_to_reference(SOURCE, REFERENCE, None, True, 4, None)

    assert [m.to("mol").magnitude for m in plain.group_moles] != [
        m.to("mol").magnitude for m in split.group_moles
    ], f"{backend}: the split decides which group each sub-fraction falls in"
    # `SOURCE`'s own moles: 1.0 + 1.5 + 2.0 + 1.0 + 0.5.
    for label, result in (("plain", plain), ("delumped", split)):
        assert sum(m.to("mol").magnitude for m in result.group_moles) == pytest.approx(
            6.0, rel=1e-12
        ), f"{backend}: {label} sums to the source's own total"


def test_the_shared_boundaries_refuse_a_slate_with_no_mole_amounts(backend: str) -> None:
    """The equal-mass cuts are placed on masses, so a slate stating none is refused.

    The mole amounts are read by that variant and not by the other - which is why the field is
    optional at all - and a silent zero would place every cut at the same place rather than saying
    it could not.
    """
    assert characterization.characterize_to_reference(SOURCE, REFERENCE, None, None, 4, None)
    with pytest.raises(InvalidInputError, match="moles"):
        characterization.characterize_to_reference(SOURCE, REFERENCE, None, None, 4, True)


def test_the_two_boundary_rules_disagree_only_over_a_cut_between_them(backend: str) -> None:
    """A source cut between the equal-mass cut and the midpoint bins differently; another does not.

    The shared cut is clamped inside the same gap the midpoint sits in, so the two rules agree
    unless a source row lies between them - which is what makes this the only slate on which the
    flag is visible at all.
    """
    reference = [cut(0.12, 412.832, 2.0), cut(0.20, 527.8, 2.0), cut(0.30, 630.83, 1.0)]
    # 460 K is above the equal-mass cut (448.397 K) and below the midpoint (470.316 K).
    straddling = [
        row(1.0, 0.10, 700.0, 378.77),
        row(1.0, 0.185, 800.0, 460.0),
        row(1.0, 0.32, 880.0, 646.972),
    ]

    plain = characterization.characterize_to_reference(straddling, reference)
    shared = characterization.characterize_to_reference(straddling, reference, None, None, 4, True)

    # The default rule bins it with the light cut, so one group holds two rows; the equal-mass rule
    # gives it a group of its own.
    assert len(plain.group_moles) == 2, f"{backend}: the midpoint rule"
    assert len(shared.group_moles) == 3, f"{backend}: the equal-mass rule"
    assert sum(m.to("mol").magnitude for m in shared.group_moles) == pytest.approx(3.0, rel=1e-12)
