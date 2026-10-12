"""`characterization.delumping`: the two conservations, and what an absent key does.

The cases themselves are run by `python/tests/models/test_model_cases.py`, which compares every
expected row against the selected backend. So this file is the model's own physics rather than a
second restatement of the capture: the two conservations the method is built from - and they are
different mechanisms, which is the reason asserting them together says something - and the sort a
lump with no boiling point gets.
"""

from __future__ import annotations

import pytest

from azoth import characterization
from azoth.core.units import ureg

q = ureg.Quantity


def lump(
    name: str,
    moles: float,
    molar_mass: float,
    density: float | None = None,
    boiling_point: float | None = None,
) -> dict[str, object]:
    """One lump, with only the fields the test is about stated.

    A field left out is the record's own absence, which is what the input kind exists for, and
    each field a caller does state is a quantity - the conversion to SI is the boundary's.
    """
    row: dict[str, object] = {
        "name": name,
        "moles": q(moles, "mol"),
        "molar_mass": q(molar_mass, "kg/mol"),
    }
    if density is not None:
        row["density"] = q(density, "kg/m**3")
    if boiling_point is not None:
        row["boiling_point"] = q(boiling_point, "K")
    return row


def test_the_moles_and_the_mass_are_each_conserved_by_a_different_mechanism() -> None:
    """**The mole weights conserve the moles; one rescale conserves the mass.**

    A port that scaled the moles instead would hold the mass and break the mole count, so
    asserting both on one answer is a claim neither mechanism can satisfy alone.
    """
    result = characterization.delumping([lump("C10_PC", 2.0, 0.18, 820.0, 500.0)], 4)

    moles = sum(f.moles.to("mol").magnitude for f in result.fractions)
    mass = sum((f.moles * f.molar_mass).to("kg").magnitude for f in result.fractions)
    assert moles == pytest.approx(2.0, rel=1e-12)
    assert mass == pytest.approx(2.0 * 0.18, rel=1e-12)


def test_a_lump_without_a_boiling_point_sorts_by_its_molar_mass() -> None:
    """Its sub-fractions state none, and their key falls back - so they come out *first*.

    The heavy lump's boiling points are finite and its key is a temperature, which is larger
    than any molar mass here. A port that carried a sentinel instead of an absence would sort
    the unset lump last, and one that filled the absence in would not state none at all.
    """
    result = characterization.delumping(
        [lump("unset", 1.0, 0.10, 750.0), lump("set", 1.0, 0.20, 840.0, 527.8)], 3
    )

    assert [f.name for f in result.fractions[:3]] == ["unset_d0", "unset_d1", "unset_d2"]
    assert [f.boiling_point for f in result.fractions[:3]] == [None, None, None]
    assert result.fractions[3].name == "set_d0"
    assert result.fractions[3].boiling_point is not None


def test_a_resolution_of_one_returns_the_lumps_in_the_order_given() -> None:
    """The class's own identity, and the order is part of it: a split sort, this does not."""
    result = characterization.delumping([lump("heavy", 2.0, 0.30), lump("light", 1.0, 0.10)], 1)

    assert [f.name for f in result.fractions] == ["heavy", "light"]
    assert result.fractions[0].moles.to("mol").magnitude == 2.0


def test_a_lump_the_method_cannot_split_is_passed_through_rather_than_refused() -> None:
    """`parentMoles > 0`, the class's own guard: it defines what it does with one."""
    result = characterization.delumping([lump("empty", 0.0, 0.20), lump("real", 1.0, 0.20)], 2)

    assert sorted(f.name for f in result.fractions) == ["empty", "real_d0", "real_d1"]


def test_the_answer_names_the_id_it_came_from() -> None:
    """A result carries its own id, so a caller can tell which calculation answered."""
    result = characterization.delumping([lump("a", 1.0, 0.10, 800.0, 400.0)], 2)

    assert result.CALC_ID == "characterization.delumping"
    assert result.is_clean
