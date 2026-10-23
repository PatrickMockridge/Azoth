"""`characterization.characterize_to_common_slate`: one slate, each fluid's own moles.

The cases themselves are run by `python/tests/models/test_model_cases.py` and compared across the
two kernels by `test_cross_impl.py`, so this file is the model's own behaviour: that the fluids
share the slate's properties and keep their mole amounts, that a zero weight removes a fluid from
the averages without removing its rows, and what the refusals name.
"""

from __future__ import annotations

import pytest

from azoth import characterization
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.units import ureg

q = ureg.Quantity


def pseudo(
    name: str,
    moles: float,
    molar_mass: float,
    *,
    fluid: int = 0,
    boiling_point: float | None = None,
) -> dict[str, object]:
    row: dict[str, object] = {
        "fluid": fluid,
        "name": name,
        "pseudo_component": True,
        "moles": q(moles, "mol"),
        "molar_mass": q(molar_mass, "kg/mol"),
        "density": q(750.0, "kg/m**3"),
        "critical_temperature": q(500.0, "K"),
        "critical_pressure": q(3.0e6, "Pa"),
        "acentric_factor": 0.4,
    }
    if boiling_point is not None:
        row["normal_boiling_point"] = q(boiling_point, "K")
    return row


def base(name: str, moles: float, *, fluid: int = 0) -> dict[str, object]:
    return {
        "fluid": fluid,
        "name": name,
        "pseudo_component": False,
        "moles": q(moles, "mol"),
    }


def test_the_fluids_share_the_slate_and_keep_their_moles() -> None:
    result = characterization.characterize_to_common_slate(
        [
            base("methane", 10.0),
            pseudo("C7_PC", 3.0, 0.095, boiling_point=370.0),
            base("methane", 4.0, fluid=1),
            pseudo("C10_PC", 4.0, 0.134, fluid=1, boiling_point=435.0),
        ],
        target_pseudo_components=1,
    )

    names = [row.name for row in result.components]
    assert names == ["methane", "PC1_PC", "methane", "PC1_PC"]
    assert list(result.fluid_index) == [0.0, 0.0, 1.0, 1.0]
    first, second = result.components[1], result.components[3]
    # The slate's properties are shared; the mole amounts are each fluid's own.
    assert first.molar_mass == second.molar_mass
    assert first.density == second.density
    assert first.moles.to("mol").magnitude == pytest.approx(3.0)
    assert second.moles.to("mol").magnitude == pytest.approx(4.0)
    assert [moles.to("mol").magnitude for moles in result.fluid_moles] == pytest.approx([13.0, 8.0])
    assert result.resolved_pseudo_components == 1


def test_a_zero_weight_fluid_keeps_its_rows_and_leaves_the_average() -> None:
    rows = [
        base("methane", 10.0),
        pseudo("C7_PC", 1.0, 0.100, boiling_point=370.0),
        pseudo("C10_PC", 1.0, 0.200, fluid=1, boiling_point=440.0),
    ]
    weighted = characterization.characterize_to_common_slate(
        rows, weights=[1.0, 0.0], target_pseudo_components=1
    )
    alone = characterization.characterize_to_common_slate(
        rows, weights=[1.0, 1.0], target_pseudo_components=1
    )

    assert [row.name for row in weighted.components] == ["methane", "PC1_PC", "PC1_PC"]
    # The slate is the first fluid's own cut, because the second fluid is weighted out of it.
    weighted_mass = weighted.components[1].molar_mass
    alone_mass = alone.components[1].molar_mass
    assert weighted_mass is not None and alone_mass is not None
    assert weighted_mass.to("kg/mol").magnitude == pytest.approx(0.100, rel=1e-12)
    assert alone_mass.to("kg/mol").magnitude != pytest.approx(0.100, rel=1e-12)
    # The excluded fluid still comes back with its own component.
    assert weighted.components[2].moles.to("mol").magnitude == pytest.approx(1.0)


def test_no_pseudo_component_returns_the_fluids_as_they_came() -> None:
    result = characterization.characterize_to_common_slate(
        [base("methane", 2.0), base("ethane", 1.0), base("methane", 3.0, fluid=1)]
    )

    assert [row.name for row in result.components] == ["methane", "ethane", "methane"]
    assert list(result.fluid_index) == [0.0, 0.0, 1.0]
    assert [moles.to("mol").magnitude for moles in result.fluid_moles] == pytest.approx([3.0, 3.0])
    assert result.resolved_pseudo_components == 1


def test_the_weights_are_one_per_fluid() -> None:
    rows = [base("methane", 1.0), pseudo("C7_PC", 1.0, 0.1, fluid=1)]
    for weights in ([1.0, 2.0, 3.0], [-1.0, 2.0], [0.0, 0.0]):
        with pytest.raises(InvalidInputError, match="weights"):
            characterization.characterize_to_common_slate(rows, weights=weights)


def test_a_target_that_is_not_positive_is_refused() -> None:
    with pytest.raises(OutOfRangeError, match="target_pseudo_components"):
        characterization.characterize_to_common_slate(
            [pseudo("C7_PC", 1.0, 0.1)], target_pseudo_components=0
        )
