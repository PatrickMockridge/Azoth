"""`characterization.combine_reservoir_fluids`, and the row-to-mixture step beside it.

The cases themselves are run by `python/tests/models/test_model_cases.py` and compared across
the two kernels by `test_cross_impl.py`, so this file is the model's own behaviour and the
assembler's: what a rebuilt cut's getter states (`cp_d`'s offset, `m`'s rewrite of the acentric
factor), and that a table assembled into a mixture is the databank's own fluid.
"""

from __future__ import annotations

import math

import pytest

from azoth import characterization, eos
from azoth.characterization.assembly import mixture_from_characterization
from azoth.core.errors import InvalidInputError
from azoth.core.units import ureg

q = ureg.Quantity


def pseudo(
    name: str,
    moles: float,
    molar_mass: float,
    *,
    fluid: int = 0,
    boiling_point: float | None = None,
    acentric_factor: float | None = None,
    attractive_m: float | None = None,
    cp_d: float | None = None,
) -> dict[str, object]:
    """One pseudo-component, with only the fields the test is about stated."""
    row: dict[str, object] = {
        "fluid": fluid,
        "name": name,
        "pseudo_component": True,
        "moles": q(moles, "mol"),
        "molar_mass": q(molar_mass, "kg/mol"),
        "density": q(750.0, "kg/m**3"),
    }
    if boiling_point is not None:
        row["normal_boiling_point"] = q(boiling_point, "K")
    if acentric_factor is not None:
        row["acentric_factor"] = acentric_factor
    if attractive_m is not None:
        row["attractive_m"] = attractive_m
    if cp_d is not None:
        row["cp_d"] = q(cp_d, "J/(mol*K**4)")
    return row


def base(name: str, moles: float, *, fluid: int = 0) -> dict[str, object]:
    return {
        "fluid": fluid,
        "name": name,
        "pseudo_component": False,
        "moles": q(moles, "mol"),
    }


def test_the_base_components_merge_and_the_cuts_are_renamed() -> None:
    result = characterization.combine_reservoir_fluids(
        2,
        [
            base("methane", 2.0),
            pseudo("C7_PC", 1.0, 0.100, boiling_point=370.0),
            base("methane", 3.0, fluid=1),
            pseudo("C10_PC", 1.0, 0.100, fluid=1, boiling_point=440.0),
        ],
    )

    assert [row.name for row in result.components] == ["methane", "PC1_PC", "PC2_PC"]
    assert result.components[0].moles.to("mol").magnitude == pytest.approx(5.0)
    assert result.components[0].pseudo_component is False
    assert result.components[0].molar_mass is None
    assert result.total_moles.to("mol").magnitude == pytest.approx(7.0)


def test_a_pool_with_no_pseudo_component_is_the_base_components() -> None:
    """The class's early return: the requested slate never appears in the answer."""
    result = characterization.combine_reservoir_fluids(
        4, [base("methane", 2.0), base("ethane", 1.0), base("methane", 3.0, fluid=1)]
    )

    assert [(row.name, row.moles.to("mol").magnitude) for row in result.components] == [
        ("methane", 5.0),
        ("ethane", 1.0),
    ]


def test_the_heat_capacity_coefficient_is_the_getter_value() -> None:
    """`getCpD` adds `1e-10` on read and `setCpD` stores what it was given.

    A one-row cut averages the row's number unchanged, so the rebuilt component's getter states
    that number plus the offset again. Reproduced rather than corrected: it is the class's own.
    """
    result = characterization.combine_reservoir_fluids(
        1,
        [
            pseudo(
                "C7_PC",
                1.0,
                0.100,
                boiling_point=370.0,
                acentric_factor=0.4,
                attractive_m=1.0,
                cp_d=1.0e-10,
            )
        ],
    )
    row = result.components[0]
    assert row.cp_d is not None
    assert row.cp_d.to("J/(mol*K**4)").magnitude == pytest.approx(2.0e-10, rel=1e-12)


def test_the_fitted_m_rewrites_the_acentric_factor() -> None:
    """`setm` solves Soave's polynomial behind the setter, and `applyTo` calls it last."""
    m = 1.217132355601029
    result = characterization.combine_reservoir_fluids(
        1,
        [
            pseudo(
                "C7_PC",
                1.0,
                0.100,
                boiling_point=370.0,
                acentric_factor=0.9,
                attractive_m=m,
                cp_d=0.0,
            )
        ],
    )
    root = (1.574 - math.sqrt(1.574**2 - 4.0 * 0.176 * (m - 0.48))) / (2.0 * 0.176)
    row = result.components[0]
    assert row.acentric_factor == pytest.approx(root, rel=1e-12)
    assert row.acentric_factor != pytest.approx(0.9)


def test_a_pseudo_row_without_a_molar_mass_is_refused() -> None:
    row = {
        "fluid": 0,
        "name": "C7_PC",
        "pseudo_component": True,
        "moles": q(1.0, "mol"),
        "density": q(750.0, "kg/m**3"),
    }
    with pytest.raises(InvalidInputError, match="components\\[0\\]\\.molar_mass"):
        characterization.combine_reservoir_fluids(2, [row])


def test_a_target_of_zero_is_refused() -> None:
    with pytest.raises(Exception, match="target_pseudo_components"):
        characterization.combine_reservoir_fluids(0, [base("methane", 1.0)])


def test_a_base_only_table_rebuilds_the_databank_mixture() -> None:
    """The same fluid the databank builds, on a flash - not merely the same constants."""
    rows = [base("methane", 1.0), base("ethane", 1.0), base("propane", 1.0)]
    assembled, ideal_gas = mixture_from_characterization(rows, eos="srk")
    databank, databank_gas = eos.components.mixture_of(["methane", "ethane", "propane"], eos="srk")

    for left, right in zip(assembled.components, databank.components, strict=True):
        assert left.Tc == right.Tc
        assert left.Pc == right.Pc
        assert left.omega == right.omega
        assert left.molar_mass == right.molar_mass
        assert left.volume_shift == right.volume_shift
        assert left.volume_shift_srk == right.volume_shift_srk
    assert ideal_gas.cp_a == databank_gas.cp_a
    assert ideal_gas.cp_b == databank_gas.cp_b
    assert ideal_gas.cp_c == databank_gas.cp_c
    assert ideal_gas.cp_d == databank_gas.cp_d
    assert ideal_gas.cp_e == databank_gas.cp_e

    z = [0.5, 0.3, 0.2]
    ours = eos.pt_flash(assembled, q(250.0, "K"), q(2.0e6, "Pa"), z)
    theirs = eos.pt_flash(databank, q(250.0, "K"), q(2.0e6, "Pa"), z)
    assert ours.phase == theirs.phase
    assert ours.vapour_fraction == pytest.approx(theirs.vapour_fraction, rel=1e-12)
    assert ours.x == pytest.approx(theirs.x, rel=1e-12)
    assert ours.y == pytest.approx(theirs.y, rel=1e-12)


def test_a_pseudo_row_carries_the_databank_derived_shifts() -> None:
    """A cut carrying a substance's own numbers is that substance's component."""
    record = eos.components.entry("n-heptane")
    cp = record.cp
    assert cp is not None
    row = {
        "fluid": 0,
        "name": "n-heptane_PC",
        "pseudo_component": True,
        "moles": q(1.0, "mol"),
        "molar_mass": record.molar_mass,
        "critical_temperature": record.Tc,
        "critical_pressure": record.Pc,
        "acentric_factor": record.omega,
        "racket_z": record.rackett_z,
        "parachor": record.parachor,
        "cp_a": q(cp[0], "J/(mol*K)"),
        "cp_b": q(cp[1], "J/(mol*K**2)"),
        "cp_c": q(cp[2], "J/(mol*K**3)"),
        "cp_d": q(cp[3] + 1.0e-10, "J/(mol*K**4)"),
    }
    assembled, ideal_gas = mixture_from_characterization([row], eos="srk")
    databank, databank_gas = eos.components.mixture_of(["n-heptane"], eos="srk")

    left, right = assembled.components[0], databank.components[0]
    assert left.volume_shift == pytest.approx(right.volume_shift, rel=1e-12)
    assert left.volume_shift_srk == pytest.approx(right.volume_shift_srk, rel=1e-12)
    assert ideal_gas.cp_d == pytest.approx(databank_gas.cp_d, rel=1e-12)


def test_a_pseudo_row_without_its_polynomial_is_refused() -> None:
    row = {
        "fluid": 0,
        "name": "C7_PC",
        "pseudo_component": True,
        "moles": q(1.0, "mol"),
        "molar_mass": q(0.1, "kg/mol"),
        "critical_temperature": q(500.0, "K"),
        "critical_pressure": q(3.0e6, "Pa"),
        "acentric_factor": 0.4,
        "cp_a": q(1.0, "J/(mol*K)"),
        "cp_b": q(0.0, "J/(mol*K**2)"),
        "cp_c": q(0.0, "J/(mol*K**3)"),
    }
    with pytest.raises(InvalidInputError, match="components\\[0\\]\\.cp_d"):
        mixture_from_characterization([row], eos="srk")
