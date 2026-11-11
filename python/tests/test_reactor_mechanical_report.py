"""``reactor_mechanical_design``, against ``process_reactor_mechanical_design.tsv``.

The Rust half holds the same numbers in ``crates/azoth-process/src/reactor/mechanical.rs``'s
``tests`` module: each row's inlet is rebuilt and every quantity ``calcDesign`` publishes is
compared with what the class printed. This is the Python half of that pair, and the two halves
are mirrors.

**The tolerance is ``1e-6`` and not tighter** because the volumetric-flow identity is not exact
between the libraries: azoth's system volume and NeqSim's ``getFlowRate("m3/hr")`` differ at about
``2e-8`` relative on these rows, and the catalyst volume and the bed drop carry it.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import InvalidInputError
from azoth.process.reference.reactor_mechanical import (
    ReactorMechanicalGeometry,
    ReactorMechanicalReport,
    reactor_mechanical_design,
)

#: These are three small designs each; the file runs on one xdist worker.
pytestmark = pytest.mark.xdist_group("reactor-mechanical")


def _relative(got: float, want: float, what: str) -> None:
    """The Rust test's own ``1e-6``, for the reason the module docstring gives."""
    assert abs(got - want) / abs(want) < 1e-6, f"{what}: {got} against {want}"


def _z(amounts: list[float]) -> list[float]:
    total = sum(amounts)
    return [amount / total for amount in amounts]


def _defaults_row() -> ReactorMechanicalReport:
    """The ammonia feed at 300 bara, with every class default."""
    return reactor_mechanical_design(
        ["hydrogen", "nitrogen", "ammonia"],
        450.0,
        300.0e5,
        _z([1.5, 0.5, 0.0]),
        1.0,
        ReactorMechanicalGeometry(),
    )


def test_the_defaults_row_floors_the_catalyst_and_sizes_the_shell() -> None:
    """**``gibbs_fixed_bed_defaults``**: the ``0.01`` m3 catalyst floor and a thick shell."""
    report = _defaults_row()

    assert report.reactor_type == "FIXED_BED"
    assert report.number_of_beds == 1
    _relative(report.vessel_diameter_m, 0.5, "vessel_diameter_m")
    _relative(report.vessel_length_m, 2.0, "vessel_length_m")
    _relative(report.shell_thickness_mm, 69.85605419136324, "shell_thickness_mm")
    _relative(report.head_thickness_mm, 69.85605419136324, "head_thickness_mm")
    _relative(report.design_pressure_bara, 330.0, "design_pressure_bara")
    _relative(report.design_temperature_c, 201.85000000000002, "design_temperature_c")
    _relative(report.catalyst_volume_m3, 0.01, "catalyst_volume_m3")
    _relative(report.catalyst_mass_kg, 8.0, "catalyst_mass_kg")
    _relative(report.bed_pressure_drop_bar, 6.101079671327595e-7, "bed_pressure_drop_bar")
    _relative(report.empty_vessel_weight_kg, 1956.1885787131052, "empty_vessel_weight_kg")
    _relative(report.total_equipped_weight_kg, 2196.7996900804355, "total_equipped_weight_kg")


def test_the_cstr_row_takes_the_lhsv_branch_and_leaves_the_diameter_floor() -> None:
    """**``gibbs_cstr``**: the LHSV branch, and the one row above the ``0.5`` m diameter floor."""
    report = reactor_mechanical_design(
        ["hydrogen", "nitrogen", "ammonia"],
        450.0,
        300.0e5,
        _z([1.5, 0.5, 0.0]),
        1.0,
        ReactorMechanicalGeometry(reactor_type="CSTR", lhsv_per_hour=1.5),
    )

    assert report.reactor_type == "CSTR"
    _relative(report.vessel_diameter_m, 0.5385511691956425, "vessel_diameter_m")
    _relative(report.vessel_length_m, 2.15420467678257, "vessel_length_m")
    _relative(report.shell_thickness_mm, 75.24211932030568, "shell_thickness_mm")
    _relative(report.catalyst_volume_m3, 0.3435016535055568, "catalyst_volume_m3")
    _relative(report.catalyst_mass_kg, 274.8013228044454, "catalyst_mass_kg")
    _relative(report.bed_pressure_drop_bar, 1.5260489476579843e-5, "bed_pressure_drop_bar")
    _relative(report.empty_vessel_weight_kg, 2444.4526271505442, "empty_vessel_weight_kg")
    _relative(report.total_equipped_weight_kg, 3006.615755120348, "total_equipped_weight_kg")


def test_the_high_pressure_row_takes_the_max_branch() -> None:
    """**``gibbs_high_pressure``**: a 1 bara feed against a stated ``150`` bara operation limit."""
    report = reactor_mechanical_design(
        ["CO2", "hydrogen", "CO", "water"],
        1000.0,
        1.0e5,
        _z([1.0, 1.0, 0.0, 0.0]),
        1.0,
        ReactorMechanicalGeometry(max_operation_pressure_bara=150.0),
    )

    _relative(report.vessel_diameter_m, 0.5, "vessel_diameter_m")
    _relative(report.vessel_length_m, 2.0, "vessel_length_m")
    _relative(report.design_pressure_bara, 165.0, "design_pressure_bara")
    _relative(report.design_temperature_c, 751.85, "design_temperature_c")
    _relative(report.shell_thickness_mm, 32.2265625, "shell_thickness_mm")
    _relative(report.catalyst_volume_m3, 0.059878658517044936, "catalyst_volume_m3")
    _relative(report.catalyst_mass_kg, 47.90292681363595, "catalyst_mass_kg")
    _relative(report.bed_pressure_drop_bar, 0.0026443223487221906, "bed_pressure_drop_bar")
    _relative(report.empty_vessel_weight_kg, 902.4448091641318, "empty_vessel_weight_kg")
    _relative(report.total_equipped_weight_kg, 1077.5844703902007, "total_equipped_weight_kg")


def test_the_head_is_always_the_shell() -> None:
    """**The branch that cannot answer**, on every captured row and on both halves of the pair:
    the head formula's denominator is the larger of the two at any positive design pressure, so
    ``max(head, shell)`` reports the shell."""
    for report in (
        _defaults_row(),
        reactor_mechanical_design(
            ["hydrogen", "nitrogen", "ammonia"],
            450.0,
            300.0e5,
            _z([1.5, 0.5, 0.0]),
            1.0,
            ReactorMechanicalGeometry(reactor_type="CSTR", lhsv_per_hour=1.5),
        ),
        reactor_mechanical_design(
            ["CO2", "hydrogen", "CO", "water"],
            1000.0,
            1.0e5,
            _z([1.0, 1.0, 0.0, 0.0]),
            1.0,
            ReactorMechanicalGeometry(max_operation_pressure_bara=150.0),
        ),
    ):
        assert report.head_thickness_mm == report.shell_thickness_mm


@pytest.mark.parametrize(
    ("geometry", "field"),
    [
        (ReactorMechanicalGeometry(catalyst_bulk_density=0.0), "catalyst_bulk_density"),
        (
            ReactorMechanicalGeometry(catalyst_particle_diameter_mm=0.0),
            "catalyst_particle_diameter",
        ),
        (ReactorMechanicalGeometry(bed_void_fraction=0.0), "bed_void_fraction"),
        (ReactorMechanicalGeometry(bed_void_fraction=1.0), "bed_void_fraction"),
        (ReactorMechanicalGeometry(number_of_beds=0), "number_of_beds"),
        (ReactorMechanicalGeometry(joint_efficiency=0.0), "joint_efficiency"),
        (ReactorMechanicalGeometry(corrosion_allowance_mm=-1.0), "corrosion_allowance"),
        (ReactorMechanicalGeometry(max_operation_pressure_bara=0.0), "max_operation_pressure"),
        (
            ReactorMechanicalGeometry(max_operation_temperature_k=0.0),
            "max_operation_temperature",
        ),
    ],
)
def test_the_undefined_configurations_are_refused_by_name(
    geometry: ReactorMechanicalGeometry, field: str
) -> None:
    with pytest.raises(InvalidInputError) as caught:
        reactor_mechanical_design(["methane"], 300.0, 10.0e5, [1.0], 1.0, geometry)
    assert caught.value.field() == field
