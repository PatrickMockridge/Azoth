"""``mechanical_design``, against ``process_column_mechanical_design.tsv``'s three binary rows.

The Rust half holds the same numbers in ``crates/azoth-process/src/column/mechanical.rs``'s own
``report`` test module: the column is solved from the row's inputs and every quantity
``calcDesign`` publishes is compared with what the class printed. This is the Python half of that
pair, and the two halves are mirrors.

**The capture's fourth row, ``absorber_mechanical``, is the Rust half's alone for now**: it is
held by ``crates/azoth-process/tests/absorber.rs``, whose column is the lean-oil absorber solved
unpinned, and the Python twin of that state is not written. The wiring this module's caller will
do puts all four rows on ``process.distillation_column``'s cases, where both implementations are
pinned to them at once.
"""

from __future__ import annotations

import math

import pytest

from azoth import ureg
from azoth.core.errors import InvalidInputError
from azoth.process.reference.distillation_column import _States
from azoth.process.reference.mechanical import (
    MechanicalGeometry,
    MechanicalReport,
    mechanical_design,
    required_wall_thickness,
    resolve_internals_type,
)
from test_designer_report import _binary_column

Q = ureg.Quantity

COMPONENTS = ["methane", "n-butane"]


def _geometry() -> MechanicalGeometry:
    """The class's own constructor defaults, which every capture row of the design states."""
    return MechanicalGeometry()


def _relative(got: float, want: float, what: str) -> None:
    """The Rust test's own ``1e-5``, which is the column port's and not this module's."""
    assert abs(got - want) / abs(want) < 1e-5, f"{what}: {got} against {want}"


@pytest.fixture(scope="module")
def binary() -> _States:
    """The capture's binary column, solved once for the four rows that read it."""
    states, _warnings = _binary_column()
    return states


def test_the_report_reproduces_the_defaults_row(binary: _States) -> None:
    """**`binary_mechanical_defaults`, and the two duties are the column's own.**"""
    report = mechanical_design(binary, COMPONENTS, _geometry())

    assert report.actual_trays == 10
    assert report.material_grade == "SA-516-70"
    _relative(report.vessel_diameter_m, 0.5, "vessel_diameter")
    _relative(report.vessel_height_m, 10.0, "vessel_height")
    _relative(report.vessel_wall_thickness_mm, 216.19496855345912, "vessel_wall_thickness_mm")
    _relative(report.flooding_factor, 0.029981966003039984, "flooding_factor")
    _relative(report.weir_loading, 4.992709047038033, "weir_loading")
    _relative(report.tray_pressure_drop_mbar, 7.137331084047494, "tray_pressure_drop_mbar")
    _relative(report.total_pressure_drop_bar, 0.017801361694102846, "total_pressure_drop_bar")

    # **The duties are pass-throughs, so their comparison is the *column* port's.** The class
    # reads `getReboiler()`/`getCondenser()`; azoth hands over the column's own duties, which
    # carry the two libraries' enthalpy offset - `test_column` holds the same numbers to `5` W.
    assert abs(report.reboiler_duty_kw - 47.78658389381015) < 5.0e-3
    assert abs(report.condenser_duty_kw - 21.323042789303567) < 5.0e-3


def test_a_stated_diameter_moves_the_thickness_and_not_the_flooding(binary: _States) -> None:
    """**`binary_mechanical_override`, which is the whole of the finding.**

    A stated `2.0` m diameter replaces the rating diameter the designer is driven at, and it
    moves **one** published quantity - the wall thickness, because it is the one recomputed after
    the swap. The flooding factor and the weir loading are the defaults row's to the last digit,
    because the class computed them one statement earlier at the `0.5` m it then threw away.
    """
    stated = MechanicalGeometry(column_diameter_override_m=2.0)
    report = mechanical_design(binary, COMPONENTS, stated)

    _relative(report.vessel_diameter_m, 2.0, "vessel_diameter")
    _relative(report.vessel_wall_thickness_mm, 864.7798742138365, "vessel_wall_thickness_mm")
    _relative(report.total_pressure_drop_bar, 0.016680181876092014, "total_pressure_drop_bar")
    # **The two that did not move**, which is the row's reason for existing.
    _relative(report.flooding_factor, 0.029981966003039984, "flooding_factor")
    _relative(report.weir_loading, 4.992709047038033, "weir_loading")
    _relative(report.vessel_height_m, 10.0, "vessel_height")


def test_the_valve_row_moves_the_souders_brown_factor_and_the_internals_type(
    binary: _States,
) -> None:
    """**`binary_mechanical_valve`**: one tray type, two answers."""
    valved = MechanicalGeometry(tray_type="valve")
    report = mechanical_design(binary, COMPONENTS, valved)

    _relative(report.flooding_factor, 0.02498497166919999, "flooding_factor")
    _relative(report.total_pressure_drop_bar, 0.017755928607031072, "total_pressure_drop_bar")
    _relative(report.vessel_diameter_m, 0.5, "vessel_diameter")
    # The weir loading is untouched by the tray type, as it is untouched by the override: it
    # reads the liquid flow, the diameter and nothing else.
    _relative(report.weir_loading, 4.992709047038033, "weir_loading")


def test_a_packed_contactor_is_refused(binary: _States) -> None:
    """**The packed contactor is refused by name rather than half-carried.**"""
    packed = MechanicalGeometry(contactor_internals_type="packed")
    with pytest.raises(InvalidInputError) as caught:
        mechanical_design(binary, COMPONENTS, packed)
    assert "contactor_internals_type" in str(caught.value)


def test_the_two_partial_inputs_are_refused(binary: _States) -> None:
    """**A tray efficiency the class divides by without a guard** is refused here instead, and a
    column with no tray is refused rather than answered with zeros.
    """
    zero = MechanicalGeometry(tray_efficiency=0.0)
    with pytest.raises(InvalidInputError) as caught:
        mechanical_design(binary, COMPONENTS, zero)
    assert "tray_efficiency" in str(caught.value)

    with pytest.raises(InvalidInputError) as empty:
        mechanical_design(_no_tray(binary), COMPONENTS, _geometry())
    assert "trays" in str(empty.value)


def _no_tray(states: _States) -> _States:
    """The solved column with its profile blanked - a column the class returns early on."""
    return states._replace(
        tray_temperature=(),
        tray_pressure=(),
        tray_gas_n=(),
        tray_liquid_n=(),
        tray_gas_z=(),
        tray_liquid_z=(),
    )


def test_the_two_helpers_are_the_class_s_own() -> None:
    """**The two functions a case can reach without a solve**, held here because the report's
    rows above do not exercise either branch that matters: `auto` and the class's own wall
    thickness, which is `bara` over `MPa` and answers millimetres.
    """
    assert resolve_internals_type(MechanicalGeometry(contactor_internals_type="auto")) == "sieve"
    inherited = MechanicalGeometry(contactor_internals_type="auto", tray_type="valve")
    assert resolve_internals_type(inherited) == "valve"
    assert resolve_internals_type(MechanicalGeometry(contactor_internals_type="valve")) == "valve"

    # The `1.1` design margin over the `0.4` allowable stress, which is the class's own mixing of
    # a bara pressure with an MPa strength - and the answer is the capture's, to the last digit.
    assert math.isclose(required_wall_thickness(0.5, 100.0), 216.19496855345912, rel_tol=1e-12)
    # The floor, which no diameter this small clears.
    assert required_wall_thickness(0.01, 1.0) == 6.0


def test_the_report_is_a_value_no_caller_can_mutate() -> None:
    """The two records are frozen, so a report handed to a caller is not a handle on the state."""
    assert MechanicalReport.__dataclass_params__.frozen
    assert MechanicalGeometry.__dataclass_params__.frozen
