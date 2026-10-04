"""``designer_report``, against ``process_internals_designer.tsv``'s first row.

The Rust half holds the same numbers in
``crates/azoth-process/src/column/designer.rs``'s own test: the column is solved from the row's
inputs and every quantity the designer publishes is compared with what the class printed. This is
the Python half of that pair, and the two halves are mirrors.
"""

from __future__ import annotations

from azoth import ureg
from azoth.core.warnings import Warning
from azoth.process.reference.designer import DesignerGeometry, designer_report
from azoth.process.reference.distillation_column import _distillation_column_states, _States

Q = ureg.Quantity

COMPONENTS = ["methane", "n-butane"]


def _binary_column() -> tuple[_States, list[Warning]]:
    """The capture's own binary column: four stages between the ends, the feed on stage 2."""
    return _distillation_column_states(
        COMPONENTS,
        Q(7.490704036290964, "mol/s"),
        [0.5, 0.5],
        Q(2.0e6, "Pa"),
        Q(300.0, "K"),
        4,
        2,
        True,
        True,
        Q(1.9e6, "Pa"),
        Q(2.0e6, "Pa"),
        # **The tolerance and the cap come before the two pins**, which is the reference's own
        # order and the reverse of the Rust outcome's.
        1.0e-6,
        200,
        Q(373.15, "K"),
        Q(253.15, "K"),
        None,
        None,
        "direct_substitution",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )


def _close(got: float, want: float) -> bool:
    return abs(got - want) / abs(want) < 1e-5


def test_the_report_reproduces_the_captures_row() -> None:
    """**The eight scalars and all six trays**, at the designer's own constructor defaults.

    The controlling tray is what makes this a measurement rather than an echo: the class picks
    the largest vapour mass flow, and this row's is tray 2 - the feed stage.
    """
    states, _warnings = _binary_column()
    report = designer_report(states, COMPONENTS, DesignerGeometry())

    assert report.controlling_tray_index == 2
    assert report.design_ok is False
    assert report.required_diameter_m == 0.5
    assert _close(report.max_percent_flood, 11.98238455201303)
    assert _close(report.min_percent_flood, 4.438168765587891)
    assert _close(report.average_tray_efficiency, 0.4025182319688279)
    assert _close(report.total_pressure_drop_pa, 1780.1361694102848)

    floods = [
        11.98238455201303,
        11.691032790829645,
        11.526376182416941,
        6.819195914181229,
        6.4771121005782195,
        4.438168765587891,
    ]
    drops = [
        262.37075147939254,
        302.46391775151983,
        320.1925564189236,
        291.20876396625795,
        293.7351343751446,
        310.16504541904646,
    ]
    efficiencies = [
        0.5763466334373019,
        0.4174415777006383,
        0.3741659842116674,
        0.3727128992183818,
        0.365627774305828,
        0.30881452293915046,
    ]
    assert len(report.trays) == 6
    for index, tray in enumerate(report.trays):
        assert _close(tray.percent_flood, floods[index]), f"tray {index} flood"
        assert _close(tray.total_pressure_drop_pa, drops[index]), f"tray {index} drop"
        assert _close(tray.tray_efficiency, efficiencies[index]), f"tray {index} efficiency"
        assert tray.design_ok is False, f"tray {index} is outside the window"


def test_an_override_at_the_sized_value_is_the_sized_row() -> None:
    """The capture's `diameter_override_0_5` row reproduces its sized row to the last digit."""
    states, _warnings = _binary_column()
    sized = designer_report(states, COMPONENTS, DesignerGeometry())
    overridden = designer_report(
        states, COMPONENTS, DesignerGeometry(column_diameter_override_m=0.5)
    )
    assert overridden.required_diameter_m == 0.5
    assert _close(overridden.total_pressure_drop_pa, sized.total_pressure_drop_pa)
    assert overridden.controlling_tray_index == sized.controlling_tray_index
