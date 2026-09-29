"""``size_column_diameter``, the tray calculation's companion, against its capture.

The sizing is a **companion of** ``hydraulics.tray_hydraulics`` rather than one of its outputs:
it mutates the object it is read from - the class writes the trial ``1.0`` m and leaves the sized
value behind - so the probe prints it last, after every other line, and the spec declares no
output for it. The eight values below are
``validation/neqsim/captures/tray_hydraulics_probe.tsv``'s own ``sized_column_diameter_m`` lines,
in the capture's order, and the Rust half holds the same eight in
``crates/azoth-hydraulics/tests/tray_hydraulics.rs``.
"""

from __future__ import annotations

import _helpers as h
from azoth import ureg
from azoth.hydraulics.reference.tray_hydraulics import size_column_diameter

Q = ureg.Quantity

CALC_ID = "hydraulics.tray_hydraulics"

SPEC = h.spec(CALC_ID)

#: The capture's own rows: label and the standard diameter the class sized to.
SIZED_DIAMETERS_M: tuple[tuple[str, float], ...] = (
    ("sieve_default", 0.8),
    ("valve_tray", 0.8),
    ("bubble_cap_tray", 0.9),
    ("sieve_high_vapor", 1.1),
    ("sieve_low_vapor", 0.5),
    ("valve_low_vapor", 0.5),
    ("sieve_stated_weir", 0.8),
    ("sieve_wide", 1.4),
    # The worked example is the first captured row, so it sizes with it.
    ("sieve_default_worked_example", 0.8),
)

ROW_INPUTS = {case["id"]: case["inputs"] for case in h.all_tests(SPEC)}


def test_the_sized_diameter_is_the_standard_table_over_the_captures_rows() -> None:
    """**Every value is a table entry, and two of the eight are the same row's load.**

    The two ``0.5`` rows are the two low-vapour ones, whose flooding velocity is a fraction of
    the default's, and the two ``0.8`` rows include the stated-weir row - because a stated weir
    length moves nothing, which the capture says in the same breath.
    """
    for label, expected in SIZED_DIAMETERS_M:
        inputs = ROW_INPUTS[label]
        sized = size_column_diameter(
            inputs["tray_type"],
            Q(inputs["column_diameter"], "m"),
            Q(inputs["tray_spacing"], "m"),
            Q(inputs["weir_height"], "m"),
            Q(inputs["weir_length"], "m"),
            inputs["downcommer_area_fraction"],
            Q(inputs["hole_diameter"], "mm"),
            inputs["hole_area_fraction"],
            inputs["design_flood_fraction"],
            Q(inputs["vapor_mass_flow"], "kg/s"),
            Q(inputs["liquid_mass_flow"], "kg/s"),
            Q(inputs["vapor_density"], "kg/m**3"),
            Q(inputs["liquid_density"], "kg/m**3"),
            Q(inputs["liquid_viscosity"], "Pa*s"),
            Q(inputs["surface_tension"], "N/m"),
            inputs["relative_volatility"],
        )
        assert sized.to("m").magnitude == expected, f"{label}: {sized} against {expected}"
