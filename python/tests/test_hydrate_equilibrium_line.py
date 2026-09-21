"""The hydrate equilibrium line, pinned to NeqSim's own capture rather than to Rust.

``test_cross_impl`` says the two kernels agree with each other, which is not the same as
either being right. This pins the Python one to ``captures/hydrate_equilibrium_probe.tsv``
directly, and the Rust test pins the Rust one to the same file.

**The agreement across the grid is not uniform**: measured against NeqSim it runs from
``7.6e-10`` at 111.6 bara to ``1.6e-7`` at 177.9, so the bar is the case's own ``1e-6``
rather than the tightest the best point manages.
"""

from __future__ import annotations

from azoth.core.units import quantity
from azoth.eos import hydrate_equilibrium_line

FEED = ["methane", "ethane", "propane", "water"]

Z = [0.7810182896688087, 0.09985170538803756, 0.020266930301532374, 0.09886307464162135]

#: The capture's ten temperatures, in grid order, at 1 to 200 bara.
TEMPERATURES = [
    258.099582916960,
    282.316876654403,
    287.762981872390,
    290.688987372908,
    292.552200127104,
    293.861737270018,
    294.863816418626,
    295.690077826984,
    296.410250673584,
    297.061002920514,
]


def test_the_line_reproduces_the_capture() -> None:
    line = hydrate_equilibrium_line(
        FEED, P_min=quantity(1.0e5, "Pa"), P_max=quantity(200.0e5, "Pa"), z=Z, eos="srk"
    )
    assert len(line.temperature) == 10
    assert len(line.pressure) == 10
    # The grid is NeqSim's own stepping, which lands on the maximum at the last point.
    for point, pressure in enumerate(line.pressure):
        assert pressure == (1.0e5 + (200.0e5 - 1.0e5) / 9.0 * point)
    for point, want in enumerate(TEMPERATURES):
        got = line.temperature[point]
        assert abs(got / want - 1.0) < 1.0e-6, f"grid point {point}: {got} against {want}"


def test_a_grid_that_is_not_a_grid_is_refused() -> None:
    import pytest

    from azoth.core.errors import InvalidInputError, OutOfRangeError

    for low, high in [(200.0, 1.0), (1.0, 1.0), (0.0, 200.0)]:
        with pytest.raises((InvalidInputError, OutOfRangeError)):
            hydrate_equilibrium_line(
                FEED,
                P_min=quantity(low * 1.0e5, "Pa"),
                P_max=quantity(high * 1.0e5, "Pa"),
                z=Z,
                eos="srk",
            )
