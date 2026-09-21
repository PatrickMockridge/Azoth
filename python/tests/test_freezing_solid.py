"""The tabulated freezing-point route, pinned to NeqSim's own capture.

``test_cross_impl`` says the two kernels agree with each other, which is not the same as either
being right. This pins the Python one to ``captures/freezing_solid_probe.tsv`` directly, and
the Rust test pins the Rust one to the same file.

**NeqSim's own test only bounds this answer.** `FreezingPointTemperatureFlashTest` asserts
`90 < T < 220` K on this feed; the capture prints `184.710072436643` at 5 bara and
`194.118909895905` at 50.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import SolverNotConvergedError
from azoth.core.units import quantity
from azoth.eos import freezing_point

FEED = ["CO2", "nitrogen", "methane", "ethane", "propane"]

Z = [
    0.0894843679470231,
    0.579634022102985,
    0.170546677734326,
    0.144227745985202,
    0.0161071862304642,
]


@pytest.mark.parametrize(("p_bar", "want"), [(5.0, 184.710072436643), (50.0, 194.118909895905)])
def test_the_tabulated_route_reproduces_the_capture(p_bar: float, want: float) -> None:
    result = freezing_point(FEED, Z, "CO2", quantity(p_bar * 1.0e5, "Pa"))
    assert result.component == "CO2"
    got = float(result.temperature.to("K").magnitude)
    assert abs(got / want - 1.0) < 1.0e-7, f"{p_bar} bar: {got} against NeqSim's {want}"


def test_methane_is_refused() -> None:
    """NeqSim's own guard: `ComponentSolid.fugcoef` returns `1e30` for methane."""
    with pytest.raises(SolverNotConvergedError):
        freezing_point(["methane"], [1.0], "methane", quantity(1.0e5, "Pa"))


def test_a_candidate_the_fluid_does_not_have_is_refused() -> None:
    from azoth.core.errors import InvalidInputError

    with pytest.raises(InvalidInputError):
        freezing_point(["methane"], [1.0], "CO2", quantity(1.0e5, "Pa"))
