"""The Guo-Finch hydrate, pinned to NeqSim's own capture rather than to Rust.

Both kernels carry this route, so `test_cross_impl` already says they agree with each other.
What that check cannot say is whether either is right, and this is the same move
`test_association_agreement.py` makes: both are pinned to the same oracle instead.

The oracle is `validation/neqsim/captures/hydrate_gf_probe.tsv`, which runs NeqSim's two
component models over one plain SRK fluid at three pressures, feeding the guests' reference
fugacities in the way `HydrateFormationTemperatureFlash` does. **The two models give
different answers**: at 50 bara the PVTsim coefficient is `6.74441265462830e-05` and the
Guo-Finch one `1.02598375004935e-04`.
"""

from __future__ import annotations

import pytest

from azoth.eos import components as databank
from azoth.eos.reference import _hydrate

NAMES = ["methane", "ethane", "propane", "water"]

#: 273.15 K, which is where the probe's fluid starts and stays: it is flashed at fixed T.
T = 273.15

#: `(P in bara, [f in bar], pvtsim coefficient, guo_finch coefficient)`, from the capture.
STATES = [
    (
        100.0,
        [70.0988594113204, 4.51610630309347, 0.532415310445719, 0.00447696002836766],
        3.34542273841746e-05,
        5.02541866879314e-05,
    ),
    (
        50.0,
        [38.7358695470026, 3.59682291143079, 0.563592533913349, 0.00425210026944710],
        6.74441265462830e-05,
        1.02598375004935e-04,
    ),
    (
        200.0,
        [120.860624796036, 5.06615849751147, 0.443338479751108, 0.00496253957920541],
        1.77990759845883e-05,
        2.60994368043168e-05,
    ),
]


def guests() -> list[_hydrate.HydrateGuest]:
    records = []
    for name in NAMES:
        entry = databank.entry(name)
        records.append(
            _hydrate.HydrateGuest(
                name=entry.name,
                langmuir_a=entry.hydrate_langmuir_a,
                langmuir_b=entry.hydrate_langmuir_b,
                guo_finch_a=entry.hydrate_guo_finch_a,
                guo_finch_b=entry.hydrate_guo_finch_b,
                former=entry.hydrate_former,
            )
        )
    return records


@pytest.mark.parametrize(("p_bara", "fugacities", "_pvtsim", "guo_finch"), STATES)
def test_the_guo_finch_coefficient_reproduces_the_capture(
    p_bara: float, fugacities: list[float], _pvtsim: float, guo_finch: float
) -> None:
    refs = [value * 1.0e5 for value in fugacities]
    p = p_bara * 1.0e5
    # The Guo-Finch route's reference is the empty lattice's, so this argument is not read:
    # passing the fluid's own water fugacity would be a PVTsim claim.
    _structure, coefficient = _hydrate.stable_structure(
        guests(), refs, _hydrate.GUO_FINCH, T, p, 0.0
    )
    assert coefficient == pytest.approx(guo_finch, rel=1.0e-9), (
        f"P = {p_bara} bara: the coefficient is {coefficient}, the capture says {guo_finch}"
    )


@pytest.mark.parametrize(("p_bara", "fugacities", "pvtsim", "guo_finch"), STATES)
def test_the_two_models_disagree_on_the_same_state(
    p_bara: float, fugacities: list[float], pvtsim: float, guo_finch: float
) -> None:
    """The capture's two columns are two numbers, and this route is the Guo-Finch one.

    The PVTsim column is compared against, not recomputed: its route needs the reference
    water fugacity, which is a cubic solve, and `test_cross_impl` pins it against the same
    capture through the public call.
    """
    refs = [value * 1.0e5 for value in fugacities]
    _structure, coefficient = _hydrate.stable_structure(
        guests(), refs, _hydrate.GUO_FINCH, T, p_bara * 1.0e5, 0.0
    )
    assert pvtsim != pytest.approx(guo_finch, rel=0.1), (
        f"P = {p_bara}: the capture's two columns are too close for this to be a test"
    )
    assert coefficient != pytest.approx(pvtsim, rel=0.1), (
        f"P = {p_bara}: this kernel gave {coefficient}, the capture's *PVTsim* column"
    )
