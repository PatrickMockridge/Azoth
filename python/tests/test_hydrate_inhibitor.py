"""The hydrate inhibitor dose, pinned to NeqSim's own capture.

``test_cross_impl`` says the two kernels agree with each other, which is not the same as either
being right. This pins the Python one to ``captures/hydrate_inhibitor_probe.tsv``, and the Rust
test pins the Rust one to the same file.

**The capture holds two fluids and this pins one.** NeqSim's own ``main`` for this class runs a
``SystemSrkCPAstatoil``; the case reproduces the same composition on a plain ``SystemSrkEos``,
which is the fluid this library can build. The two are far apart - ``0.326`` mol of MEG against
``1.663`` at 270.9 K - because the equilibrium is on water's fugacity and the association moves
it. The CPA column is recorded in the case's ``source`` rather than reproduced.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import InvalidInputError
from azoth.core.units import quantity
from azoth.eos import hydrate_inhibitor_concentration, hydrate_inhibitor_wt

NAMES = ["methane", "ethane", "propane", "i-butane", "MEG", "water"]

MOLES = [quantity(m, "mol") for m in (1.0, 0.10, 0.050, 0.0050, 0.1, 1.0)]

#: `(target in K, the capture's plain-SRK moles, its mass fraction)`.
STATES = [
    (270.9, 1.66321547941976, 0.851421604021146),
    (265.0, 2.16391196885298, 0.881734574278325),
    (275.0, 1.34741017270367, 0.822769695596269),
]


@pytest.mark.parametrize(("target", "want_moles", "want_fraction"), STATES)
def test_the_plain_cubic_column_reproduces_the_capture(
    target: float, want_moles: float, want_fraction: float
) -> None:
    result = hydrate_inhibitor_concentration(
        NAMES, MOLES, "MEG", quantity(target, "K"), quantity(100.0e5, "Pa"), eos="srk"
    )
    # The mole number's bar is looser than the fraction's, and the reason is the secant's own
    # stopping rule: it stops inside 1e-3 K, so the moles it lands on are only fixed to that
    # band divided by dT/dC.
    assert abs(result.inhibitor_moles / want_moles - 1.0) < 1.0e-6
    assert abs(result.weight_fraction / want_fraction - 1.0) < 1.0e-8
    assert abs(float(result.hydrate_temperature.to("K").magnitude) - target) <= 1.0e-3


def test_an_inhibitor_the_feed_does_not_have_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        hydrate_inhibitor_concentration(
            NAMES, MOLES, "methanol", quantity(270.9, "K"), quantity(100.0e5, "Pa"), eos="srk"
        )


#: `(target mass fraction, the capture's plain-SRK inhibitor moles)`. The wt flash's two fluids
#: agree to `3e-5`, unlike the concentration flash's factor of five: its inner step is an
#: ordinary cubic flash, so the association is not what its answer turns on.
WT_DOSES = [
    (0.30, 0.124375874492093),
    (0.50, 0.290210522059004),
    (0.70, 0.677156064897633),
]


@pytest.mark.parametrize(("target", "want"), WT_DOSES)
def test_the_weight_fraction_dose_reproduces_the_capture(target: float, want: float) -> None:
    result = hydrate_inhibitor_wt(
        NAMES, MOLES, "MEG", target, quantity(273.15, "K"), quantity(100.0e5, "Pa"), eos="srk"
    )
    assert abs(result.inhibitor_moles / want - 1.0) < 1.0e-12
    # **The target is met in the aqueous phase**, which is what the secant's residual is on.
    assert abs(result.weight_fraction - target) <= 1.0e-5
    assert result.phases == 2
