"""The Leachman dense root, which no case in the registry can state.

`eos.hydrogen_phase` declares `T` and `P` and carries the spin-isomer as a boundary-only
string, so a case cannot ask for para-hydrogen - and the dense root this file is about is the
one NeqSim's freezing operation sits on, which is para-hydrogen at the triple point. The
states below are `validation/neqsim/captures/freezing_probe.tsv`'s, which are NeqSim's own
`FreezingPointTemperatureFlashTest`'s, and they are checked here rather than in a case file
for that reason.

**The four energies and the entropy are the oracle's; `Z` is not.** The two implementations
agree in `g`, `u`, `h` and `s` to 1e-9 and differ in `Z` by `6.53e-5` relative, by the same
`6.53e-5` at all three states - a constant pressure difference at a fixed density, on an
isotherm whose `dP/drho` is about `2.3e3` against an ideal `R T` of `115`. A wrong term would
move the three states by different amounts and would move the energies with them.
"""

from __future__ import annotations

from typing import Any, Literal

import pytest

from azoth import ureg, use_backend
from azoth.eos import hydrogen_phase

Q = ureg.Quantity

#: (T in K, P in Pa, u, h, s, g) - the fluid phase at each state the freezing flash solves
#: for, from the capture.
STATES = (
    (
        13.803_299_999_737_7,
        7042.0,
        -108.521_707_866_789,
        -108.337_292_151_236,
        -6.217_028_686_431_96,
        -22.521_780_085_440_2,
    ),
    (
        13.915_093_271_702_5,
        351_270.690_962_515_2,
        -108.242_040_229_909,
        -99.066_744_688_923_8,
        -6.197_722_428_295_50,
        -12.824_859_027_069_1,
    ),
    (
        14.398_292_302_249_6,
        1_876_432.785_899_884,
        -106.685_006_050_752,
        -58.210_180_100_061_4,
        -6.110_125_058_512_88,
        29.765_186_495_707_1,
    ),
)

#: The output fields compared, and the unit each is read in.
FIELDS = (("u", "J/mol"), ("h", "J/mol"), ("s", "J/(mol*K)"), ("g", "J/mol"))


def call(t: float, p: float, backend: Literal["python", "rust"], root: str) -> Any:
    """The state at a temperature and pressure, on the root named, in one language."""
    with use_backend(backend):
        return hydrogen_phase(
            T=Q(t, "K"), P=Q(p, "Pa"), hydrogen_type="para", compressed_phase=root
        )


@pytest.mark.parametrize("backend", ["python", "rust"])
def test_the_dense_root_reproduces_the_freezing_states(
    backend: Literal["python", "rust"],
) -> None:
    for state in STATES:
        t, p = state[0], state[1]
        dense = call(t, p, backend, "liquid")

        for (field, unit), expected in zip(FIELDS, state[2:], strict=True):
            actual = getattr(dense, field).to(unit).magnitude
            assert abs(actual / expected - 1.0) <= 1e-9, (
                f"{backend} at {t} K: {field} is {actual}, the capture says {expected}"
            )

        # **And it is not the dilute root**, which is a genuine solution of the same
        # `P(rho) = p` and the one this model returned before the root was a choice - so a
        # solve that ignored the input would pass everything above but this.
        dilute = call(t, p, backend, "vapour").z_factor
        assert abs(dilute - dense.z_factor) > 0.1, (
            f"either root reproduces the other at {t} K: {dilute} against {dense.z_factor}"
        )


def test_the_two_backends_agree_on_the_dense_root() -> None:
    for state in STATES:
        t, p = state[0], state[1]
        python_side = call(t, p, "python", "liquid")
        rust_side = call(t, p, "rust", "liquid")
        assert python_side.z_factor == rust_side.z_factor, f"Z differs at {t} K"
        for field, unit in FIELDS:
            assert (
                getattr(python_side, field).to(unit).magnitude
                == getattr(rust_side, field).to(unit).magnitude
            ), f"{field} differs at {t} K"
