"""Spec-driven tests for the ``eos.rachford_rice`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import RachfordRiceResult
from azoth.eos import rachford_rice

MODEL_ID = "eos.rachford_rice"

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

#: NeqSim's ion threshold: a K-value below it takes no part in the split.
ION_THRESHOLD = 1e-30


def call(case: dict[str, Any]) -> RachfordRiceResult:
    """Run one case declared in the model spec."""
    return rachford_rice(**h.model_kwargs(SPEC, case["inputs"]))


def residual(z: list[float], k: list[float], beta: float) -> float:
    """``g(beta)`` recomputed from the inputs, ions excluded."""
    return sum(
        zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0))
        for zi, ki in zip(z, k, strict=True)
        if ki >= ION_THRESHOLD
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the model spec."""
    result = call(case)
    h.assert_close(result.beta, case["expected"]["beta"], case["tolerance"], f"{case['id']} (beta)")
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    """Guard against a spec edit that silently removes every case."""
    assert len(CASES) >= 5, f"expected several cases, found {len(CASES)}"


def test_the_reported_root_solves_the_equation() -> None:
    """The check that does not depend on knowing the answer.

    Every case whose K-values straddle one is a root of `g`, recomputed here from the
    case's own inputs. It is the one property the model claims, and a case whose
    expected value was written from a bad run would fail it.
    """
    checked = 0
    for case in CASES:
        z, k = case["inputs"]["z"], case["inputs"]["K"]
        if not (any(v > 1.0 for v in k) and any(v < 1.0 for v in k)):
            continue
        beta = call(case).beta
        assert abs(residual(z, k, beta)) <= SPEC["algorithm"]["tolerance"], (
            f"{case['id']}: g({beta}) = {residual(z, k, beta):e}, beyond the declared tolerance"
        )
        checked += 1
    assert checked >= 4, f"only {checked} cases had a root to check"


def test_a_root_outside_the_unit_interval_is_reported_not_clamped() -> None:
    """NeqSim clamps; this reports the root, and the difference is deliberate.

    NeqSim's ``calcBeta`` answers `g(0) < 0` and `g(1) > 0` with
    ``phaseFractionMinimumLimit``, measured at both of its solvers by
    ``validation/neqsim/RachfordRiceProbe.java``. ``eos.pt_flash`` reports the negative
    flash instead and explains it, so a clamp here would make that unreachable - which
    is why restoring NeqSim's two early returns has to fail a test.
    """
    subcooled_k = [5.799172708809655, 0.14913889826410245]
    negative = rachford_rice([0.1, 0.9], subcooled_k)
    assert negative.beta < -0.06, f"the subcooled state's root is below zero; got {negative.beta}"
    assert abs(residual([0.1, 0.9], subcooled_k, negative.beta)) <= 1e-10

    superheated = rachford_rice([0.5, 0.5], [1.5, 0.9])
    h.assert_close(superheated.beta, 4.0, 1e-10, "the superheated root")


@pytest.mark.parametrize(
    ("z", "k", "expected"),
    [
        ([0.5, 0.5], [0.2, 0.3], 1e-12),
        ([0.5, 0.5], [3.0, 5.0], 1.0 - 1e-12),
    ],
    ids=["all_liquid", "all_vapour"],
)
def test_a_feed_that_cannot_split_comes_back_clamped(
    z: list[float], k: list[float], expected: float
) -> None:
    """A feed with no root is NeqSim's clamp, not an error."""
    result = rachford_rice(z, k)
    h.assert_close(result.beta, expected, 1e-15, "the clamp")


def test_an_ion_is_skipped() -> None:
    """A K-value below NeqSim's threshold takes no part in the sum."""
    with_ion = rachford_rice([0.1, 0.45, 0.45], [1e-40, 7.304244305324782, 0.33749596785762953])
    h.assert_close(with_ion.beta, 0.6754007405823642, 1e-12, "the ion case")


@pytest.mark.parametrize("bad", [0.0, -0.5])
def test_a_k_value_that_is_not_positive_is_refused(bad: float) -> None:
    with pytest.raises(OutOfRangeError):
        rachford_rice([0.5, 0.5], [bad, 2.0])


def test_a_feed_and_a_k_vector_of_different_lengths_are_refused() -> None:
    with pytest.raises(InvalidInputError):
        rachford_rice([0.5, 0.5], [2.0])


def test_an_empty_feed_is_refused() -> None:
    with pytest.raises(InvalidInputError):
        rachford_rice([], [])
