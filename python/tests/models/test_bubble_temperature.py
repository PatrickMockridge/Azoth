"""Spec-driven tests for the ``eos.bubble_temperature`` model."""

from __future__ import annotations

import math
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import BubbleTemperatureResult
from azoth.eos import Mixture, bubble_temperature, component, from_names, mixture, pt_flash

MODEL_ID = "eos.bubble_temperature"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table.
METHANE = component("methane")
BUTANE = component("n-butane")
PROPANE = component("propane")


def methane_butane() -> Mixture:
    """The methane/n-butane pair, with the interaction parameter NeqSim fits for it."""
    return from_names(["methane", "n-butane"])


def ternary() -> Mixture:
    return mixture([METHANE, PROPANE, BUTANE])


def call(case: dict[str, Any]) -> BubbleTemperatureResult:
    """Run one case declared in the model spec."""
    return bubble_temperature(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    tolerance = case["tolerance"]
    expected = case["expected"]

    h.assert_close(
        result.temperature.to("K").magnitude,
        expected["temperature"],
        tolerance,
        f"{case['id']} (T)",
    )
    for name in ("z_liquid", "z_vapour", "min_t_over_tc"):
        h.assert_close(getattr(result, name), expected[name], tolerance, f"{case['id']} ({name})")
    assert result.iterations == expected["iterations"], f"{case['id']}: iteration count"
    # `residual` is not pinned: it is a weighted sum whose terms cancel, so it
    # differs between implementations in its last digits. See the spec.
    assert result.residual <= SPEC["algorithm"]["tolerance"], f"{case['id']}: residual"

    for name in ("incipient", "k"):
        actual = getattr(result, name)
        assert len(actual) == len(expected[name]), f"{case['id']}: {name} length"
        for i, (got, want) in enumerate(zip(actual, expected[name], strict=True)):
            h.assert_close(got, want, tolerance, f"{case['id']} ({name}[{i}])")

    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    assert len(CASES) >= 2, f"expected several cases, found {len(CASES)}"


def test_the_flash_agrees_at_the_returned_temperature() -> None:
    """The bubble point is where a flash's vapour fraction vanishes.

    The strongest check available without external data, and the two models are
    computed by entirely different procedures: one searches for a state where a
    phase split exists, the other for the temperature at which it stops existing.
    """
    for fluid, p, held in (
        (methane_butane(), 3950960.4937437344, [0.2, 0.8]),
        (ternary(), 4664798.1363369655, [0.2, 0.3, 0.5]),
    ):
        boundary = bubble_temperature(fluid, P=Q(p, "Pa"), x=held)
        flash = pt_flash(
            fluid,
            T=Q(boundary.temperature.to("K").magnitude, "K"),
            P=Q(p, "Pa"),
            z=held,
        )
        assert flash.beta is not None, f"P={p}: the flash found no split"
        assert abs(flash.beta) < 1e-8, f"P={p}: beta is {flash.beta}, not zero"
        for i, xi in enumerate(held):
            h.assert_close(flash.x[i], xi, 1e-8, f"P={p}: liquid {i}")
            # The incipient composition is a successive-substitution estimate from
            # `y_i = x_i K_i / S`, accurate to ~1e-8 *absolutely* rather than to the
            # flash's exact split, so the comparison is absolute rather than relative.
            assert abs(flash.y[i] - boundary.incipient[i]) < 1e-7, f"P={p}: vapour {i}"


def test_a_mixture_with_no_bubble_point_is_refused() -> None:
    """The residual cannot detect the trivial solution; the K-values can.

    `sum_i x_i K_i = 1` is satisfied by `K_i = 1` at *every* temperature, and
    `S - 1 = sum_i x_i (K_i - 1)` is a weighted sum whose terms cancel. A convergence
    test on `S` alone would accept it silently.
    """
    for p in (6.0e6, 8.0e6):
        with pytest.raises(OutOfRangeError) as excinfo:
            bubble_temperature(methane_butane(), P=Q(p, "Pa"), x=[0.2, 0.8])
        assert excinfo.value.field() == "min_t_over_tc", f"P={p}"


def test_a_genuine_bubble_point_is_never_refused() -> None:
    """The guard must not fire on a real bubble point, however slow the iteration."""
    for p, held in (
        (3950960.4937437344, [0.2, 0.8]),
        (4664798.1363369655, [0.2, 0.3, 0.5]),
    ):
        fluid = ternary() if len(held) == 3 else methane_butane()
        result = bubble_temperature(fluid, P=Q(p, "Pa"), x=held)
        max_ln_k = max(abs(math.log(k)) for k in result.k)
        assert max_ln_k > 1e-3, (
            f"P={p} {held}: max |ln K| is {max_ln_k:e}, which the guard would reject"
        )


def test_a_single_component_is_refused_and_points_at_the_right_calc() -> None:
    propane = mixture([PROPANE])
    with pytest.raises(InvalidInputError) as excinfo:
        bubble_temperature(propane, P=Q(1.0e6, "Pa"), x=[1.0])
    assert "pure_saturation" in str(excinfo.value), (
        "the refusal should name the calc that does answer it"
    )


def test_a_malformed_composition_is_refused() -> None:
    fluid = methane_butane()
    for x in ([0.2, 0.9], [0.2, 0.7], [0.6, -0.6], [0.2, 0.8, 0.0]):
        with pytest.raises(InvalidInputError):
            bubble_temperature(fluid, P=Q(3.0e6, "Pa"), x=x)


def test_the_result_is_self_consistent() -> None:
    result = bubble_temperature(methane_butane(), P=Q(3950960.4937437344, "Pa"), x=[0.2, 0.8])
    assert abs(sum(result.incipient) - 1.0) < 1e-12
    # `K_i = y_i / x_i` up to the normalisation, so relatively rather than
    # absolutely: an absolute bound here would bound how large K may be.
    for i, xi in enumerate([0.2, 0.8]):
        h.assert_close(result.incipient[i] / xi / result.k[i], 1.0, 1e-11, f"K[{i}]")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_the_two_backends_agree_on_every_spec_case() -> None:
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        assert py.iterations == rs.iterations, f"{case['id']}: iteration count"
        h.assert_close(
            py.temperature.to("K").magnitude,
            rs.temperature.to("K").magnitude,
            1e-12,
            f"{case['id']} (T)",
        )
        for name in ("incipient", "k"):
            for i, (left, right) in enumerate(
                zip(getattr(py, name), getattr(rs, name), strict=True)
            ):
                h.assert_close(left, right, 1e-12, f"{case['id']} ({name}[{i}])")


def test_the_cases_are_announced_in_the_model_docs() -> None:
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "bubble_temperature.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/bubble_temperature.md" in summary
