"""Spec-driven tests for the ``eos.dew_pressure`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import DewPressureResult
from azoth.eos import Mixture, component, dew_pressure, from_names, mixture, pt_flash

MODEL_ID = "eos.dew_pressure"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table, and this one had drifted: its methane was 0.01142 and
# 4 599 200 Pa, where COMP.csv says 0.0115 and 4 599 000.
METHANE = component("methane")
BUTANE = component("n-butane")
PROPANE = component("propane")


def methane_butane() -> Mixture:
    """The methane/n-butane pair, with the interaction parameter NeqSim fits for it.

    Resolved by name rather than assembled here, so the pair a sweep runs and the pair
    a case runs are the same fluid. The `kij` used to be an "illustrative" 0.05 - a
    number in a test file that described no fluid, and that disagreed with `INTER.csv`.
    """
    return from_names(["methane", "n-butane"])


def ternary() -> Mixture:
    return mixture([METHANE, PROPANE, BUTANE])


def call(case: dict[str, Any]) -> DewPressureResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return dew_pressure(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    tolerance = case["tolerance"]
    expected = case["expected"]

    h.assert_close(
        result.pressure.to("Pa").magnitude,
        expected["pressure"],
        tolerance,
        f"{case['id']} (pressure)",
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


def test_the_flash_agrees_at_the_returned_pressure() -> None:
    """The dew point is where a flash's vapour fraction reaches one.

    The strongest check available without external data, and the two models are
    computed by entirely different procedures: one searches for a state where a
    phase split exists, the other for the pressure at which it stops existing.
    """
    for fluid, t_c, held in (
        (methane_butane(), 300.0, [0.8, 0.2]),
        (methane_butane(), 280.0, [0.8, 0.2]),
        (ternary(), 300.0, [0.5, 0.3, 0.2]),
    ):
        boundary = dew_pressure(fluid, T=Q(t_c, "K"), y=held)
        flash = pt_flash(
            fluid, T=Q(t_c, "K"), P=Q(boundary.pressure.to("Pa").magnitude, "Pa"), z=held
        )
        assert flash.beta is not None, f"T={t_c}: the flash found no split"
        assert abs(flash.beta - 1.0) < 1e-9, f"T={t_c}: beta is {flash.beta}, not one"
        for i, yi in enumerate(held):
            h.assert_close(flash.y[i], yi, 1e-9, f"T={t_c}: vapour {i}")
            h.assert_close(flash.x[i], boundary.incipient[i], 1e-9, f"T={t_c}: liquid {i}")


def test_a_vapour_with_no_dew_point_is_refused() -> None:
    """The residual cannot detect the trivial solution; the K-values can.

    `sum_i x_i K_i = 1` is satisfied by `K_i = 1` at *every* pressure, and
    `S - 1 = sum_i x_i (K_i - 1)` is a weighted sum whose terms cancel - measured
    at 1.5e-13 while the individual K-values were 1e-7 from one. A convergence test
    on `S` alone would accept it silently and return a plausible pressure.
    """
    # Methane-rich vapours well above butane's critical temperature: no liquid can
    # condense out of these at any pressure.
    for t_c, held in ((350.0, [0.8, 0.2]), (450.0, [0.8, 0.2]), (400.0, [0.2, 0.3, 0.5])):
        fluid = ternary() if len(held) == 3 else methane_butane()
        with pytest.raises(OutOfRangeError) as excinfo:
            dew_pressure(fluid, T=Q(t_c, "K"), y=held)
        assert excinfo.value.field() == "min_t_over_tc", f"T={t_c}"


def test_a_genuine_dew_point_is_never_refused() -> None:
    """The guard must not fire on a real bubble point, however slow the iteration.

    The threshold's justification is the gap: genuine boundaries keep
    `max |ln K|` at 1.30 or above at *every* step, and states with no boundary drive
    it below 1e-6.
    """
    for t_c, held in ((300.0, [0.8, 0.2]), (280.0, [0.8, 0.2]), (300.0, [0.5, 0.3, 0.2])):
        fluid = ternary() if len(held) == 3 else methane_butane()
        result = dew_pressure(fluid, T=Q(t_c, "K"), y=held)
        max_ln_k = max(abs(__import__("math").log(k)) for k in result.k)
        assert max_ln_k > 1e-3, (
            f"T={t_c} {held}: max |ln K| is {max_ln_k:e}, which the guard would reject"
        )


def test_a_single_component_is_refused_and_points_at_the_right_calc() -> None:
    propane = mixture([PROPANE])
    with pytest.raises(InvalidInputError) as excinfo:
        dew_pressure(propane, T=Q(300.0, "K"), y=[1.0])
    assert "pure_saturation" in str(excinfo.value), (
        "the refusal should name the calc that does answer it"
    )


def test_a_malformed_composition_is_refused() -> None:
    fluid = methane_butane()
    for y in ([0.8, 0.1], [0.8, 0.3], [1.6, -0.6], [0.8, 0.2, 0.0]):
        with pytest.raises(InvalidInputError):
            dew_pressure(fluid, T=Q(300.0, "K"), y=y)


def test_the_result_is_self_consistent() -> None:
    result = dew_pressure(methane_butane(), T=Q(300.0, "K"), y=[0.8, 0.2])
    assert abs(sum(result.incipient) - 1.0) < 1e-12
    # `K_i = y_i / x_i` up to the normalisation, so relatively rather than
    # absolutely: an absolute bound here would bound how large K may be.
    for i, yi in enumerate([0.8, 0.2]):
        h.assert_close(yi / result.incipient[i] / result.k[i], 1.0, 1e-11, f"K[{i}]")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_the_two_backends_agree_on_every_spec_case() -> None:
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        assert py.iterations == rs.iterations, f"{case['id']}: iteration count"
        h.assert_close(
            py.pressure.to("Pa").magnitude,
            rs.pressure.to("Pa").magnitude,
            1e-12,
            f"{case['id']} (pressure)",
        )
        for name in ("incipient", "k"):
            for i, (left, right) in enumerate(
                zip(getattr(py, name), getattr(rs, name), strict=True)
            ):
                h.assert_close(left, right, 1e-12, f"{case['id']} ({name}[{i}])")


def test_the_cases_are_announced_in_the_model_docs() -> None:
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "dew_pressure.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/dew_pressure.md" in summary
