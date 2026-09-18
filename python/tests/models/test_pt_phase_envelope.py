"""Spec-driven tests for the ``eos.pt_phase_envelope`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError
from azoth.core.result import PtPhaseEnvelopeResult
from azoth.eos import Mixture, component, from_names, mixture, pt_phase_envelope

MODEL_ID = "eos.pt_phase_envelope"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

METHANE = component("methane")
BUTANE = component("n-butane")


def methane_butane() -> Mixture:
    return from_names(["methane", "n-butane"])


def call(case: dict[str, Any]) -> PtPhaseEnvelopeResult:
    return pt_phase_envelope(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    tolerance = case["tolerance"]
    expected = case["expected"]

    h.assert_close(
        result.critical_temperature.to("K").magnitude,
        expected["critical_temperature"],
        tolerance,
        f"{case['id']} (critical_temperature)",
    )
    h.assert_close(
        result.critical_pressure.to("Pa").magnitude,
        expected["critical_pressure"],
        tolerance,
        f"{case['id']} (critical_pressure)",
    )
    for name in (
        "cricondenbar_temperature",
        "cricondenbar_pressure",
        "cricondentherm_temperature",
        "cricondentherm_pressure",
    ):
        h.assert_close(
            getattr(result, name).to("K").magnitude
            if name.endswith("temperature")
            else getattr(result, name).to("Pa").magnitude,
            expected[name],
            tolerance,
            f"{case['id']} ({name})",
        )
    assert result.iterations == expected["iterations"], f"{case['id']}: iteration count"
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    assert len(CASES) >= 1, f"expected at least one case, found {len(CASES)}"


def test_both_branches_trace_from_the_low_pressure() -> None:
    result = pt_phase_envelope(methane_butane(), P=Q(1.0e5, "Pa"), z=[0.5, 0.5])
    assert len(result.bubble_temperature) > 10, "bubble branch too short"
    assert len(result.dew_temperature) > 10, "dew branch too short"
    for t, p in (
        (result.bubble_temperature, result.bubble_pressure),
        (result.dew_temperature, result.dew_pressure),
    ):
        assert t[-1] > t[0], "a branch should climb in temperature"
        assert p[-1] > p[0], "a branch should climb in pressure"
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_a_single_component_is_refused() -> None:
    propane = mixture([component("propane")])
    with pytest.raises(InvalidInputError) as excinfo:
        pt_phase_envelope(propane, P=Q(1.0e5, "Pa"), z=[1.0])
    assert "pure_saturation" in str(excinfo.value)


def test_a_malformed_composition_is_refused() -> None:
    fluid = methane_butane()
    for z in ([0.5, 0.6], [0.6, -0.1], [0.5, 0.4, 0.1]):
        with pytest.raises(InvalidInputError):
            pt_phase_envelope(fluid, P=Q(1.0e5, "Pa"), z=z)


def test_the_two_backends_agree_on_every_spec_case() -> None:
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        assert py.iterations == rs.iterations, f"{case['id']}: iteration count"
        # The case's *own* tolerance, not a number chosen here. For this model it is looser
        # than the others' because the critical point is the solution of a singular system:
        # the refinement stops at its best iterate rather than converging, so two
        # implementations of the same arithmetic return points about `7.8e-6` relative apart.
        # The case file carries the measurement. A second, tighter number here would be a
        # second source of truth for one declaration.
        h.assert_close(
            py.critical_temperature.to("K").magnitude,
            rs.critical_temperature.to("K").magnitude,
            case["tolerance"],
            f"{case['id']} (critical_temperature)",
        )
        for name in ("bubble_temperature", "bubble_pressure", "dew_temperature", "dew_pressure"):
            assert len(getattr(py, name)) == len(getattr(rs, name)), f"{case['id']}: {name} length"


def test_the_cases_are_announced_in_the_model_docs() -> None:
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "pt_phase_envelope.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/pt_phase_envelope.md" in summary
