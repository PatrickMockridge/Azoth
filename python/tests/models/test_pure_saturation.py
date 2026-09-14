"""Spec-driven tests for the ``eos.pure_saturation`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PureSaturationResult
from azoth.eos import pr_alpha_ab, pr_departure, pr_kappa, pr_z_factor, pure_saturation

MODEL_ID = "eos.pure_saturation"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> PureSaturationResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return pure_saturation(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the model spec."""
    result = call(case)
    h.assert_close(
        result.p_sat.to("Pa").magnitude,
        case["expected"]["p_sat"],
        case["tolerance"],
        f"{case['id']} (p_sat)",
    )
    h.assert_consistent(result, case["id"])
    assert result.iterations > 0, "the search should have taken at least one step"
    assert result.residual <= SPEC["algorithm"]["tolerance"], (
        f"{case['id']}: residual {result.residual:e} exceeds the declared tolerance"
    )


def test_every_case_ran() -> None:
    """Guard against a spec edit that silently removes every case."""
    assert len(CASES) >= 2, f"expected several cases, found {len(CASES)}"


def test_the_fugacities_agree_at_the_returned_pressure() -> None:
    """The residual vanishes at the answer, recomputed from the kernels.

    The check that does not depend on knowing the answer, and the one that matters
    for a model: `p_sat` is only meaningful if the two fugacities actually agree
    there. Recomputed here through the public kernels rather than read from the
    result, so a result that reported a residual it had not achieved fails.
    """
    for tc, pc, w, t in (
        (369.83, 4_248_000.0, 0.1523, 300.0),
        (304.13, 7_377_000.0, 0.2239, 280.0),
        (425.12, 3_796_000.0, 0.2002, 350.0),
    ):
        result = pure_saturation(Q(tc, "K"), Q(pc, "Pa"), w, Q(t, "K"))
        tr = t / tc
        pr = result.p_sat.to("Pa").magnitude / pc
        kappa = pr_kappa(w).kappa
        ab = pr_alpha_ab(kappa, tr, pr)
        z = pr_z_factor(ab.a_reduced, ab.b_reduced)
        liquid = pr_departure(ab.a_reduced, ab.b_reduced, z.z_min, kappa, tr)
        vapour = pr_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, tr)
        h.assert_close(liquid.ln_phi, vapour.ln_phi, 1e-9, f"fugacities at T={t}")
        h.assert_close(result.ln_phi, liquid.ln_phi, 1e-12, "reported ln_phi")


@pytest.mark.parametrize("t", [369.83, 370.0, 500.0])
def test_a_temperature_at_or_above_the_critical_is_refused(t: float) -> None:
    """The bound is on the ratio, and it is exclusive.

    At exactly `Tc` the two roots have merged, the residual is zero for every
    pressure, and the search would converge on nothing - so the boundary is refused
    as well as past it.
    """
    with pytest.raises(OutOfRangeError) as excinfo:
        pure_saturation(Q(369.83, "K"), Q(4_248_000.0, "Pa"), 0.1523, Q(t, "K"))
    assert excinfo.value.field() == "t_over_tc"


def test_the_result_is_clean_at_an_ordinary_state() -> None:
    result = pure_saturation(Q(369.83, "K"), Q(4_248_000.0, "Pa"), 0.1523, Q(300.0, "K"))
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_the_cases_are_announced_in_the_model_docs() -> None:
    """Every model has a generated page, and the book names it.

    A model rendered nowhere would be a spec only its authors read, which is the
    drift the generated docs exist to prevent - so this asserts the page exists and
    is in the table of contents.
    """
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "pure_saturation.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/pure_saturation.md" in summary, "the model is not in the book's contents"
