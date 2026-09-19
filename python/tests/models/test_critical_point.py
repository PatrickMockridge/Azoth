"""Spec-driven tests for the ``eos.critical_point`` model.

The spec's cases pin the two implementations to each other. These tests are for what a
case cannot say, and one of them matters more than the rest: **a pure component's
critical point is known in closed form** for Peng-Robinson, so the whole construction -
the constant-volume Hessian, the ideal part, the scaling, the cubic form and the
nesting - is checked against an analytic answer rather than against a second reading of
the same method.
"""

from __future__ import annotations

from typing import Any, cast

import pytest

import _helpers as h
from azoth import _models_gen, ureg, use_backend
from azoth.core.errors import InvalidInputError
from azoth.core.result import CriticalPointResult
from azoth.eos import Component, Mixture, component, critical_point, from_names, mixture
from azoth.eos.reference.pr_alpha_ab import OMEGA_B

MODEL_ID = "eos.critical_point"
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

#: `(1 - omega_b)/3`, which every pure Peng-Robinson fluid has at `Tr = Pr = 1`.
PURE_COMPRESSIBILITY = (1.0 - OMEGA_B) / 3.0


def methane_butane() -> Mixture:
    """The methane/n-butane pair, with the interaction parameter NeqSim fits for it.

    Resolved by name rather than assembled here, so the pair a sweep runs and the pair
    a case runs are the same fluid. The `kij` used to be an "illustrative" 0.05 - a
    number in a test file that described no fluid, and that disagreed with `INTER.csv`.
    """
    return from_names(["methane", "n-butane"])


def pure(component: Component) -> Any:
    return mixture([component])


def call(case: dict[str, Any]) -> CriticalPointResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return critical_point(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    for attribute, unit in (("tc", "K"), ("pc", "Pa"), ("vc", "m**3/mol")):
        h.assert_close(
            getattr(result, attribute).to(unit).magnitude,
            case["expected"][attribute],
            case["tolerance"],
            f"{case['id']} ({attribute})",
        )
    h.assert_close(result.z_c, case["expected"]["z_c"], case["tolerance"], f"{case['id']} (z_c)")
    h.assert_consistent(result, case["id"])


def test_every_case_ran() -> None:
    assert len(CASES) >= 2, f"expected several cases, found {len(CASES)}"


def test_a_pure_component_reproduces_the_analytic_critical_point() -> None:
    """The one check that does not come from either implementation.

    At ``Tr = Pr = 1`` a pure Peng-Robinson fluid's critical compressibility is
    ``(1 - omega_b)/3`` in closed form, so the whole construction has an independent
    answer to be compared against.

    **The tolerances are 2e-04 and 5e-06, not the 1e-10 and 1e-08 they used to be,
    and that is the price of the port.** ``(1 - omega_b)/3`` is the critical point of
    a cubic whose ``Omega_a`` and ``Omega_b`` satisfy the triple-root condition.
    NeqSim's do not - they depart by p = -5.6e-6 - so its cubic is at a triple root
    at ``Tr = 1 + 4.8e-05`` rather than at ``Tr = 1``, and that is exactly the
    residue this test now measures. It is four or five figures for every component,
    which is what a systematic departure looks like rather than an iteration that
    has not finished. See ``eos.pr_alpha_ab``'s assumptions and
    ``validation/eos/methane_butane_flash_against_neqsim.json``.
    """
    for name, substance in (
        ("propane", PROPANE),
        ("methane", METHANE),
        ("n-butane", BUTANE),
        ("carbon dioxide", component("co2")),
    ):
        t_c = substance.Tc.to_base_units().magnitude
        p_c = substance.Pc.to_base_units().magnitude
        result = critical_point(pure(substance), [1.0])

        h.assert_close(result.tc.to("K").magnitude, t_c, 2e-4, f"{name}: Tc")
        h.assert_close(result.pc.to("Pa").magnitude, p_c, 2e-4, f"{name}: Pc")
        # Measured at 2.4e-6 for all four, against an analytic 0.307398889.
        assert abs(result.z_c - PURE_COMPRESSIBILITY) < 5e-6, (
            f"{name}: Z_c is {result.z_c!r} against the analytic {PURE_COMPRESSIBILITY!r}"
        )


def test_a_mixtures_critical_compressibility_varies_with_composition() -> None:
    """The discriminating test, and the reason a pure-component check is not enough.

    Solving ``dP/dV = d2P/dV2 = 0`` at fixed composition reproduces a pure component's
    critical point exactly and returns ``(1 - omega_b)/3`` for **every** mixture,
    because in reduced variables those two conditions have a single universal root. A
    pure-component check cannot tell the two routes apart; this can, because a real
    mixture critical compressibility moves with composition and a constant does not.
    """
    seen = [
        critical_point(methane_butane(), [fraction, 1.0 - fraction]).z_c
        for fraction in (0.2, 0.4, 0.6, 0.8)
    ]
    spread = max(seen) - min(seen)
    assert spread > 0.1, (
        f"Z_c varied by only {spread} across the composition range: {seen}. A spread near "
        f"zero is what the mechanical conditions give, and is not a mixture critical point."
    )


def test_a_binarys_critical_locus_falls_between_its_pure_endpoints() -> None:
    """A shape check, and the one that catches an iteration on the wrong root.

    The critical locus of a binary is a continuous curve from one pure component to the
    other, so a point outside the bracket or out of order is wrong whatever value it is
    near.
    """
    fluid = methane_butane()
    previous = BUTANE.Tc.to_base_units().magnitude
    for fraction in (0.2, 0.4, 0.6, 0.8):
        t_c = critical_point(fluid, [fraction, 1.0 - fraction]).tc.to("K").magnitude
        assert METHANE.Tc.to_base_units().magnitude < t_c < previous, (
            f"Tc = {t_c} at z = {fraction} is outside the pure endpoints or out of order"
        )
        previous = t_c


def test_a_malformed_composition_is_refused() -> None:
    """Checked, not renormalised: a caller's error must not be invisible downstream."""
    fluid = methane_butane()
    for label, z in (
        ("too short", [0.5]),
        ("too long", [0.5, 0.25, 0.25]),
        ("negative", [1.5, -0.5]),
        ("does not sum to one", [0.5, 0.6]),
    ):
        with pytest.raises(InvalidInputError):
            critical_point(fluid, z)
            pytest.fail(f"{label} was accepted")


def test_the_two_backends_agree_on_every_spec_case() -> None:
    """Including the iteration counts, which must match exactly.

    A ULP difference that flips one iteration shows up here as ``26 vs 27`` rather than
    as a mystifying drift in the fifth digit, which is what makes this the sharpest
    cheap check the two languages have.
    """
    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        for name, attribute in (("Tc", "tc"), ("Pc", "pc"), ("Vc", "vc")):
            unit = {"Tc": "K", "Pc": "Pa", "Vc": "m**3/mol"}[name]
            h.assert_close(
                getattr(py, attribute).to(unit).magnitude,
                getattr(rs, attribute).to(unit).magnitude,
                1e-9,
                f"{case['id']} ({name})",
            )
        h.assert_close(py.z_c, rs.z_c, 1e-9, f"{case['id']} (Z_c)")
        assert py.iterations == rs.iterations, (
            f"{case['id']}: the two backends took {py.iterations} and {rs.iterations} "
            f"iterations, so they did not run the same procedure"
        )


def test_the_cases_are_announced_in_the_model_docs() -> None:
    from pathlib import Path

    root = Path(__file__).resolve().parents[3]
    page = root / "docs" / "src" / "eos" / "critical_point.md"
    assert page.is_file(), f"{page} was not generated"
    summary = (root / "docs" / "src" / "SUMMARY.md").read_text(encoding="utf-8")
    assert "eos/critical_point.md" in summary


def test_an_associating_mixture_is_refused_rather_than_answered_from_the_cubic() -> None:
    """The critical point is the criticality matrix, and this mixture has no matrix.

    `Mixture::criticality_matrix` is the Helmholtz Hessian reshaped, and the Hessian
    carries the cubic's second derivatives alone: the association's second composition
    derivative is not assembled in either kernel. So for an associating mixture the
    matrix would be a *different fluid's*, and the critical point it places would be a
    point this mixture does not have.

    The rule is the one the heat-capacity departure follows - report the association's
    absence rather than the cubic's value - and it is the same in both kernels because
    the refusal is written in both. Before the association crossed the Python-to-Rust
    boundary neither kernel could be handed an associating mixture here at all, so this
    is a case the transport fix created as much as it exposed.
    """
    from azoth.eos import components as databank

    associating = databank.from_names(["water", "methanol"], eos="srk", associating=True)
    for backend in ("python", "rust"):
        expected = pytest.raises(InvalidInputError, match="associates")
        with use_backend(cast("Any", backend)), expected:
            critical_point(associating, [0.6, 0.4])

    # The control: a mixture that does not associate is answered, in both kernels.
    plain = databank.from_names(["methane", "n-butane"])
    answers = []
    for backend in ("python", "rust"):
        with use_backend(cast("Any", backend)):
            answers.append(critical_point(plain, [0.6, 0.4]))
    h.assert_close(answers[1].tc.to("K").magnitude, answers[0].tc.to("K").magnitude, 1e-9, "tc")
