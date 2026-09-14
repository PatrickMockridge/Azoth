"""Spec-driven tests for ``eos.pr_z_factor``."""

from __future__ import annotations

import re
import tomllib
from itertools import pairwise
from pathlib import Path
from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrZFactorResult, RootStructure
from azoth.eos import pr_z_factor

CALC_ID = "eos.pr_z_factor"
REPO_ROOT = Path(__file__).resolve().parents[3]

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "consistency_with": "test_consistency_with",
    "monotonic": "test_monotonic",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def all_text(value: Any) -> list[str]:
    """Every string in a parsed document, in document order."""
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        return [text for item in value.values() for text in all_text(item)]
    if isinstance(value, list):
        return [text for item in value for text in all_text(item)]
    return []


def cubic(a: float, b: float, z: float) -> float:
    """The cubic the calc solves, written out here rather than imported.

    Deliberately a second transcription of the published equation: an imported
    helper would make a wrong coefficient in the implementation wrong in the test
    too, which is the failure this file exists to catch.
    """
    c2 = -(1.0 - b)
    c1 = a - 3.0 * b * b - 2.0 * b
    c0 = -(a * b - b * b - b * b * b)
    return ((z + c2) * z + c1) * z + c0


def call(case: dict[str, Any]) -> PrZFactorResult:
    return pr_z_factor(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    a = float(case["inputs"]["a_reduced"])
    b = float(case["inputs"]["b_reduced"])

    for field in ("z_min", "z_max"):
        actual = getattr(result, field)
        h.assert_close(
            actual,
            h.expected(case, field),
            case.get("tolerance", 1e-12),
            f"{case['id']} ({field})",
        )
        # The residual check comes free with every case, and it is the one check
        # that does not depend on the expected value being right - which matters
        # most at the critical point, where the spec asserts loosely because the
        # answer genuinely is not determined.
        residual = cubic(a, b, actual)
        assert abs(residual) < 1e-12, (
            f"{case['id']}: {field} = {actual} does not satisfy the cubic: f = {residual:e}"
        )

    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity in ("z_min", "z_max"):
            return float(getattr(result, quantity))
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_the_derivations_coefficients_are_what_the_cubic_produces() -> None:
    """Every coefficient quoted in the spec is one the equation actually produces.

    This exists because the prose in a `derivation` is otherwise unchecked. The
    spec prints the cubic's coefficients for each case - a reader is invited to type
    them into a polynomial solver - and nothing verified them against the equation
    until this test. When the spec was first written they were wrong in all four
    cases, which is exactly the failure mode: coefficients that look plausible,
    describe a slightly different cubic, and are wrong in a way no test was
    comparing.

    Checking for presence rather than parsing the prose: the coefficient triple is
    computed from the case's own inputs and each value's `repr` must appear
    somewhere in the spec's text. A mis-typed coefficient does not match.

    The spec is read from its file rather than from the generated registry, because
    the registry flattens `verification` to its status string and drops the notes -
    and the notes are half of what this checks.
    """
    spec_path = REPO_ROOT / "specs" / "calcs" / "eos" / "pr_z_factor.toml"
    document = tomllib.loads(spec_path.read_text(encoding="utf-8"))

    # The coefficient lines the spec prints, with their signs. Matched on the whole
    # `z**3 c2 z**2 + c1 z - c0 = 0` shape rather than by collecting loose numbers,
    # because a coefficient's sign appears in the prose as a separate token - `- 0.97`
    # rather than `-0.97` when the line wraps - so harvesting decimals alone loses it,
    # and a sign error is exactly the kind of slip this is meant to catch.
    #
    # Whitespace is flattened first, because several of those lines wrap. The prose is
    # read from the parsed document rather than from the file's text: how a value is
    # spelled in the file is the file format's business, not the arithmetic's.
    flat = re.sub(r"\s+", " ", " ".join(all_text(document)))
    printed: set[tuple[float, float, float]] = set()
    for terms in re.findall(
        r"z\*\*3 ([+-]) ([\d.]+) z\*\*2 ([+-]) ([\d.]+) z ([+-]) ([\d.]+) = 0", flat
    ):
        sign2, mag2, sign1, mag1, sign0, mag0 = terms
        printed.add(
            (
                float(f"{sign2}{mag2}"),
                float(f"{sign1}{mag1}"),
                float(f"{sign0}{mag0}"),
            )
        )

    assert printed, "the spec prints no cubic in the expected `z**3 ... = 0` form"

    checked = 0
    cases = [("worked_example", document["worked_example"])] + [
        (c["id"], c) for c in document["tests"] if c.get("status") == "active" and "inputs" in c
    ]
    for case_id, case in cases:
        if "expected" not in case:
            continue
        a = float(case["inputs"]["a_reduced"])
        b = float(case["inputs"]["b_reduced"])
        coefficients = (
            -(1.0 - b),
            a - 3.0 * b * b - 2.0 * b,
            -(a * b - b * b - b * b * b),
        )
        assert coefficients in printed, (
            f"{case_id}: the spec does not print the cubic {coefficients}, which is what the "
            f"equation produces for A = {a}, B = {b}. It prints {sorted(printed)}. The "
            f"derivation's arithmetic needs correcting, not the test."
        )
        checked += 1

    assert checked >= 4, f"expected a printed cubic for each of several cases, got {checked}"


def test_consistency_with() -> None:
    """The two `consistency_with` properties the spec declares.

    The first is the residual check, swept over parameters rather than only over
    the named cases: every returned root satisfies the cubic. The second is that
    the admissible count is one or three and never two, which is the theorem in
    ``RootStructure`` and depends on the cubic's constant term being exactly what
    the spec says it is.
    """
    # Every returned root satisfies the cubic.
    for a in (0.05, 0.1, 0.20206500174625697, 0.3, 0.4572355289213822):
        for b in (0.005, 0.012, 0.02431127309496514, 0.05, 0.07779607390388846):
            result = pr_z_factor(a, b)
            for z in (result.z_min, result.z_max):
                assert abs(cubic(a, b, z)) < 1e-12, f"A={a}, B={b}, z={z}"

    # The admissible count is never two.
    seen: set[RootStructure] = set()
    for a_step in range(0, 25):
        for b_step in range(1, 20):
            a = 0.6 * a_step / 24.0
            b = 0.002 + 0.098 * b_step / 19.0
            result = pr_z_factor(a, b)
            seen.add(result.root_structure)
            # The theorem stated on the polynomial as well: f(B) is -2*B**2 < 0,
            # and that is what forces the count.
            assert cubic(a, b, b) < 0.0, f"f(B) should be negative for A={a}, B={b}"
    assert seen == {RootStructure.ONE_ROOT, RootStructure.THREE_ROOTS}, seen


def test_the_middle_root_is_not_returned() -> None:
    """The returned pair is the outermost two, not the two smallest.

    Counted by sign change on a grid, because the endpoints are themselves roots
    and a sign test at them would be comparing rounding noise to zero. Exactly one
    root lies strictly between, so exactly one sign change.
    """
    a, b = 0.20206500174625697, 0.02431127309496514
    result = pr_z_factor(a, b)
    assert result.root_structure is RootStructure.THREE_ROOTS

    steps = 4000
    changes = 0
    previous = cubic(a, b, result.z_min + (result.z_max - result.z_min) * 0.5 / steps)
    for i in range(2, steps):
        z = result.z_min + (result.z_max - result.z_min) * i / steps
        value = cubic(a, b, z)
        if previous * value < 0.0:
            changes += 1
        previous = value
    assert changes == 1, f"exactly one root should lie between, found {changes}"


def test_monotonic() -> None:
    """Both extreme roots fall as ``a_reduced`` rises and rise as ``b_reduced`` rises."""
    for low, high in pairwise((0.10, 0.15, 0.25, 0.30)):
        assert (
            pr_z_factor(high, 0.02431127309496514).z_max
            < pr_z_factor(low, 0.02431127309496514).z_max
        ), f"z_max must fall as A rises: {low} -> {high}"

    for low, high in pairwise((0.005, 0.012, 0.024, 0.06)):
        assert (
            pr_z_factor(0.20206500174625697, high).z_max
            > pr_z_factor(0.20206500174625697, low).z_max
        ), f"z_max must rise as B rises: {low} -> {high}"


def test_root_structure_matches_the_outputs() -> None:
    """One admissible root means the two outputs are equal, and three means they differ.

    The enum output cannot ride through the spec's test cases - `$defs.value` is a
    number or a list of strings, so an enum's expected value has nowhere to live -
    which is why it is asserted here instead.
    """
    cases = (
        (0.20206500174625697, 0.02431127309496514, RootStructure.THREE_ROOTS),
        (0.25, 0.03125, RootStructure.THREE_ROOTS),
        (0.08448417260831159, 0.025932024634629486, RootStructure.ONE_ROOT),
        (0.4572355289213822, 0.07779607390388846, RootStructure.ONE_ROOT),
    )
    for a, b, expected in cases:
        result = pr_z_factor(a, b)
        assert result.root_structure is expected, f"for A={a}, B={b}"
        if expected is RootStructure.ONE_ROOT:
            assert result.z_min == result.z_max
        else:
            assert result.z_max > result.z_min


def test_zero_a_reduced_is_the_hard_sphere_limit_and_is_allowed() -> None:
    """The bound is inclusive at zero on purpose, and two roots are discarded there.

    With no attraction the cubic has roots -0.0587, 0.0101 and 1.0243; the first two
    are at or below B, so they are unphysical and only the third survives. Worth
    asserting because it is counter-intuitive: three real roots does not mean three
    admissible ones.
    """
    b = 0.02431127309496514
    result = pr_z_factor(0.0, b)
    assert result.root_structure is RootStructure.ONE_ROOT
    assert result.z_max > b
    assert abs(cubic(0.0, b, result.z_max)) < 1e-12
    assert cubic(0.0, b, 0.0) > 0.0, "a root lies between 0 and B"
    assert cubic(0.0, b, b) < 0.0, "and B is not one of them"


@pytest.mark.parametrize("b", [0.0, -0.01])
def test_a_non_positive_b_reduced_is_an_error(b: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pr_z_factor(0.2, b)
    assert excinfo.value.field() == "b_reduced"


def test_a_negative_a_reduced_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pr_z_factor(-0.01, 0.024)
    assert excinfo.value.field() == "a_reduced"


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run."""
    declared = {str(c["property"]) for c in PROPERTIES if c.get("property")}
    missing = declared - set(_PROPERTY_IMPLEMENTATIONS)
    assert not missing, (
        f"{CALC_ID}: spec declares property tests {sorted(missing)} that no test implements"
    )
    for property_name in sorted(declared):
        function = globals().get(_PROPERTY_IMPLEMENTATIONS[property_name])
        assert callable(function), f"{CALC_ID}: property {property_name!r} has no implementation"
        function()


def test_the_skipped_case_says_why_it_is_skipped() -> None:
    h.assert_skips_are_explained(SPEC)
    declared = {str(c["property"]) for c in SKIPPED if c.get("property")}
    assert "unit_round_trip" in declared, "this spec skips the unit round trip on purpose"
