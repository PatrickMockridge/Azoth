"""Spec-driven tests for ``hydraulics.crane_k_factors``.

# What these tests can and cannot check

The coefficients in ``data/fittings/crane_k_factors.csv`` are estimated dummy
values, not engineering data. There is therefore no correct answer for this calc
to be checked against, and nothing here can detect a wrong coefficient.

These tests validate the arithmetic, the registry lookup and the error handling.
They are honest about the gap: the numerical assertions compare against values
derived from the same dummy coefficients, so they would all still pass if the
coefficients were nonsense.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError, UnknownFittingError
from azoth.core.result import KFactorsResult
from azoth.hydraulics import crane_k_factors
from azoth.hydraulics.reference.fittings import find_fitting, registry

CALC_ID = "hydraulics.crane_k_factors"

SPEC = h.spec(CALC_ID)
#: Property test name -> the test function implementing it. Name-based rather
#: than a direct reference because the functions are defined further down.
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {"symmetry": "test_symmetry"}

ALL_CASES = h.all_tests(SPEC)
# Property tests assert invariants rather than comparing against a spec value, so
# they have no `inputs` and are exercised by their own named test functions rather
# than by the parametrized runner below.
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> KFactorsResult:
    return crane_k_factors(h.list_input(case, "fittings"), h.input_(case, "f_t"))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.k_total,
        h.expected(case, "k_total"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (k_total)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "f_t":
            return h.input_(_case, "f_t")
        if quantity == "n_fittings":
            return float(len(result.components))
        return None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 3, f"expected several active cases, found {total}"


def test_symmetry() -> None:
    """K is linear in ``f_t`` and additive over fittings.

    These are the invariants that would break if the summation or the data lookup
    were subtly wrong.
    """
    fittings = ["90_elbow", "gate_valve_open"]
    single = crane_k_factors(fittings, 0.018).k_total
    double = crane_k_factors(fittings, 0.036).k_total
    h.assert_close(double, 2.0 * single, 1e-15, "K is linear in f_t")

    reversed_k = crane_k_factors(list(reversed(fittings)), 0.018).k_total
    h.assert_close(reversed_k, single, 1e-15, "K is order-independent")

    part_a = crane_k_factors(["90_elbow"], 0.018).k_total
    part_b = crane_k_factors(["gate_valve_open"], 0.018).k_total
    h.assert_close(single, part_a + part_b, 1e-15, "K is additive over fittings")


def test_unknown_fitting_is_an_error_not_a_zero() -> None:
    """Silently treating an unknown fitting as zero loss would under-report
    pressure drop, which is the dangerous direction to be wrong in."""
    with pytest.raises(UnknownFittingError) as excinfo:
        crane_k_factors(["90_elbow", "no_such_fitting"], 0.018)
    assert "no_such_fitting" in str(excinfo.value)


def test_empty_fitting_list_is_an_error() -> None:
    """An empty list would return K = 0, which reads as "no fittings" when the
    caller may have meant to supply some."""
    with pytest.raises(OutOfRangeError) as excinfo:
        crane_k_factors([], 0.018)
    assert excinfo.value.field() == "n_fittings"


def test_non_positive_friction_factor_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        crane_k_factors(["90_elbow"], 0.0)
    assert excinfo.value.field() == "f_t"


def test_components_echo_the_registry_coefficients() -> None:
    """The per-fitting breakdown must come from the registry rather than being
    recomputed, so a caller can see which coefficients were used."""
    result = crane_k_factors(["90_elbow", "gate_valve_open"], 0.018)
    for component in result.components:
        row = find_fitting(component.fitting_id)
        assert component.n_ld == row.n_ld, "n_ld must come from the registry"
        h.assert_close(component.k, 0.018 * row.n_ld, 1e-15, f"k for {component.fitting_id}")
    h.assert_close(result.k_total, sum(c.k for c in result.components), 1e-15, "k_total is the sum")


def test_the_registry_is_entirely_estimated_dummy_data_right_now() -> None:
    """Fail loudly the day someone populates the registry from a real source.

    At that point this assertion breaks, and whoever did the work is forced to
    update the spec's verification notes and the docs rather than leaving them
    claiming the data is placeholder. Deleting this test is the correct response
    to that failure, not updating the expected number.
    """
    rows = registry()
    estimated = [row.id for row in rows if row.is_estimated]
    assert len(estimated) == len(rows), (
        f"some rows are no longer estimated_dummy ({estimated} of {len(rows)}). If the "
        f"coefficients have now been verified against the primary standard, delete this "
        f"test and update the spec's verification block."
    )


def test_an_estimated_row_says_so_in_its_citation() -> None:
    """The marker and the prose must agree, or a reader skimming the citation
    would not realise the number is a placeholder."""
    for row in registry():
        if row.is_estimated:
            assert "DUMMY" in row.citation.upper(), (
                f"{row.id} is marked estimated but its citation does not say so"
            )


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run.

    Not a name check: for each property the spec asks for, the function that
    implements it is looked up and called here. A spec asking for an invariant
    that nothing checks is the failure mode this exists to catch, and a mapping
    that pointed at a renamed or deleted function would otherwise pass silently.
    """
    declared = {str(c["property"]) for c in PROPERTIES if c.get("property")}
    missing = declared - set(_PROPERTY_IMPLEMENTATIONS)
    assert not missing, (
        f"{CALC_ID}: spec declares property tests {sorted(missing)} that no test implements"
    )
    for property_name in sorted(declared):
        function_name = _PROPERTY_IMPLEMENTATIONS[property_name]
        function = globals().get(function_name)
        assert callable(function), (
            f"{CALC_ID}: property {property_name!r} claims to be implemented by "
            f"{function_name}, which does not exist"
        )
        function()
