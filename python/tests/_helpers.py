"""Shared helpers for the spec-driven tests.

Mirrors ``crates/chemeng-hydraulics/tests/common/mod.rs``. The tests in this
directory are not hand-written per calc: they walk the spec's ``tests`` list and
execute whatever it declares, so a test added to a spec YAML runs in both
languages with no new test code. That is the whole point of generating the
registry on both sides.

``python/tests`` is on ``pythonpath`` (see ``pyproject.toml``) so this module is
importable as a top-level ``_helpers``.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable, Iterable, Mapping
from typing import Any

import pint

from chemeng._registry_gen import BY_ID, CALCS
from chemeng.core.result import FlowRegime
from chemeng.core.units import to_si
from chemeng.core.warnings import Warning, WarningCode

__all__ = [
    "CALCS",
    "all_tests",
    "assert_close",
    "assert_consistent",
    "assert_results_equal",
    "assert_skips_are_explained",
    "assert_warnings_agree_with_spec",
    "expected",
    "input_",
    "list_input",
    "spec",
]


def spec(calc_id: str) -> dict[str, Any]:
    """Fetch a spec, failing loudly if the id is wrong."""
    assert calc_id in BY_ID, f"no spec for {calc_id!r}; registry has {sorted(BY_ID)}"
    return BY_ID[calc_id]


def all_tests(spec_: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The worked example plus every other test, in spec order.

    Mirrors ``CalcSpec::all_tests`` in Rust. The worked example is a separate
    block in the spec rather than an entry in ``tests``, so its inputs, expected
    values and tolerance have to be folded into a runnable case. Returning it
    alongside the rest means a test file can walk one list instead of special
    casing the example everywhere.
    """
    example = spec_["worked_example"]
    declared = next(
        (t for t in spec_["tests"] if t["type"] == "worked_example"),
        {"id": f"{str(spec_['id']).rpartition('.')[2]}_worked_example", "status": "active"},
    )
    worked = {
        "id": declared["id"],
        "type": "worked_example",
        "status": declared["status"],
        "inputs": example["inputs"],
        "expected": example["expected"],
        "tolerance": example["tolerance"],
    }
    if "skip_reason" in declared:
        worked["skip_reason"] = declared["skip_reason"]
    others = [t for t in spec_["tests"] if t["type"] != "worked_example"]
    return [worked, *others]


def assert_close(actual: float, expected_value: float, tolerance: float, context: str) -> None:
    """Relative comparison with a readable failure message.

    Relative rather than absolute, because the quantities here span nine orders of
    magnitude (Re ~ 1e5, f ~ 1e-2) and a single absolute tolerance would be
    meaningless at both ends.
    """
    assert actual == actual, f"{context}: got NaN, expected {expected_value}"
    assert actual not in (float("inf"), float("-inf")), (
        f"{context}: got {actual}, expected {expected_value}"
    )
    scale = max(abs(expected_value), float.fromhex("0x0.0000000000001p-1022"))
    relative = abs(actual - expected_value) / scale
    assert relative <= tolerance, (
        f"{context}: got {actual}, expected {expected_value}\n"
        f"  relative error {relative:.3e} exceeds tolerance {tolerance:.3e}"
    )


def input_(case: Mapping[str, Any], name: str) -> float:
    """A scalar input from a test case."""
    inputs = case["inputs"]
    assert name in inputs, (
        f"test {case['id']!r} does not supply input {name!r}; it has {sorted(inputs)}"
    )
    return float(inputs[name])


def expected(case: Mapping[str, Any], name: str) -> float:
    """An expected output from a test case."""
    outputs = case["expected"]
    assert name in outputs, (
        f"test {case['id']!r} does not assert {name!r}; it asserts {sorted(outputs)}"
    )
    return float(outputs[name])


def list_input(case: Mapping[str, Any], name: str) -> list[str]:
    """A list-valued input from a test case."""
    inputs = case["inputs"]
    assert name in inputs, f"test {case['id']!r} does not supply list {name!r}"
    value = inputs[name]
    assert isinstance(value, list), f"{name!r} is not a list in test {case['id']!r}"
    return value


def assert_skips_are_explained(spec_: Mapping[str, Any]) -> None:
    """Assert that a skipped test says why it is skipped.

    A skip with no reason is indistinguishable from an oversight, and the whole
    convention of marking a source ``TODO: source needed`` depends on the reason
    travelling with the test.
    """
    for case in spec_["tests"]:
        if case["status"] != "active":
            reason = case.get("skip_reason", "").strip()
            assert reason, f"{spec_['id']}::{case['id']} is skipped with no reason"


def assert_warnings_agree_with_spec(
    spec_: Mapping[str, Any],
    warnings: Iterable[Warning],
    resolve: Callable[[str], float | None],
    context: str,
) -> None:
    """Assert the warnings a result carries agree exactly with the spec's own
    range checks.

    This is the strongest available check on the warning mechanism, and it is
    fully generic: it reads the bounds out of the spec and confirms the
    implementation warned precisely when a warning-severity bound was violated -
    no more, and no fewer.

    It catches the two failure modes that matter. Warning when nothing is wrong is
    warning fatigue, which teaches callers to ignore the field. Failing to warn
    when something is wrong is the dangerous one: an out-of-range result
    indistinguishable from a validated one.
    """
    warnings = tuple(warnings)

    def has(code: str, field: str) -> bool:
        return any(w.code == code and w.field == field for w in warnings)

    # Grouped by (quantity, code), not examined one at a time. A spec may
    # legitimately carry several bounds on one quantity sharing a code, so "bound
    # A was satisfied" only implies "no warning" if no other bound in the group
    # fired.
    fired: dict[tuple[str, str], bool] = {}
    unresolvable: set[str] = set()

    for raw in spec_["valid_range"]:
        quantity = raw["quantity"]
        code = raw.get("code", "OUT_OF_VALID_RANGE")
        value = resolve(quantity)
        if value is None:
            unresolvable.add(quantity)
            continue
        from chemeng.core.range import RangeCheck

        check = RangeCheck.from_spec(raw)
        should_warn = check.severity.value == "warning" and check.violated(value)
        key = (quantity, code)
        fired[key] = fired.get(key, False) or should_warn

    for quantity in sorted(unresolvable):
        assert has(WarningCode.RANGE_CHECK_SKIPPED, quantity), (
            f"{context}: `{quantity}` could not be resolved, so its check did not run, "
            f"but no RANGE_CHECK_SKIPPED warning was emitted - an unchecked value must "
            f"not look like a checked one"
        )

    for (quantity, code), should_warn in fired.items():
        warned = has(WarningCode(code), quantity)
        if should_warn:
            assert warned, (
                f"{context}: `{quantity}` violates a warning-severity spec bound yet no "
                f"{code} warning was emitted"
            )
        else:
            assert not warned, (
                f"{context}: `{quantity}` satisfies every warning-severity spec bound "
                f"for {code}, yet a {code} warning was emitted"
            )


def assert_consistent(result: Any, context: str) -> None:
    """Assert a result is internally consistent for reporting purposes."""
    assert result.is_clean == (not result.warnings), (
        f"{context}: is_clean disagrees with the warning list"
    )
    for warning in result.warnings:
        assert warning.message.strip(), f"{context}: warning {warning.code} has an empty message"


def assert_results_equal(left: Any, right: Any, tolerance: float, context: str) -> None:
    """Assert two results agree, field by field, within a tolerance.

    Used by the cross-implementation tests. Compares in three passes, because the
    fields are of three kinds: scalars within a tolerance, quantities converted to
    base SI before comparing, and the non-numeric fields - regime and warnings -
    which must match exactly.

    Comparing warnings *structurally, including their order*, is deliberate. A
    Rust implementation that forgot to emit RANGE_CHECK_SKIPPED would produce
    identical numbers and a different warning list, and a numeric-only comparison
    would call that a pass.
    """
    assert type(left).__name__ == type(right).__name__, (
        f"{context}: different result types: {type(left).__name__} vs {type(right).__name__}"
    )

    # `dataclasses.fields`, not `vars()`: the result dataclasses use slots,
    # so they have no `__dict__` for `vars()` to read.
    for spec_field in dataclasses.fields(left):
        field = spec_field.name
        a = getattr(left, field)
        b = getattr(right, field)

        if isinstance(a, pint.Quantity) or isinstance(b, pint.Quantity):
            assert isinstance(a, pint.Quantity) and isinstance(b, pint.Quantity), (
                f"{context}.{field}: one side is a bare number and the other a quantity"
            )
            a_base = a.to_base_units().magnitude
            b_base = b.to_base_units().magnitude
            assert_close(float(a_base), float(b_base), tolerance, f"{context}.{field}")
        elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
            assert_close(float(a), float(b), tolerance, f"{context}.{field}")
        elif isinstance(a, tuple) and a and isinstance(a[0], Warning):
            assert _warnings_equal(a, b), (
                f"{context}.{field}: warnings differ\n  python: {a}\n  rust:   {b}"
            )
        elif isinstance(a, tuple) and a and hasattr(a[0], "k"):
            for index, (ca, cb) in enumerate(zip(a, b, strict=True)):
                assert ca.fitting_id == cb.fitting_id, f"{context}.{field}[{index}] id"
                assert_close(ca.n_ld, cb.n_ld, tolerance, f"{context}.{field}[{index}].n_ld")
                assert_close(ca.k, cb.k, tolerance, f"{context}.{field}[{index}].k")
        else:
            assert a == b, f"{context}.{field}: {a!r} != {b!r}"


def _warnings_equal(left: tuple[Warning, ...], right: Any) -> bool:
    """Compare warnings by code, field and order - not by message text.

    The message is prose for a human, and the two implementations format the
    numbers inside it differently: Python renders ``0.0`` and ``1e-06`` where
    Rust renders ``0`` and ``0.000001``. Neither is wrong, and forcing them to
    agree would mean reimplementing IEEE-754 shortest-representation formatting
    identically in both languages - real work, for a difference no caller can
    branch on.

    What a caller *does* branch on is the code, and that is compared exactly,
    along with the field and the order. Message text is required to be non-empty
    on both sides so a missing explanation is still caught.
    """
    right = tuple(right)
    if len(left) != len(right):
        return False
    for a, b in zip(left, right, strict=True):
        if a.code != b.code or a.field != b.field:
            return False
        if not a.message.strip() or not str(b.message).strip():
            return False
    return True


def regime_of(name: str) -> FlowRegime:
    """Look up a flow regime by its string value, for assertions."""
    return FlowRegime(name)


def to_si_(value: Any, spec_unit: str, field: str) -> float:
    """Re-exported for tests that need to convert an expected value."""
    return to_si(value, spec_unit, field)
