"""Shared helpers for the spec-driven tests.

Mirrors ``crates/azoth-hydraulics/tests/common/mod.rs``. The tests in this
directory are not hand-written per calc: they walk the spec's ``tests`` list and
execute whatever it declares, so a test added to a spec file runs in both
languages with no new test code. That is the whole point of generating the
registry on both sides.

``python/tests`` is on ``pythonpath`` (see ``pyproject.toml``) so this module is
importable as a top-level ``_helpers``.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

import pint

from azoth._dispatch import fluid_inputs
from azoth._registry_gen import BY_ID, CALCS
from azoth.core.result import FlowRegime
from azoth.core.units import quantity, to_si
from azoth.core.warnings import Warning, WarningCode

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
    "kwargs_for",
    "list_input",
    "model_kwargs",
    "range_checks_that_may_skip",
    "spec",
]


#: Result fields that report a solver's own convergence residual rather than a physical
#: answer. They are the one kind of field whose cross-language agreement is bounded by the
#: model's declared ``algorithm.tolerance`` and not by the case's answer tolerance. The name
#: is the one this spec system uses for that quantity throughout ``specs/models/``.
#:
#: ``tm`` is here for the same reason and is the clearer case of it: a tangent-plane distance
#: at a *trivial* stationary point - a trial that converged to the feed - is zero in exact
#: arithmetic and lands within an ulp of it, at ``~1e-16``, with whatever sign the last
#: floating-point operation left. Measured, `eos.tp_multiflash`'s CO2/methane/nc10 case has
#: one implementation at ``6.66e-16`` and the other at ``-1.78e-15``: a relative comparison
#: calls that a factor-of-three divergence, and the quantity it is measuring is not there.
#:
#: ``balance_error`` is the third and the same thing again: it is `sum_p beta_p x_ip - z_i`,
#: zero **by construction** wherever a model assembles its phases from the same `z` it was
#: given, so what crosses the boundary is the last ulp of the arithmetic that reached it -
#: measured, `eos.hydrate_fraction`'s first case is ``0.0`` in Rust and ``1.11e-16`` in
#: Python. The model's case records the zero because that is the claim; comparing the two
#: kernels to each other at that magnitude compares rounding, not the claim.
_DIAGNOSTIC_FIELDS: frozenset[str] = frozenset(
    {
        "balance_error",
        "residual",
        "tm",
        # `reactions.reactive_phase_equilibrium`'s certificate. Each is a difference of
        # `O(100)` logarithms - `sum(nu_i (ln x_i + ln gamma_i)) - ln K` - whose own value is
        # of order `1e-8` on a converged solve, so a rounding path difference of `1e-12`
        # between the two kernels is a *hundred percent* of what is being compared. They are
        # diagnostics and are compared against the bound the solver declared, which is what
        # makes the comparison a measurement rather than an impossible strictness.
        "max_reaction_log_residual",
        "net_charge_moles",
        # `reactive_hybrid_eos_ge_flash`'s, and the same quantity one level out: a net charge
        # left over after the reactions have conserved it exactly. Measured, the two kernels
        # leave `1.65e-24` and `4.96e-24` - two roundings of a number that is zero, whose
        # *ratio* is a factor of six. The class's own gate is `1e-8` and the case records it.
        "charge_residual",
        "max_element_residual",
        # `eos.hybrid_eos_ge_flash`'s two. Both are of order `1e-13`: one is a material
        # balance that is zero **by construction** and carries only the last ulp of the two
        # kernels' arithmetic, the other a cross-role fugacity spread the same size. The
        # contract's `1e-7` and `1e-5` are what the model claims and the case records them
        # as such; between the two kernels what is being compared is rounding.
        "max_material_balance_residual",
        "max_log_fugacity_residual",
        # `reactions.reactive_hybrid_eos_ge_flash`'s. It is the quantity the coupled loop
        # *stops on* - the sum of `|x_old - x_new|` over the brine at the last pass - so the
        # two kernels land either side of `1e-13` for the same reason an `error` does, and a
        # relative comparison of two numbers that small measures the rounding that decided
        # where each stopped. The loop's other stopping quantity, `residual`, is here already.
        "chemical_deviation",
        # `reactions.reactive_tp_flash`'s two. The element residual is a scaled
        # root-mean-square of `A n - b`, so on a converged solve it is of order `1e-6` and a
        # difference in the last few digits of the arithmetic is a large fraction of it -
        # measured, the water-gas shift comes out `9.7e-7` in Rust and `1.7e-6` in Python
        # while the compositions agree to `1e-7`. `residual` is `max` of that and the worst
        # potential error, and is the same kind of quantity one level along.
        "element_residual",
    }
)

#: Fields that are **not compared across implementations at all**, because their value is a
#: property of the arithmetic rather than of the model. An iteration count is the clearest
#: case: two implementations of the same loop and the same stopping rule reach the same
#: answer in a different number of steps whenever one of them rounds differently on the way,
#: and `eos.tp_flash_saft`'s do - seven steps against fourteen, on the same state, to the same
#: K-values. The tolerance the loop stops at is asserted through `residual`, which carries
#: the bound the solver declared.
#:
#: `error` is the second, and for the same reason one level along: it is the quantity an
#: iteration *stops on*, so it reports where the loop is when the two codes part rather
#: than what they were computing. `reactions.chemical_equilibrium` is the case that
#: needed it - the two kernels converge to the same composition, one pass apart, and
#: report errors five times apart because that is where each stopped. **The value is still
#: pinned against NeqSim**, by the model's own case; this only drops the comparison of the
#: two kernels to each other.
#:
#: `converged` is the third and is the verdict *on* `error`: it is true when that quantity
#: came in under a bound the loop relaxes by 1.5x after its fifteenth pass, so where the
#: two codes stop decides it and nothing else does. On `reactions.reactive_phase_equilibrium`'s
#: first case the reference reports `2.27e-8` and the extension `3.72e-8` against a relaxed
#: bound of `3.375e-8`, and the flag crosses with it. **It is still asserted** - by the
#: model's own crate test, per case, and by `reactions.chemical_equilibrium`'s.
#:
#: `reactive_tp_flash`'s `total_iterations` is the fourth and is the same quantity under
#: NeqSim's name for it: the passes every inner solve took, summed over the outer loop. The
#: two kernels take 35 and 37 of them on the water-gas shift and agree on the composition to
#: `1e-7`, which is what a path quantity looks like - the loop it counts is not the answer.
#:
#: `reactive_ph_flash`'s two are the fifth and sixth, and they are the same quantity one level
#: out: the secant's temperature steps and every inner flash's passes summed. On the water-gas
#: shift the two kernels take 21 and 11 outer passes to temperatures that agree to `0.06 K` on
#: `600 K`, because a secant's path depends on both enthalpy curves and those are built by
#: different implementations. The spec says so and the case reports the counts rather than
#: pinning them.
_UNCOMPARED_FIELDS: frozenset[str] = frozenset(
    {
        "iterations",
        "error",
        "converged",
        "total_iterations",
        "outer_iterations",
        "total_inner_iterations",
        # **`distillation_column`'s three closure measures, which are where each solve stopped
        # rather than what it converged to.** The first is a mean tray-temperature change under
        # the substitution core and the scaled MESH norm under the mesh solve; the other two are
        # the products' imbalances, and all three are numbers of order `1e-9` against a feed. On
        # the mesh case's row the two implementations report `9.2816e-9` and `9.2872e-9` for the
        # first and `8.385e-11` and `8.409e-11` for the second, so the relative agreement is the
        # *seventh* digit of a quantity that is already a residual. The crate's own column tests
        # assert all three per case against the capture, which is where the claim belongs; the
        # case's tolerance is set by the profile and the duties, which agree to `2e-5`.
        "temperature_residual",
        "mass_residual",
        "energy_residual",
        # `reactive_hybrid_eos_ge_flash`'s pass count, which is the clearest case of the rule
        # above: **the two kernels take three passes and four on the three-phase fluid**, and
        # the fourth moves the brine's bicarbonate `2e-7`. The loops stop on a composition
        # deviation that is the difference of two nearly equal numbers, so a trajectory
        # difference of `1e-8` decides it. Each kernel's own test asserts the count against the
        # capture where the two agree, which is where the claim belongs.
        "passes",
    }
)


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
        from azoth.core.range import RangeCheck

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


def assert_results_equal(
    left: Any,
    right: Any,
    tolerance: float,
    context: str,
    diagnostic_bound: float | None = None,
) -> None:
    """Assert two results agree, field by field, within a tolerance.

    Used by the cross-implementation tests. Compares in three passes, because the
    fields are of three kinds: scalars within a tolerance, quantities converted to
    base SI before comparing, and the non-numeric fields - regime and warnings -
    which must match exactly.

    Comparing warnings *structurally, including their order*, is deliberate. A
    Rust implementation that forgot to emit RANGE_CHECK_SKIPPED would produce
    identical numbers and a different warning list, and a numeric-only comparison
    would call that a pass.

    ``diagnostic_bound`` is the absolute bound applied to a solver diagnostic (see
    ``_DIAGNOSTIC_FIELDS``). It is the model's own declared convergence tolerance,
    not the case's answer tolerance: a residual is meaningful only to the precision
    at which the solver declared it converged, and comparing it relatively is not a
    stricter test but an impossible one.
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

        if field in _UNCOMPARED_FIELDS:
            # One side has to have produced *something*; that the loop ran is the claim.
            assert a is not None and b is not None, f"{context}.{field}: absent on one side"
            continue
        if field in _MAY_AGREE_ON_NAN and _both_nan(a, b):
            # NeqSim's own readers return NaN where there is nothing to read, and the two
            # kernels agreeing on that is the claim rather than a failure to make one.
            continue
        # **A diagnostic is checked before the kind of value it is.** These are the fields
        # whose value is a *residue*, and one of them carries a unit: the relative
        # comparison the quantity branch applies would be reading a cancellation's last
        # digits as if they were the answer.
        if field in _DIAGNOSTIC_FIELDS and diagnostic_bound is not None:
            # A diagnostic may carry a unit - a residual in moles - and the bound it is
            # compared against is a bare magnitude in the same unit, so the quantities are
            # stripped here rather than after the conversion the other branch does.
            if isinstance(a, pint.Quantity) or isinstance(b, pint.Quantity):
                assert isinstance(a, pint.Quantity) and isinstance(b, pint.Quantity), (
                    f"{context}.{field}: one side is a bare number and the other a quantity"
                )
                a = a.to_base_units().magnitude
                b = b.to_base_units().magnitude
            _assert_diagnostic(a, b, diagnostic_bound, f"{context}.{field}")
        elif isinstance(a, pint.Quantity) or isinstance(b, pint.Quantity):
            assert isinstance(a, pint.Quantity) and isinstance(b, pint.Quantity), (
                f"{context}.{field}: one side is a bare number and the other a quantity"
            )
            a_base = a.to_base_units().magnitude
            b_base = b.to_base_units().magnitude
            assert_close(float(a_base), float(b_base), tolerance, f"{context}.{field}")
        elif field in _DIAGNOSTIC_FIELDS and diagnostic_bound is not None:
            _assert_diagnostic(a, b, diagnostic_bound, f"{context}.{field}")
        elif isinstance(a, tuple) and a and isinstance(a[0], Warning):
            assert _warnings_equal(a, b), (
                f"{context}.{field}: warnings differ\n  python: {a}\n  rust:   {b}"
            )
        elif isinstance(a, tuple) and a and hasattr(a[0], "k"):
            for index, (ca, cb) in enumerate(zip(a, b, strict=True)):
                assert ca.fitting_id == cb.fitting_id, f"{context}.{field}[{index}] id"
                assert_close(ca.n_ld, cb.n_ld, tolerance, f"{context}.{field}[{index}].n_ld")
                assert_close(ca.k, cb.k, tolerance, f"{context}.{field}[{index}].k")
        elif isinstance(a, (tuple, list)) and a and isinstance(a[0], pint.Quantity):
            # A *vector* output whose entries carry units - a list of saturation
            # pressures, a surface tension per component. Compared entry by entry in base
            # SI, for the reason the scalar branch above converts: the two sides may state
            # the same pressure in different units.
            assert len(a) == len(b), f"{context}.{field}: {len(a)} entries against {len(b)}"
            for index, (qa, qb) in enumerate(zip(a, b, strict=True)):
                assert isinstance(qb, pint.Quantity), (
                    f"{context}.{field}[{index}]: one side is a bare number and the other "
                    f"a quantity"
                )
                assert_close(
                    float(qa.to_base_units().magnitude),
                    float(qb.to_base_units().magnitude),
                    tolerance,
                    f"{context}.{field}[{index}]",
                )
        elif isinstance(a, (tuple, list)) and _is_numeric_nested(a):
            # A vector or matrix output - a composition, a set of K-values, the stationary
            # compositions of a stability trial - compared entry by entry within the
            # tolerance, at any nesting depth. Two implementations agree to the last bit
            # only by luck, so exact equality here reports a rounding difference as a
            # divergence in the physics.
            _assert_nested(a, b, tolerance, f"{context}.{field}")
        elif isinstance(a, (tuple, list)) and _is_quantity_nested(a):
            # A **matrix of quantities** - `reactions.reactive_tp_flash`'s per-phase mole
            # numbers are the first - which is neither of the two branches above: the
            # vector branch expects quantities one level down and the numeric one expects
            # plain numbers at every level. Each entry is converted to base SI and
            # compared within the tolerance, as the vector branch does, because the two
            # sides may state the same amount in different units.
            assert isinstance(b, (tuple, list)) and len(a) == len(b), (
                f"{context}.{field}: {len(a)} rows against "
                f"{len(b) if isinstance(b, (tuple, list)) else type(b).__name__}"
            )
            for index, (row_a, row_b) in enumerate(zip(a, b, strict=True)):
                assert isinstance(row_b, (tuple, list)) and len(row_a) == len(row_b), (
                    f"{context}.{field}[{index}]: one side is not a row"
                )
                for column, (qa, qb) in enumerate(zip(row_a, row_b, strict=True)):
                    assert isinstance(qb, pint.Quantity), (
                        f"{context}.{field}[{index}][{column}]: one side is a bare number "
                        f"and the other a quantity"
                    )
                    assert_close(
                        float(qa.to_base_units().magnitude),
                        float(qb.to_base_units().magnitude),
                        tolerance,
                        f"{context}.{field}[{index}][{column}]",
                    )
        elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
            assert_close(float(a), float(b), tolerance, f"{context}.{field}")
        else:
            assert a == b, f"{context}.{field}: {a!r} != {b!r}"


#: Fields where two NaNs are an agreement rather than a failure.
#:
#: `reactions.reactive_phase_equilibrium`'s three residuals are NeqSim's own readers, and each
#: returns `NaN` when there is no reactive phase to read: a skipped phase has nothing to
#: certify. One NaN against a number is still a failure - that is one kernel computing what
#: the other refused - and a case pinning a number catches a NaN wherever both should have one.
_MAY_AGREE_ON_NAN = frozenset(
    {"max_reaction_log_residual", "net_charge_moles", "max_element_residual"}
)


def _both_nan(a: Any, b: Any) -> bool:
    """Whether both sides are NaN, however they are wrapped.

    A quantity, a bare number or a vector of either: the three fields this is used for are
    scalars in both languages, and the vector case is here so that the rule cannot be
    quietly wrong if one of them ever carries a vector.
    """

    def one(value: Any) -> bool:
        if isinstance(value, pint.Quantity):
            value = value.to_base_units().magnitude
        if isinstance(value, (int, float)):
            return math.isnan(float(value))
        if isinstance(value, (tuple, list)):
            return all(one(entry) for entry in value)
        return False

    return one(a) and one(b)


def _is_quantity_nested(value: Any) -> bool:
    """Whether a sequence is sequences of quantities, at any depth below one row.

    A matrix of amounts: rows of quantities. Kept apart from `_is_numeric_nested`, which
    would answer ``False`` for every quantity - a pint quantity is not an ``int`` or a
    ``float`` - so without this branch a matrix of units falls through to exact equality
    and a rounding difference between the two kernels reads as a divergence.
    """
    if not isinstance(value, (tuple, list)) or not value:
        return False
    for row in value:
        if not isinstance(row, (tuple, list)):
            return False
        if not all(isinstance(item, pint.Quantity) for item in row):
            return False
    return True


def _is_numeric_nested(value: Any) -> bool:
    """Whether a sequence is numbers all the way down, at any depth.

    A vector is one level; a matrix is two. Both are compared entry by entry rather
    than by equality, so both have to be recognised here.
    """
    if isinstance(value, (tuple, list)):
        return len(value) == 0 or all(_is_numeric_nested(item) for item in value)
    return isinstance(value, (int, float))


def _assert_nested(a: Any, b: Any, tolerance: float, context: str) -> None:
    """Compare two nested numeric sequences entry by entry."""
    assert len(a) == len(b), f"{context}: {len(a)} entries against {len(b)}"
    for index, (ca, cb) in enumerate(zip(a, b, strict=True)):
        if _is_numeric_nested(ca) and isinstance(ca, (tuple, list)):
            _assert_nested(ca, cb, tolerance, f"{context}[{index}]")
        else:
            assert_close(float(ca), float(cb), tolerance, f"{context}[{index}]")


def _assert_diagnostic(a: Any, b: Any, bound: float, context: str) -> None:
    """Compare a solver diagnostic on an absolute bound.

    A residual is not a result, it is the solver's evidence that it converged, and it
    is meaningful only down to the tolerance the solver declares. It is near zero by
    construction, so a relative comparison divides by the thing that is vanishing:
    for ``residual = |S - 1|`` with ``S`` near 1, the subtraction pins the absolute
    error at the rounding of ``S`` (~1e-16) however small the residual becomes, and
    the attainable relative accuracy ``1e-16 / residual`` diverges as the residual
    falls. No relative tolerance can be met here, at any value.

    A diagnostic can also be a *vector* - a tangent-plane distance per trial - and
    then the bound applies to each entry, because each entry is one solver's evidence
    about one trial and the two trials are independent iterations.

    Both sides failing to converge is agreement, so NaN equals NaN.
    """
    if isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        assert len(a) == len(b), (
            f"{context}: {len(a)} diagnostic(s) against {len(b)} ({a!r} against {b!r})"
        )
        for index, (ca, cb) in enumerate(zip(a, b, strict=True)):
            _assert_diagnostic(ca, cb, bound, f"{context}[{index}]")
        return
    a_nan = a != a
    b_nan = b != b
    if a_nan or b_nan:
        assert a_nan and b_nan, (
            f"{context}: one implementation did not converge and the other did "
            f"({a!r} against {b!r})"
        )
        return
    difference = abs(float(a) - float(b))
    assert difference <= bound, (
        f"{context}: residuals differ by {difference:.3e}, beyond the model's own "
        f"declared convergence tolerance {bound:.3e} ({a!r} against {b!r})"
    )


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


def kwargs_for(calc: Mapping[str, Any], inputs: Mapping[str, Any]) -> dict[str, Any]:
    """Turn a spec's declared inputs into keyword arguments for the real call.

    Both the type and the unit come from the spec's own ``inputs`` block, so a
    validation case or a test case needs no unit handling of its own, and adding
    a calc needs no change here.

    Note that ``type`` is absent for a plain quantity - it defaults to
    ``quantity``, and the specs only spell it out for the exceptions. Treating a
    missing ``type`` as anything other than the default passes bare floats where
    quantities are required, which is exactly the mistake this library exists to
    make impossible, and which the reference implementation will reject.
    """
    declared = calc["inputs"]
    kwargs: dict[str, Any] = {}
    for name, value in inputs.items():
        declaration = declared[name]
        kind = declaration.get("type", "quantity")
        if kind == "fitting_list":
            kwargs[name] = list(value)
        elif kind in ("enum", "string"):
            kwargs[name] = value
        elif kind == "boolean":
            # A flag, not a number in any unit. `bool(value)` rather than `float(value)`:
            # `reactions.chemical_equilibrium` is the first calc whose spec declares one,
            # and a `0.0` reaching a `bool` parameter is refused by the extension.
            kwargs[name] = bool(value)
        elif kind == "quantity" and declaration.get("unit") != "dimensionless":
            kwargs[name] = quantity(float(value), declaration["unit"])
        else:
            kwargs[name] = float(value)
    return kwargs


#: The declared input a model's ``mixture`` and ``ideal_gas`` arguments are assembled
#: from. It is not a parameter of the model, it is what those two parameters are made of.
#:
#: One input, where there used to be nine. `Tc`, `Pc`, `omega`, `kij` and `cp_a`..`cp_e`
#: reached a model through `_MODEL_MIXTURE_INPUTS` and `_MODEL_IDEAL_GAS_INPUTS`, two
#: tuples that had to be kept in step with the specs by hand and with nothing checking
#: they still matched. `components` is a list of names, and `mixture_of` resolves it
#: against the databank - the same call, reading the same file, in both languages.
_MODEL_COMPONENTS_INPUT = "components"
#: The key a case uses to say its fluid runs the Wertheim association. A boundary key:
#: the runner builds the fluid with it, and the model is handed the mixture whole.
_ASSOCIATING_KEY = "associating"


def range_checks_that_may_skip(spec_: Mapping[str, Any]) -> set[str]:
    """Quantities whose bounds a case may legitimately leave unevaluated.

    A check that depends on an *optional* input cannot always run, and the spec says so by
    naming that input in ``computed_from``. That is the only licence the two
    range-check tests grant: everything else is a bound that reads as validation and
    performs none.
    """
    optional = {name for name, declared in spec_["inputs"].items() if declared.get("optional")}
    allowed: set[str] = set()
    for check in spec_.get("valid_range", []):
        dependencies = set(check.get("computed_from", [])) | {check["quantity"]}
        if dependencies & optional:
            allowed.add(check["quantity"])
    return allowed


def _scalar(unit: str, value: Any) -> Any:
    """One spec number as the argument the implementation takes."""
    return float(value) if unit == "dimensionless" else quantity(float(value), unit)


def _declared(declaration: Mapping[str, Any], value: Any) -> Any:
    """One declared input, at the shape and unit the spec gives it."""
    unit = declaration.get("unit", "dimensionless")
    kind = declaration.get("type", "quantity")
    if kind == "vector":
        return [_scalar(unit, item) for item in value]
    if kind == "matrix":
        return [[_scalar(unit, item) for item in row] for row in value]
    if kind in ("enum", "string"):
        # A symbolic input is the name itself, not a number in any unit. Only
        # `type = "enum"` reaches a function this way; `string` is what a
        # calculation's `source` field is spelled with.
        return str(value)
    if kind == "boolean":
        # A flag, for the reason `kwargs_for` gives.
        return bool(value)
    return _scalar(unit, value)


def model_kwargs(model: Mapping[str, Any], inputs: Mapping[str, Any]) -> dict[str, Any]:
    """Turn a model's declared inputs into keyword arguments for the real call.

    Distinct from :func:`kwargs_for`, because a model's arguments are objects: the spec
    declares `components` and the function takes a `mixture` - and, where an enthalpy is
    involved, an `ideal_gas` as well. This rebuilds them, which is the mirror of what
    `_rust_bridge` does in the other direction.

    **The resolution is `components.mixture_of`, not a local construction.** Given names,
    it reads the databank and applies any keycard override, so a case and a caller take
    the same path to their fluid. Building a `Mixture` here from numbers the case carried
    is what made a case's fluid a second fluid, described by a spec file rather than by
    NeqSim's tables.
    """
    from azoth.eos import components as databank

    declared = model["inputs"]
    kwargs: dict[str, Any] = {}
    # Read once, and before the branch below, because the final pass asks the same
    # question of every model whether or not it names a fluid.
    takes = parameters_of(model)

    for prefix, declared_name in fluid_inputs(model):
        names = inputs.get(declared_name)
        if names is None:
            raise AssertionError(
                f"{model['id']} declares {declared_name!r} and the case does not state it; "
                f"a fluid a case leaves out is a fluid this would build as something else"
            )
        mixture_key = f"{prefix}mixture"
        if mixture_key in takes:
            # A model may declare which cubic it is evaluated under - the gamma-phi
            # flash does, because its vapour is half of what it is. The declaration is
            # what builds the mixture, so a case states its own vapour rather than
            # inheriting whatever `mixture_of` defaults to.
            #
            # **`associating` is the same kind of statement**, and without it no case in
            # this registry could exercise an associating model at all: the runner would
            # build every fluid as a classical one, and a model whose whole behaviour
            # changes with the Wertheim term would be tested against the wrong fluid while
            # answering a plausible one. It is the model's decision - the same methanol
            # and water are a classical mixture under an SRK model and an associating one
            # under a CPA model - so a case that wants it says so.
            fluid, ideal_gas = databank.mixture_of(
                list(names),
                eos=inputs.get(f"{prefix}eos", "pr"),
                associating=bool(inputs.get(f"{prefix}associating", False)),
            )
            kwargs[mixture_key] = fluid
            ideal_gas_key = f"{prefix}ideal_gas"
            if ideal_gas_key in takes:
                kwargs[ideal_gas_key] = ideal_gas
        elif "coeffs" in takes:
            # A non-cubic reference EOS takes the coefficient sets rather than a
            # `Mixture`: `eos.bwrs_phase` is the one, and the MBWR-32 coefficients
            # resolve by name through the same databank, only the last step differs.
            kwargs["coeffs"] = [databank.bwrs_coefficients(name) for name in names]
        elif f"{prefix}components" in takes:
            # Either the EOS-CG mixture, which maps its own fixed names to indices, or a
            # process model, whose Rust side resolves the names through the same databank
            # `Stream::mixture()` does. Both take the names verbatim rather than a
            # `Mixture`, so the two languages cannot disagree about which row answered.
            kwargs[f"{prefix}components"] = list(names)
        elif "params" not in takes:
            # A pure-component model takes the constants themselves rather than a
            # `Mixture`: `eos.pure_saturation` is the one, and a saturation pressure is
            # a property of one substance. The names still go through the databank -
            # only the last step differs.
            if len(names) != 1:
                raise AssertionError(
                    f"{model['id']} takes scalar critical constants but its case names "
                    f"{len(names)} components; the case is wrong, not the model"
                )
            record = databank.entry(names[0])
            kwargs.update(Tc=record.Tc, Pc=record.Pc, omega=record.omega)
        if "params" in takes:
            # An activity model takes a resolved parameter set of its own: NRTL's
            # `alpha`/`Dij` and UNIFAC's group tables are not the critical constants a
            # `Mixture` carries, so the names resolve to the model's own record instead.
            # Resolved *beside* the mixture rather than instead of it, because the
            # gamma-phi flash takes both - one cubic vapour, one activity-model liquid.
            kwargs["params"] = parameter_set(model, names, inputs)

    fluids = fluid_inputs(model)
    # Every input the loop above consumed: the fluids' own names, and the per-fluid
    # boundary keys a case uses to describe how each one is built.
    resolved = {name for _, name in fluids}
    resolved |= {_ASSOCIATING_KEY, "eos"}
    for prefix, _ in fluids:
        resolved |= {f"{prefix}{_ASSOCIATING_KEY}", f"{prefix}eos"}
    for name, value in inputs.items():
        # Passed only what the function takes. A fluid's own inputs are resolved above -
        # into a `Mixture`, a coefficient set, or the names verbatim - and a parameter
        # set's extra inputs were consumed there too, so either would be a second, wrong
        # answer or a keyword the function does not have. `tools/gen_stub.py` asks the same
        # question of the same signature, which is what keeps the boundary and this in step.
        if name in resolved or name not in takes:
            continue
        kwargs[name] = _declared(declared[name], value)
    return kwargs


def parameters_of(model: Mapping[str, Any]) -> set[str]:
    """The parameter names a model's function actually takes.

    Read from the implementation rather than from the spec's prose. The spec declares
    *inputs* - `components`, `T`, `P` - and the function takes *objects*: a `mixture`,
    and an `ideal_gas` where an enthalpy is involved. Which objects is a fact about the
    signature, so that is where this asks, and a signature that changes fails here
    instead of silently passing the wrong argument.
    """
    import importlib
    import inspect

    namespace, _, name = model["id"].partition(".")
    module = importlib.import_module(f"azoth.{namespace}")
    function = getattr(module, name, None)
    if function is None:  # pragma: no cover - the registry contract covers this
        raise AssertionError(f"{model['id']}: azoth.{namespace} has no `{name}`")
    return set(inspect.signature(function).parameters)


#: The databank function that resolves each parameter-set record a model's ``params``
#: argument can be annotated with, keyed by the annotation's name.
#:
#: A record and the function that builds it are two names for one thing, and this is
#: where they are held together. Deriving `unifac_parameters` from `UnifacParameters`
#: by string would be a convention nothing states; this is two lines a reader can
#: check, and a model whose record is not here fails loudly rather than quietly.
#:
#: Each entry is the resolving function's name and the *other declared inputs* it
#: consumes. A record built from the component names alone needs none; UNIFAC-UMR-PRU's
#: needs the parameter set its spec declares, because NeqSim chooses that from a
#: component field this library does not carry and the caller states it instead. Those
#: inputs are consumed here and do **not** also reach the model function, which is why
#: the tuple is part of the entry rather than left to the caller to notice.
PARAMETER_RESOLVERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "GeNrtlPhaseParameters": ("ge_nrtl_phase_parameters", ()),
    "GeUnifacPhaseParameters": ("ge_unifac_phase_parameters", ()),
    "GeUniquacPhaseParameters": ("ge_uniquac_phase_parameters", ()),
    "GeVanLaarAcidPhaseParameters": ("ge_van_laar_acid_phase_parameters", ()),
    "GeWilsonPhaseParameters": ("ge_wilson_phase_parameters", ()),
    "NrtlParameters": ("nrtl_parameters", ()),
    "UnifacParameters": ("unifac_parameters", ()),
    "UnifacPsrkParameters": ("unifac_psrk_parameters", ()),
    "UnifacUmrpruParameters": ("unifac_umrpru_parameters", ("parameters",)),
    "UniquacParameters": ("uniquac_parameters", ()),
    "VanLaarAcidParameters": ("van_laar_acid_parameters", ()),
}


def parameter_set(model: Mapping[str, Any], names: Sequence[str], inputs: Mapping[str, Any]) -> Any:
    """The resolved parameters a model's ``params`` argument is annotated with.

    Which record is a fact about the annotation, so that is where this asks - the same
    move :func:`parameters_of` makes one level up, and what keeps the two from being
    able to disagree: a model whose ``params`` is annotated with a record this has
    never heard of is an error rather than a quietly wrong call.
    """
    import importlib
    import inspect

    from azoth.eos import components as databank

    namespace, _, name = model["id"].partition(".")
    module = importlib.import_module(f"azoth.{namespace}")
    annotation = inspect.signature(getattr(module, name)).parameters["params"].annotation
    try:
        resolver, extras = PARAMETER_RESOLVERS[annotation]
        resolve = getattr(databank, resolver)
    except KeyError:
        raise AssertionError(
            f"{model['id']}: `params` is annotated {annotation!r}, which is not a "
            f"parameter set this knows how to resolve; known: {sorted(PARAMETER_RESOLVERS)}"
        ) from None
    return resolve(list(names), **{extra: inputs[extra] for extra in extras})


def convergence_tolerance(spec_: dict[str, Any]) -> float | None:
    """The tolerance at which the spec says its solver converged, if it declares one.

    This is the bound a solver diagnostic is compared on. It is the spec's own number,
    so it moves with the model rather than being chosen to make a test pass.

    A *model* declares it under ``algorithm``; a *calculation* that solves an implicit
    equation declares it under ``solver``. Both name the same thing - the precision at
    which the iteration stopped - so a calc's residual is compared on its own bound
    too, rather than relatively (a residual near zero has no meaningful relative
    error, which ``_helpers._assert_diagnostic`` documents).
    """
    for block in ("algorithm", "solver"):
        value = spec_.get(block)
        if isinstance(value, dict) and value.get("tolerance") is not None:
            return float(value["tolerance"])
    return None
