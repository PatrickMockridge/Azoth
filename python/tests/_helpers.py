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
from collections.abc import Callable, Iterable, Mapping
from typing import Any

import pint

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
_DIAGNOSTIC_FIELDS: frozenset[str] = frozenset({"residual"})


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

        if isinstance(a, pint.Quantity) or isinstance(b, pint.Quantity):
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
        elif isinstance(a, (tuple, list)) and _is_numeric_nested(a):
            # A vector or matrix output - a composition, a set of K-values, the stationary
            # compositions of a stability trial - compared entry by entry within the
            # tolerance, at any nesting depth. Two implementations agree to the last bit
            # only by luck, so exact equality here reports a rounding difference as a
            # divergence in the physics.
            _assert_nested(a, b, tolerance, f"{context}.{field}")
        elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
            assert_close(float(a), float(b), tolerance, f"{context}.{field}")
        else:
            assert a == b, f"{context}.{field}: {a!r} != {b!r}"


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

    Both sides failing to converge is agreement, so NaN equals NaN.
    """
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

    names = inputs.get(_MODEL_COMPONENTS_INPUT)
    if names is not None:
        takes = parameters_of(model)
        if "mixture" in takes:
            fluid, ideal_gas = databank.mixture_of(list(names))
            kwargs["mixture"] = fluid
            if "ideal_gas" in takes:
                kwargs["ideal_gas"] = ideal_gas
        elif "coeffs" in takes:
            # A non-cubic reference EOS takes the coefficient sets rather than a
            # `Mixture`: `eos.bwrs_phase` is the one, and the MBWR-32 coefficients
            # resolve by name through the same databank, only the last step differs.
            kwargs["coeffs"] = [databank.bwrs_coefficients(name) for name in names]
        elif "components" in takes:
            # The EOS-CG mixture maps its own fixed component names to indices, so the
            # names cross the boundary verbatim rather than resolved to a `Mixture`.
            kwargs["components"] = list(names)
        else:
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

    for name, value in inputs.items():
        if name == _MODEL_COMPONENTS_INPUT:
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
