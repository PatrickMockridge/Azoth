"""The batch API, checked against the scalar API it loops over.

The batch layer's whole claim is that it is *the same calculation*, evaluated N times
with one boundary crossing. Nothing about that claim is self-evident from the code, and
a batch API that quietly disagreed with its scalar form would be the worst kind of
failure this project recognises: plausible numbers, from a path nobody cross-checked.

So the three properties below are the milestone's acceptance criteria, and each one
targets a different way the claim could be false:

* **batch-of-1 equals scalar**, field by field *and warning by warning*, on both
  backends. Catches a batch arm that reads the wrong argument, converts a unit twice, or
  drops a warning.
* **a mixed batch agrees element-wise.** Catches an arm that is right for one element and
  wrong for another - which a single-element test cannot see - and the specific case of an
  element that raises, where the policy is fail-fast.
* **cross-language batch agreement.** Catches the two Rust and Python arms disagreeing
  about a column's name, unit, kind or values.

Properties 1 and 2 run on whichever backend is selected, so under CI's
``AZOTH_REQUIRE_RUST=1`` they run on Rust; the cross-language test runs both in one
process, which is where `use_backend` earns its keep.

The cases come from the specs, not from a list here, so a test added to a spec file is
covered in batch as well as in scalar with no new code.
"""

from __future__ import annotations

import dataclasses
import importlib
import math
from collections.abc import Iterator, Mapping
from enum import Enum
from typing import Any

import pytest

import _helpers as h
from azoth._dispatch import available, resolve, use_backend
from azoth.batch._core import batchable
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import FlowRegime
from azoth.core.warnings import WarningCode

pytestmark = pytest.mark.requires_rust


def _batch_function(calc_id: str) -> Any:
    """The batch callable for a calc id, or a failure naming the missing arm.

    Looked up by the same namespace-and-function derivation the dispatcher uses, so a
    calc that resolves in scalar resolves here too. A `getattr` that misses is the
    milestone's own failure mode - a calc the registry knows and the batch layer does not
    - so the error says which one and where to add it.
    """
    namespace, _, function_name = calc_id.rpartition(".")
    module = importlib.import_module(f"azoth.batch.{namespace}")
    assert hasattr(module, function_name), (
        f"{calc_id} is batchable but azoth/batch/{namespace}.py has no {function_name}. "
        f"Add one, or explain in `_require_batchable` why it cannot have one."
    )
    return getattr(module, function_name)


def _scalar_base(value: Any) -> float | None:
    """A scalar result field as an SI base magnitude, or `None` for an absent one.

    Mirrors what the batch column holds, so the two are comparable. `bool` goes through
    `float` and lands on 1.0 or 0.0, which is what both batch arms do with a flag.
    """
    if value is None:
        return None
    if hasattr(value, "magnitude"):
        return float(value.to_base_units().magnitude)
    return float(value)


def _label(value: Any) -> str | None:
    """A scalar result field as an enum label, or `None` for an absent one.

    Any `Enum`, not `FlowRegime` by name - the same generalisation `_is_enum_hint` in
    the batch core needed once `eos.pr_z_factor`'s `root_structure` became an enum
    output outside the hydraulics namespace. `str()` on a `StrEnum` happens to give
    the value, but that is a property of `StrEnum` rather than of the rule, and a
    plain `Enum` would have come back as `RootStructure.THREE_ROOTS`.
    """
    if value is None:
        return None
    return value.value if isinstance(value, Enum) else str(value)


def _column_kinds(calc_id: str) -> dict[str, str]:
    """Reach into the batch core for the column rule, so the test cannot restate it."""
    from azoth._registry_gen import BY_ID
    from azoth.batch._core import _column_kinds as kinds

    return kinds(BY_ID[calc_id])


def assert_batch_matches_scalars(
    result: Any,
    scalars: list[Any],
    calc_id: str,
    tolerance: float,
    context: str,
) -> None:
    """Assert a batch result equals a list of scalar results, element by element.

    Compared per field rather than by dumping both to numbers, because the field names
    are the contract: a batch arm that produced the right *values* under the wrong column
    name would pass a positional comparison and break every caller.
    """
    assert len(result) == len(scalars), f"{context}: {len(result)} vs {len(scalars)} results"
    kinds = _column_kinds(calc_id)

    for field, kind in kinds.items():
        column = getattr(result, field)
        assert len(column) == len(scalars), f"{context}.{field}: column length"

        for index, scalar in enumerate(scalars):
            value = getattr(scalar, field)
            where = f"{context}.{field}[{index}]"
            if kind == "labels":
                assert column[index] == _label(value), (
                    f"{where}: got {column[index]!r}, scalar says {_label(value)!r}"
                )
                continue

            actual, want = column[index], _scalar_base(value)
            if want is None:
                # Absent in the scalar API, `NaN` in a numeric column. NaN rather than a
                # label so the column's kind cannot change with an argument, and NaN
                # rather than zero because zero is a physical claim.
                assert math.isnan(actual), f"{where}: scalar has no value, batch says {actual}"
            else:
                h.assert_close(actual, want, tolerance, where)

        if kind == "values":
            assert result.units[field] == _spec_unit(calc_id, field), (
                f"{context}.{field}: unit map says {result.units.get(field)!r}"
            )

    for index, scalar in enumerate(scalars):
        assert h._warnings_equal(result.warnings[index], scalar.warnings), (
            f"{context}.warnings[{index}]: batch and scalar disagree\n"
            f"  batch:  {result.warnings[index]}\n  scalar: {scalar.warnings}"
        )


def _spec_unit(calc_id: str, field: str) -> str:
    """The unit the spec declares for an output, or `dimensionless`."""
    from azoth._registry_gen import BY_ID

    declaration = BY_ID[calc_id]["outputs"].get(field)
    unit = None if declaration is None else declaration.get("unit")
    return "dimensionless" if unit is None else str(unit)


def _active_cases() -> Iterator[tuple[dict[str, Any], dict[str, Any]]]:
    """Every runnable case of every batchable calc, taken from the specs.

    Only cases carrying `inputs`: the specs also declare `property` cases, which assert an
    invariant rather than a value and have no inputs to array. They are covered by
    `test_properties.py` on the scalar path, and a batch call is the scalar call N times,
    so repeating them here would test the same arithmetic a third time.
    """
    covered = set(batchable())
    for calc in h.CALCS:
        if calc["id"] not in covered:
            continue
        for case in h.all_tests(calc):
            if case["status"] == "active" and "inputs" in case:
                yield calc, case


def _tolerance(case: Mapping[str, Any]) -> float:
    return float(case.get("tolerance", 1e-9))


def test_every_batchable_calc_has_a_batch_arm() -> None:
    """Coverage, as a test rather than as a promise.

    The coverage is what stops a new calc shipping scalar-only: without this, the
    extension raises `NotImplementedError` at call time and nothing in the build notices.
    Deleting any one arm makes this fail naming the calc.
    """
    for calc_id in batchable():
        _batch_function(calc_id)


def test_the_excluded_set_is_exactly_the_unbatchable_calcs() -> None:
    """The calcs with no batch form, asserted so the exclusion stays a decision.

    `crane_k_factors` takes `fittings`, a list of registry ids the spec declares with no
    unit. There is no column shape for a per-element list of names, and sharing one list
    across the batch would compute the same answer N times. `eos.antoine_vapor_pressure`
    takes `form`, a categorical enum input with no unit, which has no float-array column
    shape either. Stating both here means a future calc with a categorical input lands in
    this assertion and has to be argued for, rather than quietly joining the exclusion.
    """
    excluded = {calc["id"] for calc in h.CALCS} - set(batchable())
    assert excluded == {"hydraulics.crane_k_factors", "eos.antoine_vapor_pressure"}, sorted(
        excluded
    )


def test_batch_of_one_equals_scalar_including_warnings() -> None:
    """Property 1: one element through the batch path is the scalar path.

    Every active case of every batchable calc, on the selected backend. The cases come
    from the specs, so this covers the worked example, the reference cases, the
    out-of-range-warning cases and the boundary cases without listing any of them.
    """
    checked = 0
    for calc, case in _active_cases():
        batch = _batch_function(calc["id"])
        arrays = {name: [float(value)] for name, value in case["inputs"].items()}
        # `fittings` is the only non-numeric input and crane is not batchable, so every
        # remaining case is a plain mapping of names to numbers.
        result = batch(**arrays)
        scalar = resolve(calc["id"])(**h.kwargs_for(calc, case["inputs"]))
        assert_batch_matches_scalars(
            result, [scalar], calc["id"], _tolerance(case), f"{calc['id']}::{case['id']}"
        )
        checked += 1
    assert checked, "no active cases found; the spec walk is broken, not passing"


def test_a_mixed_batch_agrees_element_wise() -> None:
    """Property 2: a batch of deliberately unlike elements agrees with scalar per element.

    Built from three elements of the Reynolds/Darcy pair specifically, because it is the
    one calc in the tree with a *warning-severity* band (2000 < Re < 4000) and an optional
    input, so one batch can hold a clean element, a transitional one, and one whose range
    check could not run. An arm that read the wrong index, or reused the first element's
    regime, would pass a batch-of-1 test and fail here.
    """
    calc = h.spec("hydraulics.reynolds_number")
    batch = _batch_function("hydraulics.reynolds_number")
    # Re = 998 * v * 0.05 / 1e-3 = 49900 * v, so these three sit either side of the
    # 2000 and 4000 boundaries the spec declares: 1497, 2994 and 99800.
    cases = [
        {"rho": 998.0, "v": 0.03, "D": 0.05, "mu": 1e-3},  # laminar, clean
        {"rho": 998.0, "v": 0.06, "D": 0.05, "mu": 1e-3},  # transitional, warns
        {"rho": 998.0, "v": 2.0, "D": 0.05, "mu": 1e-3},  # turbulent, clean
    ]
    arrays = {name: [case[name] for case in cases] for name in cases[0]}
    result = batch(**arrays)
    scalars = [resolve("hydraulics.reynolds_number")(**h.kwargs_for(calc, case)) for case in cases]

    assert_batch_matches_scalars(
        result, scalars, "hydraulics.reynolds_number", 1e-12, "mixed reynolds"
    )

    # And the mix is what it was meant to be, so the test cannot pass by accident on three
    # elements that all behave the same way.
    assert result.regime == (FlowRegime.LAMINAR, FlowRegime.TRANSITIONAL, FlowRegime.TURBULENT)
    assert result.warning_indices(WarningCode.TRANSITIONAL_FLOW) == (1,)
    assert not result.is_clean
    assert result.has_warning(WarningCode.TRANSITIONAL_FLOW)


def test_an_error_in_one_element_fails_the_whole_batch() -> None:
    """Fail-fast, and identically on both backends.

    A partial result with silent gaps is what this project is organised against: a caller
    handed 99 values cannot tell which one is missing. So element 2 of 3 being out of
    range raises for the call, and the exception class is the one the scalar API raises -
    not the extension's - so a caller cannot tell which backend answered by what they
    caught.
    """
    batch = _batch_function("hydraulics.pump_power")
    for backend in sorted(available()):
        with use_backend(backend), pytest.raises(OutOfRangeError):
            batch(
                rho=[998.0, 998.0, 998.0],
                q=[0.01, 0.01, 0.01],
                H=[20.0, 20.0, 20.0],
                # eta = 0 makes shaft power singular at eta = 0; the spec bounds it above 0.
                eta=[0.75, 0.75, 0.0],
            )


def test_cross_language_batch_agreement() -> None:
    """Property 3: both backends, one batch, same columns and same warnings.

    Both arms in one process. The columns are compared by name and kind as well as by
    value, because the two sides build their column lists independently - Python from the
    result dataclass's type hints, Rust from a hand-written `match` - and a disagreement
    about a name or a kind would not show up as a numeric difference.
    """
    batch = _batch_function("hydraulics.darcy_weisbach")
    arrays = {
        "f": [0.02, 0.02, 0.03],
        "L": [10.0, 10.0, 25.0],
        "D": [0.05, 0.05, 0.08],
        "rho": [998.0, 998.0, 1000.0],
        "v": [0.05, 2.0, 3.0],
        "mu": [1e-3, 1e-3, 1e-3],
    }

    with use_backend("python"):
        from_python = batch(**arrays)
    with use_backend("rust"):
        from_rust = batch(**arrays)

    assert len(from_python) == len(from_rust) == 3
    kinds = _column_kinds("hydraulics.darcy_weisbach")
    assert sorted(kinds) == sorted(
        field.name
        for field in dataclasses.fields(type(from_python))
        if field.name not in ("warnings", "units")
    ), "the batch result class carries a field the column rule does not produce"

    for field, kind in kinds.items():
        assert getattr(from_python, field) == getattr(from_rust, field), (
            f"{field}: python {getattr(from_python, field)!r} != rust {getattr(from_rust, field)!r}"
        )
        if kind == "values":
            assert from_python.units[field] == from_rust.units[field], f"{field}: unit map"

    for index in range(3):
        assert h._warnings_equal(from_python.warnings[index], from_rust.warnings[index]), (
            f"warnings[{index}] differ\n  python: {from_python.warnings[index]}\n"
            f"  rust:   {from_rust.warnings[index]}"
        )


def test_an_optional_input_is_per_batch_not_per_element() -> None:
    """`mu` is present for the whole batch or for none of it.

    A partly-supplied optional input has no meaning here that a caller has asked for -
    there is no per-element absent-value encoding, and inventing one would be a design
    decision, not a convenience. Omitted, every element carries RANGE_CHECK_SKIPPED, which
    is the honest answer: a check that could not run is not a check that passed.
    """
    batch = _batch_function("hydraulics.darcy_weisbach")
    required = {"f": [0.02, 0.02], "L": [10.0, 10.0], "D": [0.05, 0.05], "rho": [998.0, 998.0]}
    result = batch(v=[2.0, 3.0], **required)

    assert result.warning_indices(WarningCode.RANGE_CHECK_SKIPPED) == (0, 1)
    assert all(math.isnan(re) for re in result.re)
    assert result.regime == (None, None)
    assert result.units["re"] == "dimensionless"


def test_inputs_are_validated_before_anything_runs() -> None:
    """The caller errors the batch layer owns, each rejected rather than guessed at.

    The public wrappers are keyword-only with required parameters, so a *missing* input is
    a `TypeError` from Python before this layer sees it - which is the better outcome and
    why the missing-input branch is exercised through `run` rather than through a wrapper.
    """
    from azoth.batch._core import run

    def _never(*_args: Any) -> None:  # pragma: no cover - the call must not get this far
        raise AssertionError("the builder ran, so validation did not happen first")

    with pytest.raises(InvalidInputError, match="missing required input"):
        run("hydraulics.reynolds_number", {"rho": [998.0], "v": [1.5], "D": [0.05]}, _never)

    with pytest.raises(InvalidInputError, match="does not declare"):
        run(
            "hydraulics.reynolds_number",
            {"rho": [998.0], "v": [1.5], "D": [0.05], "mu": [1e-3], "viscosity": [1e-3]},
            _never,
        )

    # The one branch a caller can actually reach, because no signature can express it.
    batch = _batch_function("hydraulics.reynolds_number")
    with pytest.raises(InvalidInputError, match="differing lengths"):
        batch(rho=[998.0, 998.0], v=[1.5], D=[0.05], mu=[1e-3])


def test_a_scalar_where_an_array_belongs_is_rejected() -> None:
    """No broadcasting, which is the decision that keeps a batch call honest.

    A length-one list where a length-N list was meant produces a plausible array of
    answers rather than an error, and there is no way for the call to tell that case from
    the intended one. Rejecting a bare scalar is the same instinct as the scalar API's
    refusal to take a bare float where a quantity belongs.
    """
    batch = _batch_function("hydraulics.reynolds_number")
    with pytest.raises(InvalidInputError, match="does not broadcast"):
        batch(rho=998.0, v=[1.5], D=[0.05], mu=[1e-3])


def test_quantities_converts_a_column_and_says_what_it_costs() -> None:
    """The escape hatch back to the scalar API's shape, and its refusal to invent units.

    A column is an array of SI base magnitudes, which is what makes the API fast. A caller
    who wants `pint` quantities asks for them explicitly and pays for the allocation - and
    asking for a field with no unit, or no such field, is an error rather than a silent
    dimensionless quantity.
    """
    batch = _batch_function("hydraulics.pump_power")
    result = batch(rho=[998.0], q=[0.01], H=[20.0], eta=[0.75])

    (power,) = result.quantities("power")
    assert power.check("[power]")
    assert h.to_si_(power, "W", "power") == pytest.approx(result.power[0])

    with pytest.raises(AttributeError, match="no dimensioned field"):
        result.quantities("regime")


def test_inputs_take_any_iterable_and_outputs_are_buffers() -> None:
    """The two halves of the "numpy interop without a numpy dependency" claim, tested
    without numpy installed.

    Inputs go through `sequence`, which iterates, so a generator or a tuple works as well
    as a list - and a numpy array, being iterable, is accepted by the same path rather
    than by a special case that would only be exercised when numpy happens to be present.

    Outputs are `array.array('d')`. A `memoryview` over one is the buffer-protocol read
    that `numpy.asarray` performs, so asserting on it asserts the property `numpy.asarray`
    depends on, without depending on numpy. `format` is checked too: `'d'` is what tells a
    consumer the bytes are already doubles and that no conversion is needed.
    """
    batch = _batch_function("hydraulics.pump_power")
    result = batch(
        rho=(998.0 for _ in range(2)),
        q=[0.01, 0.02],
        H=[20.0, 20.0],
        eta=[0.75, 0.75],
    )

    assert len(result) == 2
    view = memoryview(result.power)
    assert view.format == "d" and view.nbytes == 2 * 8
    assert view[0] == result.power[0]
    h.assert_close(view[1], 2 * view[0], 1e-12, "power[1] is twice power[0]")


def _an_example_call(calc_id: str) -> dict[str, list[float]]:
    """One element of inputs for a calc, taken from its first active spec case."""
    for calc, case in _active_cases():
        if calc["id"] == calc_id:
            return {name: [float(value)] for name, value in case["inputs"].items()}
    raise AssertionError(f"no active case for {calc_id}")


@pytest.mark.parametrize("calc_id", batchable())
def test_every_batch_result_uses_the_shared_repr(calc_id: str) -> None:
    """A subclass must not let the dataclass decorator write its own `__repr__`.

    `@dataclass` generates one unless it is told `repr=False`, and it *shadows* the
    inherited method silently - so the only symptom is a traceback or a log line
    containing every element of every column. There is no type error, no import error,
    and no numeric test that notices.

    Generic over every calc rather than checked on one, because that is exactly how it
    slipped: `friction_factor_haaland` had it right, and the variant that did not was in
    another module. `__len__` is checked alongside, since it is the other method the base
    class provides and a field of the wrong kind could shadow it too.
    """
    result = _batch_function(calc_id)(**_an_example_call(calc_id))
    own = type(result).__dict__
    assert "__repr__" not in own, (
        f"{type(result).__name__} defines its own __repr__; add repr=False to its @dataclass"
    )
    assert "__len__" not in own, f"{type(result).__name__} defines its own __len__"
    # Not `clean=True`: darcy_weisbach's first active case omits viscosity, so it carries
    # RANGE_CHECK_SKIPPED and is not clean. The point is the shape, not the value.
    assert repr(result) == f"{type(result).__name__}(n=1, clean={result.is_clean})"


def test_the_batch_result_reports_its_own_shape() -> None:
    """`len`, `is_clean` and `repr`, so an empty batch and a dirty one are distinguishable."""
    batch = _batch_function("hydraulics.friction_factor_haaland")
    clean = batch(re=[74850.0, 74850.0], relative_roughness=[0.0009, 0.0009])

    assert len(clean) == 2
    assert clean.is_clean
    assert clean.warnings == ((), ())
    assert clean.warning_indices(WarningCode.TRANSITIONAL_FLOW) == ()
    assert not clean.has_warning(WarningCode.TRANSITIONAL_FLOW)
    assert "n=2" in repr(clean) and "clean=True" in repr(clean)
