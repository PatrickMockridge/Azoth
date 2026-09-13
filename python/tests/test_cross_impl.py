"""Python and Rust must agree on every case the specs declare.

This is the test the project's central promise rests on: two independent
implementations, cross-checked, so that a bug has to be made twice in two
languages before it reaches a user.

Marked ``requires_rust``: it skips when the extension is not built, and *fails*
when ``AZOTH_REQUIRE_RUST=1``. That distinction is the point - a skip is fine
while developing in Python, but in CI a skip here would mean the cross-language
guarantee is verified by nothing.

Comparison is by tolerance, not bit-equality, for a documented reason: ``log10``
and ``sqrt`` are not correctly-rounded in general and libm differs between glibc,
musl and macOS. Bit-equality holds for the same platform and build, which is what
the determinism tests assert separately.
"""

from __future__ import annotations

import importlib
import inspect
from types import ModuleType
from typing import Any

import pytest

import _helpers as h
from azoth._dispatch import resolve
from azoth._registry_gen import CALCS
from azoth.core.units import quantity

pytestmark = pytest.mark.requires_rust


def _extension() -> ModuleType:
    """The compiled extension, imported by name.

    By name rather than ``from azoth import _core`` because the module does not
    exist until the bindings are built, and an attribute mypy cannot resolve is a
    worse trade than a lookup that fails clearly at runtime.
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def _kwargs(calc: dict[str, Any], inputs: dict[str, Any]) -> dict[str, Any]:
    """Turn a spec's test inputs into keyword arguments for the real signature.

    Both the type and the unit come from the spec's own ``inputs`` block, so
    adding a calc needs no change here.

    Note that `type` is absent for a plain quantity - it defaults to `quantity`,
    and the specs only spell it out for the exceptions. Treating a missing `type`
    as anything other than the default silently passes bare floats where
    quantities are required, which is exactly the mistake this library exists to
    make impossible; the reference implementation rejects it, and did.
    """
    declared = calc["inputs"]
    kwargs: dict[str, Any] = {}
    for name, value in inputs.items():
        declaration = declared[name]
        kind = declaration.get("type", "quantity")
        if kind == "fitting_list":
            kwargs[name] = list(value)
        elif kind == "quantity" and declaration.get("unit") != "dimensionless":
            kwargs[name] = quantity(float(value), declaration["unit"])
        else:
            kwargs[name] = float(value)
    return kwargs


def _active_cases() -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Every runnable case, taken from the specs rather than listed here.

    Uses the shared helper, which folds the `worked_example` block into a
    runnable case the same way the Rust codegen does. Building it separately here
    would be a second definition of "what the spec's tests are", and the two
    could disagree.
    """
    cases = []
    for calc in CALCS:
        for case in h.all_tests(calc):
            if case["status"] == "active" and case["type"] in ("worked_example", "reference"):
                cases.append((calc, case))
    return cases


CASES = _active_cases()


@pytest.mark.parametrize(
    ("calc", "case"), CASES, ids=lambda x: x["id"] if isinstance(x, dict) else ""
)
def test_python_and_rust_agree(calc: dict[str, Any], case: dict[str, Any]) -> None:
    """Run one spec case through both implementations and compare."""
    kwargs = _kwargs(calc, case["inputs"])
    tolerance = float(case.get("tolerance") or calc["worked_example"]["tolerance"])
    context = f"{calc['id']}::{case['id']}"

    # `resolve` follows the selected backend, so pin each side explicitly rather
    # than assuming which one answered.
    from azoth._dispatch import use_backend

    with use_backend("python"):
        py = resolve(calc["id"])(**kwargs)
    with use_backend("rust"):
        rs = resolve(calc["id"])(**kwargs)

    h.assert_results_equal(py, rs, tolerance, context)


def test_signatures_agree_across_languages() -> None:
    """The Rust signatures must take the same inputs the specs declare.

    Checked against the extension's own introspection, so a Rust function that
    silently gained or lost a parameter is caught even though its numbers are
    right.
    """
    core = _extension()
    for calc in CALCS:
        _, _, function_name = calc["id"].rpartition(".")
        declared = set(calc["inputs"])

        rust_parameters = set(inspect.signature(getattr(core, function_name)).parameters)
        assert rust_parameters == declared, (
            f"{calc['id']}: the Rust implementation takes {sorted(rust_parameters)} "
            f"but the spec declares {sorted(declared)}"
        )


def test_warning_codes_agree_across_languages() -> None:
    """The two warning-code sets must be identical.

    They are a cross-language contract: a caller comparing a warning from either
    implementation must not need to know which one produced it.
    """
    from azoth.core.warnings import WarningCode

    core = _extension()
    rust_codes = set(core.warning_codes())
    python_codes = {code.value for code in WarningCode}
    assert rust_codes == python_codes, (
        f"warning codes differ\n"
        f"  only in Python: {sorted(python_codes - rust_codes)}\n"
        f"  only in Rust:   {sorted(rust_codes - python_codes)}"
    )


def test_result_shapes_agree_across_languages() -> None:
    """Every result dataclass must have exactly the fields Rust reports."""
    import dataclasses

    from azoth.core.result import RESULT_TYPES

    core = _extension()
    for calc_id, result_type in RESULT_TYPES.items():
        typed: Any = result_type
        python_fields = [f.name for f in dataclasses.fields(typed)]
        rust_fields = list(core.result_fields(calc_id))
        assert python_fields == rust_fields, (
            f"{calc_id}: field mismatch\n  python: {python_fields}\n  rust:   {rust_fields}"
        )


def test_errors_are_the_same_class_object() -> None:
    """Both implementations must raise the *same* exception classes.

    Not merely classes with the same names: the same objects, so that a caller
    writing ``except OutOfRangeError`` catches errors from either backend.
    """
    from azoth.core import errors

    core = _extension()
    for name in ("OutOfRangeError", "UnknownFittingError", "SolverNotConvergedError"):
        assert getattr(errors, name) is getattr(core, name), (
            f"{name} is defined twice, so `except {name}` would catch from one backend "
            f"and not the other"
        )
