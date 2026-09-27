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
from azoth import _models_gen
from azoth._dispatch import resolve
from azoth._registry_gen import CALCS
from azoth.core.warnings import WarningCode

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


def _model_cases() -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Every case of every model, taken from the model registry.

    The sibling of :func:`_active_cases`, for the other half of the registry.
    """
    return [(model, case) for model in _models_gen.MODELS for case in model["cases"]]


MODEL_CASES = _model_cases()


def _convergence_tolerance(spec_: dict[str, Any]) -> float | None:
    """See :func:`_helpers.convergence_tolerance`; kept for this module's own callers."""
    return h.convergence_tolerance(spec_)


@pytest.mark.parametrize(
    ("model", "case"),
    MODEL_CASES,
    ids=lambda x: x["id"] if isinstance(x, dict) else "",
)
def test_python_and_rust_agree_on_models(model: dict[str, Any], case: dict[str, Any]) -> None:
    """Run one model case through both implementations and compare.

    The same check :func:`test_python_and_rust_agree` makes for calculations, over the
    ids that were not covered by it.
    """
    from azoth._dispatch import use_backend

    kwargs = h.model_kwargs(model, case["inputs"])
    tolerance = float(case["tolerance"])
    context = f"{model['id']}::{case['id']}"

    with use_backend("python"):
        py = resolve(model["id"])(**kwargs)
    with use_backend("rust"):
        rs = resolve(model["id"])(**kwargs)

    h.assert_results_equal(
        py, rs, tolerance, context, diagnostic_bound=_convergence_tolerance(model)
    )


@pytest.mark.parametrize(
    ("calc", "case"), CASES, ids=lambda x: x["id"] if isinstance(x, dict) else ""
)
def test_python_and_rust_agree(calc: dict[str, Any], case: dict[str, Any]) -> None:
    """Run one spec case through both implementations and compare."""
    kwargs = h.kwargs_for(calc, case["inputs"])
    tolerance = float(case.get("tolerance") or calc["worked_example"]["tolerance"])
    context = f"{calc['id']}::{case['id']}"

    # `resolve` follows the selected backend, so pin each side explicitly rather
    # than assuming which one answered.
    from azoth._dispatch import use_backend

    with use_backend("python"):
        py = resolve(calc["id"])(**kwargs)
    with use_backend("rust"):
        rs = resolve(calc["id"])(**kwargs)

    h.assert_results_equal(
        py, rs, tolerance, context, diagnostic_bound=_convergence_tolerance(calc)
    )


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

    from azoth._dispatch import result_types

    core = _extension()
    for calc_id, result_type in result_types().items():
        typed: Any = result_type
        python_fields = [f.name for f in dataclasses.fields(typed)]
        rust_fields = list(core.result_fields(calc_id))
        assert python_fields == rust_fields, (
            f"{calc_id}: field mismatch\n  python: {python_fields}\n  rust:   {rust_fields}"
        )


def test_every_id_agrees_across_languages() -> None:
    """Rust's `CALC_ID` const and Python's `CALC_ID` attribute must be the same id.

    The ids Rust reports come from its own associated constants - the arms of its id-to-fields
    match are `X::CALC_ID`, not retyped strings - so this compares two independent declarations
    rather than one table read twice. It matters because each side's provenance lookup is keyed
    on its own: Rust's on `CalcResult::CALC_ID`, Python's on the class attribute. An id that
    disagreed would give the two backends different blocks for the same calculation, and
    nothing else would say so.

    Both registries, because both are shipped: `calc_ids()` is the calculations and
    `model_ids()` the models, and together they are every id the extension exposes. That the
    union is the registry's own 192 is asserted rather than assumed - a list that had quietly
    lost an id would otherwise make this test pass by checking less.
    """
    from azoth._dispatch import result_types

    extension = _extension()
    types = result_types()
    rust_ids = [*extension.calc_ids(), *extension.model_ids()]

    assert len(rust_ids) == len(set(rust_ids)), "Rust lists an id twice"
    assert set(rust_ids) == set(types), (
        "the extension and the registry disagree about which ids exist: "
        f"only in Rust {sorted(set(rust_ids) - set(types))}, "
        f"only in Python {sorted(set(types) - set(rust_ids))}"
    )

    for calc_id in rust_ids:
        assert calc_id == types[calc_id].CALC_ID, (
            f"{calc_id}: Rust declares it, Python's {types[calc_id].__name__} claims "
            f"{types[calc_id].CALC_ID!r}"
        )


#: Every warning code this build declares, from the vocabulary rather than listed.
_ALL_WARNING_CODES = tuple(code.value for code in WarningCode)

#: An id to probe the merge with. Any registered id does; this one is chosen because a real
#: call to it skips a check, which the case below relies on.
_PROBE_ID = "hydraulics.darcy_weisbach"

#: Warning shapes the two merges are compared over.
#:
#: The merge is a pure function of the warnings, so the divergence worth catching is in its
#: arithmetic rather than in any one calculation's physics: a repeat that must collapse, a
#: skipped check that names no field, a code that is not a skipped check at all, and both
#: orderings that have to come out independent of the order they arrived in. Driving real
#: calculations instead would exercise one shape per case and cost a run of both backends
#: for each.
_MERGE_SHAPES: tuple[tuple[str, tuple[tuple[str, str | None], ...]], ...] = (
    ("no warnings", ()),
    ("one skipped check", (("RANGE_CHECK_SKIPPED", "re"),)),
    ("a skipped check naming no field", (("RANGE_CHECK_SKIPPED", None),)),
    ("the same field twice", (("RANGE_CHECK_SKIPPED", "re"), ("RANGE_CHECK_SKIPPED", "re"))),
    (
        "two fields, declared out of order",
        (("RANGE_CHECK_SKIPPED", "v"), ("RANGE_CHECK_SKIPPED", "D")),
    ),
    ("a code that is not a skipped check", (("OUT_OF_VALID_RANGE", "L"),)),
    (
        "both kinds, out of vocabulary order",
        (("TRANSITIONAL_FLOW", None), ("RANGE_CHECK_SKIPPED", "re")),
    ),
    ("every code at once", tuple((code, "x") for code in _ALL_WARNING_CODES)),
)


@pytest.mark.parametrize(
    ("shape", "warnings"), _MERGE_SHAPES, ids=[name for name, _ in _MERGE_SHAPES]
)
def test_the_provenance_merge_agrees_across_languages(
    shape: str, warnings: tuple[tuple[str, str | None], ...]
) -> None:
    """Rust's merge and Python's, on the same warnings, must produce the same block.

    **Both languages decide independently what a block says**: which checks were skipped and
    in what order, which warning codes were raised and in what order, whether the call was
    clean. The tables they read are held equal by id - and the static half of the block is
    compared separately - but the arithmetic over them was compared by nothing, and arithmetic
    written twice is arithmetic that can disagree. A block that disagreed would be worse than
    no block: a caller would believe the two implementations agreed about how far the answer
    can be trusted when they did not.

    Python's side is the reference's own `Provenance.of`; Rust's is its `Provenance::of`,
    reached through the extension.
    """
    from azoth.core.provenance import Provenance
    from azoth.core.warnings import Warning, WarningCode

    python_warnings = [
        Warning(code=WarningCode(code), message="a message", field=field)
        for code, field in warnings
    ]
    block = Provenance.of(static_for_probe(), python_warnings)
    rust = _extension().provenance_block(_PROBE_ID, list(warnings))

    assert rust is not None, f"{_PROBE_ID} has no block in the generated table"
    assert (rust[0], rust[1], rust[2], rust[3]) == (
        block.calc_id,
        list(block.skipped_checks),
        list(block.warning_codes),
        block.clean,
    ), shape


def static_for_probe() -> dict[str, Any]:
    """The generated static half for the probe id, so the merge has something to fill in."""
    from azoth._provenance_gen import PROVENANCE

    return PROVENANCE[_PROBE_ID]


def test_a_real_call_produces_the_same_block_on_both_backends() -> None:
    """One case that actually skips a check, so the comparison is not only synthetic.

    `darcy_weisbach` with no viscosity cannot evaluate its Reynolds-number bound, so it
    reports `RANGE_CHECK_SKIPPED` against `re` - the shape the block exists to carry.
    """
    import azoth
    from azoth._dispatch import use_backend

    calc_id = "hydraulics.darcy_weisbach"
    q = azoth.ureg.Quantity
    kwargs = dict(f=0.02, L=q(100.0, "m"), D=q(0.1, "m"), rho=q(998.0, "kg/m**3"), v=q(1.5, "m/s"))

    with use_backend("python"):
        py = resolve(calc_id)(**kwargs)
    with use_backend("rust"):
        rs = resolve(calc_id)(**kwargs)

    rust = _extension().provenance_block(calc_id, [(w.code.value, w.field) for w in rs.warnings])
    assert rust is not None
    assert (rust[1], rust[2], rust[3]) == (
        list(py.provenance.skipped_checks),
        list(py.provenance.warning_codes),
        py.provenance.clean,
    )
    assert py.provenance.skipped_checks == ("re",), "the probe no longer skips a check"


def test_an_unknown_warning_code_is_refused_rather_than_dropped() -> None:
    """A code Rust cannot name must not resolve to something plausible.

    Silently dropping it would make a block that reports fewer caveats than the call raised,
    which is the direction that matters: a caller reads a block to find out what is wrong
    with a number.
    """
    import pytest as _pytest

    with _pytest.raises(ValueError, match="not a warning code"):
        _extension().provenance_block(_PROBE_ID, [("NOT_A_CODE", None)])


def test_errors_are_the_same_class_object() -> None:
    """Both implementations must raise the *same* exception classes.

    Not merely classes with the same names: the same objects, so that a caller
    writing ``except OutOfRangeError`` catches errors from either backend.

    The names come from ``azoth.core.errors.__all__`` rather than from a list here,
    and the two sets are asserted equal in both directions. A hand-written list of
    three names is what this test used to be, and `KeycardError` was exported by the
    package and missing from the extension the whole time it passed.
    """
    from azoth.core import errors

    core = _extension()
    declared = set(errors.__all__)
    missing = sorted(name for name in declared if not hasattr(core, name))
    assert not missing, (
        f"azoth.core.errors declares {missing} but azoth._core does not export "
        f"{'them' if len(missing) > 1 else 'it'}, so `from azoth._core import "
        f"{missing[0]}` fails and the two backends do not agree on the hierarchy"
    )
    for name in sorted(declared):
        assert getattr(errors, name) is getattr(core, name), (
            f"{name} is defined twice, so `except {name}` would catch from one backend "
            f"and not the other"
        )
