"""A result, as data.

:mod:`azoth.core.serialise` writes every result in this library through one walk over
``dataclasses.fields``, and this is where that claim is held rather than asserted: **every case of
every registered calculation and model** is run and serialised. A new result is covered by being
registered, and a field of a type the walk does not know fails here rather than in a caller's
notebook.

What it cannot catch, and says so: whether the *numbers* are right. This is a shape test - the
arithmetic is `test_cross_impl.py`'s, and a serializer that wrote the wrong field's value would
pass every assertion below while `assert_results_equal` would not have moved.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import _helpers as h
import azoth
from azoth import _models_gen
from azoth._dispatch import resolve
from azoth._registry_gen import CALCS
from azoth.core import serialise
from azoth.core.units import CANONICAL_UNITS, unit_for
from azoth.core.warnings import Warning, WarningCode

REPO_ROOT = Path(__file__).resolve().parents[2]


def _cases() -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Every active worked example or reference case, as `test_cross_impl.py` takes them."""
    return [
        (calc, case)
        for calc in CALCS
        for case in h.all_tests(calc)
        if case["status"] == "active" and case["type"] in ("worked_example", "reference")
    ]


def _model_cases() -> list[tuple[dict[str, Any], dict[str, Any]]]:
    return [(model, case) for model in _models_gen.MODELS for case in model["cases"]]


def _documents(value: Any) -> list[dict[str, Any]]:
    """Every quantity-shaped mapping inside serialised data, at any depth."""
    found: list[dict[str, Any]] = []
    if isinstance(value, dict):
        if set(value) in ({"magnitude", "unit"}, {"magnitude_si", "unit"}):
            found.append(value)
        for item in value.values():
            found.extend(_documents(item))
    elif isinstance(value, list):
        for item in value:
            found.extend(_documents(item))
    return found


def _assert_serialises(result: Any, context: str) -> dict[str, Any]:
    """The shape contract, applied to one result."""
    document = result.to_json()
    loaded = json.loads(document)
    assert isinstance(loaded, dict), f"{context}: a result is an object"

    # **Every field is in the document.** A field the walk skipped would be a value a consumer
    # silently does not see, which is the failure this whole module is against.
    import dataclasses

    declared = {field.name for field in dataclasses.fields(result)}
    assert declared == set(loaded), (
        f"{context}: the document and the result disagree about the fields: "
        f"{sorted(declared ^ set(loaded))}"
    )

    # The same result serialises to the same bytes, so two runs of one calculation are comparable.
    assert document == result.to_json(), f"{context}: serialising twice differs"

    # Every quantity names a unit this library can resolve, in the vocabulary's own spelling:
    # `unit_for` is what turns it back into something `pint` reads, and that is the way back a
    # consumer takes.
    for quantity in _documents(loaded):
        name = quantity["unit"]
        assert isinstance(name, str) and name != "", f"{context}: a quantity without a unit"
        # A result's unit is one its spec declares, so it is one the vocabulary names — and the
        # name is the way back, through the same `unit_for` every input takes. The fallback to
        # `pint`'s own spelling exists for a *derived* unit no spec names, which is tested
        # separately rather than being allowed to hide here.
        assert unit_for(name), f"{context}: `{name}` is not a unit the vocabulary names"
    return loaded


@pytest.mark.parametrize(
    ("calc", "case"), _cases(), ids=lambda x: x["id"] if isinstance(x, dict) else ""
)
def test_a_calculation_result_serialises(calc: dict[str, Any], case: dict[str, Any]) -> None:
    result = resolve(calc["id"])(**h.kwargs_for(calc, case["inputs"]))
    _assert_serialises(result, f"{calc['id']}::{case['id']}")


@pytest.mark.parametrize(
    ("model", "case"), _model_cases(), ids=lambda x: x["id"] if isinstance(x, dict) else ""
)
def test_a_model_result_serialises(model: dict[str, Any], case: dict[str, Any]) -> None:
    result = resolve(model["id"])(**h.model_kwargs(model, case["inputs"]))
    _assert_serialises(result, f"{model['id']}::{case['id']}")


def test_the_unit_a_result_names_is_one_the_vocabulary_round_trips() -> None:
    """`unit_name` and `pint_name` are inverses, and the vocabulary is a bijection.

    **The bijection is the load-bearing half.** `unit_name` inverts a spec-name-to-`pint`-name
    table, and two spec names sharing one `pint` unit would make the inverse lose one of them -
    silently, and only for whichever unit came second. Asserted here so a new unit that collides
    fails in the table rather than in a document.
    """
    assert len(set(CANONICAL_UNITS.values())) == len(CANONICAL_UNITS)
    for spec_name in CANONICAL_UNITS:
        assert serialise.unit_name(serialise.pint_name(spec_name)) == spec_name
        assert azoth.ureg.Quantity(1.0, unit_for(spec_name)).units is not None


def test_a_quantity_is_the_shape_the_flowsheet_codec_writes() -> None:
    """The two layers that write a quantity are held to each other, on the key that matters.

    A stream record's five fields are SI by construction, so `executor::json` names its magnitude
    `magnitude_si` and is right to. A result's unit is whatever its spec declares, so this writer
    does not claim SI. **What both write is `unit`**, and that is the key a consumer reads without
    knowing which layer it came from - the difference is deliberate and this is where it is
    recorded rather than remembered.
    """
    fixture = REPO_ROOT / "ui/test/fixtures/envelope.json"
    if not fixture.exists():  # pragma: no cover - the fixture is committed
        pytest.skip("the envelope fixture is not in this tree")

    envelope = json.loads(fixture.read_text(encoding="utf-8"))
    stream = envelope["session"]["streams"]["p1.outlet"]["P"]
    assert set(stream) == {"magnitude_si", "unit"}

    from azoth.hydraulics import pump_power

    q = azoth.ureg.Quantity
    result = pump_power(q(998.0, "kg/m**3"), q(5.0, "m**3/s"), q(20.0, "m"), 0.75)
    (quantity,) = _documents(result.to_dict())
    assert set(quantity) == {"magnitude", "unit"}
    assert "unit" in stream, "both layers name the unit in the same key"

    # And the vocabulary's name is what a result writes, which is what the generated docs print.
    assert quantity["unit"] == "W"
    assert round(quantity["magnitude"] / 1000.0, 1) == pytest.approx(1304.9, abs=0.1)


def test_a_unit_the_vocabulary_has_no_name_for_keeps_pints_spelling() -> None:
    """The fallback, tested where it can be seen rather than relied on in the walk.

    A derived unit — a per-second-per-square-metre that no spec declares — has no vocabulary name,
    and inventing one would be worse than writing what `pint` parses.
    """
    derived = azoth.ureg.Unit("1/(s*m**2)")
    assert serialise.unit_name(derived) == str(derived)
    assert azoth.ureg.Unit(serialise.unit_name(derived)) == derived
    # And a compound unit a spec *does* name comes back in the spec's spelling, which is what the
    # generated docs print.
    assert serialise.unit_name(azoth.ureg.Unit("m**2/s")) == "m**2/s"
    assert serialise.unit_name(azoth.ureg.Unit("kg/m**3")) == "kg/m**3"


def test_a_warning_keeps_its_code_its_message_and_its_field() -> None:
    """A caveat is data, not prose: a caller switches on the code and points at the field."""
    warning = Warning(
        code=WarningCode.OUT_OF_VALID_RANGE,
        message="the friction factor correlation is not validated here",
        field="re",
    )
    assert serialise.to_dict(warning) == {
        "code": "OUT_OF_VALID_RANGE",
        "message": "the friction factor correlation is not validated here",
        "field": "re",
    }


def test_a_magnitude_that_is_not_a_number_is_written_as_no_value() -> None:
    """`NaN` means "nothing there", and a JSON document has one way to say that.

    **Run against the case that found this**, which is why the walk runs every registered case: a
    single gas phase is skipped, so there is no residual to report, and the model says so with a
    `NaN` rather than with a `None` - a distinction that did not matter until a document had to
    carry it.
    """
    for nothing in (float("nan"), float("inf"), float("-inf")):
        assert serialise.to_dict(nothing) is None
        assert serialise.to_dict({"a": [1.0, {"b": nothing}]}) == {"a": [1.0, {"b": None}]}
        assert serialise.to_dict(azoth.ureg.Quantity(nothing, "m")) == {
            "magnitude": None,
            "unit": "m",
        }
    # **The walk is what makes the text safe.** It leaves nothing non-finite behind, so the strict
    # encoder cannot raise - and `to_json` passes `allow_nan=False` precisely so that if a branch
    # here ever let one through, the write fails instead of emitting `NaN`.
    assert json.dumps(serialise.to_dict({"a": [float("nan")]}), allow_nan=False) == '{"a": [null]}'
    with pytest.raises(ValueError):
        json.dumps({"a": float("nan")}, allow_nan=False)


def test_the_skipped_case_that_found_this_serialises() -> None:
    """The case itself, kept as the regression: a skipped solve is a document, not an error."""
    model = next(m for m in _models_gen.MODELS if m["id"] == "reactions.reactive_phase_equilibrium")
    case = next(c for c in model["cases"] if c["id"] == "a_single_gas_phase_is_skipped")
    result = resolve(model["id"])(**h.model_kwargs(model, case["inputs"]))
    document = result.to_dict()
    assert document["max_reaction_log_residual"] is None
    # **The absence is explained by a declared field and not by a warning**, which is better than
    # the reverse: a consumer reading the document sees `skipped: true` beside the `null` and does
    # not have to know what a missing residual implies.
    assert document["skipped"] is True
    assert json.loads(result.to_json())["max_reaction_log_residual"] is None


def test_a_field_of_an_unknown_type_is_refused_with_its_path() -> None:
    """Refused rather than stringified: `str()` of an object is a value nothing reads back."""
    with pytest.raises(TypeError) as raised:
        serialise.to_dict({"outer": {"inner": object()}})
    assert "result.outer.inner" in str(raised.value)
    assert "object" in str(raised.value)


def test_a_batch_result_serialises_too() -> None:
    """A batch is a result with a different shape, and it gets the same document.

    **It is not a `_HasWarnings`** - its warnings are one tuple per element - so it carries the
    two methods of its own rather than inheriting them, and both delegate to the same writer. This
    is what says so: a batch computed from a registered case's own inputs, written and read back.
    """
    import importlib

    calc = next(c for c in CALCS if c["id"] == "hydraulics.reynolds_number")
    case = next(c for c in h.all_tests(calc) if c["status"] == "active" and "inputs" in c)
    namespace, _, name = calc["id"].partition(".")
    batch = getattr(importlib.import_module(f"azoth.batch.{namespace}"), name)
    result = batch(**{key: [float(value)] for key, value in case["inputs"].items()})

    document = result.to_dict()
    assert set(document) >= {"units", "warnings"}
    assert result.to_json() == serialise.to_json(result), "one writer, two doors"
    loaded = json.loads(result.to_json())
    for quantity in _documents(loaded):
        assert unit_for(quantity["unit"])


def test_both_writers_agree_that_a_non_finite_magnitude_cannot_be_a_number() -> None:
    """The session codec refuses one; this writes `null`. Both are right about their own case.

    A session that ran cannot produce a `NaN`, so refusing it there is refusing a computation gone
    wrong. A result carries one on purpose - a skipped solve has no residual - so refusing it here
    would leave a legitimate state with no document at all.
    """
    from azoth.process import run_flowsheet

    demo = (REPO_ROOT / "specs" / "flowsheets" / "demo.toml").read_text()
    # The session codec writes real magnitudes and refuses anything else; the comparison here is
    # that neither writer ever puts a bare `NaN` in a document.
    assert "NaN" not in run_flowsheet(demo).document
