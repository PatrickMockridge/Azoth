"""The model vocabulary is one set, agreed on by three separate places.

A model spec names its *scheme* in ``algorithm.scheme``, and both implementations
must run that scheme. The name is a cross-language contract in the same way
``solver.kind`` is, and it has to hold across:

* ``specs/schema/model.schema.json``, which is what a model spec is validated against;
* the Python implementation, which reads the scheme from ``azoth._models_gen``;
* ``crates/azoth-eos/src/model_gen.rs``, reachable from here through
  ``azoth._core.model_schemes()``.

The difference from the solver contract is that a model's scheme is a *name for a
procedure* rather than a pointer to a shared function - the two implementations write
the loop separately and must be told they are writing the same one. So the name is
the whole of the contract, and this file is what makes it one.

It also asserts the two registries stay separate. A model id in ``calc_ids()`` would
break ``test_registration_completeness``'s exact-equality assertion against
``specs/calcs/``, so the separation is deliberate and load-bearing rather than an
oversight.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth import _models_gen

REPO_ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "model.schema.json"


def _extension() -> ModuleType:
    """The compiled extension, imported by name."""
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def model_result_types() -> dict[str, type[Any]]:
    """The result dataclass each model produces, derived from the models themselves.

    From `azoth._dispatch.result_types`, which reads each implementation's return
    annotation, filtered to the model ids. It used to be a hand-maintained table in
    `core.result`, asserted equal to the registry by the test below - which is a list
    kept in step with a list, so the test could only tell you that you had forgotten.
    """
    from azoth._dispatch import result_types

    return {key: value for key, value in result_types().items() if key in schema_model_ids()}


def schema_model_ids() -> set[str]:
    """The model ids the spec tree declares."""
    return {str(m["id"]) for m in _models_gen.MODELS}


def test_there_are_models_to_check() -> None:
    """The checks below would pass vacuously on an empty registry."""
    assert _models_gen.MODELS, "no models in the registry - is specs/models/ populated?"


def test_every_model_spec_is_in_the_generated_registry() -> None:
    """`specs/models/**` and `_models_gen.MODELS` are the same set.

    The generated table is the source of truth the implementations read, so a spec
    that did not reach it would be a model nothing runs.
    """
    spec_ids = {
        str(raw["id"])
        for path in sorted((REPO_ROOT / "specs" / "models").rglob("*.yaml"))
        for raw in [__import__("yaml").safe_load(path.read_text(encoding="utf-8"))]
    }
    registry_ids = schema_model_ids()
    assert spec_ids == registry_ids, (
        f"specs/models and the generated registry differ\n"
        f"  only on disk: {sorted(spec_ids - registry_ids)}\n"
        f"  only generated: {sorted(registry_ids - spec_ids)}"
    )


@pytest.mark.requires_rust
def test_the_rust_registry_covers_the_same_models() -> None:
    """`_core.model_ids()` and the generated Python registry are one set."""
    extension_ids = set(_extension().model_ids())
    python_ids = schema_model_ids()
    assert extension_ids == python_ids, (
        f"azoth._core.model_ids() and _models_gen differ\n"
        f"  missing from the extension: {sorted(python_ids - extension_ids)}\n"
        f"  not in the registry: {sorted(extension_ids - python_ids)}"
    )


@pytest.mark.requires_rust
def test_every_scheme_is_reachable_from_both_languages() -> None:
    """Each model's scheme name is the one the Rust side implements.

    The whole of the contract for a procedure: a model spec naming
    `saturation_pressure_bisection` while the Rust implementation ran a different loop
    would agree with the Python side on every case only by luck, and the mismatch would
    be a name nothing compares.

    A `direct` model has no scheme, so the assertion is that the list is *empty* rather
    than that some placeholder matches. That is the difference between a contract and a
    comparison that passes vacuously.
    """
    core = _extension()
    for model in _models_gen.MODELS:
        rust_schemes = core.model_schemes(model["id"])
        expected = [model["algorithm"]["scheme"]] if "algorithm" in model else []
        assert rust_schemes == expected, (
            f"{model['id']}: the spec implies {expected} but Rust reports {rust_schemes}"
        )


@pytest.mark.requires_rust
def test_every_model_kind_is_reachable_from_both_languages() -> None:
    """Each model's `kind` is the one the Rust side reports.

    `kind` decides whether a model has an algorithm block at all, so it is load-bearing
    rather than descriptive - and a spec field that decides something and is checked by
    nothing is exactly the defect this project is organised against. Held to the Rust
    table the same way `scheme` is.
    """
    core = _extension()
    for model in _models_gen.MODELS:
        assert core.model_kind(model["id"]) == model["kind"], (
            f"{model['id']}: the spec says {model['kind']!r} but Rust reports "
            f"{core.model_kind(model['id'])!r}"
        )


def test_the_schema_requires_every_field_a_model_needs() -> None:
    """The schema's required list covers what the implementations read.

    A field the implementations depend on but the schema does not require would be
    optional in the spec and mandatory in the code - so a spec omitting it would
    validate and then fail at the generator, which is later than it should.
    """
    schema: dict[str, Any] = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    required = set(schema["required"])
    needed = {
        "id",
        "name",
        "kind",
        "inputs",
        "outputs",
        # Not `verification`: it is optional now. A status on every spec was a
        # mandatory badge on everything, which said nothing about any one of them.
        "implementations",
        "cases",
    }
    missing = needed - required
    assert not missing, f"the model schema does not require {sorted(missing)}"

    # `algorithm` is required of a *procedure* and forbidden of a direct model, which
    # is what lets a model with no iteration live in the same tree. Asserted as a
    # conditional rather than as a required key, because a schema that required it
    # outright would force a direct model to carry a vacuous block.
    assert "algorithm" not in required, (
        "`algorithm` must not be unconditionally required: a direct model has none, and "
        "the conditional in allOf is what expresses that"
    )
    assert "allOf" in schema, "the kind-conditional on `algorithm` is missing"

    algorithm_required = set(schema["$defs"]["algorithm"]["required"])
    assert {"scheme", "tolerance", "max_iterations", "convergence"} <= algorithm_required, (
        "the algorithm block must require the parameters both implementations read"
    )


def test_every_model_declares_a_kind_the_schema_allows() -> None:
    """`kind` is `procedure` or `direct`, and a procedure has an algorithm.

    The schema expresses this conditionally, so nothing else would catch a spec that
    named a third kind or a direct model carrying an algorithm - both of which the
    generator would silently accept a default for.
    """
    schema: dict[str, Any] = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    allowed = set(schema["properties"]["kind"]["enum"])
    for model in _models_gen.MODELS:
        assert model["kind"] in allowed, f"{model['id']}: kind {model['kind']!r}"
        if model["kind"] == "procedure":
            assert "algorithm" in model, f"{model['id']}: a procedure needs an algorithm"
        else:
            assert "algorithm" not in model, (
                f"{model['id']}: a direct model must not carry an algorithm block - "
                f"there is no loop for one to describe"
            )


def test_the_model_registry_is_not_the_calc_registry() -> None:
    """A model id must not leak into `calc_ids()`.

    Asserted because the alternative is silent: `test_registration_completeness`
    compares `calc_ids()` to the calc registry by exact equality, so a model id
    appearing there would fail *that* test with a confusing message rather than this
    one with a clear one.
    """
    from azoth._registry_gen import BY_ID

    overlap = schema_model_ids() & set(BY_ID)
    assert not overlap, f"ids in both registries: {sorted(overlap)}"


def test_the_bridge_carries_a_model_arm_for_every_model() -> None:
    """Every registered model has a bridge function to resolve to.

    The bridge used to keep a model table separate from the calc one, and both are
    gone: `resolve` derives the function from the id, so what can go wrong is a
    bridge function that was never written. A model without one would resolve to the
    Python reference on one backend and raise on the other, which is a cross-language
    disagreement rather than a missing function - so it is checked here rather than
    left to the first caller who switches backends.
    """
    from azoth._rust_bridge import resolve

    for model_id in sorted(schema_model_ids()):
        assert callable(resolve(model_id)), f"{model_id} resolves to nothing in the bridge"


@pytest.mark.requires_rust
def test_the_result_shapes_agree_across_languages() -> None:
    """Every model result dataclass must have exactly the fields Rust reports.

    The same check `test_cross_impl.py` runs for the calc registry, which the models
    were outside of until the flash made it matter: the flash's result has thirteen
    fields, one of them optional and five of them vectors, and a shape mismatch is
    invisible to every numerical test - the numbers agree perfectly and only the
    attribute differs.
    """
    import dataclasses

    core = _extension()
    assert model_result_types(), "no model results to check"
    for model_id, result_type in model_result_types().items():
        typed: Any = result_type
        python_fields = [f.name for f in dataclasses.fields(typed)]
        rust_fields = list(core.result_fields(model_id))
        assert python_fields == rust_fields, (
            f"{model_id}: field mismatch\n  python: {python_fields}\n  rust:   {rust_fields}"
        )


def test_the_model_result_table_covers_every_model() -> None:
    """Every model resolves to a result type, and nothing else does.

    The derivation reads a return annotation, so a model whose implementation
    annotates nothing - or annotates something that is not a dataclass - is caught
    here rather than at the first caller who asks for a field. That is what the
    hand-maintained table this replaced was for, and it caught only the case where
    somebody remembered to add a row.
    """

    assert set(model_result_types()) == schema_model_ids(), (
        f"the model result table and the registry differ\n"
        f"  missing: {sorted(schema_model_ids() - set(model_result_types()))}\n"
        f"  extra:   {sorted(set(model_result_types()) - schema_model_ids())}"
    )


def test_every_model_scheme_and_inner_scheme_is_named() -> None:
    """A nested scheme is part of the procedure, so the spec has to declare it.

    The flash runs Rachford-Rice every outer iteration, and the spec says so in
    `algorithm.inner`. A spec that named only the outer scheme would leave the
    inner loop's tolerance undocumented and free to differ between implementations -
    which is the whole failure this tree exists to prevent.
    """
    for model in _models_gen.MODELS:
        if model["kind"] != "procedure":
            # A direct model has no algorithm to name, and asserting anything about
            # one would be asserting something about a block it must not have.
            continue
        algorithm = model["algorithm"]
        assert "initialisation" in algorithm or algorithm["scheme"].endswith("bisection"), (
            f"{model['id']}: a scheme that needs a starting point must declare "
            f"`initialisation`, or two implementations take different paths"
        )
        inner = algorithm.get("inner")
        if inner is not None:
            assert inner["scheme"] != algorithm["scheme"], (
                f"{model['id']}: the inner scheme must be a different procedure"
            )
            assert inner["tolerance"] < algorithm["tolerance"], (
                f"{model['id']}: an inner solve must be tighter than the outer loop "
                f"it feeds, or the outer residual measures the inner tolerance"
            )
