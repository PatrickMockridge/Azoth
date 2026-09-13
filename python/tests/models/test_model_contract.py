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

    The whole of the contract: a model spec naming `saturation_pressure_bisection`
    while the Rust implementation ran a different loop would agree with the Python
    side on every case only by luck, and the mismatch would be a name nothing
    compares.
    """
    for model in _models_gen.MODELS:
        rust_schemes = _extension().model_schemes(model["id"])
        assert rust_schemes == [model["algorithm"]["scheme"]], (
            f"{model['id']}: the spec names {model['algorithm']['scheme']!r} but Rust "
            f"reports {rust_schemes}"
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
        "algorithm",
        "inputs",
        "outputs",
        "verification",
        "implementations",
        "cases",
    }
    missing = needed - required
    assert not missing, f"the model schema does not require {sorted(missing)}"

    algorithm_required = set(schema["$defs"]["algorithm"]["required"])
    assert {"scheme", "tolerance", "max_iterations", "convergence"} <= algorithm_required, (
        "the algorithm block must require the parameters both implementations read"
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
    """`_MODEL_IMPLEMENTATIONS` covers exactly the model registry.

    Kept in its own table rather than the calc one, because that table is asserted to
    be exactly the calc registry in both directions. This is the test that makes the
    separation hold: a model added to the registry without a bridge entry would
    otherwise resolve to the Python reference for one backend and raise for the
    other, which is a cross-language disagreement rather than a missing function.
    """
    from azoth._rust_bridge import _MODEL_IMPLEMENTATIONS

    assert set(_MODEL_IMPLEMENTATIONS) == schema_model_ids(), (
        f"the bridge's model table and the registry differ\n"
        f"  missing from the bridge: {sorted(schema_model_ids() - set(_MODEL_IMPLEMENTATIONS))}\n"
        f"  not in the registry: {sorted(set(_MODEL_IMPLEMENTATIONS) - schema_model_ids())}"
    )


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

    from azoth.core.result import MODEL_RESULT_TYPES

    core = _extension()
    assert MODEL_RESULT_TYPES, "no model results to check"
    for model_id, result_type in MODEL_RESULT_TYPES.items():
        typed: Any = result_type
        python_fields = [f.name for f in dataclasses.fields(typed)]
        rust_fields = list(core.result_fields(model_id))
        assert python_fields == rust_fields, (
            f"{model_id}: field mismatch\n  python: {python_fields}\n  rust:   {rust_fields}"
        )


def test_the_model_result_table_covers_every_model() -> None:
    """`MODEL_RESULT_TYPES` is exactly the model registry.

    Kept separate from the calc table, which is asserted to be exactly the calc
    registry - so a model result appearing there would break that contract rather
    than extend it. A model in the registry with no result type would be a shape
    nothing checks, which is what this table exists to prevent.
    """
    from azoth.core.result import MODEL_RESULT_TYPES

    assert set(MODEL_RESULT_TYPES) == schema_model_ids(), (
        f"the model result table and the registry differ\n"
        f"  missing: {sorted(schema_model_ids() - set(MODEL_RESULT_TYPES))}\n"
        f"  extra:   {sorted(set(MODEL_RESULT_TYPES) - schema_model_ids())}"
    )


def test_every_model_scheme_and_inner_scheme_is_named() -> None:
    """A nested scheme is part of the procedure, so the spec has to declare it.

    The flash runs Rachford-Rice every outer iteration, and the spec says so in
    `algorithm.inner`. A spec that named only the outer scheme would leave the
    inner loop's tolerance undocumented and free to differ between implementations -
    which is the whole failure this tree exists to prevent.
    """
    for model in _models_gen.MODELS:
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
