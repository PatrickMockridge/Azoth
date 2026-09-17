#!/usr/bin/env python3
"""Generate the model registries from `specs/models/**`.

The sibling of `gen_registry.py`, and deliberately a sibling rather than a branch of
it. A calc's spec fixes an equation; a model's fixes a *procedure* - a bracketing
rule, a tolerance, an iteration cap - and the parameters that have to reach both
implementations are therefore a different set. Running both through one generator
would mean one file emitting two unrelated table shapes.

What is *not* duplicated: the emitters for bounds and test cases are imported from
`gen_registry.py`, because a model's bounds are the same kind of bound and its cases
the same kind of case. A second implementation of `emit_range_check` would be a
second definition of what a bound is, which is the drift this project is built
against.

    python tools/gen_models.py            # write
    python tools/gen_models.py --check    # fail if the tree is not what this emits
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import schema_registry
from gen_registry import (
    ROOT,
    emit_range_check,
    emit_test_case,
    rust_f64,
    rust_str,
)

MODEL_DIR = ROOT / "specs" / "models"
CASE_DIR = ROOT / "specs" / "cases"


def crate_for(namespace: str) -> str:
    """The crate that owns a namespace, mirroring `gen_registry.crate_for`."""
    return f"azoth-{namespace.replace('_', '-')}"


def rust_out_for(namespace: str) -> Path:
    """Where a namespace's generated model table goes.

    Fails rather than falling back when the crate does not exist, for the same reason
    `gen_registry` does: emitting into a directory that is not there is a mistake
    worth being told about, and emitting into some *other* crate's directory is a
    worse one.
    """
    crate = ROOT / "crates" / crate_for(namespace)
    if not (crate / "src").is_dir():
        sys.exit(
            f"gen_models: no crate for namespace '{namespace}' at {crate}. Create it "
            f"before adding a model to that namespace."
        )
    return crate / "src" / "model_gen.rs"


PY_OUT = ROOT / "python" / "src" / "azoth" / "_models_gen.py"


def namespace_of(model_id: str) -> str:
    return model_id.split(".", 1)[0]


def ident(model_id: str) -> str:
    """`eos.pure_saturation` -> `PURE_SATURATION`."""
    return model_id.split(".", 1)[1].upper()


def load_models() -> list[dict[str, Any]]:
    """Every model spec, validated against the schema, in id order."""
    paths = sorted(MODEL_DIR.rglob("*.toml"))
    if not paths:
        return []

    try:
        from jsonschema import Draft202012Validator
    except ImportError:  # pragma: no cover
        sys.exit("gen_models requires jsonschema")

    schema = schema_registry.load("model.schema.json")
    validator = Draft202012Validator(schema, registry=schema_registry.registry())

    models: list[dict[str, Any]] = []
    for path in paths:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
        errors = sorted(validator.iter_errors(raw), key=lambda e: list(e.path))
        if errors:
            for error in errors:
                where = "/".join(str(p) for p in error.path) or "(root)"
                print(f"gen_models: {path.relative_to(ROOT)}: {where}: {error.message}")
            sys.exit(1)
        raw["_path"] = str(path.relative_to(ROOT))
        models.append(raw)

    models.sort(key=lambda m: m["id"])
    seen: set[str] = set()
    for model in models:
        if model["id"] in seen:
            sys.exit(f"gen_models: duplicate model id '{model['id']}'")
        seen.add(model["id"])

    # The instances. A model is a type; the machines to run live in their own files,
    # against their own schema, so the type is never edited to make a test pass.
    # Every schema, because a case's values are the model's - vectors and matrices
    # are inputs here, where a calc's are scalars - and the model schema `$ref`s the
    # calc schema in turn.
    case_schema = schema_registry.load("case.schema.json")
    case_validator = Draft202012Validator(case_schema, registry=schema_registry.registry())

    by_id = {model["id"]: model for model in models}
    for model in models:
        model["cases"] = []
    for path in sorted(CASE_DIR.rglob("*.toml")):
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
        errors = sorted(case_validator.iter_errors(raw), key=lambda e: list(e.path))
        if errors:
            for error in errors:
                where = "/".join(str(p) for p in error.path) or "(root)"
                print(f"gen_models: {path.relative_to(ROOT)}: {where}: {error.message}")
            sys.exit(1)
        target = by_id.get(raw["model"])
        if target is None:
            sys.exit(
                f"gen_models: {path.relative_to(ROOT)} instantiates {raw['model']!r}, "
                f"which is no model. A case file that names nothing is a file nobody runs."
            )
        if target["cases"]:
            sys.exit(f"gen_models: two case files for {raw['model']!r}")
        target["cases"] = raw["cases"]

    for model in models:
        if not model["cases"]:
            sys.exit(
                f"gen_models: {model['id']} has no cases. Every model ships at least one "
                f"instance in specs/cases/, or nothing verifies it."
            )
    return models


def emit_rust_algorithm(
    algorithm: dict[str, Any], inner_ref: str, indent: str, fallback_ref: str = "None"
) -> str:
    """One `ModelAlgorithm` literal.

    `inner_ref` and `fallback_ref` are the Rust expressions for this scheme's `inner`
    and `fallback` fields. A nested scheme is a `static` of its own rather than a
    literal inlined here, because a `ModelAlgorithm` refers to its nested ones by
    `&'static`, so they have to have somewhere to live. Passing the references in
    rather than patching placeholders afterwards keeps the callers - leaf, parent and
    fallback - from having to agree about the shape of the text.
    """
    bracket = algorithm.get("bracket")
    if bracket is None:
        bracket_literal = "None"
    else:
        bracket_literal = (
            "Some(ModelBracket {\n"
            f"{indent}    scheme: {rust_str(bracket['scheme'])},\n"
            f"{indent}    lower: {rust_f64(bracket['lower'])},\n"
            f"{indent}    upper: {rust_f64(bracket['upper'])},\n"
            f"{indent}    steps: {bracket['steps']},\n"
            f"{indent}}})"
        )
    initialisation = algorithm.get("initialisation")
    initialisation_literal = (
        "None" if initialisation is None else f"Some({rust_str(initialisation)})"
    )
    initial_temperature = algorithm.get("initial_temperature")
    initial_temperature_literal = (
        "None" if initial_temperature is None else f"Some({rust_f64(initial_temperature)})"
    )
    fallback = algorithm.get("fallback")
    if fallback is None:
        fallback_literal = "None"
    else:
        fallback_literal = (
            "Some(&ModelFallback {\n"
            f"{indent}    after: {fallback['after']},\n"
            f"{indent}    algorithm: &{fallback_ref},\n"
            f"{indent}}})"
        )
    return (
        "ModelAlgorithm {\n"
        f"{indent}    scheme: {rust_str(algorithm['scheme'])},\n"
        f"{indent}    convergence: {rust_str(algorithm['convergence'])},\n"
        f"{indent}    tolerance: {rust_f64(algorithm['tolerance'])},\n"
        f"{indent}    max_iterations: {algorithm['max_iterations']},\n"
        f"{indent}    bracket: {bracket_literal},\n"
        f"{indent}    initialisation: {initialisation_literal},\n"
        f"{indent}    initial_temperature: {initial_temperature_literal},\n"
        f"{indent}    inner: {inner_ref},\n"
        f"{indent}    fallback: {fallback_literal},\n"
        f"{indent}}}"
    )


def emit_rust_model(model: dict[str, Any]) -> str:
    """One registry entry, plus whatever static its algorithm needs.

    A `direct` model has no algorithm at all - the schema forbids one rather than
    allowing a vacuous one - so the algorithm static is emitted only for a
    `procedure`. That is the whole of the difference between the two kinds here; the
    checks and the cases are emitted the same way, because a direct model has bounds
    and cases like any other.
    """
    ident_ = ident(model["id"])
    kind = model.get("kind", "procedure")

    algorithm_static = ""
    algorithm_literal = "None"
    if kind == "procedure":
        algorithm = model["algorithm"]
        # A nested scheme becomes its own `static`, emitted before the one that points
        # at it. One level is all any spec here uses today; the recursion in
        # `emit_rust_algorithm` costs nothing and means a second level is a spec
        # change rather than a generator change.
        inner = algorithm.get("inner")
        fallback = algorithm.get("fallback")
        inner_static = ""
        fallback_static = ""
        inner_ref = "None"
        if inner is not None:
            inner_static = (
                f"static {ident_}_INNER: ModelAlgorithm = "
                f"{emit_rust_algorithm(inner, 'None', '    ')};\n\n"
            )
            inner_ref = f"Some(&{ident_}_INNER)"
        fallback_ref = "None"
        if fallback is not None:
            fallback_static = (
                f"static {ident_}_FALLBACK: ModelAlgorithm = "
                f"{emit_rust_algorithm(fallback['algorithm'], 'None', '    ')};\n\n"
            )
            fallback_ref = f"{ident_}_FALLBACK"
        outer = emit_rust_algorithm(algorithm, inner_ref, "    ", fallback_ref)
        algorithm_static = (
            inner_static
            + fallback_static
            + f"static {ident_}_ALGORITHM: ModelAlgorithm = {outer};\n\n"
        )
        algorithm_literal = f"Some(&{ident_}_ALGORITHM)"

    checks = "".join(
        f"        SpecCheck {{\n"
        f"            on_input: {str(check['quantity'] in model['inputs']).lower()},\n"
        f"{emit_range_check(check, model['id'])}"
        f"        }},\n"
        for check in model.get("valid_range", [])
    )
    cases = "".join(
        emit_test_case(
            test_id=case["id"],
            kind="case",
            status="active",
            tolerance=case["tolerance"],
            inputs=case["inputs"],
            expected=case["expected"],
            indent="        ",
        )
        for case in model["cases"]
    )
    return f"""
static {ident_}_CHECKS: &[SpecCheck] = &[
{checks}];

static {ident_}_CASES: &[TestCase] = &[
{cases}];

{algorithm_static}/// Registry entry for `{model["id"]}`.
pub static {ident_}_SPEC: ModelSpec = ModelSpec {{
    id: {rust_str(model["id"])},
    kind: {rust_str(kind)},
    algorithm: {algorithm_literal},
    checks: {ident_}_CHECKS,
    cases: {ident_}_CASES,
}};
"""


def emit_rust(models: list[dict[str, Any]], namespace: str) -> str:
    mine = [m for m in models if namespace_of(m["id"]) == namespace]
    body = "".join(emit_rust_model(m) for m in mine)

    # Only the items this namespace's models actually reference.
    #
    # The list used to be emitted whole, which is correct for `eos` - it has a bracketed
    # model, so `ModelBracket` is used - and wrong for the first namespace that has
    # none. `process` is that namespace: its `separator` names no bracket because it has
    # no scan to bracket, so the import was unused and `clippy -D warnings` failed on a
    # file nobody is allowed to edit. Emitting what the body names is the same rule the
    # rest of this generator follows, which is that the spec decides and the generator
    # reports.
    imports = [
        "Band",
        "ModelAlgorithm",
        "ModelBracket",
        "ModelFallback",
        "ModelSpec",
        "RangeCheck",
        "Severity",
        "SpecCheck",
        "TestCase",
        "WarningCode",
    ]
    banner = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n//!\n"
        "//! Generated by `tools/gen_models.py` from:\n"
        + "".join(f"//!   - {m['_path']}\n" for m in mine)
        + "//!\n"
        "//! Regenerate with `python tools/gen_models.py`; CI runs `--check` and fails\n"
        "//! on any difference.\n\n"
        "use azoth_core::{\n"
        f"    {', '.join(name for name in imports if name in body)},\n"
        "};\n"
    )
    entries = "".join(f"    &{ident(m['id'])}_SPEC,\n" for m in mine)
    lookup = f"""
static ALL_MODELS: &[&ModelSpec] = &[
{entries}];

/// Every model in this namespace, in id order.
#[must_use]
pub fn models() -> &'static [&'static ModelSpec] {{
    ALL_MODELS
}}

/// One model by id, or `None`. Used by the contract tests.
#[must_use]
pub fn model(id: &str) -> Option<&'static ModelSpec> {{
    ALL_MODELS.iter().copied().find(|m| m.id == id)
}}
"""
    return banner + body + lookup


def emit_python(models: list[dict[str, Any]]) -> str:
    header = (
        '"""GENERATED FILE - DO NOT EDIT BY HAND.\n\n'
        "Generated by `tools/gen_models.py` from:\n"
        + "".join(f"  - {m['_path']}\n" for m in models)
        + '"""\n\n'
        "from __future__ import annotations\n\n"
        "from typing import Any, Final\n\n"
        "MODELS: Final[tuple[dict[str, Any], ...]] = (\n"
    )
    entries = "".join(f"    {m!r},\n" for m in models)
    lookups = """

MODEL_BY_ID: Final[dict[str, dict[str, Any]]] = {m["id"]: m for m in MODELS}


def model(model_id: str) -> dict[str, Any]:
    \"\"\"One model's spec, by id.

    Raises:
        KeyError: for an unknown id, naming what is known - the same rule
            `_registry_gen.spec` follows, so a typo fails loudly rather than
            returning something plausible.
    \"\"\"
    try:
        return MODEL_BY_ID[model_id]
    except KeyError:
        raise KeyError(
            f"unknown model {model_id!r}; known models: {sorted(MODEL_BY_ID)}"
        ) from None
"""
    return header + entries + ")\n" + lookups


def main() -> None:
    check = "--check" in sys.argv
    models = load_models()

    outputs: dict[Path, str] = {PY_OUT: emit_python(models)}
    namespaces = sorted({namespace_of(m["id"]) for m in models})
    for namespace in namespaces:
        outputs[rust_out_for(namespace)] = emit_rust(models, namespace)

    stale: list[Path] = []
    for path, text in outputs.items():
        rendered = text
        if path.suffix == ".rs":
            proc = subprocess.run(
                ["rustfmt", "--edition", "2024"],
                input=text,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode == 0:
                rendered = proc.stdout
        existing = path.read_text(encoding="utf-8") if path.exists() else None
        if existing == rendered:
            continue
        if check:
            stale.append(path)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(rendered, encoding="utf-8")
            print(f"gen_models: wrote {path.relative_to(ROOT)}")

    if check and stale:
        for path in stale:
            print(f"gen_models: {path.relative_to(ROOT)} is out of date")
        sys.exit("gen_models: generated files are stale; run tools/gen_models.py")
    print(f"gen_models: {len(models)} model(s)" + (", generated files up to date" if check else ""))


if __name__ == "__main__":
    main()
