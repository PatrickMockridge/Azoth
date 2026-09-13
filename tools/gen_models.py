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

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gen_registry import (
    ROOT,
    emit_range_check,
    emit_test_case,
    rust_f64,
    rust_str,
)

MODEL_DIR = ROOT / "specs" / "models"
SCHEMA_PATH = ROOT / "specs" / "schema" / "model.schema.json"
CALC_SCHEMA_PATH = ROOT / "specs" / "schema" / "calc.schema.json"


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
    paths = sorted(MODEL_DIR.rglob("*.yaml"))
    if not paths:
        return []

    try:
        from jsonschema import Draft202012Validator
        from referencing import Registry, Resource
    except ImportError:  # pragma: no cover
        sys.exit("gen_models requires jsonschema and referencing")

    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    calc = json.loads(CALC_SCHEMA_PATH.read_text(encoding="utf-8"))
    registry = Registry().with_resource(calc["$id"], Resource.from_contents(calc))
    validator = Draft202012Validator(schema, registry=registry)

    models: list[dict[str, Any]] = []
    for path in paths:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
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
    return models


def emit_rust_model(model: dict[str, Any]) -> str:
    ident_ = ident(model["id"])
    algorithm = model["algorithm"]
    bracket = algorithm["bracket"]

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

/// Registry entry for `{model["id"]}`.
pub static {ident_}_SPEC: ModelSpec = ModelSpec {{
    id: {rust_str(model["id"])},
    verification: {rust_str(model["verification"]["status"])},
    algorithm: ModelAlgorithm {{
        scheme: {rust_str(algorithm["scheme"])},
        convergence: {rust_str(algorithm["convergence"])},
        tolerance: {rust_f64(algorithm["tolerance"])},
        max_iterations: {algorithm["max_iterations"]},
        bracket: ModelBracket {{
            scheme: {rust_str(bracket["scheme"])},
            lower: {rust_f64(bracket["lower"])},
            upper: {rust_f64(bracket["upper"])},
            steps: {bracket["steps"]},
        }},
    }},
    checks: {ident_}_CHECKS,
    cases: {ident_}_CASES,
}};
"""


def emit_rust(models: list[dict[str, Any]], namespace: str) -> str:
    mine = [m for m in models if namespace_of(m["id"]) == namespace]
    banner = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n//!\n"
        "//! Generated by `tools/gen_models.py` from:\n"
        + "".join(f"//!   - {m['_path']}\n" for m in mine)
        + "//!\n"
        "//! Regenerate with `python tools/gen_models.py`; CI runs `--check` and fails\n"
        "//! on any difference.\n\n"
        "use azoth_core::{\n"
        "    Band, ModelAlgorithm, ModelBracket, ModelSpec, RangeCheck, Severity, SpecCheck,\n"
        "    TestCase, WarningCode,\n"
        "};\n"
    )
    body = "".join(emit_rust_model(m) for m in mine)
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
