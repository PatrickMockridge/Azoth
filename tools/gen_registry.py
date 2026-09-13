#!/usr/bin/env python3
"""Generate language-level registries from the calc specs.

Reads every spec under specs/calcs/ and emits:

  crates/azoth-<namespace>/src/spec_gen.rs   one file per namespace, as Rust statics
  python/src/azoth/_registry_gen.py          every spec, as Python data

The Rust side is split by namespace because a crate is the unit of compilation: a
calc reads its bounds from a table in its own crate, and a crate must not have to
depend on another namespace's crate to see its own spec. The types those tables are
built from are shared, and live in `azoth_core::spec`.

The Python side is one file. It all ships as one package, and `azoth.hydraulics`
and `azoth.thermal` are subpackages of a single distribution, so there is nothing
for a split to buy.

Both outputs are committed and drift-checked in CI, so the specs are genuinely
the single source of truth rather than a document that is supposed to match the
code.

Why generate Rust rather than have Rust parse the YAML: parsing would mean
depending on a YAML crate (serde_yaml is deprecated, and serde_yml is an
unrelated low-trust fork), and it would make the spec's contents available only
at runtime, so a malformed bound would be a runtime error instead of a build
failure. Generating a `&'static` table moves that to compile time and keeps the
dependency tree free of a parser that exists for one purpose.

Design note on the generated range checks: they are emitted as full struct
literals rather than through the convenience constructors in azoth-core. The
constructors assume inclusive bounds, so a spec with, say, an exclusive lower
and inclusive upper bound would be silently generated wrong. Struct literals
cannot be wrong that way, and a reviewer reading the generated file sees the
exact bound semantics rather than having to recall what a constructor does.

Usage:
    python tools/gen_registry.py            # write the files
    python tools/gen_registry.py --check    # fail if they are out of date
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("gen_registry requires PyYAML: pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"
CRATES_DIR = ROOT / "crates"
PY_OUT = ROOT / "python" / "src" / "azoth" / "_registry_gen.py"

GENERATED_BANNER = "GENERATED FILE - DO NOT EDIT BY HAND."


def namespace_of(spec: dict[str, Any]) -> str:
    """The namespace a spec belongs to, which is its id's first segment."""
    namespace: str = spec["id"].split(".")[0]
    return namespace


def crate_for(namespace: str) -> str:
    """The crate that owns a namespace's generated tables.

    Derived rather than looked up in a hand-maintained map, so that adding a
    namespace stays "add a directory under specs/calcs/" - the same rule the rest
    of the tooling follows. A namespace whose crate does not exist is caught by
    `rust_out_for` rather than silently writing tables into some default crate.
    """
    return f"azoth-{namespace.replace('_', '-')}"


def rust_out_for(namespace: str) -> Path:
    """Where a namespace's generated Rust tables live.

    Fails rather than falling back. The fallback would be writing one namespace's
    specs into another namespace's crate, which compiles, passes the drift check
    against itself, and puts a calc's bounds in a crate that calc does not depend
    on - so the file would simply never be read.
    """
    crate = crate_for(namespace)
    source_dir = CRATES_DIR / crate / "src"
    if not source_dir.is_dir():
        sys.exit(
            f"gen_registry: specs/calcs/{namespace}/ exists but crates/{crate}/src/ "
            f"does not. A namespace needs a crate to generate into; either create "
            f"crates/{crate}/, or the namespace is misspelled."
        )
    return source_dir / "spec_gen.rs"


def rustfmt(source: str) -> str:
    """Format generated Rust with rustfmt.

    The output is committed and CI runs `cargo fmt --check` over the whole
    workspace, so generated code has to be formatted like hand-written code.
    Formatting here rather than emitting carefully aligned strings means the
    generator never has to know rustfmt's rules - and if rustfmt is unavailable
    the source is still valid, just unformatted, which `cargo fmt --check` will
    report rather than something silently drifting.
    """
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024", "--emit", "stdout"],
            input=source,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        print("gen_registry: rustfmt not found; emitting unformatted Rust", file=sys.stderr)
        return source
    if proc.returncode != 0 or not proc.stdout.strip():
        print(
            f"gen_registry: rustfmt failed, emitting unformatted Rust\n{proc.stderr}",
            file=sys.stderr,
        )
        return source
    return proc.stdout


def rust_str(value: str) -> str:
    """A Rust string literal for a bounded piece of prose.

    Whitespace runs are collapsed to single spaces. Spec rationales are YAML
    folded scalars, so they arrive with embedded newlines; those are a YAML
    formatting artefact, not meaning, and collapsing them keeps the generated
    file readable and the message rendering identical on both sides.
    """
    flat = " ".join(value.split())
    escaped = flat.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def rust_f64(value: float) -> str:
    """A Rust `f64` literal.

    The `float()` coercion is load-bearing, not defensive. YAML parses a bare
    `0` as an integer, and Python's repr then produces `0`, which is an integer
    literal in Rust and will not compile where an `f64` is expected. Forcing the
    type first means every emitted literal carries a '.' or an 'e'.

    repr() round-trips exactly, so no precision is lost in translation.
    """
    return repr(float(value))


def rust_opt_f64(value: float | None) -> str:
    return "None" if value is None else f"Some({rust_f64(value)})"


def rust_slice_str(items: list[str]) -> str:
    inner = ", ".join(rust_str(i) for i in items)
    return f"&[{inner}]"


def camel_variant(screaming_snake: str) -> str:
    """`TRANSITIONAL_FLOW` -> `TransitionalFlow`, matching the Rust enum."""
    return "".join(part.capitalize() for part in screaming_snake.split("_"))


def collect_numbers(mapping: dict[str, Any]) -> list[tuple[str, float]]:
    return [(k, v) for k, v in mapping.items() if isinstance(v, (int, float))]


def collect_lists(mapping: dict[str, Any]) -> list[tuple[str, list[str]]]:
    return [(k, v) for k, v in mapping.items() if isinstance(v, list)]


def emit_range_check(check: dict[str, Any], spec_id: str) -> str:
    severity = "Severity::Error" if check["severity"] == "error" else "Severity::Warning"
    band = "Band::Inside" if check.get("when", "outside") == "inside" else "Band::Outside"
    # Warning codes are emitted by their Rust variant name. The schema restricts
    # `code` to the known set, and spec_lint checks the same strings against the
    # Python enum, so a typo here cannot reach the generated file.
    code = "WarningCode::" + camel_variant(check.get("code", "OUT_OF_VALID_RANGE"))
    rationale = check.get("rationale", "")
    if not rationale:
        raise SystemExit(
            f"{spec_id}: range check on '{check['quantity']}' has no rationale. "
            f"Every bound must explain why it exists - a bound with no reason is a "
            f"bound nobody dares change."
        )
    return (
        "            check: RangeCheck {\n"
        f"                quantity: {rust_str(check['quantity'])},\n"
        f"                min: {rust_opt_f64(check.get('min'))},\n"
        f"                min_inclusive: {str(check.get('min_inclusive', True)).lower()},\n"
        f"                max: {rust_opt_f64(check.get('max'))},\n"
        f"                max_inclusive: {str(check.get('max_inclusive', True)).lower()},\n"
        f"                band: {band},\n"
        f"                severity: {severity},\n"
        f"                code: {code},\n"
        f"                rationale: {rust_str(rationale)},\n"
        "            },\n"
    )


def emit_test_case(
    test_id: str,
    kind: str,
    status: str,
    tolerance: float,
    inputs: dict[str, Any],
    expected: dict[str, Any],
    indent: str,
    property_name: str | None = None,
    skip_reason: str | None = None,
) -> str:
    numbers = collect_numbers(inputs)
    lists = collect_lists(inputs)
    expected_numbers = collect_numbers(expected)

    def pairs(items: list[tuple[str, float]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {rust_f64(v)})" for k, v in items)
        return f"&[{inner}]"

    def list_pairs(items: list[tuple[str, list[str]]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {rust_slice_str(v)})" for k, v in items)
        return f"&[{inner}]"

    prop = f"Some({rust_str(property_name)})" if property_name else "None"
    skip = f"Some({rust_str(skip_reason)})" if skip_reason else "None"

    return (
        f"{indent}TestCase {{\n"
        f"{indent}    id: {rust_str(test_id)},\n"
        f"{indent}    kind: {rust_str(kind)},\n"
        f"{indent}    property: {prop},\n"
        f"{indent}    status: {rust_str(status)},\n"
        f"{indent}    skip_reason: {skip},\n"
        f"{indent}    tolerance: {rust_f64(tolerance)},\n"
        f"{indent}    numbers: {pairs(numbers)},\n"
        f"{indent}    lists: {list_pairs(lists)},\n"
        f"{indent}    expected: {pairs(expected_numbers)},\n"
        f"{indent}}},\n"
    )


def emit_rust(specs: list[dict[str, Any]], source_files: list[str], namespace: str) -> str:
    out: list[str] = []
    out.append(
        f"""//! {GENERATED_BANNER}
//!
//! Generated by `tools/gen_registry.py` from:
//!
"""
    )
    for src in source_files:
        out.append(f"//!   - {src}\n")
    out.append(
        f"""//!
//! Tables for the `{namespace}` namespace. Every namespace has its own generated
//! file, because a crate is the unit of compilation and a calculation must be able
//! to read its own bounds without depending on another namespace's crate. The
//! types these tables are built from are shared - see `azoth_core::spec`.
//!
//! These tables are what make the specs authoritative at runtime rather than
//! merely descriptive. Each calc reads its own range checks from here, so a
//! bound changed in a spec file changes the code's behaviour with no second
//! edit - and `cargo test` fails if the two ever disagree.
"""
    )
    out.append(
        """
use azoth_core::{
    Band, CalcSpec, RangeCheck, Severity, SpecCheck, SolverSpec, TestCase, WarningCode,
};
"""
    )

    # Per-calc statics.
    for spec in specs:
        ident = spec["id"].split(".")[-1].upper()
        inputs = spec["inputs"]

        checks: list[str] = []
        for check in spec["valid_range"]:
            on_input = check["quantity"] in inputs
            checks.append("        SpecCheck {\n")
            checks.append(f"            on_input: {str(on_input).lower()},\n")
            checks.append(emit_range_check(check, spec["id"]))
            checks.append("        },\n")

        solver = spec.get("solver")
        solver_str = "None"
        if solver:
            solver_str = (
                "Some(SolverSpec {\n"
                f"            kind: {rust_str(solver['kind'])},\n"
                f"            tolerance: {rust_f64(solver['tolerance'])},\n"
                f"            max_iterations: {solver['max_iterations']},\n"
                f"            initial_guess: {rust_f64(solver['initial_guess'])},\n"
                f"            convergence: {rust_str(solver['convergence'])},\n"
                "        })"
            )

        example = spec["worked_example"]
        out.append(f"/// Registry entry for `{spec['id']}`.\n")
        out.append(f"static {ident}_CHECKS: &[SpecCheck] = &[\n{''.join(checks)}];\n\n")

        others = [t for t in spec["tests"] if t["type"] != "worked_example"]
        out.append(f"static {ident}_TESTS: &[TestCase] = &[\n")
        for test in others:
            out.append(
                emit_test_case(
                    test_id=test["id"],
                    kind=test["type"],
                    status=test["status"],
                    tolerance=test.get("tolerance", example["tolerance"]),
                    inputs=test.get("inputs", {}),
                    expected=test.get("expected", {}),
                    indent="    ",
                    property_name=test.get("property"),
                    skip_reason=test.get("skip_reason"),
                )
            )
        out.append("];\n\n")

        out.append(
            f"""/// Registered spec for `{spec["id"]}`.
///
/// Public and addressable directly, so a calc can hold `&{ident}_SPEC` with no
/// lookup and no failure path. A calc whose spec is missing is a build-time
/// invariant, not a runtime condition, and this shape makes it unrepresentable
/// rather than something to handle.
pub static {ident}_SPEC: CalcSpec = CalcSpec {{
    id: {rust_str(spec["id"])},
    verification: {rust_str(spec["verification"]["status"])},
    checks: {ident}_CHECKS,
    solver: {solver_str},
    worked_example: {
                emit_test_case(
                    # The spec's own worked-example test id, not a synthesized one, so the two
                    # languages report the same name for the same case.
                    test_id=next(
                        (t["id"] for t in spec["tests"] if t["type"] == "worked_example"),
                        spec["id"].split(".")[-1] + "_worked_example",
                    ),
                    kind="worked_example",
                    status=next(
                        (t["status"] for t in spec["tests"] if t["type"] == "worked_example"),
                        "active",
                    ),
                    tolerance=example["tolerance"],
                    inputs=example["inputs"],
                    expected=example["expected"],
                    indent="        ",
                    # strip() leaves the trailing comma that emit_test_case adds for list
                    # elements; the template supplies its own, and two in a row will not
                    # compile.
                )
                .strip()
                .removesuffix(",")
            },
    tests: {ident}_TESTS,
}};

"""
        )

    names = ", ".join(f"&{s['id'].split('.')[-1].upper()}_SPEC" for s in specs)
    out.append(
        f"""/// Every calculation in the registry, sorted by id.
static ALL_SPECS: &[&CalcSpec] = &[{names}];

/// All specs, in a stable order.
#[must_use]
pub fn specs() -> &'static [&'static CalcSpec] {{
    ALL_SPECS
}}

/// Look up a spec by id.
#[must_use]
pub fn spec(id: &str) -> Option<&'static CalcSpec> {{
    ALL_SPECS.iter().copied().find(|s| s.id == id)
}}
"""
    )
    return "".join(out)


def emit_python(specs: list[dict[str, Any]], source_files: list[str]) -> str:
    out: list[str] = []
    out.append(
        f'''"""{GENERATED_BANNER}

Generated by `tools/gen_registry.py` from:
'''
    )
    for src in source_files:
        out.append(f"  - {src}\n")
    out.append(
        '''"""

from __future__ import annotations

from typing import Any, Final

#: Every calc spec, as plain data, in a stable order.
CALCS: Final[tuple[dict[str, Any], ...]] = (
'''
    )
    for spec in specs:
        out.append("    {\n")
        out.append(f'        "id": {spec["id"]!r},\n')
        out.append(f'        "name": {spec["name"]!r},\n')
        out.append(f'        "equation": {spec["equation"]!r},\n')
        out.append(f'        "verification": {spec["verification"]["status"]!r},\n')
        out.append(f'        "source": {spec["source"]!r},\n')
        out.append(f'        "inputs": { {k: dict(v) for k, v in spec["inputs"].items()}!r},\n')
        out.append(f'        "outputs": { {k: dict(v) for k, v in spec["outputs"].items()}!r},\n')
        out.append(f'        "valid_range": {spec["valid_range"]!r},\n')
        out.append(f'        "assumptions": {spec["assumptions"]!r},\n')
        out.append(f'        "references": {spec["references"]!r},\n')
        out.append(f'        "solver": {spec.get("solver")!r},\n')
        out.append(f'        "data": {spec.get("data")!r},\n')
        out.append(f'        "worked_example": {spec["worked_example"]!r},\n')
        out.append(f'        "implementations": {spec["implementations"]!r},\n')
        out.append(f'        "tests": {spec["tests"]!r},\n')
        out.append("    },\n")
    out.append(
        """)


#: Lookup by calc id.
BY_ID: Final[dict[str, dict[str, Any]]] = {c["id"]: c for c in CALCS}


def spec(calc_id: str) -> dict[str, Any]:
    \"\"\"Fetch a spec by id, raising a clear error if it is not a calc we ship.\"\"\"
    try:
        return BY_ID[calc_id]
    except KeyError:
        raise KeyError(
            f"unknown calc {calc_id!r}; known ids: {sorted(BY_ID)}"
        ) from None
"""
    )
    return "".join(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; fail if the generated files are out of date",
    )
    args = parser.parse_args()

    paths = sorted(SPEC_DIR.rglob("*.yaml"))
    if not paths:
        sys.exit(f"gen_registry: no specs found under {SPEC_DIR}")

    # Paired and sorted together. The specs used to be sorted by id while the
    # source paths kept directory order - which happened to agree, because a
    # namespace's directory name is also its id prefix, so the two orders always
    # matched. The lists are consumed in parallel and are now split per namespace,
    # so an accidental disagreement would attribute one namespace's spec files to
    # another's tables and be very hard to see.
    loaded: list[tuple[str, dict[str, Any]]] = [
        (str(p.relative_to(ROOT)), yaml.safe_load(p.read_text(encoding="utf-8"))) for p in paths
    ]
    loaded.sort(key=lambda pair: pair[1]["id"])

    # Guard against two specs claiming the same id: the generated tables key on
    # it, and a duplicate would silently drop one calc from the registry.
    ids = [spec["id"] for _, spec in loaded]
    duplicates = {i for i in ids if ids.count(i) > 1}
    if duplicates:
        sys.exit(f"gen_registry: duplicate calc id(s): {sorted(duplicates)}")

    by_namespace: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for source_file, spec in loaded:
        by_namespace.setdefault(namespace_of(spec), []).append((source_file, spec))

    outputs: dict[Path, str] = {
        PY_OUT: emit_python([spec for _, spec in loaded], [src for src, _ in loaded])
    }
    for namespace in sorted(by_namespace):
        entries = by_namespace[namespace]
        outputs[rust_out_for(namespace)] = rustfmt(
            emit_rust([spec for _, spec in entries], [src for src, _ in entries], namespace)
        )

    stale: list[Path] = []
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if args.check:
            current = path.read_text(encoding="utf-8") if path.exists() else ""
            if current != content:
                stale.append(path)
        else:
            path.write_text(content, encoding="utf-8")

    if args.check:
        if stale:
            for path in stale:
                print(f"gen_registry: {path.relative_to(ROOT)} is out of date", file=sys.stderr)
            print("Run `python tools/gen_registry.py` to regenerate.", file=sys.stderr)
            return 1
        print(f"gen_registry: {len(loaded)} spec(s), generated files up to date")
        return 0

    for path in outputs:
        print(f"gen_registry: wrote {path.relative_to(ROOT)}")
    print(f"gen_registry: {len(loaded)} spec(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
